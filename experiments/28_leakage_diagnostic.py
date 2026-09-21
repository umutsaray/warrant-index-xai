"""Forensic patient-overlap diagnostic for the three frozen TCGA models.

This script is deliberately read-only with respect to data/model artifacts. It
loads the existing models, reproduces their stored test C-indices, and scores
the patient-disjoint subset. It never fits or modifies a model.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored
import torch


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results" / "leakage_diagnostic"
REPORT = ROOT / "LEAKAGE_DIAGNOSTIC.md"


def frame(name: str) -> pd.DataFrame:
    return pd.read_parquet(SOURCE / name).set_index("sample_key")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def iso_mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat().replace("+00:00", "Z")


def cindex(y: pd.DataFrame, risk: np.ndarray, mask: np.ndarray) -> float:
    return float(
        concordance_index_censored(
            y.event.to_numpy()[mask].astype(bool), y.time.to_numpy()[mask], np.asarray(risk)[mask]
        )[0]
    )


def load_deepsurv_module():
    path = SOURCE / "train_explainable_dl.py"
    spec = importlib.util.spec_from_file_location("frozen_deepsurv_training_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    x_train = frame("X_train_selected.parquet")
    x_test = frame("X_test_selected.parquet")
    y_train = frame("y_train.parquet").loc[x_train.index]
    y_test = frame("y_test.parquet").loc[x_test.index]
    assert x_train.shape == (673, 616) and x_test.shape == (169, 616)
    assert y_train.shape == (673, 2) and y_test.shape == (169, 2)
    assert x_train.index.is_unique and x_test.index.is_unique
    assert all(str(value).startswith("TCGA-") and len(str(value)) >= 12 for value in [*x_train.index, *x_test.index])

    # TCGA participant barcodes are the first 12 characters; the stored
    # sample_key is the first 15 characters and therefore preserves sample type.
    train_patients = x_train.index.astype(str).str[:12]
    test_patients = x_test.index.astype(str).str[:12]
    train_patient_set = set(train_patients)
    test_patient_set = set(test_patients)
    overlap = train_patient_set & test_patient_set
    disjoint_mask = np.array([patient not in train_patient_set for patient in test_patients], dtype=bool)
    all_mask = np.ones(len(x_test), dtype=bool)
    counts = {
        "train_rows": len(x_train),
        "test_rows": len(x_test),
        "train_unique_patients": len(train_patient_set),
        "test_unique_patients": len(test_patient_set),
        "overlap_unique_patients": len(overlap),
        "overlap_train_rows": int(sum(patient in overlap for patient in train_patients)),
        "overlap_test_rows": int(sum(patient in overlap for patient in test_patients)),
        "disjoint_test_rows": int(disjoint_mask.sum()),
        "disjoint_test_unique_patients": len(test_patient_set - train_patient_set),
        "original_test_events": int(y_test.event.sum()),
        "disjoint_test_events": int(y_test.event.to_numpy()[disjoint_mask].sum()),
    }
    assert counts == {
        "train_rows": 673,
        "test_rows": 169,
        "train_unique_patients": 628,
        "test_unique_patients": 164,
        "overlap_unique_patients": 23,
        "overlap_train_rows": 25,
        "overlap_test_rows": 23,
        "disjoint_test_rows": 146,
        "disjoint_test_unique_patients": 141,
        "original_test_events": 26,
        "disjoint_test_events": 19,
    }

    # Frozen CoxPH: same scaler/features/model used by train_baselines.py.
    cox_artifact = joblib.load(SOURCE / "coxph_baseline.joblib")
    cox_values = cox_artifact["scaler"].transform(x_test[cox_artifact["features"]].astype(float))
    cox_risk = cox_artifact["model"].predict_partial_hazard(
        pd.DataFrame(cox_values, columns=cox_artifact["features"])
    ).to_numpy().reshape(-1)

    # Frozen RSF: no refit; prediction only.
    rsf_artifact = joblib.load(SOURCE / "rsf_baseline.joblib")
    rsf_risk = rsf_artifact["model"].predict(
        x_test[rsf_artifact["features"]].to_numpy(dtype=np.float32)
    )

    # Frozen DeepSurv checkpoint: reconstruct architecture and load state_dict.
    deep_module = load_deepsurv_module()
    checkpoint = torch.load(SOURCE / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    deep_model = deep_module.OmicsDeepSurv(checkpoint["block_dims"], checkpoint["dropout"])
    deep_model.load_state_dict(checkpoint["state_dict"])
    deep_model.eval()
    deep_values = x_test[checkpoint["ordered_columns"]].to_numpy(dtype=np.float32)
    deep_values = (deep_values - checkpoint["means"]) / checkpoint["scales"]
    deep_risk = deep_module.predict_risk(
        deep_model,
        deep_module.make_blocks(deep_values, checkpoint["slices"], torch.device("cpu")),
    )

    baseline = json.loads((SOURCE / "baseline_results.json").read_text(encoding="utf-8"))
    reported = {
        "CoxPH": float(baseline["models"]["coxph"]["test_c_index"]),
        "RSF": float(baseline["models"]["random_survival_forest"]["test_c_index"]),
        "DeepSurv": float(checkpoint["test_c_index"]),
    }
    risks = {"CoxPH": cox_risk, "RSF": rsf_risk, "DeepSurv": deep_risk}
    scores = {}
    for model, risk in risks.items():
        original = cindex(y_test, risk, all_mask)
        disjoint = cindex(y_test, risk, disjoint_mask)
        assert abs(original - reported[model]) < 1e-12
        scores[model] = {
            "reported": reported[model],
            "reproduced_original_169": original,
            "disjoint_146": disjoint,
            "delta_disjoint_minus_original": disjoint - original,
        }

    # DeepSurv early-stopping split was also row-level. Diagnose without fitting.
    from sklearn.model_selection import train_test_split

    fit_idx, val_idx = train_test_split(
        np.arange(len(x_train)), test_size=0.15, random_state=42, stratify=y_train.event.to_numpy()
    )
    fit_patients = set(train_patients[fit_idx])
    val_patients = set(train_patients[val_idx])
    internal_overlap = fit_patients & val_patients
    deep_internal = {
        "fit_rows": len(fit_idx),
        "validation_rows": len(val_idx),
        "overlap_unique_patients": len(internal_overlap),
        "overlap_fit_rows": int(sum(patient in internal_overlap for patient in train_patients[fit_idx])),
        "overlap_validation_rows": int(sum(patient in internal_overlap for patient in train_patients[val_idx])),
    }

    input_names = [
        "X_train_selected.parquet", "X_test_selected.parquet", "y_train.parquet", "y_test.parquet",
        "coxph_baseline.joblib", "rsf_baseline.joblib", "explainable_dl_model.pt",
        "impute_and_split.py", "select_features.py", "train_baselines.py", "train_explainable_dl.py",
    ]
    inventory = {
        name: {
            "path": str((SOURCE / name).resolve()),
            "bytes": (SOURCE / name).stat().st_size,
            "sha256": sha256(SOURCE / name),
            "mtime_utc": iso_mtime(SOURCE / name),
        }
        for name in input_names
    }

    # Cross-check the hashes frozen into prior Warrant Index artifacts.
    phase3 = json.loads((ROOT / "results/phase3_brca/metadata.json").read_text(encoding="utf-8"))
    phase6 = json.loads((ROOT / "results/phase6_deepsurv_pilot/metadata.json").read_text(encoding="utf-8"))
    for name in ["X_train_selected.parquet", "X_test_selected.parquet", "y_train.parquet", "y_test.parquet", "coxph_baseline.joblib"]:
        assert inventory[name]["sha256"] == phase3["input_sha256"][name]
    for name in ["X_train_selected.parquet", "X_test_selected.parquet", "y_train.parquet", "y_test.parquet", "explainable_dl_model.pt"]:
        assert inventory[name]["sha256"] == phase6["input_sha256"][name]

    payload = {
        "scope": "read-only forensic diagnostic; no fitting or model modification",
        "verdict": "Models were trained and tested on a sample-disjoint but patient-overlapping (dirty) split.",
        "patient_key_rule": "first 12 characters of the 15-character stored TCGA sample_key",
        "counts": counts,
        "scores": scores,
        "deepsurv_internal_validation": deep_internal,
        "inventory": inventory,
        "hash_cross_checks_passed": True,
        "patient_or_sample_identifiers_exported": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "diagnostic.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    src = str(SOURCE.resolve()).replace("\\", "/")
    wr = str(ROOT.resolve()).replace("\\", "/")
    report = f"""# Hasta-düzeyi veri sızıntısı teşhisi

**Tek cümlelik sonuç: Modeller örnek-düzeyinde ayrık fakat hasta-düzeyinde çakışmalı (kirli) split üzerinde eğitildi ve test edildi.**

Bu teşhis salt okunur yapıldı. Hiçbir model eğitilmedi, hiçbir artefaktın ağırlığı değiştirilmedi ve WI/certification koduna dokunulmadı. Raporda gerçek barcode veya hasta/örnek kimliği yoktur.

## ADIM 0 — Dosya ve kod zinciri

### Kullanılan veri ve model dosyaları

Üç model aynı seçilmiş matris çiftini kullandı:

- Train özellikleri: `{src}/X_train_selected.parquet` — 673 satır, 616 özellik.
- Test özellikleri: `{src}/X_test_selected.parquet` — 169 satır, 616 özellik.
- Train sonuçları: `{src}/y_train.parquet` — 673 satır.
- Test sonuçları: `{src}/y_test.parquet` — 169 satır.
- Modeller: `{src}/coxph_baseline.joblib`, `{src}/rsf_baseline.joblib`, `{src}/explainable_dl_model.pt`.

`{src}/train_baselines.py:46-49` bu dört parquet dosyasını doğrudan okur. CoxPH train fit'i `{src}/train_baselines.py:64-81`, 169 satırlık test prediction/C-index hesabı `{src}/train_baselines.py:85-92` içindedir. RSF fit'i `{src}/train_baselines.py:97-108`, aynı test dosyasında prediction/C-index hesabı `{src}/train_baselines.py:109-113` içindedir. Sonuçlar `{src}/baseline_results.json:12` ve `{src}/baseline_results.json:27` satırlarına yazılmıştır.

DeepSurv aynı dört dosyayı `{src}/train_explainable_dl.py:144-147` satırlarında okur; test tensörünü `{src}/train_explainable_dl.py:175-181` arasında kurar, mevcut test riskini ve C-index'ini `{src}/train_explainable_dl.py:227-228` satırlarında hesaplar ve checkpoint'e `{src}/train_explainable_dl.py:232-246` arasında kaydeder. Checkpoint içindeki test C-index 0.7456922162804516'dır.

Dosya hashleri daha sonraki audit metadata'sıyla birebir eşleşir: `{wr}/results/phase3_brca/metadata.json` içindeki `input_sha256` CoxPH/veri dosyalarını; `{wr}/results/phase6_deepsurv_pilot/metadata.json` içindeki `input_sha256` DeepSurv/veri dosyalarını doğrular. Tam hash ve zaman damgaları `{wr}/results/leakage_diagnostic/diagnostic.json` içindedir.

### “23 hasta” kaynağı

Audit script'i aynı dosya çiftini okur: `{wr}/experiments/06_brca_real.py:225-228`. Hasta anahtarını indeksin ilk 12 karakteriyle kurar (`:238-239`), train'de bulunan hastaların test satırlarını atar (`:240-241`) ve kalan test havuzunda hasta başına tek satır tutar (`:242-243`). Yani “23” farklı veya önceden temizlenmiş bir dosya çiftinden gelmez; modellerin kullandığı aynı 673/169 dosyadan gelir. Bu filtre model eğitimine/test C-index hesabına geri uygulanmamış, yalnız sonraki XAI audit query havuzuna uygulanmıştır.

## ADIM 1 — Bağımsız hasta çakışması hesabı

Kaynak üretim kodu `{src}/prepare_dataset.py:19,29-31` ile `sample_key` değerini TCGA barcode'unun ilk **15** karakteri olarak tanımlar. `{src}/impute_and_split.py:203-213` yalnız bu 15-karakter örnek anahtarının benzersizliğini kontrol eder; `{src}/impute_and_split.py:221-230` aynı örnek anahtarlarını satır düzeyinde train/test'e böler ve yalnız tam örnek-anahtarı kesişimini yasaklar. TCGA hasta düzeyi ilk 12 karakterdir. Bu nedenle farklı örnek türlerine ait 15-karakter anahtarlar aynı hastayı iki split'e taşıyabilir.

Orijinal dosyalardan sıfırdan yapılan ilk-12-karakter karşılaştırması:

| Ölçü | Sayı |
|---|---:|
| Train satırı | 673 |
| Test satırı | 169 |
| Train benzersiz hasta | 628 |
| Test benzersiz hasta | 164 |
| Her iki split'te bulunan benzersiz hasta | **23** |
| Ortak hastalara ait train satırı | 25 |
| Ortak hastalara ait test satırı | 23 |
| Ortak hastalar çıkarılınca test satırı | 146 |
| Disjoint test benzersiz hasta | 141 |
| Orijinal / disjoint test olayı | 26 / 19 |

Bağımsız sonuç önceki 23 sayısıyla tam eşleşir. Audit'teki 141 sayısı, 23 ortak-hasta test satırının çıkarılmasından sonra kalan 146 satırda hasta başına tek örnek tutulmasıyla oluşur.

## ADIM 2 — İşlem sırası ve kirli/temiz karar

`{src}/select_features.py:173-177,215-216` sample-level split'ten gelen X_train/X_test satırlarını koruyarak 616-feature dosyalarını yazar. CoxPH ve RSF bunları doğrudan `{src}/train_baselines.py:46-49` üzerinden; DeepSurv `{src}/train_explainable_dl.py:144-147` üzerinden kullanır. Eğitim kodlarında ilk-12-karakter hasta filtresi yoktur.

Dosya zamanları bu sıralamayla uyumludur: seçilmiş train/test dosyaları 2026-08-27 21:12:32 UTC, saklı CoxPH/RSF/DeepSurv modelleri sırasıyla 21:19:26, 21:19:30 ve 21:19:40 UTC'dir. Hasta-düzeyi audit filtresi `{wr}/experiments/06_brca_real.py:238-243` içinde yer alır; script ve Faz 3 metadata'sı 2026-09-17 tarihlidir. Kardeş dizin bir Git deposu değildir; commit geçmişi yoktur, bu nedenle kanıt çalıştırma kodu, SHA256 eşleşmeleri ve dosya zamanlarıyla sınırlıdır.

Ek teşhis: DeepSurv erken-durdurma ayrımı da train satırları üzerinde `{src}/train_explainable_dl.py:166-172` ile yapılır. Aynı seed ve stratifikasyon yeniden kurulduğunda 572 fit ve 101 validation satırı arasında 12 ortak hasta vardır; bu hastalar her tarafta 12 satıra karşılık gelir. Bu da hasta-gruplu bir validation ayrımı değildir.

## ADIM 3 — Mevcut modellerde disjoint test alt-kümesi

Hiçbir model yeniden fit edilmedi. Saklı scaler/ağırlıklar aynen yüklenerek önce 169 satırlık orijinal test C-index'leri yeniden üretildi; ardından yalnız 23 ortak-hasta test satırı çıkarılarak kalan 146 satırda skorlandı.

| Model | Raporlanan C-index | Yeniden üretilen 169-satır | Disjoint 146-satır | Δ (disjoint−orijinal) |
|---|---:|---:|---:|---:|
| CoxPH | {scores['CoxPH']['reported']:.6f} | {scores['CoxPH']['reproduced_original_169']:.6f} | {scores['CoxPH']['disjoint_146']:.6f} | {scores['CoxPH']['delta_disjoint_minus_original']:+.6f} |
| RSF | {scores['RSF']['reported']:.6f} | {scores['RSF']['reproduced_original_169']:.6f} | {scores['RSF']['disjoint_146']:.6f} | {scores['RSF']['delta_disjoint_minus_original']:+.6f} |
| DeepSurv | {scores['DeepSurv']['reported']:.6f} | {scores['DeepSurv']['reproduced_original_169']:.6f} | {scores['DeepSurv']['disjoint_146']:.6f} | {scores['DeepSurv']['delta_disjoint_minus_original']:+.6f} |

CoxPH farkı 0.01'in altında; RSF ve DeepSurv mutlak farkları 0.01–0.02 aralığındadır. Disjoint alt-kümede RSF ve DeepSurv skorları düşmek yerine yükselmiştir; dolayısıyla bu tanı, çakışmanın raporlanan C-index'leri yapay olarak yükselttiğine dair ampirik işaret vermez. Bununla birlikte bağımsız test varsayımı 169-satırlık raporlanan metrik için ihlal edilmiştir. Sayısal etki mütevazı olsa da metodolojik sorun makale açısından ihmal edilebilir değildir; tamamen hasta-ayrık performansı belirlemek, bu teşhisin kapsamı dışında kalan hasta-gruplu yeniden split/eğitim gerektirir.

## Bütünlük ve gizlilik kontrolleri

- Dört veri dosyası ve üç model yalnız okundu; model `.fit` çağrısı yapılmadı.
- Orijinal üç C-index 1e-12 toleransında yeniden üretildi.
- Faz 3 ve Faz 6 metadata SHA256 kayıtları mevcut dosyalarla eşleşti.
- `diagnostic.json` ve bu rapor gerçek barcode, patient ID veya sample ID içermez.
"""
    if re.search(r"\bTCGA-[A-Z0-9]{2}-[A-Z0-9]{4}\b", report):
        raise RuntimeError("Patient identifier detected in report")
    REPORT.write_text(report, encoding="utf-8")
    print(json.dumps({"verdict": payload["verdict"], "overlap_patients": 23, "scores": scores}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
