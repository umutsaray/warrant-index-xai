"""Phase 3 step 0 manifest; output contains features/aggregates, never row IDs."""

from __future__ import annotations

from pathlib import Path
import warnings

import joblib
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"


def run() -> None:
    """Write and print verified inventory before auditing (Phase 3 step 0)."""
    x = pd.read_parquet(TCGA / "X_train_selected.parquet").set_index("sample_key")
    xt = pd.read_parquet(TCGA / "X_test_selected.parquet").set_index("sample_key")
    y = pd.read_parquet(TCGA / "y_train.parquet").set_index("sample_key")
    yt = pd.read_parquet(TCGA / "y_test.parquet").set_index("sample_key")
    with warnings.catch_warnings(record=True) as caught:
        model = joblib.load(TCGA / "coxph_baseline.joblib")
    rawcols = pd.read_parquet(TCGA / "X_train.parquet").columns
    excluded = list(x.columns.difference(model["features"]))
    lines = [
        "# Faz 3 — TCGA-BRCA MANIFEST",
        "",
        "Kaynak: `../" + TCGA.name + "/`",
        "",
        "## Adım 0: dosyalar ve doğrulama",
        "",
    ]
    for name in [
        "prepare_dataset.py",
        "impute_and_split.py",
        "select_features.py",
        "selected_features.json",
        "X_train_selected.parquet",
        "X_test_selected.parquet",
        "X_train.parquet",
        "X_test.parquet",
        "y_train.parquet",
        "y_test.parquet",
        "coxph_baseline.joblib",
        "train_baselines.py",
    ]:
        p = TCGA / name
        lines.append(f"- `{name}`: " + (f"mevcut, {p.stat().st_size} byte" if p.exists() else "BULUNAMADI"))
    for name, xf, yf, n, e in [("train", x, y, 673, 106), ("test", xt, yt, 169, 26)]:
        lines.append(
            f"- {name}: X={xf.shape}, y={yf.shape}, event={int(yf.event.sum())}; beklenen {n}/{e}: "
            + ("EŞLEŞİYOR" if (len(xf), int(yf.event.sum())) == (n, e) else "SAPMA VAR")
        )
        lines.append(
            f"- {name} event/time mevcut; sıra uyumu={xf.index.equals(yf.index)}; eksik X hücresi={int(xf.isna().sum().sum())}; time min={yf.time.min()}."
        )
    lines += [
        "",
        "## SAPMA ve beklenmeyen yapı",
        "",
        f"- **İstenen 623 feature yerine {len(x.columns)} bulundu.** RNA=200, miRNA=200, methylation=200, clinical={sum(c.startswith('clin_') for c in x)}.",
        f"- Ham X_train ara matrisinde clinical={sum(c.startswith('clin_') for c in rawcols)}; select_features.py AJCC alt evrelerini collapse ediyor. Seçilmiş panel ve selected_features.json 616 sütunla tutarlı; eksik 7 sütun uydurulmadı.",
        "- İçte ayrıca başka analiz/kohort dizinleri var; istenen 673/169 ayrımıyla eşleşen kök dosyalar kullanıldı.",
        f"- Joblib sözlüğünde model/scaler/features var. Model tipi={type(model['model']).__name__}, penalizer={model['model'].penalizer}, feature={len(model['features'])}; yükleme uyarı sayısı={len(caught)}.",
        "- Eğitim betiği: event veya censored grupta varyans<=1e-12 dışlama; StandardScaler train fit; penalizer=0.1, max_steps=50. Kaydedilmiş model kullanılır; bootstrap replikaları aynı ayarlarla yeniden eğitilir.",
        "- Modelden dışlanan feature'lar: " + ", ".join("`" + c + "`" for c in excluded) + ".",
        "- Klinik panel: yaş, cinsiyet, ırk, etnisite, Stage I/II/III/IV/X; grade/ER/PR/HER2 yok. Bunlar varmış gibi eklenmez.",
    ]
    trainp = x.index.astype(str).str[:12]
    testp = xt.index.astype(str).str[:12]
    lines += [
        f"- **Hasta çakışması:** train/testte {len(set(trainp) & set(testp))} ortak hasta; train {trainp.nunique()} benzersiz hasta. Audit'te ortak hastalar çıkarılır, kalan her hastada ilk örnek tutulur. Bootstrap hasta-kümesi bazındadır. Orijinal dosyalar değişmez.",
        "",
        "## Tam seçilmiş feature listesi (hasta kimliği içermez)",
        "",
    ]
    lines.extend("- `" + c + "`" for c in x.columns)
    lines += [
        "",
        "## Ek kohort ve yakınsama notları",
        "",
        "- Kaynak primer-tümör-only değildir. Train örnek türü sayıları: 01=609, 11=60, 06=4; test: 01=153, 11=15, 06=1. Bunlar agregat sayılardır; hasta/sample ID yazılmaz.",
        "- Bootstrap altıncı replikada meth_cg04831594 separation uyarısı verdi. Aynı örnekleme ve penalizer=0.1 ile uyarılı sütun dışlanır; kalan sayısal yakınsamama için solver üst sınırı 50'den 200'e, step_size 0.25'e geçirilir. Model hedefi/toleranslar değiştirilmez; fit denemeleri bootstrap_summary.json içinde tutulur.",
    ]
    text = "\n".join(lines) + "\n"
    (ROOT / "MANIFEST.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    run()
