#!/usr/bin/env python3
"""Classify frozen certified warrants without refitting or resampling."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from time import perf_counter

import joblib
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_HASHES = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def md(frame: pd.DataFrame) -> str:
    labels = [str(x) for x in frame.columns]
    rows = [[str(x) for x in row] for row in frame.itertuples(index=False, name=None)]
    widths = [max(len(labels[i]), *(len(row[i]) for row in rows)) for i in range(len(labels))]
    line = lambda row: "| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(row))) + " |"
    return "\n".join([line(labels), "| " + " | ".join("---" for _ in labels) + " |", *(line(r) for r in rows)])


def display_maps() -> tuple[dict[str, str], dict[str, str]]:
    names: dict[str, str] = {}
    for filename in ("PHASE73_COXPH_FULLTABLE.csv", "PHASE73_RSF_FULLTABLE.csv"):
        table = pd.read_csv(ROOT / filename)
        names.update(dict(zip(table["Feature ID"], table["Feature"], strict=False)))
    mapping = pd.read_csv(ROOT / "results/phase75_fold1_panel/feature_mapping.csv")
    names.update(dict(zip(mapping.feature, mapping.symbol, strict=False)))
    clinical = {
        "clin_numeric__age_at_index.demographic": "Age at index",
        **{
            f"clin_categorical__ajcc_pathologic_stage.diagnoses_Stage {stage}": f"AJCC Stage {stage}"
            for stage in ("I", "II", "III", "IV")
        },
    }
    names.update(clinical)
    return names, clinical


def certified(path: str) -> list[str]:
    table = pd.read_csv(ROOT / path)
    return table.loc[table.certified_BH_alpha_0_05.astype(bool), "feature"].tolist()


def classify(model: str, features: list[str], tensor: np.ndarray, tensor_features: list[str],
             scope: str, names: dict[str, str]) -> pd.DataFrame:
    rows = []
    for feature in features:
        values = tensor[..., tensor_features.index(feature)]
        mean = float(values.mean())
        rows.append({
            "Model": model,
            "Feature": names.get(feature, feature),
            "Feature ID": feature,
            "Mean signed attribution": mean,
            "Non-zero replicate frequency": int(np.count_nonzero(np.any(values != 0.0, axis=1))),
            "Non-zero query frequency": int(np.count_nonzero(np.any(values != 0.0, axis=0))),
            "Query denominator": int(values.shape[1]),
            "Averaging scope": scope,
            "Label": "informative" if mean != 0.0 else "degenerate",
        })
    return pd.DataFrame(rows)


def main() -> None:
    all_started = perf_counter()
    actual = {name: sha(ROOT / "warrantxai" / name) for name in EXPECTED_HASHES}
    if actual != EXPECTED_HASHES:
        raise RuntimeError(f"Frozen formula hashes changed: {actual}")
    names, _ = display_maps()

    step1_started = perf_counter()
    cox_cert = certified("results/phase72_coxph_certification/feature_certification.csv")
    rsf_cert = certified("results/phase72_rsf_certification/feature_certification.csv")
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as saved:
        cox = classify("CoxPH-v2", cox_cert, saved["source_phi"], saved["features"].astype(str).tolist(),
                       "153 TCGA audit patients × 50 replicas", names)
    with np.load(ROOT / "results/phase72_rsf/attributions_n50.npz", allow_pickle=False) as saved:
        rsf = classify("RSF", rsf_cert, saved["candidate"][0], saved["features"].astype(str).tolist(),
                       "fixed 12 TCGA queries × 50 replicas", names)
    counts = pd.concat([cox, rsf]).groupby(["Model", "Label"]).size().to_dict()
    expected_counts = {
        ("CoxPH-v2", "informative"): 8,
        ("RSF", "informative"): 4,
        ("RSF", "degenerate"): 8,
    }
    if counts != expected_counts:
        raise RuntimeError(f"Fold-0 classification assertion failed: observed={counts}, expected={expected_counts}")
    fold0 = pd.concat([cox, rsf], ignore_index=True)
    fold0.to_csv(ROOT / "PHASE76_WARRANT_CLASSIFICATION_FOLD0.csv", index=False)
    step1 = perf_counter() - step1_started
    report0 = f"""# Faz 7.6 — Fold-0 warrant sınıflandırması

Etiket kuralı önceden tanımlıdır: saklı signed attribution tensörünün aritmetik ortalaması exact `0.0` ise `degenerate`, değilse `informative`. Ek test veya tolerans kullanılmadı.

{md(fold0)}

Assertion geçti: CoxPH-v2 **8 informative / 0 degenerate**; RSF **4 informative / 8 degenerate**. Birleşim 14 benzersiz feature, model-feature düzeyinde 20 satırdır. Süre **{step1:.3f} sn**.
"""
    (ROOT / "PHASE76_WARRANT_CLASSIFICATION_FOLD0.md").write_text(report0, encoding="utf-8")

    step2_started = perf_counter()
    bundle = joblib.load(ROOT / "results/phase7_retrain/rsf_baseline.joblib")
    model, model_features = bundle["model"], list(bundle["features"])
    try:
        model.feature_importances_
        importance_status = "available"
    except NotImplementedError:
        importance_status = "unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError"
    split_rows = []
    label_by_feature = rsf.set_index("Feature ID")["Label"].to_dict()
    for feature in rsf_cert:
        index = model_features.index(feature)
        tree_count = int(sum(np.any(tree.tree_.feature == index) for tree in model.estimators_))
        split_count = int(sum(np.count_nonzero(tree.tree_.feature == index) for tree in model.estimators_))
        split_rows.append({
            "Feature": names.get(feature, feature), "Feature ID": feature,
            "Label": label_by_feature[feature], "feature_importances_": importance_status,
            "Trees with ≥1 split": tree_count, "Total split nodes": split_count,
        })
    split = pd.DataFrame(split_rows)
    step2 = perf_counter() - step2_started
    deg_used = split[(split.Label == "degenerate") & (split["Trees with ≥1 split"] > 0)]
    split_note = f"""# Faz 7.6 — Fold-0 RSF split-kullanım diagnostik

Ana fold-0 RSF artefaktındaki 500 fitted ağacın `tree_.feature` dizileri doğrudan okundu; model yeniden fit edilmedi. `feature_importances_` bu scikit-survival sınıfında mevcut görünse de çağrı doğrudan `NotImplementedError` üretir; bu yüzden değer uydurulmadı.

{md(split)}

Beklenen tam uyum görülmedi. Sekiz degenerate feature'ın altısı hiçbir ağaçta kullanılmadı; **OR6C1 10/500 ağaçta 10 split**, **AJCC Stage IV 25/500 ağaçta 25 split** içinde kullanıldığı halde saklı fixed-query candidate attribution'ları exact sıfırdır. Bu tablo yalnız gözlenen tutarsızlığı kaydeder. Süre **{step2:.3f} sn**.
"""
    (ROOT / "PHASE76_RSF_SPLIT_USAGE.md").write_text(split_note, encoding="utf-8")
    split.to_csv(ROOT / "results/phase76_rsf_split_usage.csv", index=False)

    step3_started = perf_counter()
    fold1_cert = certified("results/phase75_fold1_rsf_certification/feature_certification.csv")
    with np.load(ROOT / "results/phase75_fold1_rsf/attributions_n50.npz", allow_pickle=False) as saved:
        fold1 = classify("RSF", fold1_cert, saved["candidate"][0], saved["features"].astype(str).tolist(),
                         "fixed 12 TCGA queries × 50 replicas", names)
    fold1.to_csv(ROOT / "PHASE76_WARRANT_CLASSIFICATION_FOLD1_RSF.csv", index=False)
    fold1_counts = fold1.Label.value_counts().to_dict()
    step3 = perf_counter() - step3_started
    report1 = f"""# Faz 7.6 — Fold-1 RSF warrant sınıflandırması

Fold-1'in yedi BH-sertifikalı RSF feature'ı aynı exact-zero kuralıyla sınıflandırıldı.

{md(fold1)}

Sonuç: **{fold1_counts.get('informative', 0)} informative / {fold1_counts.get('degenerate', 0)} degenerate**. Süre **{step3:.3f} sn**.
"""
    (ROOT / "PHASE76_WARRANT_CLASSIFICATION_FOLD1_RSF.md").write_text(report1, encoding="utf-8")

    total = perf_counter() - all_started
    summary = f"""# Faz 7.6 — Özet

Bu faz yalnız saklı attribution tensörlerini ve fitted fold-0 RSF ağaçlarını okudu. Yeni model fit, explainer, bootstrap, kalibrasyon veya sertifikasyon çalıştırılmadı; mevcut BH/BY kararları değiştirilmedi.

| Adım | Sonuç | Süre |
|---|---|---:|
| Fold-0 sınıflandırma | CoxPH-v2 8/0; RSF 4/8 informative/degenerate; assertion geçti | {step1:.3f} sn |
| Fold-0 RSF split kullanımı | 8 degenerate feature'ın 6'sı kullanılmamış; OR6C1 10, AJCC Stage IV 25 ağaçta kullanılmış | {step2:.3f} sn |
| Fold-1 RSF sınıflandırma | {fold1_counts.get('informative', 0)}/{fold1_counts.get('degenerate', 0)} informative/degenerate | {step3:.3f} sn |

Toplam wall-clock **{total:.3f} sn**. Formül SHA256 değerleri değişmedi: `components.py` `{actual['components.py']}`, `aggregation.py` `{actual['aggregation.py']}`.
"""
    (ROOT / "PHASE76_SUMMARY.md").write_text(summary, encoding="utf-8")
    print(json.dumps({"fold0": {str(k): v for k, v in counts.items()}, "fold1": fold1_counts,
                      "degenerate_split_used": deg_used[["Feature", "Trees with ≥1 split"]].to_dict("records"),
                      "seconds": total}, indent=2))


if __name__ == "__main__":
    main()
