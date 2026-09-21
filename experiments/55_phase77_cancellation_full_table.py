#!/usr/bin/env python3
"""Verify cell-level zero attributions and build the frozen fold-0 supplementary table."""

from __future__ import annotations

import hashlib
from pathlib import Path
from time import perf_counter

import joblib
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def md(frame: pd.DataFrame) -> str:
    # Compact GitHub table; CSV remains the machine-readable primary artifact.
    cols = [str(x) for x in frame.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    return "\n".join(lines)


def names() -> dict[str, str]:
    result: dict[str, str] = {}
    for path in (ROOT / "PHASE73_COXPH_FULLTABLE.csv", ROOT / "PHASE73_RSF_FULLTABLE.csv"):
        table = pd.read_csv(path)
        result.update(dict(zip(table["Feature ID"], table["Feature"], strict=False)))
    result.update({
        "clin_numeric__age_at_index.demographic": "Age at index",
        **{f"clin_categorical__ajcc_pathologic_stage.diagnoses_Stage {x}": f"AJCC Stage {x}"
           for x in ("I", "II", "III", "IV")},
    })
    return result


def load_tensor(path: str, key: str, site: int | None = None) -> tuple[np.ndarray, list[str]]:
    with np.load(ROOT / path, allow_pickle=False) as saved:
        tensor = saved[key]
        if site is not None:
            tensor = tensor[site]
        return tensor.copy(), saved["features"].astype(str).tolist()


def cell_metrics(values: np.ndarray) -> dict:
    nonzero = values != 0.0
    count = int(np.count_nonzero(nonzero))
    return {
        "All cells exactly zero": bool(np.all(values == 0.0)),
        "M_j mean absolute attribution": float(np.abs(values).mean()),
        "Non-zero cell count": count,
        "Cell denominator": int(values.size),
        "P_j non-zero frequency": count / values.size,
        "Q_j max absolute attribution": float(np.abs(values).max()),
    }


def cancellation_rows(fold: int, model: str, selected: pd.DataFrame, tensor: np.ndarray,
                      features: list[str], display: dict[str, str]) -> list[dict]:
    rows = []
    for item in selected.itertuples(index=False):
        feature = getattr(item, "Feature_ID")
        values = tensor[:, :, features.index(feature)]
        rows.append({"Fold": fold, "Model": model, "Feature": display.get(feature, feature),
                     "Feature ID": feature, "Phase 7.6 label": item.Label, **cell_metrics(values)})
    return rows


def feature_type(feature: str) -> str:
    if feature.startswith("rna_"):
        return "RNA"
    if feature.startswith("clin_numeric__"):
        return "clinical-numeric"
    return "clinical-categorical"


def main() -> None:
    total_started = perf_counter()
    actual = {name: sha(ROOT / "warrantxai" / name) for name in EXPECTED}
    if actual != EXPECTED:
        raise RuntimeError(f"Frozen formula hashes changed: {actual}")
    display = names()

    step1_started = perf_counter()
    phase76_f0 = pd.read_csv(ROOT / "PHASE76_WARRANT_CLASSIFICATION_FOLD0.csv")
    phase76_f1 = pd.read_csv(ROOT / "PHASE76_WARRANT_CLASSIFICATION_FOLD1_RSF.csv")
    cox_tensor, cox_features = load_tensor("results/phase72_cox/attributions_n50.npz", "source_phi")
    rsf0_tensor, rsf0_features = load_tensor("results/phase72_rsf/attributions_n50.npz", "candidate", 0)
    rsf1_tensor, rsf1_features = load_tensor("results/phase75_fold1_rsf/attributions_n50.npz", "candidate", 0)
    selected_cox = phase76_f0[phase76_f0.Model == "CoxPH-v2"].rename(columns={"Feature ID": "Feature_ID"})
    selected_rsf0 = phase76_f0[phase76_f0.Model == "RSF"].rename(columns={"Feature ID": "Feature_ID"})
    selected_rsf1 = phase76_f1.rename(columns={"Feature ID": "Feature_ID"})
    rows = []
    rows += cancellation_rows(0, "CoxPH-v2", selected_cox, cox_tensor, cox_features, display)
    rows += cancellation_rows(0, "RSF", selected_rsf0, rsf0_tensor, rsf0_features, display)
    rows += cancellation_rows(1, "RSF", selected_rsf1[selected_rsf1.Label == "degenerate"],
                              rsf1_tensor, rsf1_features, display)
    check = pd.DataFrame(rows)
    degenerate = check[check["Phase 7.6 label"] == "degenerate"]
    failures = degenerate[~degenerate["All cells exactly zero"]]
    if len(failures):
        detail = []
        for row in failures.itertuples(index=False):
            tensor, feats = (rsf0_tensor, rsf0_features) if row.Fold == 0 else (rsf1_tensor, rsf1_features)
            values = tensor[:, :, feats.index(getattr(row, "Feature_ID"))]
            indices = np.argwhere(values != 0.0)
            detail.append({"fold": row.Fold, "feature": getattr(row, "Feature_ID"),
                           "nonzero": [(int(t), int(q), float(values[t, q]),
                                        "positive" if values[t, q] > 0 else "negative") for t, q in indices]})
        (ROOT / "PHASE77_CANCELLATION_FAILURE.txt").write_text(repr(detail), encoding="utf-8")
        raise RuntimeError(f"Degenerate cell-level zero assertion failed: {detail}")
    assert len(degenerate) == 10
    assert (degenerate[["M_j mean absolute attribution", "P_j non-zero frequency",
                        "Q_j max absolute attribution"]] == 0.0).all().all()
    informative = check[check["Phase 7.6 label"] == "informative"]
    assert (~informative["All cells exactly zero"]).all()
    check.to_csv(ROOT / "PHASE77_CANCELLATION_CHECK.csv", index=False)
    step1 = perf_counter() - step1_started
    report1 = f"""# Faz 7.7 — Tekil-değer cancellation kontrolü

Metrikler, hiçbir agregasyon öncesinde saklı `(replica, query)` hücrelerinden hesaplandı. RSF için her feature 50×12=600 hücre; CoxPH-v2 için 50×153=7650 hücredir. Tarihsel kolon adı “All 600 cells” olmakla birlikte CoxPH satırlarının gerçek paydası `Cell denominator=7650` olarak açıkça saklanmıştır.

{md(check)}

Ana assertion geçti: fold-0'daki sekiz ve fold-1'deki iki ayrı degenerate RSF sonucu dahil **10/10** satırda bütün 600 tekil değer exact `0.0`; dolayısıyla M_j=P_j=Q_j=0. Karşılaştırma satırlarının tamamında informative özellikler için exact-zero assertion `False`, M_j>0 ve P_j>0'dır. Süre **{step1:.3f} sn**.
"""
    (ROOT / "PHASE77_CANCELLATION_CHECK.md").write_text(report1, encoding="utf-8")

    step2_started = perf_counter()
    sources = {
        "CoxPH-v2": ("results/phase72_cox/feature_report.csv",
                     "results/phase72_coxph_certification/feature_certification.csv", cox_tensor, cox_features),
        "RSF": ("results/phase72_rsf/feature_report.csv",
                "results/phase72_rsf_certification/feature_certification.csv", rsf0_tensor, rsf0_features),
    }
    deep_tensor, deep_features = load_tensor("results/phase72_deepsurv/attributions_n50.npz", "candidate", 0)
    sources["DeepSurv"] = ("results/phase72_deepsurv/feature_report.csv",
                           "results/phase72_deepsurv_certification/feature_certification.csv",
                           deep_tensor, deep_features)
    rsf_bundle = joblib.load(ROOT / "results/phase7_retrain/rsf_baseline.joblib")
    rsf_model, rsf_model_features = rsf_bundle["model"], list(rsf_bundle["features"])
    full_rows = []
    for model, (report_path, cert_path, tensor, tensor_features) in sources.items():
        report = pd.read_csv(ROOT / report_path).set_index("feature")
        cert = pd.read_csv(ROOT / cert_path).set_index("feature")
        for feature in tensor_features if model != "CoxPH-v2" else report.index.tolist():
            # Cox tensor has all 610 fitted features; the report index is the common 27-feature family.
            if feature not in report.index:
                continue
            values = tensor[:, :, tensor_features.index(feature)]
            metrics = cell_metrics(values)
            bh = bool(cert.loc[feature, "certified_BH_alpha_0_05"])
            split_trees = "N/A"
            if model == "RSF":
                feature_index = rsf_model_features.index(feature)
                split_trees = int(sum(np.any(tree.tree_.feature == feature_index) for tree in rsf_model.estimators_))
            full_rows.append({
                "Feature": display.get(feature, feature), "Feature ID": feature,
                "Feature type": feature_type(feature), "Model": model,
                "Mean signed attribution": float(values.mean()),
                "M_j mean absolute attribution": metrics["M_j mean absolute attribution"],
                "P_j non-zero frequency": metrics["P_j non-zero frequency"],
                "Q_j max absolute attribution": metrics["Q_j max absolute attribution"],
                "F": report.loc[feature, "F"], "S": report.loc[feature, "S"],
                "U": report.loc[feature, "U"], "G": report.loc[feature, "G_magnitude"],
                "G_signed": report.loc[feature, "G_signed"], "R": report.loc[feature, "R"],
                "WI": report.loc[feature, "WI"],
                "Point-adequate": "Evet" if bool(report.loc[feature, "point_adequate"]) else "Hayır",
                "p_joint": cert.loc[feature, "p_joint"], "BH-adjusted p": cert.loc[feature, "q_BH"],
                "BY-adjusted p": cert.loc[feature, "q_BY"],
                "BH certified": "Evet" if bh else "Hayır",
                "BY certified": "Evet" if bool(cert.loc[feature, "certified_BY_alpha_0_05"]) else "Hayır",
                "Informative/Degenerate label": ("informative" if values.mean() != 0.0 else "degenerate")
                                                  if bh else "N/A (not certified)",
                "RSF trees with split": split_trees,
            })
    full = pd.DataFrame(full_rows)
    model_order = {"CoxPH-v2": 0, "RSF": 1, "DeepSurv": 2}
    feature_order = {f: i for i, f in enumerate(rsf0_features)}
    full["_m"] = full.Model.map(model_order)
    full["_f"] = full["Feature ID"].map(feature_order)
    full = full.sort_values(["_f", "_m"]).drop(columns=["_m", "_f"]).reset_index(drop=True)
    assert len(full) == 81 and full.groupby("Model").size().eq(27).all()
    full.to_csv(ROOT / "PHASE77_FULL_FEATURE_TABLE_FOLD0.csv", index=False)
    step2 = perf_counter() - step2_started
    report2 = f"""# Faz 7.7 — Fold-0 tam feature supplementary tablosu

Bu tablo 27 ortak feature × 3 model = **81 satır** içerir. Attribution kapsamı CoxPH-v2 için 153 TCGA audit hastası; RSF ve DeepSurv için saklı deterministic 12 TCGA query'dir. RSF split sütunu ana 500-ağaç fold-0 modelinden doğrudan okunmuştur. Sertifikasyon kararları ve component skorları değiştirilmeden Faz 7.2 artefaktlarından alınmıştır.

{md(full)}

Üretim süresi **{step2:.3f} sn**.
"""
    (ROOT / "PHASE77_FULL_FEATURE_TABLE_FOLD0.md").write_text(report2, encoding="utf-8")

    total = perf_counter() - total_started
    summary = f"""# Faz 7.7 — Özet

| Adım | Sonuç | Süre |
|---|---|---:|
| Cancellation kontrolü | 10/10 degenerate RSF fold-feature satırında 600/600 exact sıfır; M_j=P_j=Q_j=0 | {step1:.3f} sn |
| Tam supplementary tablo | 27 feature × 3 model = 81 satır; bütün component, inferans, attribution ve RSF split kolonları mevcut | {step2:.3f} sn |

Toplam wall-clock **{total:.3f} sn**. Yeni bootstrap, model fit, explainer, kalibrasyon veya sertifikasyon yapılmadı. Formül SHA256 değişmedi: `components.py` `{actual['components.py']}`, `aggregation.py` `{actual['aggregation.py']}`.
"""
    (ROOT / "PHASE77_SUMMARY.md").write_text(summary, encoding="utf-8")
    print({"degenerate_exact_zero": int(degenerate["All cells exactly zero"].sum()),
           "cancellation_rows": len(check), "full_rows": len(full), "seconds": total})


if __name__ == "__main__":
    main()
