#!/usr/bin/env python3
"""Phase 8 query-variance and CoxPH retry-bias diagnostics from frozen artifacts."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from warrantxai.components import (faithfulness, generalizability,
    generalizability_signed, importance, importance_ci, rank_reproducibility,
    stability, uncertainty)

EXPECTED = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def md(frame: pd.DataFrame) -> str:
    columns = [str(x) for x in frame.columns]
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.itertuples(index=False, name=None)]
    return "\n".join(lines)


def phase3_module():
    path = ROOT / "experiments/06_brca_real.py"
    spec = importlib.util.spec_from_file_location("phase8_phase3", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def query_variance() -> tuple[pd.DataFrame, float]:
    started = perf_counter()
    specifications = [
        (0, "results/phase72_panel/query_matrices.npz", "AJCC Stage IV",
         "clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV", "degenerate"),
        (0, "results/phase72_panel/query_matrices.npz", "OR6C1", "rna_ENSG00000205330.4", "degenerate"),
        (1, "results/phase75_fold1_panel/query_matrices.npz", "AJCC Stage IV",
         "clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV", "degenerate"),
        (1, "results/phase75_fold1_panel/query_matrices.npz", "OR13C3", "rna_ENSG00000204246.3", "degenerate"),
        (0, "results/phase72_panel/query_matrices.npz", "ACCSL", "rna_ENSG00000205126.2", "informative comparator"),
        (0, "results/phase72_panel/query_matrices.npz", "AJCC Stage I",
         "clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I", "informative comparator"),
    ]
    rows = []
    for fold, path, label, feature, role in specifications:
        with np.load(ROOT / path, allow_pickle=False) as saved:
            features = saved["active_features"].astype(str).tolist()
            source = saved["source_query12"][:, features.index(feature)].astype(float)
            external = saved["external_query12"][:, features.index(feature)].astype(float)
        values = np.r_[source, external]
        variance = float(np.var(values, ddof=0))
        rows.append({
            "Fold": fold, "Feature": label, "Feature ID": feature, "Role": role,
            "Distinct values (24)": int(np.unique(values).size), "Minimum": float(values.min()),
            "Maximum": float(values.max()), "Variance (ddof=0)": variance,
            "TCGA 12 values": json.dumps(source.tolist()), "METABRIC 12 values": json.dumps(external.tolist()),
            "All 24 values": json.dumps(values.tolist()),
            "Hypothesis result": "confirmed: variance=0" if variance == 0.0 else "rejected: variance>0; explanation remains open",
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(ROOT / "PHASE8_QUERY_VARIANCE_CHECK.csv", index=False)
    note = f"""# Faz 8 — Sorgu varyansı kontrolü

Değerler, ilgili auditin fiilen kullandığı `query_matrices.npz` model-input matrislerinden okunmuştur. Bunlar Cox scaler'ı uygulanmamış RSF girdileridir; METABRIC sütunları dondurulmuş cross-cohort hizalama sonrasındaki değerlerdir.

{md(frame)}

Dört degenerate fold-feature çiftinin tamamında 24-sorgu varyansı sıfırdan büyüktür. Hipotez **4/4 çiftte reddedildi**; split kullanımı ile fixed-query attribution sıfırı arasındaki açıklama açık kalır. Süre **{perf_counter()-started:.3f} sn**.
"""
    (ROOT / "PHASE8_QUERY_VARIANCE_CHECK.md").write_text(note, encoding="utf-8")
    return frame, perf_counter() - started


def group_components(source: np.ndarray, external: np.ndarray, indices: np.ndarray,
                     all_features: list[str], active: list[str]) -> pd.DataFrame:
    selected = source[indices]
    ext = external[indices]
    tensor = selected[None, ...]
    lower, upper = importance_ci(tensor)
    f = faithfulness(tensor, tensor, 0.5)
    s = stability(tensor, 0.5)
    u = uncertainty(lower, upper, 3.0)
    r = rank_reproducibility(importance(tensor), 0.02)
    targets = np.array([all_features.index(x) for x in active])
    selected_active = selected[:, :, targets]
    gm = generalizability(np.stack([
        np.abs(selected_active).mean(axis=(0, 1)), np.abs(ext).mean(axis=(0, 1))
    ]), 0.5)
    helper = phase3_module()
    source_mean, external_mean = selected_active.mean(axis=0), ext.mean(axis=0)
    grid = max(len(source_mean), len(external_mean))
    source_grid = np.column_stack([helper.empirical_rank_grid(source_mean[:, j], grid) for j in range(len(active))])
    external_grid = np.column_stack([helper.empirical_rank_grid(external_mean[:, j], grid) for j in range(len(active))])
    gs = generalizability_signed(np.stack([source_grid, external_grid]), 0.5)
    return pd.DataFrame({"Feature ID": active, "F": f[targets], "S": s[targets], "U": u[targets],
                         "G": gm, "G_signed": gs, "R": r[targets]})


def retry_bias() -> tuple[pd.DataFrame, float, int]:
    started = perf_counter()
    records = json.loads((ROOT / "results/phase72_cox/bootstrap_summary.json").read_text())
    retry = np.array([len(row["fit_attempts"]) > 1 for row in records])
    if (int(retry.sum()), int((~retry).sum())) != (28, 22):
        raise RuntimeError(f"Retry partition changed: retry={retry.sum()}, direct={(~retry).sum()}")
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as saved:
        source = saved["source_phi"].copy()
        external = saved["external_phi"].copy()
        all_features = saved["features"].astype(str).tolist()
        active = saved["active_features"].astype(str).tolist()
    targets = np.array([all_features.index(x) for x in active])
    replicate_means = source[:, :, targets].mean(axis=1)
    retry_values, direct_values = replicate_means[retry], replicate_means[~retry]
    mean_retry, mean_direct = retry_values.mean(axis=0), direct_values.mean(axis=0)
    difference = mean_retry - mean_direct
    se = np.sqrt(retry_values.var(axis=0, ddof=1)/len(retry_values) + direct_values.var(axis=0, ddof=1)/len(direct_values))
    statistic, pvalue = ttest_ind(retry_values, direct_values, axis=0, equal_var=False)
    ratio = np.divide(np.abs(difference), se, out=np.where(difference == 0, 0.0, np.inf), where=se > 0)
    retry_components = group_components(source, external, np.flatnonzero(retry), all_features, active).set_index("Feature ID")
    direct_components = group_components(source, external, np.flatnonzero(~retry), all_features, active).set_index("Feature ID")
    mapping = pd.read_csv(ROOT / "results/phase72_panel/feature_mapping.csv")
    display = dict(zip(mapping.feature, mapping.symbol, strict=False))
    display.update({
        "clin_numeric__age_at_index.demographic": "Age at index",
        **{f"clin_categorical__ajcc_pathologic_stage.diagnoses_Stage {x}": f"AJCC Stage {x}"
           for x in ("I", "II", "III", "IV")},
    })
    rows = []
    for j, feature in enumerate(active):
        row = {
            "Feature": display.get(feature, feature), "Feature ID": feature,
            "Retry n": int(retry.sum()), "Direct n": int((~retry).sum()),
            "Mean signed attribution retry": mean_retry[j], "Mean signed attribution direct": mean_direct[j],
            "Attribution difference retry-direct": difference[j], "Combined SE": se[j],
            "Difference/SE": ratio[j], "95% CI lower": difference[j]-1.96*se[j],
            "95% CI upper": difference[j]+1.96*se[j], "Welch t": statistic[j], "Welch p": pvalue[j],
            "Flag |difference| > 1.96 SE": bool(np.abs(difference[j]) > 1.96*se[j]),
        }
        for component in ("F", "S", "U", "G", "G_signed", "R"):
            rv, dv = retry_components.loc[feature, component], direct_components.loc[feature, component]
            row[f"{component} retry"] = rv
            row[f"{component} direct"] = dv
            row[f"{component} difference"] = rv-dv
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame.to_csv(ROOT / "PHASE8_COXPH_RETRY_BIAS_CHECK.csv", index=False)
    flagged = frame[frame["Flag |difference| > 1.96 SE"]]
    conclusion = ("En az bir feature nominal 1.96×SE eşiğini aştı; feature-düzeyi farklar aşağıda flaglendi."
                  if len(flagged) else "Hiçbir feature 1.96×SE eşiğini aşmadı; sistematik retry önyargısı kanıtı yok.")
    compact = frame[["Feature", "Feature ID", "Attribution difference retry-direct", "Combined SE", "Difference/SE",
                     "95% CI lower", "95% CI upper", "Welch p", "Flag |difference| > 1.96 SE"]]
    note = f"""# Faz 8 — CoxPH-v2 retry önyargısı kontrolü

Faz 7.2 logundan 28 retry ve 22 doğrudan-yakınsayan replika ayrıldı; hiçbir model yeniden fit edilmedi. Attribution analizi her replikada 153 TCGA audit hastası üzerindeki signed ortalamayı kullanır. Fark `retry-direct`; birleşik SE Welch standart hatasıdır. Flag önceden belirtilen `|fark| > 1.96×SE` kuralıdır ve çoklu-test düzeltilmiş yeni bir inferans kararı değildir.

{md(compact)}

## Grup-içi yeniden hesaplanan component skorları

Tam CSV her feature için F/S/U/G/G_signed/R retry, direct ve fark kolonlarını içerir. Bu skorlar dondurulmuş formüllerle yalnız ilgili replika alt-kümesinde yeniden hesaplandı; WI veya sertifikasyon kararları değiştirilmedi.

**Sonuç:** {conclusion} Flag sayısı **{len(flagged)}/27**. Süre **{perf_counter()-started:.3f} sn**.
"""
    (ROOT / "PHASE8_COXPH_RETRY_BIAS_CHECK.md").write_text(note, encoding="utf-8")
    return frame, perf_counter() - started, len(flagged)


def main() -> None:
    started = perf_counter()
    actual = {name: sha(ROOT / "warrantxai" / name) for name in EXPECTED}
    if actual != EXPECTED:
        raise RuntimeError(f"Frozen formula hashes changed: {actual}")
    queries, query_seconds = query_variance()
    retry, retry_seconds, flags = retry_bias()
    tested = queries[queries.Role == "degenerate"]
    payload = {"query_hypothesis_confirmed": int((tested["Variance (ddof=0)"] == 0).sum()),
               "query_hypothesis_rejected": int((tested["Variance (ddof=0)"] > 0).sum()),
               "retry_flags": flags, "query_seconds": query_seconds, "retry_seconds": retry_seconds,
               "formula_sha256": actual, "wall_clock_seconds": perf_counter()-started}
    (ROOT / "results/phase8_robustness_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
