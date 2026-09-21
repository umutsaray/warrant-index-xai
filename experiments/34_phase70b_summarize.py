#!/usr/bin/env python3
"""Summarize Phase 7.0b fold, feature-selection, CoxPH, and DeepSurv diagnostics."""

from __future__ import annotations

import hashlib
import itertools
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
SIBLING = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results" / "phase70b_collapse_diagnostic"
PHASE7 = ROOT / "results" / "phase7_retrain"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    started = time.perf_counter()
    folds = json.loads((OUT / "fold_scores.json").read_text(encoding="utf-8"))["folds"]
    retries = json.loads((OUT / "cox_fold_retries.json").read_text(encoding="utf-8"))
    fold_table = pd.DataFrame(folds)
    fold_table["coxph_converged_c_index"] = [
        retries.get(str(fold), {}).get("retry_c_index", score)
        for fold, score in zip(fold_table.fold, fold_table.coxph_test_c_index)
    ]
    fold_table.to_csv(OUT / "fold_scores.csv", index=False)

    old = json.loads((SIBLING / "selected_features.json").read_text(encoding="utf-8"))
    new = json.loads((PHASE7 / "selected_features.json").read_text(encoding="utf-8"))
    rows = []
    for block in ("rna", "mirna", "methylation"):
        old_set, new_set = set(old["selected_by_block"][block]), set(new["selected_by_block"][block])
        rows.append(
            {
                "block": block,
                "old_count": len(old_set),
                "new_count": len(new_set),
                "common_count": len(old_set & new_set),
                "common_fraction_of_200": len(old_set & new_set) / 200,
                "jaccard": len(old_set & new_set) / len(old_set | new_set),
            }
        )
    old_clinical = set(old["clinical_features_retained_without_univariate_filter"])
    new_clinical = set(new["clinical_features_retained_without_univariate_filter"])
    rows.append(
        {
            "block": "clinical",
            "old_count": len(old_clinical),
            "new_count": len(new_clinical),
            "common_count": len(old_clinical & new_clinical),
            "common_fraction_of_200": len(old_clinical & new_clinical) / len(old_clinical),
            "jaccard": len(old_clinical & new_clinical) / len(old_clinical | new_clinical),
        }
    )
    old_all, new_all = set(old["final_features"]), set(new["final_features"])
    rows.append(
        {
            "block": "total",
            "old_count": len(old_all),
            "new_count": len(new_all),
            "common_count": len(old_all & new_all),
            "common_fraction_of_200": len(old_all & new_all) / len(old_all),
            "jaccard": len(old_all & new_all) / len(old_all | new_all),
        }
    )
    stability = pd.DataFrame(rows).rename(columns={"common_fraction_of_200": "common_fraction"})
    stability.to_csv(OUT / "old_new_feature_stability.csv", index=False)

    fold_selections = [new] + [
        json.loads((OUT / f"fold_{fold}" / "selected_features.json").read_text(encoding="utf-8"))
        for fold in range(1, 5)
    ]
    grouped_pairwise = {}
    for block in ("rna", "mirna", "methylation"):
        counts = []
        for left, right in itertools.combinations(fold_selections, 2):
            counts.append(len(set(left["selected_by_block"][block]) & set(right["selected_by_block"][block])))
        grouped_pairwise[block] = {"mean_common": float(np.mean(counts)), "minimum": min(counts), "maximum": max(counts)}
    counts = [
        len(set(left["final_features"]) & set(right["final_features"]))
        for left, right in itertools.combinations(fold_selections, 2)
    ]
    grouped_pairwise["total"] = {"mean_common": float(np.mean(counts)), "minimum": min(counts), "maximum": max(counts)}

    deep = json.loads((OUT / "deepsurv_diagnostic.json").read_text(encoding="utf-8"))
    deep_summary = {}
    selected_rows = []
    for key in ("row_level_validation", "patient_grouped_validation"):
        regime = deep[key]
        history = regime["history"]
        best_test = max(history, key=lambda row: row["test_c_index"])
        correlation = spearmanr(
            [row["validation_c_index"] for row in history],
            [row["test_c_index"] for row in history],
        )
        deep_summary[key] = {
            "fit_validation_patient_overlap": regime["fit_validation_patient_overlap"],
            "best_validation_epoch": regime["best_epoch"],
            "best_validation_c_index": regime["best_validation_c_index"],
            "best_checkpoint_test_c_index": regime["best_checkpoint_test_c_index"],
            "maximum_observed_test_epoch": best_test["epoch"],
            "maximum_observed_test_c_index": best_test["test_c_index"],
            "validation_test_spearman_rho": float(correlation.statistic),
            "wall_clock_seconds": regime["wall_clock_seconds"],
        }
        selected_rows.extend({"regime": key, **row} for row in regime["selected_epoch_metrics"])
    pd.DataFrame(selected_rows).to_csv(OUT / "deepsurv_selected_epochs.csv", index=False)

    formula_paths = [ROOT / "warrantxai/components.py", ROOT / "warrantxai/aggregation.py"]
    payload = {
        "scope": "diagnostic only; no audit pipeline or formula change",
        "fold_summary": {
            "rsf_mean": float(fold_table.rsf_test_c_index.mean()),
            "rsf_minimum": float(fold_table.rsf_test_c_index.min()),
            "rsf_maximum": float(fold_table.rsf_test_c_index.max()),
            "rsf_fold0_rank_of_5": int(fold_table.rsf_test_c_index.rank(method="min").iloc[0]),
            "converged_cox_mean": float(fold_table.coxph_converged_c_index.mean()),
            "converged_cox_minimum": float(fold_table.coxph_converged_c_index.min()),
            "converged_cox_maximum": float(fold_table.coxph_converged_c_index.max()),
        },
        "old_new_feature_stability": rows,
        "grouped_fold_pairwise_feature_stability": grouped_pairwise,
        "cox_retries": retries,
        "deepsurv": deep_summary,
        "formula_sha256": {path.name: sha256(path) for path in formula_paths},
        "audit_pipeline_run": False,
        "patient_or_sample_identifiers_exported": False,
        "summary_wall_clock_seconds": time.perf_counter() - started,
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
