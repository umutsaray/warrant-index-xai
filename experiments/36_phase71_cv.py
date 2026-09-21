#!/usr/bin/env python3
"""Complete frozen three-model five-fold grouped CV and summarize scores."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import torch


ROOT = Path(__file__).resolve().parents[1]
SIBLING = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
PHASE7 = ROOT / "results" / "phase7_retrain"
DIAGNOSTIC = ROOT / "results" / "phase70b_collapse_diagnostic"
OUT = ROOT / "results" / "phase71_cv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize(values: list[float]) -> dict:
    return {
        "scores": values,
        "mean": sum(values) / len(values),
        "minimum": min(values),
        "maximum": max(values),
        "range_width": max(values) - min(values),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    protected_names = [
        "coxph_baseline.joblib",
        "coxph_baseline_v2.joblib",
        "rsf_baseline.joblib",
        "explainable_dl_model.pt",
    ]
    protected_before = {name: sha256(PHASE7 / name) for name in protected_names}

    fold_table = pd.read_csv(DIAGNOSTIC / "fold_scores.csv").sort_values("fold")
    deep_records = []
    fold0 = torch.load(PHASE7 / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    deep_records.append(
        {
            "fold": 0,
            "test_c_index": float(fold0["test_c_index"]),
            "best_validation_c_index": float(fold0["best_validation_c_index"]),
            "best_epoch": int(fold0["best_epoch"]),
            "stopped_epoch": int(fold0["stopped_epoch"]),
            "elapsed_seconds": float(fold0["elapsed_seconds"]),
            "reference_model_reused": True,
        }
    )
    for fold in range(1, 5):
        directory = DIAGNOSTIC / f"fold_{fold}"
        command = [sys.executable, str(SIBLING / "train_explainable_dl.py"), "--input-dir", str(directory)]
        print(f"DeepSurv fold {fold}", flush=True)
        run_started = time.perf_counter()
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (OUT / f"deepsurv_fold_{fold}.log").write_text(completed.stdout, encoding="utf-8")
        if completed.returncode:
            raise RuntimeError(f"DeepSurv fold {fold} failed")
        checkpoint = torch.load(directory / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
        if checkpoint["block_dims"] != {"rna": 200, "mirna": 200, "methylation": 200, "clinical": 16}:
            raise AssertionError(f"Unexpected block dimensions in fold {fold}")
        if checkpoint["dropout"] != 0.3 or checkpoint["random_state"] != 42:
            raise AssertionError(f"Frozen DeepSurv settings changed in fold {fold}")
        deep_records.append(
            {
                "fold": fold,
                "test_c_index": float(checkpoint["test_c_index"]),
                "best_validation_c_index": float(checkpoint["best_validation_c_index"]),
                "best_epoch": int(checkpoint["best_epoch"]),
                "stopped_epoch": int(checkpoint["stopped_epoch"]),
                "elapsed_seconds": float(checkpoint["elapsed_seconds"]),
                "process_wall_clock_seconds": time.perf_counter() - run_started,
                "reference_model_reused": False,
            }
        )

    deep_table = pd.DataFrame(deep_records).sort_values("fold")
    combined = fold_table[
        [
            "fold",
            "train_rows",
            "test_rows",
            "train_event_rate",
            "test_event_rate",
            "coxph_converged_c_index",
            "rsf_test_c_index",
        ]
    ].merge(
        deep_table[["fold", "test_c_index", "best_validation_c_index", "best_epoch", "stopped_epoch"]],
        on="fold",
        validate="one_to_one",
    ).rename(columns={"test_c_index": "deepsurv_test_c_index"})
    combined.to_csv(OUT / "cv_fold_scores.csv", index=False)
    deep_table.to_csv(OUT / "deepsurv_fold_details.csv", index=False)

    protected_after = {name: sha256(PHASE7 / name) for name in protected_names}
    if protected_after != protected_before:
        raise AssertionError("A fold-0 reference artifact changed during CV")
    payload = {
        "design": "five-fold patient-grouped CV; train-only preprocessing/selection within every fold",
        "models": {
            "CoxPH": summarize(combined.coxph_converged_c_index.tolist()),
            "RSF": summarize(combined.rsf_test_c_index.tolist()),
            "DeepSurv": summarize(combined.deepsurv_test_c_index.tolist()),
        },
        "deepsurv_fold_details": deep_records,
        "coxph_retry_rule": {
            "trigger": "any convergence warning",
            "penalizer": 0.1,
            "max_steps": 200,
            "step_size": 0.25,
            "retried_folds": [0, 2],
            "all_final_cv_fits_warning_free": True,
        },
        "deep_surv_frozen_settings": {
            "architecture": "RNA 200-128-32; miRNA 200-64-16; methylation 200-128-32; clinical 16-32-8; concat 88-64-16-1",
            "dropout": 0.3,
            "optimizer": "AdamW",
            "learning_rate": 0.001,
            "weight_decay": 0.0001,
            "max_epochs": 300,
            "patience": 30,
            "seed": 42,
        },
        "source_hashes": {
            "train_explainable_dl.py": sha256(SIBLING / "train_explainable_dl.py"),
            "train_baselines.py": sha256(SIBLING / "train_baselines.py"),
            "select_features.py": sha256(SIBLING / "select_features.py"),
            "impute_and_split.py": sha256(SIBLING / "impute_and_split.py"),
        },
        "protected_reference_hashes_before": protected_before,
        "protected_reference_hashes_after": protected_after,
        "protected_reference_artifacts_unchanged": True,
        "inherited_phase70b_fold_sweep_seconds": json.loads(
            (DIAGNOSTIC / "fold_scores.json").read_text(encoding="utf-8")
        )["wall_clock_seconds"],
        "phase71_incremental_wall_clock_seconds": time.perf_counter() - started,
        "audit_pipeline_run": False,
        "patient_or_sample_identifiers_exported": False,
    }
    (OUT / "cv_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    main()
