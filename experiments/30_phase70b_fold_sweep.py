#!/usr/bin/env python3
"""Run the Phase 7.0 preprocessing/selection/CoxPH/RSF chain on all grouped folds."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SIBLING = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results" / "phase70b_collapse_diagnostic"
PHASE7 = ROOT / "results" / "phase7_retrain"


def run_logged(command: list[str], log_path: Path) -> None:
    started = time.perf_counter()
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(f"$ {' '.join(command)}\n")
        handle.write(completed.stdout)
        handle.write(f"\nexit={completed.returncode}; elapsed={time.perf_counter() - started:.3f}s\n")
    if completed.returncode:
        raise RuntimeError(f"Command failed; see {log_path}: {command}")


def score_record(fold: int, directory: Path, *, reused: bool) -> dict:
    baseline = json.loads((directory / "baseline_results.json").read_text(encoding="utf-8"))
    selection = json.loads((directory / "selected_features.json").read_text(encoding="utf-8"))
    split_rows = (directory / "patient_grouped_fold_candidates.csv")
    import pandas as pd

    candidate = pd.read_csv(split_rows).set_index("fold").loc[fold]
    return {
        "fold": fold,
        "reused_phase7_artifacts": reused,
        "train_rows": int(candidate["train_rows"]),
        "test_rows": int(candidate["test_rows"]),
        "train_event_rate": float(candidate["train_event_rate"]),
        "test_event_rate": float(candidate["test_event_rate"]),
        "selected_feature_count": int(selection["final_feature_count"]),
        "coxph_test_c_index": baseline["models"]["coxph"]["test_c_index"],
        "coxph_warning_count": len(baseline["models"]["coxph"]["warnings"]),
        "coxph_excluded_columns": baseline["models"]["coxph"]["complete_separation_columns_excluded"],
        "rsf_test_c_index": baseline["models"]["random_survival_forest"]["test_c_index"],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    records = [score_record(0, PHASE7, reused=True)]
    for fold in range(1, 5):
        print(f"fold {fold}: preprocessing", flush=True)
        directory = OUT / f"fold_{fold}"
        directory.mkdir(parents=True, exist_ok=True)
        log_path = directory / "pipeline.log"
        log_path.write_text("", encoding="utf-8")
        run_logged(
            [
                sys.executable,
                str(SIBLING / "impute_and_split.py"),
                "--input-dir",
                str(SIBLING),
                "--output-dir",
                str(directory),
                "--test-fold",
                str(fold),
            ],
            log_path,
        )
        print(f"fold {fold}: train-only feature selection", flush=True)
        run_logged(
            [sys.executable, str(SIBLING / "select_features.py"), "--input-dir", str(directory)],
            log_path,
        )
        # The full 72k-column intermediates are no longer needed after selection.
        for name in ("X_train.parquet", "X_test.parquet"):
            (directory / name).unlink()
        print(f"fold {fold}: CoxPH + RSF", flush=True)
        run_logged(
            [sys.executable, str(SIBLING / "train_baselines.py"), "--input-dir", str(directory)],
            log_path,
        )
        records.append(score_record(fold, directory, reused=False))
        # Diagnostic folds retain scores/data but not duplicate large fitted models.
        for name in ("coxph_baseline.joblib", "rsf_baseline.joblib"):
            (directory / name).unlink()
    payload = {
        "scope": "five-fold grouped split diagnostic; CoxPH and RSF only",
        "folds": records,
        "wall_clock_seconds": time.perf_counter() - started,
        "audit_pipeline_run": False,
    }
    (OUT / "fold_scores.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    main()
