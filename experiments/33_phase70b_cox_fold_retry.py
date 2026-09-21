#!/usr/bin/env python3
"""Retry warning-bearing grouped-fold CoxPH fits with the frozen remediation protocol."""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from lifelines import CoxPHFitter
from sklearn.preprocessing import StandardScaler
from sksurv.metrics import concordance_index_censored


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase70b_collapse_diagnostic"
PHASE7 = ROOT / "results" / "phase7_retrain"


def read_frame(directory: Path, name: str) -> pd.DataFrame:
    return pq.read_table(directory / name).to_pandas().set_index("sample_key")


def retry(directory: Path) -> dict:
    baseline = json.loads((directory / "baseline_results.json").read_text(encoding="utf-8"))
    excluded = baseline["models"]["coxph"]["complete_separation_columns_excluded"]
    x_train, x_test = read_frame(directory, "X_train_selected.parquet"), read_frame(directory, "X_test_selected.parquet")
    y_train, y_test = read_frame(directory, "y_train.parquet"), read_frame(directory, "y_test.parquet")
    features = [column for column in x_train if column not in excluded]
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(x_train[features]).astype("float64")
    test_scaled = scaler.transform(x_test[features]).astype("float64")
    frame = pd.DataFrame(train_scaled, columns=features)
    frame["time"] = y_train["time"].to_numpy()
    frame["event"] = y_train["event"].to_numpy()
    model = CoxPHFitter(penalizer=0.1)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model.fit(
            frame,
            duration_col="time",
            event_col="event",
            show_progress=False,
            fit_options={"max_steps": 200, "step_size": 0.25},
        )
    risk = model.predict_partial_hazard(pd.DataFrame(test_scaled, columns=features)).to_numpy()
    score = float(
        concordance_index_censored(
            y_test["event"].astype(bool).to_numpy(), y_test["time"].to_numpy(), risk
        )[0]
    )
    coefficients = model.params_.to_numpy(float)
    return {
        "original_c_index": baseline["models"]["coxph"]["test_c_index"],
        "original_warnings": baseline["models"]["coxph"]["warnings"],
        "retry_c_index": score,
        "retry_warnings": [str(item.message) for item in caught],
        "max_abs_coefficient": float(np.abs(coefficients).max()),
        "hazard_ratio_range": [float(np.exp(coefficients).min()), float(np.exp(coefficients).max())],
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    fold_scores = json.loads((OUT / "fold_scores.json").read_text(encoding="utf-8"))
    results = {}
    for record in fold_scores["folds"]:
        if record["coxph_warning_count"]:
            fold = record["fold"]
            directory = PHASE7 if fold == 0 else OUT / f"fold_{fold}"
            print(f"retry fold {fold}", flush=True)
            results[str(fold)] = retry(directory)
    (OUT / "cox_fold_retries.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
