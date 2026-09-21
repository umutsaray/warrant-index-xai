#!/usr/bin/env python3
"""Inspect Phase 7.0 CoxPH coefficients and retry optimizer convergence."""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from lifelines import CoxPHFitter
from sksurv.metrics import concordance_index_censored


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "phase7_retrain"
OUT = ROOT / "results" / "phase70b_collapse_diagnostic"


def read_frame(name: str) -> pd.DataFrame:
    return pq.read_table(DATA / name).to_pandas().set_index("sample_key")


def coefficient_summary(values: np.ndarray) -> dict:
    absolute = np.abs(values)
    return {
        "minimum": float(values.min()),
        "maximum": float(values.max()),
        "maximum_absolute": float(absolute.max()),
        "median_absolute": float(np.median(absolute)),
        "count_abs_gt_1": int((absolute > 1).sum()),
        "count_abs_gt_2": int((absolute > 2).sum()),
        "count_abs_gt_5": int((absolute > 5).sum()),
        "hazard_ratio_minimum": float(np.exp(values).min()),
        "hazard_ratio_maximum": float(np.exp(values).max()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    primary_results = json.loads((DATA / "baseline_results.json").read_text(encoding="utf-8"))
    stored = joblib.load(DATA / "coxph_baseline.joblib")
    x_train, x_test = read_frame("X_train_selected.parquet"), read_frame("X_test_selected.parquet")
    y_train, y_test = read_frame("y_train.parquet"), read_frame("y_test.parquet")
    features = stored["features"]
    scaler = stored["scaler"]
    train_scaled = scaler.transform(x_train[features]).astype("float64")
    test_scaled = scaler.transform(x_test[features]).astype("float64")
    frame = pd.DataFrame(train_scaled, columns=features)
    frame["time"] = y_train["time"].to_numpy()
    frame["event"] = y_train["event"].to_numpy()

    started = time.perf_counter()
    retry = CoxPHFitter(penalizer=0.1)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        retry.fit(
            frame,
            duration_col="time",
            event_col="event",
            show_progress=False,
            fit_options={"max_steps": 200, "step_size": 0.25},
        )
    retry_warnings = [str(item.message) for item in caught]
    retry_risk = retry.predict_partial_hazard(pd.DataFrame(test_scaled, columns=features)).to_numpy()
    retry_cindex = float(
        concordance_index_censored(
            y_test["event"].astype(bool).to_numpy(), y_test["time"].to_numpy(), retry_risk
        )[0]
    )
    original_coefficients = stored["model"].params_.to_numpy(dtype=float)
    retry_coefficients = retry.params_.to_numpy(dtype=float)
    largest = (
        pd.DataFrame(
            {
                "feature": features,
                "original_coefficient": original_coefficients,
                "retry_coefficient": retry_coefficients,
            }
        )
        .assign(original_absolute=lambda x: x.original_coefficient.abs())
        .sort_values("original_absolute", ascending=False)
        .head(20)
    )
    largest.to_csv(OUT / "cox_largest_coefficients.csv", index=False)
    payload = {
        "primary": {
            "test_c_index": primary_results["models"]["coxph"]["test_c_index"],
            "warnings": primary_results["models"]["coxph"]["warnings"],
            "excluded_columns": primary_results["models"]["coxph"]["complete_separation_columns_excluded"],
            "feature_count": len(features),
            "coefficients": coefficient_summary(original_coefficients),
        },
        "retry": {
            "penalizer": 0.1,
            "max_steps": 200,
            "step_size": 0.25,
            "test_c_index": retry_cindex,
            "warnings": retry_warnings,
            "coefficients": coefficient_summary(retry_coefficients),
            "elapsed_seconds": time.perf_counter() - started,
        },
        "maximum_absolute_coefficient_change": float(
            np.max(np.abs(original_coefficients - retry_coefficients))
        ),
    }
    (OUT / "cox_diagnostic.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
