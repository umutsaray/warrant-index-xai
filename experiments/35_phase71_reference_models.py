#!/usr/bin/env python3
"""Freeze the clean fold-0 reference model set and create CoxPH v2 via optimizer retry."""

from __future__ import annotations

import hashlib
import json
import time
import warnings
from pathlib import Path

import joblib
import pandas as pd
import pyarrow.parquet as pq
import torch
from lifelines import CoxPHFitter
from sksurv.metrics import concordance_index_censored


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "phase7_retrain"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_frame(name: str) -> pd.DataFrame:
    return pq.read_table(DATA / name).to_pandas().set_index("sample_key")


def main() -> None:
    started = time.perf_counter()
    protected = ["coxph_baseline.joblib", "rsf_baseline.joblib", "explainable_dl_model.pt"]
    before = {name: sha256(DATA / name) for name in protected}
    v1 = joblib.load(DATA / "coxph_baseline.joblib")
    features, scaler = v1["features"], v1["scaler"]
    x_train, x_test = read_frame("X_train_selected.parquet"), read_frame("X_test_selected.parquet")
    y_train, y_test = read_frame("y_train.parquet"), read_frame("y_test.parquet")
    train_scaled = scaler.transform(x_train[features]).astype("float64")
    test_scaled = scaler.transform(x_test[features]).astype("float64")
    frame = pd.DataFrame(train_scaled, columns=features)
    frame["time"] = y_train["time"].to_numpy()
    frame["event"] = y_train["event"].to_numpy()
    model = CoxPHFitter(penalizer=0.1)
    fit_started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model.fit(
            frame,
            duration_col="time",
            event_col="event",
            show_progress=False,
            fit_options={"max_steps": 200, "step_size": 0.25},
        )
    fit_seconds = time.perf_counter() - fit_started
    warning_messages = [str(item.message) for item in caught]
    risk = model.predict_partial_hazard(pd.DataFrame(test_scaled, columns=features)).to_numpy()
    score = float(
        concordance_index_censored(
            y_test["event"].astype(bool).to_numpy(), y_test["time"].to_numpy(), risk
        )[0]
    )
    if warning_messages:
        raise RuntimeError(f"CoxPH v2 retry still warned: {warning_messages}")
    if abs(score - 0.667704280155642) > 1e-12:
        raise AssertionError(f"Unexpected CoxPH v2 score: {score}")
    v2_path = DATA / "coxph_baseline_v2.joblib"
    joblib.dump(
        {
            "model": model,
            "scaler": scaler,
            "features": features,
            "optimization": {"penalizer": 0.1, "max_steps": 200, "step_size": 0.25},
            "test_c_index": score,
        },
        v2_path,
    )
    after = {name: sha256(DATA / name) for name in protected}
    if after != before:
        raise AssertionError("A protected fold-0 v1/RSF/DeepSurv artifact changed")
    baseline = json.loads((DATA / "baseline_results.json").read_text(encoding="utf-8"))
    deep = torch.load(DATA / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    payload = {
        "selection_rule": "pre-specified grouped fold 0; no performance-based model/fold selection",
        "models": {
            "CoxPH": {
                "path": str(v2_path),
                "sha256": sha256(v2_path),
                "test_c_index": score,
                "version": "v2 optimizer-retry official reference",
                "penalizer": 0.1,
                "max_steps": 200,
                "step_size": 0.25,
                "warnings": warning_messages,
                "fit_seconds": fit_seconds,
                "v1_preserved_path": str(DATA / "coxph_baseline.joblib"),
                "v1_preserved_sha256": before["coxph_baseline.joblib"],
            },
            "RSF": {
                "path": str(DATA / "rsf_baseline.joblib"),
                "sha256": before["rsf_baseline.joblib"],
                "test_c_index": baseline["models"]["random_survival_forest"]["test_c_index"],
                "unchanged_from_phase70": True,
            },
            "DeepSurv": {
                "path": str(DATA / "explainable_dl_model.pt"),
                "sha256": before["explainable_dl_model.pt"],
                "test_c_index": float(deep["test_c_index"]),
                "best_epoch": int(deep["best_epoch"]),
                "unchanged_from_phase70": True,
            },
        },
        "protected_artifact_hashes_before": before,
        "protected_artifact_hashes_after": after,
        "protected_artifacts_unchanged": before == after,
        "audit_pipeline_run": False,
        "patient_or_sample_identifiers_exported": False,
        "wall_clock_seconds": time.perf_counter() - started,
    }
    (DATA / "reference_models.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
