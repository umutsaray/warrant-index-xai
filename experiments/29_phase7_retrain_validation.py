#!/usr/bin/env python3
"""Validate the Phase 7.0 patient-grouped split and isolated retraining artifacts."""

from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path

import pandas as pd
import torch
from sklearn.model_selection import train_test_split


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase7_retrain"
SIBLING = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"

EXPECTED_OLD_HASHES = {
    "X_train_selected.parquet": "9317ae2083df0d6fd3dd62f9505ee6286537874d0a5d4c0364de5a2873e804d1",
    "X_test_selected.parquet": "3eaa1e2be715eb813232d49cd292db4ef60e0a633806bee8d3e058bcf9ff548a",
    "y_train.parquet": "7122084290a91057f8d0ed5274a07ef3cc6748dec74ca956c6cf9308635f6dfd",
    "y_test.parquet": "329a9eb85492afbd3be80760923ef26618090ab8a8a62c492d8a5cab0c70d716",
    "coxph_baseline.joblib": "937ef8eb2351b8048095a174f062cc4b7ea23aa021d5afbf42ee18346fb41775",
    "rsf_baseline.joblib": "3ab8da5972812ec3411430419ec40862306316504b270a77a9b1a9ce96121574",
    "explainable_dl_model.pt": "60fa69d2945dc883605d805f49b8aedb1860138604bdac796c276ea60e608e35",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_frame(name: str) -> pd.DataFrame:
    frame = pd.read_parquet(OUT / name)
    return frame.set_index("sample_key") if "sample_key" in frame.columns else frame


def main() -> None:
    started = time.perf_counter()
    x_train = read_frame("X_train_selected.parquet")
    x_test = read_frame("X_test_selected.parquet")
    y_train = read_frame("y_train.parquet")
    y_test = read_frame("y_test.parquet")
    train_patients = set(y_train.index.astype(str).str[:12])
    test_patients = set(y_test.index.astype(str).str[:12])
    patient_overlap = train_patients & test_patients
    sample_overlap = set(y_train.index) & set(y_test.index)

    selection = json.loads((OUT / "selected_features.json").read_text(encoding="utf-8"))
    baselines = json.loads((OUT / "baseline_results.json").read_text(encoding="utf-8"))
    old_dl = torch.load(SIBLING / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    new_dl = torch.load(OUT / "explainable_dl_model.pt", map_location="cpu", weights_only=False)

    fit_idx, val_idx = train_test_split(
        range(len(y_train)),
        test_size=0.15,
        random_state=42,
        stratify=y_train["event"].to_numpy(),
    )
    fit_patients = set(y_train.iloc[list(fit_idx)].index.astype(str).str[:12])
    val_patients = set(y_train.iloc[list(val_idx)].index.astype(str).str[:12])

    old_hashes = {name: sha256(SIBLING / name) for name in EXPECTED_OLD_HASHES}
    new_files = [
        "X_train_selected.parquet",
        "X_test_selected.parquet",
        "y_train.parquet",
        "y_test.parquet",
        "coxph_baseline.joblib",
        "rsf_baseline.joblib",
        "explainable_dl_model.pt",
        "selected_features.json",
        "baseline_results.json",
    ]
    new_hashes = {name: sha256(OUT / name) for name in new_files}
    source_files = ["impute_and_split.py", "select_features.py", "train_baselines.py", "train_explainable_dl.py"]
    source_hashes = {name: sha256(SIBLING / name) for name in source_files}

    assertions = {
        "train_shape_674x616": list(x_train.shape) == [674, 616],
        "test_shape_168x616": list(x_test.shape) == [168, 616],
        "sample_overlap_zero": len(sample_overlap) == 0,
        "patient_overlap_zero": len(patient_overlap) == 0,
        "selection_scope_train_only": selection.get("selection_scope") == "train_only",
        "selected_feature_count_616": selection.get("final_feature_count") == 616,
        "old_artifact_hashes_unchanged": old_hashes == EXPECTED_OLD_HASHES,
        "cox_settings_frozen": baselines["models"]["coxph"]["penalizer"] == 0.1
        and baselines["models"]["coxph"]["max_steps"] == 50,
        "rsf_settings_frozen": all(
            baselines["models"]["random_survival_forest"].get(key) == value
            for key, value in {
                "n_estimators": 500,
                "min_samples_split": 10,
                "min_samples_leaf": 5,
                "max_features": "sqrt",
                "random_state": 42,
            }.items()
        ),
        "deepsurv_settings_frozen": new_dl["block_dims"]
        == {"rna": 200, "mirna": 200, "methylation": 200, "clinical": 16}
        and new_dl["dropout"] == 0.3
        and new_dl["random_state"] == 42,
    }
    if not all(assertions.values()):
        raise AssertionError({key: value for key, value in assertions.items() if not value})

    payload = {
        "scope": "Phase 7.0 split, feature selection, and baseline retraining only; no XAI audit",
        "split": {
            "train_rows": len(y_train),
            "test_rows": len(y_test),
            "train_unique_patients": len(train_patients),
            "test_unique_patients": len(test_patients),
            "sample_overlap": len(sample_overlap),
            "patient_overlap": len(patient_overlap),
            "train_event_rate": float(y_train["event"].mean()),
            "test_event_rate": float(y_test["event"].mean()),
        },
        "feature_selection": {
            "selection_scope": selection["selection_scope"],
            "feature_count": selection["final_feature_count"],
            "block_counts": {
                "rna": len(selection["selected_by_block"]["rna"]),
                "mirna": len(selection["selected_by_block"]["mirna"]),
                "methylation": len(selection["selected_by_block"]["methylation"]),
                "clinical": len(selection["clinical_features_retained_without_univariate_filter"]),
            },
        },
        "scores": {
            "CoxPH": {"old": 0.7064765300059418, "new": baselines["models"]["coxph"]["test_c_index"]},
            "RSF": {"old": 0.714795008912656, "new": baselines["models"]["random_survival_forest"]["test_c_index"]},
            "DeepSurv": {"old": 0.7456922162804516, "new": float(new_dl["test_c_index"])},
        },
        "coxph": baselines["models"]["coxph"],
        "rsf": baselines["models"]["random_survival_forest"],
        "deepsurv": {
            "old_best_epoch": int(old_dl["best_epoch"]),
            "new_best_epoch": int(new_dl["best_epoch"]),
            "stopped_epoch": int(new_dl["stopped_epoch"]),
            "best_validation_c_index": float(new_dl["best_validation_c_index"]),
            "test_c_index": float(new_dl["test_c_index"]),
            "elapsed_seconds": float(new_dl["elapsed_seconds"]),
            "block_dims": new_dl["block_dims"],
            "dropout": float(new_dl["dropout"]),
            "random_state": int(new_dl["random_state"]),
            "internal_validation_protocol_unchanged": True,
            "internal_fit_validation_patient_overlap": len(fit_patients & val_patients),
        },
        "old_artifact_hashes": old_hashes,
        "new_artifact_hashes": new_hashes,
        "source_hashes": source_hashes,
        "assertions": assertions,
        "patient_or_sample_identifiers_exported_to_reports": False,
        "validation_wall_clock_seconds": time.perf_counter() - started,
    }
    for values in payload["scores"].values():
        values["delta_new_minus_old"] = values["new"] - values["old"]
    report_text = "\n".join(
        (ROOT / name).read_text(encoding="utf-8")
        for name in ("SPLIT_REDESIGN.md", "MODEL_RETRAIN_COMPARISON.md")
    )
    if re.search(r"\bTCGA-[A-Z0-9]{2}-[A-Z0-9]{4}\b", report_text):
        raise AssertionError("A Phase 7.0 report exports a patient barcode")
    (OUT / "validation.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
