#!/usr/bin/env python3
"""Diagnose DeepSurv checkpoint selection and internal patient overlap."""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import torch
from sklearn.model_selection import StratifiedGroupKFold, train_test_split


ROOT = Path(__file__).resolve().parents[1]
SIBLING = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
DATA = ROOT / "results" / "phase7_retrain"
OUT = ROOT / "results" / "phase70b_collapse_diagnostic"


def load_training_module():
    spec = importlib.util.spec_from_file_location("phase7_train_dl", SIBLING / "train_explainable_dl.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_frame(name):
    return pq.read_table(DATA / name).to_pandas().set_index("sample_key")


def grouped_split(y):
    groups = y.index.astype(str).str[:12].to_numpy()
    splitter = StratifiedGroupKFold(n_splits=7, shuffle=True, random_state=42)
    candidates = list(splitter.split(np.arange(len(y)), y["event"].to_numpy(), groups))
    overall = float(y["event"].mean())
    return min(
        candidates,
        key=lambda split: (
            abs(len(split[1]) / len(y) - 0.15),
            abs(float(y.iloc[split[1]]["event"].mean()) - overall),
        ),
    )


def run_training(module, x_train, x_test, y_train, y_test, fit_idx, val_idx, label):
    module.set_seed(42)
    columns_by_block = {
        block: [c for c in x_train.columns if c.startswith(prefix)]
        for block, prefix in module.PREFIXES.items()
    }
    ordered = [c for block in module.PREFIXES for c in columns_by_block[block]]
    x_train, x_test = x_train[ordered], x_test[ordered]
    slices, cursor = {}, 0
    for block in module.PREFIXES:
        end = cursor + len(columns_by_block[block])
        slices[block] = (cursor, end)
        cursor = end
    fit_values = x_train.iloc[fit_idx].to_numpy(dtype=np.float32)
    val_values = x_train.iloc[val_idx].to_numpy(dtype=np.float32)
    test_values = x_test.to_numpy(dtype=np.float32)
    means = fit_values.mean(axis=0, dtype=np.float64).astype(np.float32)
    scales = fit_values.std(axis=0, dtype=np.float64).astype(np.float32)
    scales[scales == 0] = 1.0
    fit_values = (fit_values - means) / scales
    val_values = (val_values - means) / scales
    test_values = (test_values - means) / scales
    device = torch.device("cpu")
    fit_blocks = module.make_blocks(fit_values, slices, device)
    val_blocks = module.make_blocks(val_values, slices, device)
    test_blocks = module.make_blocks(test_values, slices, device)
    fit_time = torch.as_tensor(y_train.iloc[fit_idx]["time"].to_numpy(), dtype=torch.float32)
    fit_event = torch.as_tensor(y_train.iloc[fit_idx]["event"].to_numpy(), dtype=torch.float32)
    model = module.OmicsDeepSurv({b: len(columns_by_block[b]) for b in module.PREFIXES}, dropout=0.3)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=0.0001)
    best_val, best_epoch, best_state, wait = -np.inf, 0, None, 0
    history = []
    targets = {1, 10, 26, 38, 50}
    started = time.perf_counter()
    for epoch in range(1, 301):
        model.train()
        optimizer.zero_grad()
        risk = model(fit_blocks)
        loss = module.cox_negative_log_likelihood(risk, fit_time, fit_event)
        loss.backward()
        optimizer.step()
        fit_risk = module.predict_risk(model, fit_blocks)
        val_risk = module.predict_risk(model, val_blocks)
        test_risk = module.predict_risk(model, test_blocks)
        fit_c = module.cindex(
            y_train.iloc[fit_idx]["event"].to_numpy(), y_train.iloc[fit_idx]["time"].to_numpy(), fit_risk
        )
        val_c = module.cindex(
            y_train.iloc[val_idx]["event"].to_numpy(), y_train.iloc[val_idx]["time"].to_numpy(), val_risk
        )
        test_c = module.cindex(y_test["event"].to_numpy(), y_test["time"].to_numpy(), test_risk)
        history.append(
            {"epoch": epoch, "loss": float(loss.item()), "fit_c_index": fit_c, "validation_c_index": val_c, "test_c_index": test_c}
        )
        if val_c > best_val + 1e-6:
            best_val, best_epoch = val_c, epoch
            best_state = {key: value.detach().clone() for key, value in model.state_dict().items()}
            wait = 0
        else:
            wait += 1
        if wait >= 30:
            break
    model.load_state_dict(best_state)
    best_test = module.cindex(
        y_test["event"].to_numpy(), y_test["time"].to_numpy(), module.predict_risk(model, test_blocks)
    )
    fit_patients = set(y_train.iloc[fit_idx].index.astype(str).str[:12])
    val_patients = set(y_train.iloc[val_idx].index.astype(str).str[:12])
    return {
        "label": label,
        "fit_rows": len(fit_idx),
        "validation_rows": len(val_idx),
        "fit_validation_patient_overlap": len(fit_patients & val_patients),
        "best_epoch": best_epoch,
        "stopped_epoch": history[-1]["epoch"],
        "best_validation_c_index": best_val,
        "best_checkpoint_test_c_index": best_test,
        "wall_clock_seconds": time.perf_counter() - started,
        "selected_epoch_metrics": [row for row in history if row["epoch"] in targets or row["epoch"] == best_epoch],
        "final_epoch_metrics": history[-1],
        "history": history,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    module = load_training_module()
    x_train, x_test = read_frame("X_train_selected.parquet"), read_frame("X_test_selected.parquet")
    y_train, y_test = read_frame("y_train.parquet"), read_frame("y_test.parquet")
    indices = np.arange(len(x_train))
    row_fit, row_val = train_test_split(
        indices, test_size=0.15, random_state=42, stratify=y_train["event"].to_numpy()
    )
    group_fit, group_val = grouped_split(y_train)
    payload = {
        "stored_intermediate_checkpoints_found": False,
        "diagnostic_retraining_reason": "No intermediate checkpoints were stored; deterministic CPU retraining was measured below one second in Phase 7.0.",
        "row_level_validation": run_training(
            module, x_train, x_test, y_train, y_test, row_fit, row_val, "frozen row-level validation"
        ),
        "patient_grouped_validation": run_training(
            module, x_train, x_test, y_train, y_test, group_fit, group_val, "diagnostic patient-grouped validation"
        ),
    }
    (OUT / "deepsurv_diagnostic.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
