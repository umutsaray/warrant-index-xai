"""Integrity checks for the Phase 7.4 tau/eta sensitivity exports."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_tau_reference_reproduces_phase72_and_point_counts_are_monotone() -> None:
    table = pd.read_csv(ROOT / "PHASE74_TAU_SENSITIVITY.csv")
    assert table.tau.tolist() == [0.70, 0.75, 0.80, 0.85, 0.90]
    reference = table.loc[table.tau == 0.80].iloc[0]
    assert tuple(reference[["CoxPH_point", "CoxPH_BH", "CoxPH_BY"]].astype(int)) == (9, 8, 8)
    assert tuple(reference[["RSF_point", "RSF_BH", "RSF_BY"]].astype(int)) == (17, 12, 12)
    assert tuple(reference[["DeepSurv_point", "DeepSurv_BH", "DeepSurv_BY"]].astype(int)) == (0, 0, 0)
    for column in ("CoxPH_point", "RSF_point", "DeepSurv_point"):
        assert table[column].is_monotonic_decreasing
    low = table.loc[table.tau == 0.70].iloc[0]
    assert tuple(low[["DeepSurv_point", "DeepSurv_BH", "DeepSurv_BY"]].astype(int)) == (21, 6, 6)


def test_eta_reference_reproduces_rank_bottlenecks_without_refits() -> None:
    table = pd.read_csv(ROOT / "PHASE74_ETA_SENSITIVITY.csv")
    reference = table.loc[table.eta == 0.02].set_index("model")
    assert reference.R_bottleneck_count.astype(int).to_dict() == {
        "CoxPH": 24, "RSF": 24, "DeepSurv": 27
    }
    wide = table.loc[table.eta == 0.10].set_index("model")
    assert wide.R_bottleneck_count.astype(int).to_dict() == {
        "CoxPH": 4, "RSF": 25, "DeepSurv": 25
    }
    feature = pd.read_csv(ROOT / "PHASE74_ETA_FEATURE_VALUES.csv")
    assert len(feature) == 4 * 3 * 27
    metadata = json.loads((ROOT / "results/phase74_metadata.json").read_text())
    assert metadata["feasibility"]["feasible_without_patient_bootstrap_or_model_explainer"]
    assert not metadata["model_or_explainer_calls"]
    assert not metadata["new_patient_cluster_bootstrap_fits"]
