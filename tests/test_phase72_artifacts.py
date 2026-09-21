"""Integrity checks for the clean patient-grouped Phase 7.2 audit artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}


def test_phase72_panel_and_point_audits_are_aligned() -> None:
    panel = json.loads((ROOT / "results/phase72_panel/metadata.json").read_text())
    m = panel["new"]["common_features"]
    assert panel["new"]["tcga_audit_patients_after_one_sample_per_patient"] == 153
    with np.load(ROOT / "results/phase72_panel/query_matrices.npz", allow_pickle=False) as saved:
        assert saved["source_query12"].shape == (12, m)
        assert saved["external_query12"].shape == (12, m)
    for model in ("cox", "rsf", "deepsurv"):
        report = pd.read_csv(ROOT / f"results/phase72_{model}/feature_report.csv")
        assert len(report) == m
        assert report.feature.tolist() == panel["active_features"]
        assert np.isfinite(report[["F", "S", "U", "G_claim", "R", "WI"]]).all().all()
        assert not any("patient" in column.lower() or "sample" in column.lower() for column in report.columns)


def test_phase72_frozen_hashes_and_deepsurv_additivity() -> None:
    actual = {
        name: hashlib.sha256((ROOT / "warrantxai" / name).read_bytes()).hexdigest()
        for name in EXPECTED
    }
    assert actual == EXPECTED
    deep = json.loads((ROOT / "results/phase72_deepsurv/metadata.json").read_text())
    assert deep["additivity"]["max_candidate_relative_error"] <= 0.05
    assert deep["additivity"]["max_reference_relative_error"] <= 0.05


def test_phase72_certification_obeys_joint_gate() -> None:
    calibration = json.loads((ROOT / "results/phase72_calibration/metadata.json").read_text())
    for model in ("coxph", "rsf", "deepsurv"):
        certification = ROOT / f"results/phase72_{model}_certification/feature_certification.csv"
        assert certification.exists() == calibration["all_models_gate_passed"]
        if certification.exists():
            table = pd.read_csv(certification)
            component_p = table[[f"p_{name}" for name in ("F", "S", "U", "G_claim", "R")]]
            np.testing.assert_allclose(table.p_joint, component_p.max(axis=1))
