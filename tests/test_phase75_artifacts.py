import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_fold1_panel_and_privacy_contract():
    meta = json.loads((ROOT / "results/phase75_fold1_panel/metadata.json").read_text())
    assert meta["new"]["common_features"] == 80
    assert meta["new"]["tcga_audit_patients_after_one_sample_per_patient"] == 153
    assert meta["new"]["metabric_complete_patients"] == 1466
    assert meta["privacy"]["identifiers_exported"] is False
    with np.load(ROOT / "results/phase75_fold1_panel/query_matrices.npz") as panel:
        assert panel["source_query12"].shape == (12, 80)
        assert panel["external_query12"].shape == (12, 80)


def test_available_model_gate_and_certification_counts():
    gate = json.loads((ROOT / "results/phase75_fold1_calibration/metadata.json").read_text())
    assert gate["all_available_models_gate_passed"] is True
    assert set(gate["nuisance"]) == {"RSF_N50", "DeepSurv_N50"}
    summary = json.loads(
        (ROOT / "results/phase75_fold1_calibration/certification_summary.json").read_text()
    )["counts"]
    assert summary["RSF"] == {"family_m": 80, "point": 10, "BH": 7, "BY": 7}
    assert summary["DeepSurv"] == {"family_m": 80, "point": 7, "BH": 0, "BY": 0}
    assert summary["CoxPH"]["status"].startswith("NA:")


def test_cox_failure_is_preserved_and_not_certified():
    failure = json.loads(
        (ROOT / "results/phase75_fold1_cox/failed_replicate_15.json").read_text()
    )
    assert failure["replicate"] == 15
    assert [item["fit_options"] for item in failure["fit_attempts"]] == [
        {"max_steps": 50},
        {"max_steps": 200, "step_size": 0.25},
    ]
    assert "failed to converge sufficiently" in failure["final_warning_messages"][0]
    assert not (ROOT / "results/phase75_fold1_cox_certification").exists()


def test_formula_hashes_and_comparison_rows():
    expected = {
        "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
        "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    }
    actual = {
        name: hashlib.sha256((ROOT / "warrantxai" / name).read_bytes()).hexdigest()
        for name in expected
    }
    assert actual == expected
    comparison = pd.read_csv(ROOT / "PHASE75_MULTIFOLD_COMPARISON.csv")
    assert len(comparison) == 6
    fold1 = comparison[comparison.Fold == 1].set_index("Model")
    assert fold1.loc["RSF", "BH_certified"] == 7
    assert fold1.loc["DeepSurv", "BH_certified"] == 0
    assert pd.isna(fold1.loc["CoxPH-v2", "BH_certified"])
