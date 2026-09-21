import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_all_ten_degenerate_fold_feature_rows_are_cellwise_zero():
    table = pd.read_csv(ROOT / "PHASE77_CANCELLATION_CHECK.csv")
    degenerate = table[table["Phase 7.6 label"] == "degenerate"]
    assert len(degenerate) == 10
    assert degenerate["All cells exactly zero"].all()
    assert (degenerate["Non-zero cell count"] == 0).all()
    assert (degenerate["M_j mean absolute attribution"] == 0.0).all()
    assert (degenerate["P_j non-zero frequency"] == 0.0).all()
    assert (degenerate["Q_j max absolute attribution"] == 0.0).all()


def test_informative_comparison_rows_are_not_cellwise_zero():
    table = pd.read_csv(ROOT / "PHASE77_CANCELLATION_CHECK.csv")
    informative = table[table["Phase 7.6 label"] == "informative"]
    assert len(informative) == 12
    assert (~informative["All cells exactly zero"]).all()
    assert (informative["M_j mean absolute attribution"] > 0).all()
    assert (informative["P_j non-zero frequency"] > 0).all()


def test_full_fold0_table_shape_decisions_and_labels():
    table = pd.read_csv(ROOT / "PHASE77_FULL_FEATURE_TABLE_FOLD0.csv")
    assert len(table) == 81
    assert table["Feature ID"].nunique() == 27
    assert table.groupby("Model").size().to_dict() == {
        "CoxPH-v2": 27, "DeepSurv": 27, "RSF": 27,
    }
    certified = table[table["BH certified"] == "Evet"]
    assert certified.groupby(["Model", "Informative/Degenerate label"]).size().to_dict() == {
        ("CoxPH-v2", "informative"): 8,
        ("RSF", "degenerate"): 8,
        ("RSF", "informative"): 4,
    }
    assert (table.loc[table["BH certified"] != "Evet", "Informative/Degenerate label"]
            == "N/A (not certified)").all()
    assert table.loc[table.Model == "RSF", "RSF trees with split"].notna().all()
    assert table.loc[table.Model != "RSF", "RSF trees with split"].isna().all()


def test_formula_hashes_unchanged():
    expected = {
        "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
        "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    }
    assert {
        name: hashlib.sha256((ROOT / "warrantxai" / name).read_bytes()).hexdigest()
        for name in expected
    } == expected
