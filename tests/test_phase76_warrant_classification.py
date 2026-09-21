import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_fold0_predeclared_counts_and_exact_zero_rule():
    table = pd.read_csv(ROOT / "PHASE76_WARRANT_CLASSIFICATION_FOLD0.csv")
    counts = table.groupby(["Model", "Label"]).size().to_dict()
    assert counts == {
        ("CoxPH-v2", "informative"): 8,
        ("RSF", "degenerate"): 8,
        ("RSF", "informative"): 4,
    }
    assert (table.loc[table.Label == "degenerate", "Mean signed attribution"] == 0.0).all()
    assert (table.loc[table.Label == "informative", "Mean signed attribution"] != 0.0).all()
    assert len(table) == 20
    assert table["Feature ID"].nunique() == 14


def test_rsf_split_usage_records_observed_inconsistency():
    split = pd.read_csv(ROOT / "results/phase76_rsf_split_usage.csv")
    degenerate = split[split.Label == "degenerate"].set_index("Feature")
    assert len(degenerate) == 8
    assert int((degenerate["Trees with ≥1 split"] == 0).sum()) == 6
    assert degenerate.loc["OR6C1", "Trees with ≥1 split"] == 10
    assert degenerate.loc["AJCC Stage IV", "Trees with ≥1 split"] == 25
    assert split["feature_importances_"].str.contains("NotImplementedError").all()


def test_fold1_rsf_counts():
    table = pd.read_csv(ROOT / "PHASE76_WARRANT_CLASSIFICATION_FOLD1_RSF.csv")
    assert table.Label.value_counts().to_dict() == {"informative": 5, "degenerate": 2}
    assert len(table) == 7


def test_formula_hashes_unchanged():
    expected = {
        "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
        "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    }
    assert {
        name: hashlib.sha256((ROOT / "warrantxai" / name).read_bytes()).hexdigest()
        for name in expected
    } == expected
