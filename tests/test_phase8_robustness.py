import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_query_variance_hypothesis_rejected_for_all_four_null_pairs():
    table = pd.read_csv(ROOT / "PHASE8_QUERY_VARIANCE_CHECK.csv")
    tested = table[table.Role == "degenerate"]
    assert len(tested) == 4
    assert (tested["Variance (ddof=0)"] > 0).all()
    assert tested["Hypothesis result"].str.startswith("rejected").all()
    assert table["All 24 values"].map(lambda value: len(json.loads(value)) == 24).all()


def test_retry_partition_and_nominal_flags():
    table = pd.read_csv(ROOT / "PHASE8_COXPH_RETRY_BIAS_CHECK.csv")
    assert len(table) == 27
    assert (table["Retry n"] == 28).all()
    assert (table["Direct n"] == 22).all()
    assert int(table["Flag |difference| > 1.96 SE"].sum()) == 5
    for component in ("F", "S", "U", "G", "G_signed", "R"):
        assert {f"{component} retry", f"{component} direct", f"{component} difference"} <= set(table)


def test_formula_hashes_and_summary():
    expected = {
        "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
        "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    }
    actual = {
        name: hashlib.sha256((ROOT / "warrantxai" / name).read_bytes()).hexdigest()
        for name in expected
    }
    assert actual == expected
    summary = json.loads((ROOT / "results/phase8_robustness_summary.json").read_text())
    assert summary["query_hypothesis_confirmed"] == 0
    assert summary["query_hypothesis_rejected"] == 4
    assert summary["retry_flags"] == 5
