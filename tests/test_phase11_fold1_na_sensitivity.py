import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_phase11_directional_counts_and_decisions():
    data = pd.read_csv(ROOT / "PHASE11_FOLD1_NA_SENSITIVITY.csv")
    assert len(data) == 160
    assert {"F", "S", "U", "G_magnitude", "G_signed", "G_claim", "R", "WI"} <= set(data)
    a = data[data.scenario == "A_exclude_rep15_N49"]
    b = data[data.scenario == "B_rescue_rep15_penalizer_0.15_N50"]
    assert (len(a), len(b)) == (80, 80)
    assert (int(a.point_adequate.sum()), int(a.BH.sum()), int(a.BY.sum())) == (19, 13, 11)
    assert (int(b.point_adequate.sum()), int(b.BH.sum()), int(b.BY.sum())) == (20, 13, 11)
    assert set(a.loc[a.BH, "feature"]) == set(b.loc[b.BH, "feature"])
    assert set(b.loc[b.point_adequate, "feature"]) - set(a.loc[a.point_adequate, "feature"]) == {
        "clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II"
    }


def test_phase11_is_diagnostic_and_formula_frozen():
    metadata = json.loads(
        (ROOT / "results/phase11_fold1_na_sensitivity/metadata.json").read_text()
    )
    assert metadata["diagnostic_only"] is True
    assert metadata["official_NA_preserved"] is True
    assert metadata["rescue_record"]["penalizer"] == 0.15
    assert metadata["rescue_record"]["converged"] is True
    expected = metadata["formula_sha256"]
    for filename in ("components.py", "aggregation.py"):
        actual = hashlib.sha256((ROOT / "warrantxai" / filename).read_bytes()).hexdigest()
        assert actual == expected[filename]
