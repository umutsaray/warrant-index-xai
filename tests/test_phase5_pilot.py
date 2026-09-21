import importlib.util
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def pilot_module():
    spec = importlib.util.spec_from_file_location("phase5_pilot", ROOT / "experiments/08_certification_pilot.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_centered_bootstrap_pvalue_hand_calculation():
    module = pilot_module()
    point = np.array([0.75, 0.5])
    draws = np.array([[1.0, 0.4], [1.0, 0.5], [0.75, 0.6], [0.5, 0.55]])
    # Margins are .25 and 0.  Counts of delta >= margin are 2 and 3;
    # the finite Monte Carlo correction adds one to numerator/denominator.
    np.testing.assert_allclose(module.centered_bootstrap_pvalue(point, draws, 0.5), [3 / 5, 4 / 5])


def test_phase5_pilot_is_calibrated_and_stops_before_fdr():
    result = pd.read_csv(ROOT / "results/phase5_certification_pilot/pilot_pvalues.csv")
    assert len(result) == 10
    assert set(result.groupby("model").size()) == {5}
    high = result[result.selection_reason == "high_margin"]
    assert (high.p_F == 1 / 5000).all()
    near = result[result.selection_reason == "near_threshold"]
    assert near.p_R.between(0.4, 0.75).all()
    low_rsf = result[(result.model == "RSF_N20") & (result.selection_reason == "phase4.1_low_WI")].iloc[0]
    assert low_rsf.G_claim == 0 and low_rsf.p_G_claim == 1
    assert not any(column.startswith(("q_", "BH_", "BY_")) for column in result.columns)
    metadata = json.loads((ROOT / "results/phase5_certification_pilot/metadata.json").read_text())
    assert metadata["method"]["full_BH_BY_applied"] is False
    assert metadata["formula_module_sha256_unchanged_from_phase4_2"] is True
    assert metadata["privacy_check"]["passed"] is True


def test_phase5_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    tcga = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(tcga / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [
        ROOT / "PHASE5_DESIGN_NOTE.md",
        *list((ROOT / "results/phase5_certification_pilot").glob("*")),
    ]
    for path in paths:
        if path.is_file() and pattern.search(path.read_text(encoding="utf-8")) is not None:
            raise AssertionError(f"Identifier leakage in {path.name}")
