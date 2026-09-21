import importlib.util
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase51_null_calibration"


def calibration_module():
    spec = importlib.util.spec_from_file_location("phase51", ROOT / "experiments/09_null_calibration.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_known_null_structures_have_population_target_point_eight():
    # Analytic identities used by the simulation: E|U-V|=1/3 and
    # P(Bernoulli(.5) differs)=.5.
    assert 1 - 0.6 / 3 == 0.8
    assert 1 - 0.4 * 0.5 == 0.8
    module = calibration_module()
    rng = np.random.default_rng(31)
    uniform = rng.uniform(size=(50, 2, 5))
    indices = rng.integers(0, 50, size=(20, 50))
    point, sampled = module.frozen_shape_statistics(uniform, np.full(5, 0.05), indices, "single_active")
    assert point.shape == (5,) and sampled.shape == (20, 5)
    assert np.isfinite(sampled).all()


def test_calibration_gate_fails_and_full_fdr_stays_stopped():
    rates = pd.read_csv(OUT / "type1_rates.csv")
    assert len(rates) == 60
    assert set(rates.model) == {"CoxPH_N50", "RSF_N20"}
    assert set(rates.null_scenario) == {"single_active", "tie"}
    assert set(rates.component) == {"F", "S", "U", "G_claim", "R"}
    assert set(rates.nominal_alpha) == {0.01, 0.05, 0.1}
    flagged = rates[rates.anti_conservative_gate]
    assert set(flagged.component) == {"F", "U", "G_claim"}
    assert flagged.model.eq("RSF_N20").all() and flagged.null_scenario.eq("single_active").all()
    assert rates.loc[rates.model == "RSF_N20", "inference_scope"].str.contains(r"12\+12", regex=True).all()
    metadata = json.loads((OUT / "metadata.json").read_text())
    assert metadata["calibration_gate_passed"] is False
    assert metadata["full_BH_BY_applied"] is False
    assert metadata["formula_module_sha256_unchanged_from_phase4_2"] is True


def test_phase51_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    tcga = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(tcga / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [ROOT / "PHASE5_1_NULL_CALIBRATION_REPORT.md", *list(OUT.glob("*"))]
    for path in paths:
        if path.is_file() and pattern.search(path.read_text(encoding="utf-8")) is not None:
            raise AssertionError(f"Identifier leakage in {path.name}")
