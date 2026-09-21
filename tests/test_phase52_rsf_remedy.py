import json
from pathlib import Path
import re

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase52_rsf_remedy"


def test_remedy_candidates_and_gate_decision():
    summary = pd.read_csv(OUT / "candidate_gate_summary.csv").set_index("candidate")
    assert bool(summary.loc["recursive_studentized", "gate_passed"])
    assert not bool(summary.loc["recursive_basic", "gate_passed"])
    assert not bool(summary.loc["m_out_of_n_root", "gate_passed"])
    rates = pd.read_csv(OUT / "candidate_type1_rates.csv")
    assert len(rates) == 90
    passing = rates[rates.candidate == "recursive_studentized"]
    assert len(passing) == 30 and not passing.anti_conservative_gate.any()
    assert rates.inference_scope.str.contains(r"12\+12", regex=True).all()
    metadata = json.loads((OUT / "metadata.json").read_text())
    assert metadata["recommended_candidate"] == "recursive_studentized"
    assert metadata["rsf_adim2_run"] is False
    assert metadata["formula_module_sha256_unchanged_from_phase4_2"] is True


def test_actual_five_feature_remedy_pilot():
    pilot = pd.read_csv(OUT / "actual_5feature_pilot.csv")
    assert len(pilot) == 5
    assert pilot.inference_scope.str.contains(r"12\+12", regex=True).all()
    high = pilot[pilot.selection_reason == "high_margin"].iloc[0]
    assert high.p_joint_recursive_studentized == 0.0002
    near = pilot[pilot.selection_reason == "near_threshold"].iloc[0]
    assert 0.4 <= near.p_R_recursive_studentized <= 0.8
    low = pilot[pilot.selection_reason == "phase4.1_low_WI"].iloc[0]
    assert low.G_claim == 0 and low.p_G_claim_recursive_studentized == 1


def test_phase52_rsf_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    tcga = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(tcga / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [ROOT / "PHASE5_2_RSF_REMEDY_REPORT.md", *list(OUT.glob("*"))]
    for path in paths:
        if path.is_file() and pattern.search(path.read_text(encoding="utf-8")) is not None:
            raise AssertionError(f"Identifier leakage in {path.name}")
