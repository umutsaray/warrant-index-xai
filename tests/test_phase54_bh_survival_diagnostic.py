"""Regression checks for the read-only Phase 5.4 diagnostic."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase54_bh_survival_diagnostic"


def test_phase54_reconstructs_existing_certification_counts_and_margins():
    data = pd.read_csv(OUT / "point_pass_feature_diagnostics.csv")
    counts = data.groupby("model").agg(
        point_pass=("feature", "size"),
        certified=("BH_status", lambda values: (values == "BH certified").sum()),
        point_only=("BH_status", lambda values: (values == "point only").sum()),
    )
    assert counts.loc["CoxPH N=50"].to_dict() == {"point_pass": 22, "certified": 8, "point_only": 14}
    assert counts.loc["RSF N=20"].to_dict() == {"point_pass": 23, "certified": 15, "point_only": 8}
    np.testing.assert_allclose(data.WI_margin, data.WI - 0.8)
    np.testing.assert_allclose(data.bottleneck_margin, data.bottleneck_value - 0.8)
    np.testing.assert_allclose(data.WI_margin, data.bottleneck_margin)


def test_phase54_key_diagnostic_findings_are_data_backed():
    data = pd.read_csv(OUT / "point_pass_feature_diagnostics.csv")
    missed = data.loc[data.BH_status == "point only"]
    assert set(missed.loc[missed.model == "CoxPH N=50", "bottleneck"]) == {"R"}
    assert set(missed.loc[missed.model == "RSF N=20", "bottleneck"]) == {"R"}
    assert (data.loc[(data.model == "CoxPH N=50") & (data.WI_margin > 0.03)].shape[0]) == 9
    assert (data.loc[(data.model == "RSF N=20") & (data.WI_margin > 0.03)].shape[0]) == 16
    low = pd.read_csv(OUT / "margin_0_003_comparison.csv").set_index("model")
    assert low.loc["CoxPH N=50", "mean_p_bottleneck"] < low.loc["RSF N=20", "mean_p_bottleneck"]
    assert not (data.loc[data.WI_margin < 0.03, "BH_status"] == "BH certified").any()


def test_phase54_is_read_only_and_formula_hashes_are_frozen():
    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["new_bootstrap_or_model_calls"] is False
    assert metadata["formula_module_sha256_unchanged"] is True
    for relative, digest in metadata["source_files"].items():
        path = ROOT / Path(relative)
        assert hashlib.file_digest(path.open("rb"), "sha256").hexdigest() == digest
    assert metadata["privacy_check"]["passed"] is True


def test_phase54_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    source = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(source / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [ROOT / "PHASE5_4_BH_SURVIVAL_DIAGNOSTIC.md", *OUT.glob("*")]
    for path in paths:
        if path.is_file() and path.suffix != ".png":
            assert pattern.search(path.read_text(encoding="utf-8")) is None, path.name
