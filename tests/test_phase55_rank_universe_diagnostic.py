"""Regression checks for the diagnostic-only Phase 5.5 calculation."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase55_rank_universe_diagnostic"


def test_phase55_reproduces_frozen_r_and_changes_only_diagnostic_r():
    data = pd.read_csv(OUT / "all_77_feature_comparison.csv")
    frozen = pd.read_csv(ROOT / "results/phase53_cox_recursive/feature_certification.csv")
    merged = data.merge(frozen[["feature", "R", "WI"]], on="feature", validate="one_to_one")
    np.testing.assert_allclose(merged.R_614, merged.R, atol=0, rtol=0)
    np.testing.assert_allclose(merged.WI_614, merged.WI, atol=0, rtol=0)
    np.testing.assert_allclose(data.delta_R_77_minus_614, data.R_77 - data.R_614)
    np.testing.assert_allclose(data.WI_77, data[["F", "S", "U", "G_claim", "R_77"]].min(axis=1))
    assert not (OUT / "feature_certification.csv").exists()


def test_phase55_hypothesis_verdict_is_supported_by_distances():
    data = pd.read_csv(OUT / "all_77_feature_comparison.csv")
    same_22 = data.loc[data.cox_original_point_pass.astype(bool)]
    assert len(same_22) == 22
    assert np.isclose(same_22.delta_R_77_minus_614.mean(), 0.0028304648349971443)
    assert (same_22.WI_margin_614 > 0.03).sum() == 9
    assert (same_22.WI_margin_77 > 0.03).sum() == 9
    assert data.cox_diagnostic_point_pass.sum() == 23
    changed = pd.read_csv(OUT / "diagnostic_status_changes.csv")
    assert changed.feature.tolist() == ["rna_ENSG00000139874.6"]
    distances = pd.read_csv(OUT / "distance_to_rsf.csv")
    assert not distances.closer_after_restriction.astype(bool).any()
    assert (distances.R77_distance_to_RSF > distances.original_distance_to_RSF).all()


def test_phase55_metadata_proves_read_only_scope_and_frozen_hashes():
    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["diagnostic_only"] is True
    assert metadata["official_certification_outputs_modified"] is False
    assert metadata["new_bootstrap_resampling"] is False
    assert metadata["model_fit_prediction_or_explainer_calls"] is False
    assert metadata["validation"]["recomputed_R614_max_abs_error"] == 0
    assert metadata["formula_module_sha256_unchanged"] is True
    assert metadata["verdict"] == "hypothesis_not_supported"
    for relative, digest in metadata["source_sha256"].items():
        path = ROOT / Path(relative)
        assert hashlib.file_digest(path.open("rb"), "sha256").hexdigest() == digest


def test_phase55_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    source = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(source / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [ROOT / "PHASE5_5_RANK_UNIVERSE_DIAGNOSTIC.md", *OUT.glob("*")]
    for path in paths:
        if path.is_file() and path.suffix != ".png":
            assert pattern.search(path.read_text(encoding="utf-8")) is None, path.name
