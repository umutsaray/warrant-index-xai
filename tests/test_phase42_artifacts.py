"""Artifact and privacy checks for the approved Phase 4.2 N=20 RSF run."""

from pathlib import Path
import hashlib
import json
import re

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase4_rsf"
pytestmark = pytest.mark.skipif(
    not (OUT / "metadata.json").exists(), reason="Run the optional Phase 4.2 experiment first"
)


def test_phase42_n20_scope_scores_and_provenance():
    report = pd.read_csv(OUT / "feature_report.csv")
    metadata = json.loads((OUT / "metadata.json").read_text())
    summary = json.loads((OUT / "bootstrap_summary.json").read_text())
    assert len(report) == 77
    assert report.feature.str.startswith("rna_").sum() == 72
    assert report.feature.str.startswith("clin_").sum() == 5
    components = report[["F", "S", "U", "G_claim", "R"]]
    assert np.isfinite(components).all().all()
    assert ((components >= 0) & (components <= 1)).all().all()
    np.testing.assert_allclose(report.WI, components.min(axis=1))
    assert metadata["arguments"] == {
        "replicas": 20,
        "anchors_per_site": 12,
        "background_patients": 16,
        "candidate_permutations": 8,
        "reference_permutations": 128,
        "seed": 20260917,
    }
    assert metadata["model"]["main_model_retrained"] is False
    assert metadata["formula_module_sha256_unchanged_from_phase21"] is True
    assert metadata["results"]["point_adequate_n20_sampled_background"] == 23
    assert [row["replicate"] for row in summary] == list(range(20))
    assert all(not row["warning_categories"] for row in summary)
    assert all(len(row["site_records"]) == 2 for row in summary)
    source = ROOT / "experiments/07_rsf_pilot.py"
    assert metadata["source_sha256"]["experiments/07_rsf_pilot.py"] == hashlib.sha256(
        source.read_bytes()
    ).hexdigest()
    with np.load(OUT / "attributions_n20.npz", allow_pickle=False) as values:
        assert values["candidate"].shape == values["reference"].shape == (2, 20, 12, 77)


def test_phase42_nested_n2_comparison():
    metadata = json.loads((OUT / "metadata.json").read_text())
    result = metadata["results"]
    assert result["point_adequate_n2_fixed_background_pilot"] == 43
    assert result["point_adequate_n2_sampled_background"] == 35
    assert result["n2_sampled_passing_still_passing_n20"] == 22
    assert result["n2_sampled_passing_now_failing_n20"] == 13
    assert result["n20_newly_passing_vs_n2_sampled"] == 1
    comparison = pd.read_csv(OUT / "n2_n20_comparison.csv")
    assert len(comparison) == 77


def test_phase42_outputs_contain_no_patient_or_sample_identifiers():
    raw = ROOT / "data/external/metabric"
    samples = json.loads((raw / "samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    pytest.importorskip("pyarrow")
    tcga = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(tcga / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [
        ROOT / "MANIFEST.md",
        ROOT / "PHASE4_DESIGN_NOTE.md",
        ROOT / "PHASE4_FULL_RUN_REPORT.md",
        *OUT.rglob("*.csv"),
        *OUT.rglob("*.md"),
        *OUT.rglob("*.json"),
        *OUT.rglob("*.txt"),
    ]
    for path in paths:
        if pattern.search(path.read_text(encoding="utf-8")) is not None:
            pytest.fail(f"Identifier leakage in {path.name}", pytrace=False)
