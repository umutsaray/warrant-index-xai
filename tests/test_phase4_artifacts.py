"""Optional checks for the intentionally bounded Phase 4 RSF pilot."""

from pathlib import Path
import hashlib
import json
import re

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase4_rsf_pilot"
pytestmark = pytest.mark.skipif(
    not (OUT / "metadata.json").exists(), reason="Run the optional Phase 4 pilot first"
)


def test_phase4_pilot_scope_and_formula_identity():
    report = pd.read_csv(OUT / "feature_report.csv")
    metadata = json.loads((OUT / "metadata.json").read_text())
    assert len(report) == 77
    assert report.feature.str.startswith("rna_").sum() == 72
    assert report.feature.str.startswith("clin_").sum() == 5
    components = report[["F", "S", "U", "G_claim", "R"]]
    assert np.isfinite(components).all().all()
    assert ((components >= 0) & (components <= 1)).all().all()
    np.testing.assert_allclose(report.WI, components.min(axis=1))
    assert metadata["arguments"] == {
        "replicas": 2,
        "anchors_per_site": 12,
        "candidate_permutations": 8,
        "reference_permutations": 128,
        "seed": 20260917,
    }
    assert metadata["explainer"]["tree_shap_supported"] is False
    assert metadata["explainer"]["candidate_reference_budget_ratio"] == 16
    assert metadata["explainer"]["max_additivity_error"] < 1e-10
    assert all(not row["warning_categories"] for row in metadata["training"])
    model_path = Path(metadata["model"]["source"])
    assert metadata["model"]["sha256"] == hashlib.file_digest(model_path.open("rb"), "sha256").hexdigest()


def test_phase4_outputs_contain_no_patient_or_sample_identifiers():
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
        *OUT.rglob("*.csv"),
        *OUT.rglob("*.md"),
        *OUT.rglob("*.json"),
        *OUT.rglob("*.txt"),
    ]
    for path in paths:
        if pattern.search(path.read_text(encoding="utf-8")) is not None:
            pytest.fail(f"Identifier leakage in {path.name}", pytrace=False)


def test_phase41_sensitivity_decision_and_scope():
    report = pd.read_csv(OUT / "sensitivity_report.csv")
    metadata = json.loads((OUT / "sensitivity_metadata.json").read_text())
    assert len(report) == 8
    assert report.selection_reason.value_counts().to_dict() == {
        "low_WI": 3,
        "high_WI": 3,
        "seeded_random": 2,
    }
    assert metadata["arguments"]["replicas"] == 2
    assert metadata["arguments"]["query_rows_per_site"] == 12
    assert metadata["arguments"]["background_patients"] == 16
    assert metadata["background_small"] is False
    assert metadata["reference_budget_small"] is True
    assert report.F_budget_delta_256_minus_128.abs().max() < 0.01
    assert report.G_claim_budget_delta_256_minus_128.eq(0).all()
    assert report.G_claim_background_delta.abs().max() > 0.01
    assert metadata["recommendation"] == "Use sampled background (k=16) and 8/128 for N=50"
