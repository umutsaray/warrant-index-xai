"""Checks for the approved RSF Phase 5.3A full-family run."""

import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from warrantxai.certification import fdr_correction


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase53_rsf_certification"
COMPONENTS = ("F", "S", "U", "G_claim", "R")


def test_rsf_recursive_studentized_full_family_iut_and_fdr():
    result = pd.read_csv(OUT / "feature_certification.csv")
    assert len(result) == 77 and result.feature.nunique() == 77
    np.testing.assert_allclose(result.p_joint, result[[f"p_{x}" for x in COMPONENTS]].max(axis=1))
    rejected_bh, q_bh = fdr_correction(result.p_joint, alpha=0.05, method="bh")
    rejected_by, q_by = fdr_correction(result.p_joint, alpha=0.05, method="by")
    np.testing.assert_allclose(result.q_BH, q_bh)
    np.testing.assert_allclose(result.q_BY, q_by)
    np.testing.assert_array_equal(result.certified_BH_alpha_0_05, rejected_bh)
    np.testing.assert_array_equal(result.certified_BY_alpha_0_05, rejected_by)
    assert result.inference_scope.str.contains(r"12\+12", regex=True).all()
    assert result.inference_scope.str.contains("per-permutation-path uncertainty is excluded", regex=False).all()


def test_rsf_certification_expected_counts_and_metadata():
    source = pd.read_csv(ROOT / "results/phase4_rsf/feature_report.csv")
    result = pd.read_csv(OUT / "feature_certification.csv")
    assert set(result.feature) == set(source.feature)
    assert result.point_WI_ge_0_8.sum() == 23
    assert result.certified_BH_alpha_0_05.sum() == 15
    assert result.certified_BY_alpha_0_05.sum() == 12
    assert not (result.certified_BH_alpha_0_05 & ~result.point_WI_ge_0_8).any()
    summary = pd.read_csv(OUT / "summary.csv")
    assert summary.inference_scope.str.contains(r"12\+12", regex=True).all()
    metadata = json.loads((OUT / "metadata.json").read_text())
    assert metadata["family"]["m"] == 77
    assert metadata["method"]["name"] == "recursive_studentized"
    assert metadata["method"]["full_BH_BY_applied"] is True
    assert metadata["formula_module_sha256_unchanged_from_phase5_2"] is True
    assert metadata["privacy_check"]["passed"] is True
    assert metadata["model_or_explainer_calls"] is False
    assert "12+12" in metadata["rsf_conditionality_warning"]


def test_phase53_rsf_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    tcga = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(tcga / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [
        ROOT / "PHASE5_3_RSF_CERTIFICATION_REPORT.md",
        ROOT / "PHASE5_3_REPORT.md",
        *list(OUT.glob("*")),
    ]
    for path in paths:
        if path.is_file() and pattern.search(path.read_text(encoding="utf-8")) is not None:
            raise AssertionError(f"Identifier leakage in {path.name}")
