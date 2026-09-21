"""Regression checks for Phase 5.3B CoxPH method consistency."""

import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from warrantxai.certification import fdr_correction


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase53_cox_recursive"
COMPONENTS = ("F", "S", "U", "G_claim", "R")


def test_recursive_cox_full_family_iut_and_fdr_outputs():
    result = pd.read_csv(OUT / "feature_certification.csv")
    assert len(result) == 77 and result.feature.nunique() == 77
    np.testing.assert_allclose(result.p_G_claim, result[["p_G_magnitude_leaf", "p_G_signed_leaf"]].max(axis=1))
    np.testing.assert_allclose(result.p_joint, result[[f"p_{name}" for name in COMPONENTS]].max(axis=1))
    rejected_bh, q_bh = fdr_correction(result.p_joint, alpha=0.05, method="bh")
    rejected_by, q_by = fdr_correction(result.p_joint, alpha=0.05, method="by")
    np.testing.assert_allclose(result.q_BH, q_bh)
    np.testing.assert_allclose(result.q_BY, q_by)
    np.testing.assert_array_equal(result.certified_BH_alpha_0_05, rejected_bh)
    np.testing.assert_array_equal(result.certified_BY_alpha_0_05, rejected_by)


def test_recursive_cox_matches_basic_certification_set_and_isolates_method_difference():
    comparison = pd.read_csv(OUT / "method_comparison.csv")
    assert comparison.certified_BH_alpha_0_05_basic_root.sum() == 8
    assert comparison.certified_BH_alpha_0_05_recursive_studentized.sum() == 8
    np.testing.assert_array_equal(
        comparison.certified_BH_alpha_0_05_basic_root,
        comparison.certified_BH_alpha_0_05_recursive_studentized,
    )
    # The identical bootstrap stream isolates recursive G_claim leaf testing;
    # single-leaf basic S/R and degenerate/high-margin F/U remain unchanged.
    for component in ("F", "S", "U", "R"):
        np.testing.assert_allclose(
            comparison[f"p_{component}_basic_root"],
            comparison[f"p_{component}_recursive_studentized"],
        )
    assert (comparison.p_G_claim_recursive_studentized >= comparison.p_G_claim_basic_root).all()
    metadata = json.loads((OUT / "metadata.json").read_text())
    assert metadata["method"]["name"] == "recursive_studentized"
    assert metadata["method"]["full_BH_BY_applied"] is True
    assert metadata["counts"]["BH_intersection"] == 8
    assert metadata["counts"]["basic_only"] == metadata["counts"]["recursive_only"] == 0
    assert metadata["formula_module_sha256_unchanged_from_phase5_pilot"] is True


def test_phase53_cox_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    tcga = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(tcga / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [ROOT / "PHASE5_3_COX_METHOD_COMPARISON_REPORT.md", *list(OUT.glob("*"))]
    for path in paths:
        if path.is_file() and pattern.search(path.read_text(encoding="utf-8")) is not None:
            raise AssertionError(f"Identifier leakage in {path.name}")
