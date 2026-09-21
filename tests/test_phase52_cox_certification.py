"""Checks for the approved CoxPH-only Phase 5.2A full-family run."""

import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from warrantxai.certification import fdr_correction


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase52_cox_certification"
COMPONENTS = ("F", "S", "U", "G_claim", "R")


def test_cox_full_family_iut_and_fdr_outputs():
    result = pd.read_csv(OUT / "feature_certification.csv")
    assert len(result) == 77 and result.feature.nunique() == 77
    np.testing.assert_allclose(result.p_joint, result[[f"p_{x}" for x in COMPONENTS]].max(axis=1))
    rejected_bh, q_bh = fdr_correction(result.p_joint, alpha=0.05, method="bh")
    rejected_by, q_by = fdr_correction(result.p_joint, alpha=0.05, method="by")
    np.testing.assert_allclose(result.q_BH, q_bh)
    np.testing.assert_allclose(result.q_BY, q_by)
    np.testing.assert_array_equal(result.certified_BH_alpha_0_05, rejected_bh)
    np.testing.assert_array_equal(result.certified_BY_alpha_0_05, rejected_by)


def test_cox_certification_uses_all_assessed_features_and_expected_counts():
    source = pd.read_csv(ROOT / "results/phase3_brca/feature_report.csv")
    source = source.loc[np.isfinite(source[list(COMPONENTS)]).all(axis=1)]
    result = pd.read_csv(OUT / "feature_certification.csv")
    assert set(result.feature) == set(source.feature)
    assert result.point_WI_ge_0_8.sum() == 22
    assert result.certified_BH_alpha_0_05.sum() == 8
    assert result.certified_BY_alpha_0_05.sum() == 8
    assert not (result.certified_BH_alpha_0_05 & ~result.point_WI_ge_0_8).any()
    metadata = json.loads((OUT / "metadata.json").read_text())
    assert metadata["family"]["m"] == 77
    assert "no WI preselection" in metadata["family"]["definition"]
    assert metadata["method"]["full_BH_BY_applied"] is True
    assert metadata["formula_module_sha256_unchanged_from_phase5_pilot"] is True
    assert metadata["privacy_check"]["passed"] is True


def test_phase52_cox_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    tcga = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(tcga / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [ROOT / "PHASE5_2_COX_CERTIFICATION_REPORT.md", *list(OUT.glob("*"))]
    for path in paths:
        if path.is_file() and pattern.search(path.read_text(encoding="utf-8")) is not None:
            raise AssertionError(f"Identifier leakage in {path.name}")
