"""Optional real-data artifact checks; never display patient identifiers."""

from pathlib import Path
import json
import re

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase3_brca"
pytestmark = pytest.mark.skipif(
    not (OUT / "metadata.json").exists(), reason="Run the optional Phase 3 experiment first"
)


def test_real_feature_scope_and_frozen_common_q_identity():
    report = pd.read_csv(OUT / "feature_report.csv")
    assert len(report) == 616
    assert report.feature.is_unique
    missing = report.feature.str.startswith(("mirna_", "meth_"))
    assert missing.sum() == 400
    assert report.loc[missing, ["G_magnitude", "G_signed", "G_claim", "WI"]].isna().all().all()
    retained = missing & ~report.status.str.contains("excluded from Cox")
    assert report.loc[retained, ["F", "S", "U", "R"]].notna().all().all()
    assessed = report.WI.notna()
    assert assessed.sum() == 77
    np.testing.assert_allclose(report.loc[assessed, "G_claim_range_overlap_diagnostic"], 1)
    assert (report.loc[assessed, ["G_magnitude", "G_signed", "G_claim"]].to_numpy() < 1).any()
    np.testing.assert_allclose(
        report.loc[assessed, "WI"], report.loc[assessed, ["F", "S", "U", "G_claim", "R"]].min(axis=1)
    )
    assert report.loc[report.status.str.contains("excluded from Cox"), ["F", "S", "U", "R"]].isna().all().all()


def test_real_patient_separation_and_oracle_verification():
    metadata = json.loads((OUT / "metadata.json").read_text())
    stats = metadata["stats"]
    assert stats["patient_overlap"] == 23
    assert stats["audit_unique_patients"] < stats["original_test"]
    assert stats["oracle_max_abs_error"] < 1e-8
    assert stats["external"]["survival_used_for_fitting_or_tuning"] is False
    assert metadata["arguments"]["generalizability_mode"] == "claim"
    assert metadata["arguments"]["phase"] == "3.2"
    assert stats["formula_module_sha256_unchanged_from_phase21"] is True
    assert metadata["source_sha256"] and metadata["input_sha256"]


def test_phase31_uses_different_real_cohort_query_sets():
    diagnostics = pd.read_csv(OUT / "external_support_diagnostics.csv")
    assessed = diagnostics.status.eq("assessed (different real cohort query distributions)")
    rows = diagnostics.loc[assessed]
    assert len(rows) == 77
    assert rows.query_source_site0.eq("TCGA clean disjoint audit patients").all()
    assert rows.query_source_site1.eq("METABRIC actual patient vectors").all()
    assert rows.site0_query_n.eq(141).all()
    assert (rows.site1_query_n > rows.site0_query_n).all()
    assert (~rows.query_value_distributions_identical.astype(bool)).all()
    assert (rows.site0_query_sha256 != rows.site1_query_sha256).all()
    comparison = pd.read_csv(OUT / "phase31_previous_pass_comparison.csv")
    assert len(comparison) == 27


def test_phase32_uses_50_replicates_and_preserves_n6_results():
    metadata = json.loads((OUT / "metadata.json").read_text())
    assert metadata["arguments"]["replicates"] == 50
    report = pd.read_csv(OUT / "feature_report.csv")
    assert {"R_n6", "WI_n6", "point_adequate_n6", "R_n50"} <= set(report)
    uncertainty = pd.read_csv(OUT / "phase32_R_uncertainty.csv")
    comparison = pd.read_csv(OUT / "phase32_previous_pass_comparison.csv")
    assert len(uncertainty) == len(comparison) == 27
    lookup = report.set_index("feature")
    np.testing.assert_allclose(uncertainty.R_n6, lookup.loc[uncertainty.feature, "R_n6"])
    np.testing.assert_allclose(uncertainty.R_n50, lookup.loc[uncertainty.feature, "R"])
    assert uncertainty[["IQR_width_n6", "IQR_width_n50"]].ge(0).all().all()
    training = json.loads((OUT / "bootstrap_summary.json").read_text())
    assert [row["replicate"] for row in training] == list(range(50))
    assert all(not row["warning_categories"] for row in training)
    assert all(row["fit_attempts"] for row in training)


def test_reports_contain_no_patient_or_sample_identifiers():
    raw = ROOT / "data/external/metabric"
    samples = json.loads((raw / "samples.json").read_bytes())
    identifiers = {s[k] for s in samples for k in ["sampleId", "patientId"]}
    # Include actual TCGA keys and patient-prefix keys, only in memory.
    pytest.importorskip("pyarrow")
    tcga = (
        ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    )
    for split in ["train", "test"]:
        keys = pd.read_parquet(tcga / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(k[:12] for k in keys)
    pattern = re.compile("|".join(re.escape(v) for v in sorted(identifiers, key=len, reverse=True)))
    paths = [
        ROOT / "MANIFEST.md",
        ROOT / "PHASE3_REPORT_TR.md",
        *OUT.rglob("*.csv"),
        *OUT.rglob("*.md"),
        *OUT.rglob("*.json"),
        *OUT.rglob("*.txt"),
    ]
    for path in paths:
        # Report filenames only on failure, never matched identity text.
        if pattern.search(path.read_text(encoding="utf-8")) is not None:
            pytest.fail(f"Identifier leakage in {path.name}", pytrace=False)
