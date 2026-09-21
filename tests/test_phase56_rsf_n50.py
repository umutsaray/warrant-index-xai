"""End-to-end artifact checks for Phase 5.6 RSF N=50."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from warrantxai.certification import fdr_correction


ROOT = Path(__file__).resolve().parents[1]
N50 = ROOT / "results/phase56_rsf_n50"
GATE = ROOT / "results/phase56_rsf_n50_null_gate"
CERT = ROOT / "results/phase56_rsf_n50_certification"
COMPARE = ROOT / "results/phase56_comparison"
COMPONENTS = ("F", "S", "U", "G_claim", "R")


def test_phase56_n50_full_run_preserves_n20_and_has_all_replicas():
    report = pd.read_csv(N50 / "feature_report.csv")
    old = pd.read_csv(ROOT / "results/phase4_rsf/feature_report.csv")
    assert len(report) == 77 and report.feature.nunique() == 77
    assert report.inference_scope.str.contains(r"12\+12", regex=True).all()
    assert report.point_adequate.sum() == 20
    assert report.bottleneck.value_counts().to_dict() == {"R": 72, "G_claim": 5}
    for component in (*COMPONENTS, "WI"):
        np.testing.assert_allclose(report[f"{component}_n20"], old[component], atol=1e-12)
    summary = json.loads((N50 / "bootstrap_summary.json").read_text(encoding="utf-8"))
    assert len(summary) == 50
    assert [row["replicate"] for row in summary] == list(range(50))
    assert all(not row["warning_categories"] and row["fit_attempts"] for row in summary)
    assert all(len(row["site_records"]) == 2 for row in summary)
    assert all("12+12" in row["inference_scope"] for row in summary)
    metadata = json.loads((N50 / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["model"]["main_model_retrained"] is False
    assert metadata["model"]["first_20_reused_from_phase4_2"] is True
    assert metadata["model"]["new_models_fit_this_phase"] == 30
    assert metadata["runtime"]["new_replicates_20_to_49_seconds"] < 4 * 3600


def test_phase56_n50_null_gate_passes_all_holm_controls():
    rates = pd.read_csv(GATE / "type1_rates.csv")
    assert len(rates) == 30
    assert set(rates.null_scenario) == {"single_active", "tie"}
    assert set(rates.nominal_alpha) == {0.01, 0.05, 0.10}
    assert set(rates.component) == set(COMPONENTS)
    assert not rates.anti_conservative_gate.any()
    assert rates.inference_scope.str.contains(r"12\+12", regex=True).all()
    metadata = json.loads((GATE / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["gate_passed"] is True
    assert metadata["arguments"]["outer"] == 2000
    assert metadata["arguments"]["inner"] == 999
    assert metadata["model_or_explainer_calls"] is False


def test_phase56_n50_certification_iut_fdr_and_counts():
    result = pd.read_csv(CERT / "feature_certification.csv")
    assert len(result) == 77
    np.testing.assert_allclose(result.p_joint, result[[f"p_{c}" for c in COMPONENTS]].max(axis=1))
    rejected_bh, q_bh = fdr_correction(result.p_joint, alpha=0.05, method="bh")
    rejected_by, q_by = fdr_correction(result.p_joint, alpha=0.05, method="by")
    np.testing.assert_array_equal(result.certified_BH_alpha_0_05, rejected_bh)
    np.testing.assert_array_equal(result.certified_BY_alpha_0_05, rejected_by)
    np.testing.assert_allclose(result.q_BH, q_bh)
    np.testing.assert_allclose(result.q_BY, q_by)
    assert result.point_WI_ge_0_8.sum() == 20
    assert result.certified_BH_alpha_0_05.sum() == 10
    assert result.certified_BY_alpha_0_05.sum() == 10
    assert result.inference_scope.str.contains(r"12\+12", regex=True).all()


def test_phase56_maturation_verdict_and_frozen_hashes():
    transition = pd.read_csv(COMPARE / "feature_transition.csv")
    assert (transition.point_transition == "retained").sum() == 16
    assert (transition.point_transition == "lost").sum() == 7
    assert (transition.point_transition == "new").sum() == 4
    selected = transition.loc[transition.point_n20.astype(bool)]
    assert (selected.delta_R < 0).sum() == 19
    r_summary = pd.read_csv(COMPARE / "R_delta_summary.csv").set_index("scope")
    assert abs(r_summary.loc["all 77", "mean_delta_R"]) < 0.001
    assert r_summary.loc["all 77", "wilcoxon_p"] > 0.05
    metadata = json.loads((COMPARE / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["verdict"] == "partially_supported"
    assert metadata["formula_module_sha256_unchanged"] is True
    for relative, digest in metadata["source_sha256"].items():
        assert hashlib.file_digest((ROOT / relative).open("rb"), "sha256").hexdigest() == digest


def test_phase56_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    source = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(source / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [
        ROOT / "PHASE5_6_STEP1_RSF_N50_REPORT.md", ROOT / "PHASE5_6_STEP2_NULL_GATE_REPORT.md",
        ROOT / "PHASE5_6_STEP3_CERTIFICATION_REPORT.md", ROOT / "PHASE5_6_REPORT.md",
        *N50.glob("*"), *GATE.glob("*"), *CERT.glob("*"), *COMPARE.glob("*"),
    ]
    for path in paths:
        if path.is_file() and path.suffix not in {".png", ".npz"}:
            assert pattern.search(path.read_text(encoding="utf-8")) is None, path.name
