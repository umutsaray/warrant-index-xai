"""End-to-end artifact checks for Phase 6.2 DeepSurv N=50."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from warrantxai.certification import fdr_correction


ROOT = Path(__file__).resolve().parents[1]
N50 = ROOT / "results/phase62_deepsurv_n50"
GATE = ROOT / "results/phase62_deepsurv_n50_null_gate"
CERT = ROOT / "results/phase62_deepsurv_n50_certification"
COMPARE = ROOT / "results/phase62_three_model_comparison"
COMPONENTS = ("F", "S", "U", "G_claim", "R")
WARNING = (
    "DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested "
    "DeepExplainer k=16/128 TCGA backgrounds."
)


def test_phase62_n50_is_complete_and_preserves_phase61_n2():
    report = pd.read_csv(N50 / "feature_report.csv")
    old = pd.read_csv(ROOT / "results/phase61_deepsurv_explainer_revision/feature_report.csv")
    assert len(report) == 77 and report.feature.nunique() == 77
    assert report.inference_scope.eq(WARNING).all()
    assert report.point_adequate.sum() == 8
    assert report.bottleneck.value_counts().to_dict() == {"R": 77}
    for component in (*COMPONENTS, "WI"):
        np.testing.assert_allclose(report[f"{component}_n2"], old[component], atol=2e-12)
    boot = json.loads((N50 / "bootstrap_summary.json").read_text(encoding="utf-8"))
    assert len(boot) == 50
    assert [row["replicate"] for row in boot] == list(range(50))
    assert all(len(row["sites"]) == 2 for row in boot)
    metadata = json.loads((N50 / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["model"]["main_model_retrained"] is False
    assert metadata["model"]["new_models_fit_this_phase"] == 48
    assert metadata["additivity"]["shap_internal_absolute_assertion_fallbacks"] == 2
    assert metadata["additivity"]["max_new_reference_relative_error"] < 0.05
    assert metadata["runtime"]["observed_total_wall_seconds_across_attempts"] < 30 * 60


def test_phase62_fresh_null_gate_passes_all_holm_controls():
    rates = pd.read_csv(GATE / "type1_rates.csv")
    assert len(rates) == 30
    assert set(rates.null_scenario) == {"single_active", "tie"}
    assert set(rates.nominal_alpha) == {0.01, 0.05, 0.10}
    assert set(rates.component) == set(COMPONENTS)
    assert rates.inference_scope.eq(WARNING).all()
    assert not rates.anti_conservative_gate.any()
    metadata = json.loads((GATE / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["gate_passed"] is True
    assert metadata["arguments"]["outer"] == 2000
    assert metadata["arguments"]["inner"] == 999
    assert metadata["model_or_explainer_calls"] is False


def test_phase62_full_family_certification_has_zero_bh_by_discoveries():
    result = pd.read_csv(CERT / "feature_certification.csv")
    assert len(result) == 77
    assert result.inference_scope.eq(WARNING).all()
    np.testing.assert_allclose(result.p_joint, result[[f"p_{c}" for c in COMPONENTS]].max(axis=1))
    rejected_bh, q_bh = fdr_correction(result.p_joint, alpha=0.05, method="bh")
    rejected_by, q_by = fdr_correction(result.p_joint, alpha=0.05, method="by")
    np.testing.assert_array_equal(result.certified_BH_alpha_0_05, rejected_bh)
    np.testing.assert_array_equal(result.certified_BY_alpha_0_05, rejected_by)
    np.testing.assert_allclose(result.q_BH, q_bh)
    np.testing.assert_allclose(result.q_BY, q_by)
    assert result.point_WI_ge_0_8.sum() == 8
    assert result.certified_BH_alpha_0_05.sum() == 0
    assert result.certified_BY_alpha_0_05.sum() == 0


def test_phase62_three_model_comparison_and_winners_curse_diagnosis():
    summary = pd.read_csv(COMPARE / "three_model_summary.csv").set_index("model")
    assert summary[["point_WI_ge_0_8", "BH", "BY"]].astype(int).to_dict("index") == {
        "CoxPH N=50": {"point_WI_ge_0_8": 22, "BH": 8, "BY": 8},
        "RSF N=50": {"point_WI_ge_0_8": 20, "BH": 10, "BY": 10},
        "DeepSurv N=50": {"point_WI_ge_0_8": 8, "BH": 0, "BY": 0},
    }
    margins = pd.read_csv(COMPARE / "margin_distribution.csv").set_index("model")
    assert margins.loc["DeepSurv N=50", "share_0.05_0.10"] == 0
    distances = pd.read_csv(COMPARE / "margin_distances.csv").set_index("comparison")
    assert distances.loc["DeepSurv_vs_CoxPH", "wasserstein_margin_distance"] < distances.loc[
        "DeepSurv_vs_RSF", "wasserstein_margin_distance"
    ]
    transition = pd.read_csv(COMPARE / "deepsurv_n2_n50_transition.csv")
    assert transition.transition.value_counts().to_dict() == {
        "neither": 59, "lost": 10, "new": 7, "retained": 1
    }


def test_phase62_hashes_and_privacy_are_clean():
    for directory in [N50, GATE, CERT, COMPARE]:
        metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
        assert metadata["formula_module_sha256_unchanged"] is True
        for relative, digest in metadata["source_sha256"].items():
            path = ROOT / relative
            if path.exists():
                with path.open("rb") as handle:
                    assert hashlib.file_digest(handle, "sha256").hexdigest() == digest

    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    source = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(source / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [
        ROOT / "PHASE6_2_STEP1_N50_REPORT.md", ROOT / "PHASE6_2_STEP2_NULL_GATE_REPORT.md",
        ROOT / "PHASE6_2_STEP3_CERTIFICATION_REPORT.md", ROOT / "PHASE6_2_REPORT.md",
        *N50.glob("*"), *GATE.glob("*"), *CERT.glob("*"), *COMPARE.glob("*"),
    ]
    for path in paths:
        if path.is_file() and path.suffix not in {".npz"}:
            assert pattern.search(path.read_text(encoding="utf-8")) is None, path.name
