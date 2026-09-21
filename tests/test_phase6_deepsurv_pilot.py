"""Artifact and honesty checks for the Phase 6 DeepSurv N=2 design pilot."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase6_deepsurv_pilot"
WARNING = (
    "DeepSurv pilot conditions on fixed 12+12 queries and stochastic "
    "GradientExplainer 8/128 attribution draws; N=2 is a design pilot, "
    "not a certification result."
)


def test_phase6_pilot_assesses_complete_77_feature_panel():
    report = pd.read_csv(OUT / "feature_report.csv")
    assert len(report) == 77 and report.feature.nunique() == 77
    assert report.inference_scope.eq(WARNING).all()
    assert report.point_adequate.sum() == 15
    assert report.bottleneck.value_counts().to_dict() == {"R": 75, "G_claim": 2}
    assert report.block.value_counts().to_dict() == {"rna": 72, "clin": 5}
    np.testing.assert_allclose(
        report.WI,
        report[["F", "S", "U", "G_claim", "R"]].min(axis=1),
        atol=1e-12,
    )


def test_phase6_sampled_background_and_reference_health_checks_are_small():
    sensitivity = pd.read_csv(OUT / "sensitivity_report.csv")
    assert len(sensitivity) == 8
    assert sensitivity.inference_scope.eq(WARNING).all()
    assert sensitivity.delta_F_k16_minus_k8.abs().max() < 0.001
    assert sensitivity.delta_G_k16_minus_k8.abs().max() < 0.005
    assert sensitivity.delta_F_ref256_minus_128.abs().max() < 0.001


def test_phase6_explainers_are_really_supported_but_candidate_gate_stops_scale_up():
    probe = json.loads((OUT / "explainer_probe.json").read_text(encoding="utf-8"))
    assert probe["DeepExplainer"]["supported"] is True
    assert probe["GradientExplainer"]["supported"] is True
    assert probe["DeepExplainer"]["shape"] == [3, 616]
    assert probe["GradientExplainer"]["shape"] == [3, 616]
    assert probe["DeepExplainer"]["max_abs_additivity_error"] < 1e-5

    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["runtime"]["max_candidate_additivity_error"] > 0.4
    assert metadata["runtime"]["max_reference_additivity_error"] < 0.1
    assert metadata["adim3"]["design_gate"] == "revision_required_before_N20"
    assert metadata["adim3"]["full_scale_started"] is False


def test_phase6_two_bootstraps_are_complete_and_formula_sources_are_frozen():
    boot = json.loads((OUT / "bootstrap_summary.json").read_text(encoding="utf-8"))
    assert [row["replicate"] for row in boot] == [0, 1]
    assert all(not row["warning_categories"] for row in boot)
    assert all(len(row["sites"]) == 2 for row in boot)
    assert all(row["inference_scope"] == WARNING for row in boot)

    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["formula_module_sha256_unchanged"] is True
    for relative, digest in metadata["formula_module_sha256"].items():
        with (ROOT / relative).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == digest
    for relative, digest in metadata["source_sha256"].items():
        with (ROOT / relative).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == digest


def test_phase6_outputs_contain_no_patient_or_sample_identifiers():
    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    source = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(source / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    paths = [ROOT / "PHASE6_DESIGN_NOTE.md", *OUT.glob("*")]
    for path in paths:
        if path.is_file() and path.suffix not in {".npz"}:
            assert pattern.search(path.read_text(encoding="utf-8")) is None, path.name
