"""Artifact checks for the Phase 6.1 DeepSurv explainer revision."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase61_deepsurv_explainer_revision"
WARNING = (
    "DeepSurv Phase 6.1 conditions on fixed 12+12 queries and deterministic nested "
    "DeepExplainer k=16/128 TCGA backgrounds; N=2 is a design result, not certification."
)


def test_phase61_gradient_budget_curve_and_deep_background_design():
    comparison = pd.read_csv(OUT / "option_comparison.csv")
    gradient = comparison[comparison.option == "A_GradientExplainer"].sort_values("candidate_budget")
    deep = comparison[comparison.option == "B_DeepExplainer"].set_index("design")
    assert gradient[["candidate_budget", "reference_budget"]].to_records(index=False).tolist() == [
        (8, 128), (32, 256), (64, 512), (128, 1024), (256, 2048), (512, 4096)
    ]
    assert gradient.loc[gradient.additivity_pass, "candidate_budget"].min() == 256
    assert not gradient.loc[gradient.candidate_budget <= 128, "additivity_pass"].any()
    assert bool(deep.loc["16/64", "additivity_pass"])
    assert bool(deep.loc["16/128", "additivity_pass"])
    assert deep.loc["16/128", "max_relative_additivity_error_candidate"] < 1e-6
    assert deep.loc["16/128", "max_relative_additivity_error_reference"] < 1e-6


def test_phase61_frozen_f_detects_controlled_corruption_for_both_options():
    curves = pd.read_csv(OUT / "controlled_corruption.csv")
    assert curves.design.nunique() == 2
    for _, curve in curves.groupby("design"):
        curve = curve.sort_values("added_candidate_error")
        assert curve.added_candidate_error.tolist() == [0.0, 0.1, 0.25, 0.5]
        assert np.all(np.diff(curve.selected_mean_F) < 0)
        assert curve.iloc[0].selected_mean_F - curve.iloc[2].selected_mean_F > 0.20
        assert curve.iloc[2].selected_below_tau >= 7


def test_phase61_recommended_n2_pilot_is_complete_and_stops_before_scale_up():
    report = pd.read_csv(OUT / "feature_report.csv")
    assert len(report) == 77 and report.feature.nunique() == 77
    assert report.inference_scope.eq(WARNING).all()
    assert report.point_adequate.sum() == 11
    assert report.bottleneck.value_counts().to_dict() == {"R": 75, "G_claim": 2}
    np.testing.assert_allclose(
        report.WI, report[["F", "S", "U", "G_claim", "R"]].min(axis=1), atol=1e-12
    )
    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["recommendation"] == "DeepExplainer_k16_k128"
    assert metadata["n2"]["full_scale_started"] is False
    assert all(metadata["controlled_corruption_pass"].values())


def test_phase61_hashes_replicas_and_privacy_are_clean():
    metadata = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["formula_module_sha256_unchanged"] is True
    for relative, digest in metadata["formula_module_sha256"].items():
        with (ROOT / relative).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == digest
    for relative, digest in metadata["source_sha256"].items():
        with (ROOT / relative).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == digest
    boot = json.loads((OUT / "bootstrap_summary.json").read_text(encoding="utf-8"))
    assert [row["replicate"] for row in boot] == [0, 1]
    assert all(len(row["sites"]) == 2 for row in boot)

    samples = json.loads((ROOT / "data/external/metabric/samples.json").read_bytes())
    identifiers = {row[key] for row in samples for key in ["sampleId", "patientId"]}
    source = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
    for split in ["train", "test"]:
        keys = pd.read_parquet(source / f"y_{split}.parquet")["sample_key"].astype(str)
        identifiers.update(keys)
        identifiers.update(key[:12] for key in keys)
    pattern = re.compile("|".join(re.escape(value) for value in sorted(identifiers, key=len, reverse=True)))
    for path in [ROOT / "PHASE6_1_EXPLAINER_REVISION.md", *OUT.glob("*")]:
        if path.is_file() and path.suffix != ".npz":
            assert pattern.search(path.read_text(encoding="utf-8")) is None, path.name
