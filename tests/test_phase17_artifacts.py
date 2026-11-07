"""Integrity and frozen-result checks for the Phase 17 exploratory package."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase17"
HASHES = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    "certification.py": "76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487",
}


def test_v1_reproduction_gate_is_exact() -> None:
    frame = pd.read_csv(OUT / "PHASE17_ADIM0_V1_REPRODUCTION.csv")
    assert len(frame) == 81
    # The producing gate checks all five components and p_joint with
    # atol=1e-12, then checks both decisions exactly before setting this flag.
    assert frame["official_exact_match"].all()


def test_grid_primary_and_certification_counts() -> None:
    grid = pd.read_csv(OUT / "PHASE17_ADIM2_POINT_GRID.csv")
    assert len(grid) == 3 * 45 * 27
    primary_grid = grid.query("k == 1 and rho == 1 and kappa == 0.25")
    assert primary_grid.groupby("model")["point_adequate6"].sum().to_dict() == {
        "CoxPH-v2": 0,
        "DeepSurv": 0,
        "RSF": 0,
    }
    primary = pd.read_csv(OUT / "PHASE17_ADIM3_V2_CERT_PRIMARY.csv")
    assert len(primary) == 81
    assert int(primary["certified_BH"].sum()) == 0
    assert int(primary["certified_BY"].sum()) == 0
    neighbors = pd.read_csv(OUT / "PHASE17_ADIM3_NEIGHBORS.csv")
    assert len(neighbors) == 324
    assert int(neighbors["certified_BH"].sum()) == 0
    assert int(neighbors["certified_BY"].sum()) == 0


def test_synthetic_output_is_complete_20_repeat_gate_result() -> None:
    frame = pd.read_csv(OUT / "PHASE17_ADIM5_SYNTH.csv")
    assert len(frame) == 20 * 30
    assert frame["repeat"].nunique() == 20
    assert frame["feature"].nunique() == 30
    assert frame.groupby("truth_class")["feature"].nunique().to_dict() == {
        "A": 8,
        "B": 6,
        "C": 6,
        "D": 5,
        "E": 5,
    }
    assert int(frame.query("truth_class == 'B'")["v1_BH"].sum()) == 120
    assert int(frame.query("truth_class == 'B'")["v2_BH"].sum()) == 0


def test_n_calibration_and_frozen_formula_hashes() -> None:
    calibration = pd.read_csv(OUT / "PHASE17_ADIM6_N_CALIBRATION.csv")
    assert len(calibration) == 18
    assert not calibration["anti_conservative_holm"].any()
    for filename, expected in HASHES.items():
        actual = hashlib.sha256((ROOT / "warrantxai" / filename).read_bytes()).hexdigest()
        assert actual == expected


def test_report_keeps_required_frozen_analysis_statement() -> None:
    report = (ROOT / "PHASE17_V2.md").read_text(encoding="utf-8")
    required = (
        "Bu analiz v1'in dondurulmuş ayarlarını, Tablo 3–9'u, resmi BH/BY "
        "kararlarını, saklı atıf tensörlerini, öznitelik panelini ve CoxPH-v2 "
        "retry protokolünü DEĞİŞTİRMEMEKTEDİR; v2 sonuçları ayrı, post hoc "
        "etiketli, keşifsel bir analizdir."
    )
    assert report.rstrip().endswith(required)
