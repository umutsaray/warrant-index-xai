"""Integrity checks for the post hoc Phase 18 tolerance-profile package."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase18"
HASHES = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    "certification.py": "76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487",
}


def test_phase17_reproduction_gate_passed() -> None:
    gate = pd.read_csv(OUT / "PHASE18_ADIM0_GATE.csv")
    assert len(gate) == 21
    assert gate["exact_match"].all()
    assert gate.query("gate == 'phase17_primary_81'").iloc[0]["successes"] == 81


def test_point_grid_matched_scale_and_bootstrap_shapes() -> None:
    grid = pd.read_csv(OUT / "PHASE18_ADIM2_POINT_GRID.csv")
    matched = pd.read_csv(OUT / "PHASE18_ADIM2_MATCHED_SCALE.csv")
    boot = pd.read_csv(OUT / "PHASE18_ADIM2_BOOT_CELLS.csv")
    assert len(grid) == 52_200
    assert len(matched) == 1_800
    assert len(boot) == 1_800
    assert (matched["v1_WI5"] == matched["matched_v2_WI5"]).all()
    assert not boot.query("truth_class == 'B'")["certified_BH"].any()
    assert not boot.query("truth_class == 'B'")["certified_BY"].any()


def test_auc_kstar_and_real_profile_results() -> None:
    auc = pd.read_csv(OUT / "PHASE18_ADIM3_AUC.csv")
    assert len(auc) == 116
    primary = auc.query("variant == 'v2_WI6' and k == 1 and kappa == 0.25").set_index("positive_class")
    assert primary.loc["A", "auc_feature_repeat"] > 0.95
    assert primary.loc["D", "auc_feature_repeat"] > 0.94
    kstar = pd.read_csv(OUT / "PHASE18_ADIM4_KSTAR.csv")
    assert len(kstar) == 2_488
    assert set(kstar["k_star"].dropna()) == {0.25}
    held_out = kstar.query("dataset == 'EVAL_boot_feature'")
    assert not held_out["certified_BH"].any()
    assert not held_out["certified_BY"].any()
    profile = pd.read_csv(OUT / "PHASE18_ADIM5_KPROFILE.csv")
    entries = pd.read_csv(OUT / "PHASE18_ADIM5_KENTER_81.csv")
    assert len(profile) == 2_187
    assert len(entries) == 81
    assert entries.groupby("model")["right_censored_128"].apply(lambda x: (~x).sum()).to_dict() == {
        "CoxPH-v2": 8,
        "DeepSurv": 0,
        "RSF": 2,
    }


def test_frozen_hashes_and_required_report_ending() -> None:
    for filename, expected in HASHES.items():
        actual = hashlib.sha256((ROOT / "warrantxai" / filename).read_bytes()).hexdigest()
        assert actual == expected
    required = (
        "Bu analiz v1 ve v2 birincil ayarlarını, Tablo 3–9'u, BH/BY kararlarını, "
        "saklı atıf tensörlerini, öznitelik panelini ve CoxPH-v2 retry protokolünü "
        "DEĞİŞTİRMEMEKTEDİR."
    )
    report = (ROOT / "PHASE18_KAPANIS.md").read_text(encoding="utf-8")
    assert report.rstrip().endswith(required)
