from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_phase12_artifacts_are_complete_and_diagnostic():
    scale = pd.read_csv(ROOT / "PHASE12_ADIM1_SCALE_GRID.csv")
    tie = pd.read_csv(ROOT / "PHASE12_ADIM2_TIE_R.csv")
    query = pd.read_csv(ROOT / "PHASE12_ADIM3_QUERY_COHORT.csv")
    agg = pd.read_csv(ROOT / "PHASE12_ADIM4_AGGREGATORS.csv")
    assert len(scale) == 3 * 27 * 5
    assert len(tie) == 3 * 27 * 8
    assert set(query.fold) == {0, 1}
    assert len(agg) == 3 * 27
    assert scale.groupby("model").scale_s_m.nunique().eq(1).all()
    assert tie[tie.grid.eq("absolute")].groupby("model").eta.nunique().eq(4).all()
    assert agg.groupby("model").min_gt_0_8.sum().to_dict() == {
        "CoxPH-v2": 9, "DeepSurv": 0, "RSF": 17
    }
    report = (ROOT / "PHASE12_SCALE_TIE_QUERY_AUDIT.md").read_text(encoding="utf-8")
    assert "diagnostic/exploratory" in report
    assert "BH/BY kararlarını ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR" in report


def test_phase12_primary_eta_reproduces_phase74():
    tie = pd.read_csv(ROOT / "PHASE12_ADIM2_TIE_R.csv")
    old = pd.read_csv(ROOT / "PHASE74_ETA_FEATURE_VALUES.csv")
    names = {"CoxPH-v2": "CoxPH", "RSF": "RSF", "DeepSurv": "DeepSurv"}
    for new_name, old_name in names.items():
        new = tie[(tie.model == new_name) & (tie.grid == "absolute") & (tie.eta == .02)].set_index("feature").R
        reference = old[(old.model == old_name) & (old.eta == .02)].set_index("feature").R
        pd.testing.assert_series_equal(new.loc[reference.index], reference, check_names=False, atol=1e-12)
