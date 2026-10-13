from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_phase13_outputs_cover_fold0_and_saved_fold1_analogs():
    split = pd.read_csv(ROOT / "PHASE13_ADIM1_SPLIT_SIDES.csv")
    paths = pd.read_csv(ROOT / "PHASE13_ADIM2_PATH_VISITS.csv")
    effects = pd.read_csv(ROOT / "PHASE13_ADIM3_EXTREME_COALITION.csv")
    summary = pd.read_csv(ROOT / "PHASE13_FEATURE_SUMMARY.csv")
    assert len(summary) == 27
    assert set(summary.fold) == {0}
    assert set(split.fold) == {0, 1}
    assert set(paths.state) == {"all_active_query", "all_active_base"}
    assert {"delta_last", "delta_first"}.issubset(effects.columns)
    assert set(split.loc[split.fold.eq(1), "feature"]) == {"AJCC Stage IV", "OR13C3"}


def test_phase13_h1_counts_and_split_reference():
    summary = pd.read_csv(ROOT / "PHASE13_FEATURE_SUMMARY.csv")
    zero = summary[summary.TCGA_M_j.eq(0)]
    assert len(zero) == 11
    assert zero.TCGA_separating_max.eq(0).all()
    assert zero.TCGA_nonzero_delta_last.eq(0).all()
    assert zero.TCGA_nonzero_delta_first.eq(0).all()
    expected = pd.read_csv(ROOT / "results/phase76_rsf_split_usage.csv").set_index("Feature ID")
    actual = summary.set_index("feature_id")
    pd.testing.assert_series_equal(
        actual.loc[expected.index, "split_trees"].astype(int),
        expected["Trees with ≥1 split"].astype(int),
        check_names=False,
    )
    report = (ROOT / "PHASE13_RSF_TRAVERSAL.md").read_text(encoding="utf-8")
    assert "H1'in aynı-taraf veya hiç-ziyaret koşulu **11/11**" in report
    assert "DEĞİŞTİRMEMEKTEDİR" in report
