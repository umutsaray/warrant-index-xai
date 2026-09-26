from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_phase10_fold1_harmonization_is_explicitly_partial():
    data = pd.read_csv(ROOT / "PHASE10_FOLD1_HARMONIZED_27FEATURE.csv")
    assert len(data) == 81
    assert set(data["model"]) == {"COXPH", "RSF", "DeepSurv"}
    availability = data.groupby("feature")["available_in_fold1_model_panel"].first()
    assert int(availability.sum()) == 10
    assert int((~availability).sum()) == 17
    missing = data.loc[~data["available_in_fold1_model_panel"]]
    assert missing["fold1_WI"].isna().all()


def test_phase10_retry_depth_decision_sets():
    data = pd.read_csv(ROOT / "PHASE10_COXPH_RETRY_BIAS_DEPTH.csv")
    assert len(data) == 27
    assert int(data["point_direct22"].sum()) == 9
    assert int(data["BH_full50"].sum()) == 8
    assert int(data["BY_full50"].sum()) == 8
    assert int(data["BH_direct22"].sum()) == 6
    assert int(data["BY_direct22"].sum()) == 6
    dropped = set(data.loc[data["BH_full50"] & ~data["BH_direct22"], "symbol"])
    added = set(data.loc[~data["BH_full50"] & data["BH_direct22"], "symbol"])
    assert dropped == {"CA1", "OR13C3"}
    assert added == set()


def test_phase10_query_resampling_grid_and_original_reproduction():
    data = pd.read_csv(ROOT / "PHASE10_QUERY_RESAMPLING_SENSITIVITY.csv")
    assert len(data) == 2 * 27 * 7
    assert set(data["model"]) == {"RSF", "DeepSurv"}
    assert set(data["query_per_cohort"]) == {12, 24, 48}
    assert (data.groupby(["model", "feature"]).size() == 7).all()
    assert (data.query("query_per_cohort == 12").groupby("model").size() == 27).all()
    for model, path in {
        "RSF": ROOT / "results" / "phase72_rsf" / "feature_report.csv",
        "DeepSurv": ROOT / "results" / "phase72_deepsurv" / "feature_report.csv",
    }.items():
        original = data[(data.model == model) & (data.query_per_cohort == 12)].set_index("feature")
        frozen = pd.read_csv(path).set_index("feature")
        for column in ["F", "S", "U", "G_claim", "R"]:
            assert (original[column] - frozen[column]).abs().max() < 1e-12
