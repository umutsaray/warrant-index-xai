import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase70b_collapse_diagnostic"


def test_five_fold_scores_and_cox_retries() -> None:
    folds = json.loads((OUT / "fold_scores.json").read_text(encoding="utf-8"))["folds"]
    assert len(folds) == 5
    assert [row["fold"] for row in folds] == [0, 1, 2, 3, 4]
    assert all(row["selected_feature_count"] == 616 for row in folds)
    assert [round(row["rsf_test_c_index"], 6) for row in folds] == [
        0.538521,
        0.789609,
        0.636232,
        0.700471,
        0.650748,
    ]
    retries = json.loads((OUT / "cox_fold_retries.json").read_text(encoding="utf-8"))
    assert set(retries) == {"0", "2"}
    assert retries["0"]["retry_warnings"] == []
    assert abs(retries["0"]["retry_c_index"] - 0.667704280155642) < 1e-12
    assert abs(retries["2"]["retry_c_index"] - 0.6811594202898551) < 1e-12


def test_feature_selection_and_deepsurv_diagnostics() -> None:
    summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
    stability = {row["block"]: row for row in summary["old_new_feature_stability"]}
    assert stability["rna"]["common_count"] == 43
    assert stability["mirna"]["common_count"] == 151
    assert stability["methylation"]["common_count"] == 123
    assert stability["clinical"]["common_count"] == 16
    assert stability["total"]["common_count"] == 333
    assert summary["grouped_fold_pairwise_feature_stability"]["total"] == {
        "mean_common": 304.9,
        "minimum": 277,
        "maximum": 328,
    }
    row = summary["deepsurv"]["row_level_validation"]
    grouped = summary["deepsurv"]["patient_grouped_validation"]
    assert row["fit_validation_patient_overlap"] == 15
    assert row["best_validation_epoch"] == 38
    assert row["maximum_observed_test_epoch"] == 8
    assert grouped["fit_validation_patient_overlap"] == 0
    assert grouped["best_checkpoint_test_c_index"] < 0.5
    assert row["validation_test_spearman_rho"] < 0
    assert grouped["validation_test_spearman_rho"] < 0


def test_phase70b_scope_and_report_privacy() -> None:
    summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
    assert summary["audit_pipeline_run"] is False
    assert summary["patient_or_sample_identifiers_exported"] is False
    report = (ROOT / "PHASE7_COLLAPSE_DIAGNOSTIC.md").read_text(encoding="utf-8")
    assert "tek bir leakage-düzeltilmiş" in report
    assert "F/S/U/G/G_signed/R" in report
    assert not re.search(r"\bTCGA-[A-Z0-9]{2}-[A-Z0-9]{4}\b", report)
