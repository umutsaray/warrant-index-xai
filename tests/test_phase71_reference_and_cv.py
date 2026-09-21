import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_reference_models_are_frozen_and_versioned() -> None:
    payload = json.loads(
        (ROOT / "results/phase7_retrain/reference_models.json").read_text(encoding="utf-8")
    )
    assert payload["protected_artifacts_unchanged"] is True
    assert payload["models"]["CoxPH"]["warnings"] == []
    assert abs(payload["models"]["CoxPH"]["test_c_index"] - 0.667704280155642) < 1e-12
    assert payload["models"]["CoxPH"]["path"].endswith("coxph_baseline_v2.joblib")
    assert payload["models"]["RSF"]["unchanged_from_phase70"] is True
    assert payload["models"]["DeepSurv"]["unchanged_from_phase70"] is True
    assert (ROOT / "results/phase7_retrain/coxph_baseline.joblib").is_file()
    assert (ROOT / "results/phase7_retrain/coxph_baseline_v2.joblib").is_file()


def test_five_fold_three_model_cv_values() -> None:
    payload = json.loads((ROOT / "results/phase71_cv/cv_summary.json").read_text(encoding="utf-8"))
    expected = {
        "CoxPH": [0.667704280155642, 0.7058823529411765, 0.6811594202898551, 0.6551133932391956, 0.6501726121979287],
        "RSF": [0.5385214007782101, 0.7896092743666809, 0.636231884057971, 0.7004706889174155, 0.6507479861910241],
        "DeepSurv": [0.4964980544747082, 0.7191927866036926, 0.65, 0.5318784766795036, 0.6628308400460299],
    }
    for model, scores in expected.items():
        assert payload["models"][model]["scores"] == scores
        assert abs(payload["models"][model]["mean"] - sum(scores) / 5) < 1e-12
        assert payload["models"][model]["minimum"] == min(scores)
        assert payload["models"][model]["maximum"] == max(scores)
    assert payload["coxph_retry_rule"]["retried_folds"] == [0, 2]
    assert payload["coxph_retry_rule"]["all_final_cv_fits_warning_free"] is True
    assert payload["protected_reference_artifacts_unchanged"] is True
    assert payload["audit_pipeline_run"] is False


def test_phase71_reports_distinguish_three_score_types_and_hide_ids() -> None:
    names = ["PHASE7_REFERENCE_MODELS.md", "PHASE7_CV_SUMMARY.md", "RESULTS_LEDGER_PHASE7.md"]
    text = "\n".join((ROOT / name).read_text(encoding="utf-8") for name in names)
    assert "Eski, kirli split C-index" in text
    assert "Audit referansı: temiz fold 0" in text
    assert "5-fold CV genelleme tahmini" in text
    assert "coxph_baseline_v2.joblib" in text
    assert "F/S/U/G/G_signed/R" in text
    assert not re.search(r"\bTCGA-[A-Z0-9]{2}-[A-Z0-9]{4}\b", text)
