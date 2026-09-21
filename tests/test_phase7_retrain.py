import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_validation() -> dict:
    return json.loads(
        (ROOT / "results/phase7_retrain/validation.json").read_text(encoding="utf-8")
    )


def test_phase7_patient_grouped_split_and_frozen_settings() -> None:
    payload = load_validation()
    assert payload["split"] == {
        "train_rows": 674,
        "test_rows": 168,
        "train_unique_patients": 616,
        "test_unique_patients": 153,
        "sample_overlap": 0,
        "patient_overlap": 0,
        "train_event_rate": 0.1572700296735905,
        "test_event_rate": 0.15476190476190477,
    }
    assert payload["feature_selection"]["selection_scope"] == "train_only"
    assert payload["feature_selection"]["feature_count"] == 616
    assert payload["feature_selection"]["block_counts"] == {
        "rna": 200,
        "mirna": 200,
        "methylation": 200,
        "clinical": 16,
    }
    assert all(payload["assertions"].values())


def test_phase7_scores_and_artifact_isolation() -> None:
    payload = load_validation()
    expected = {
        "CoxPH": (0.7064765300059418, 0.6311284046692607),
        "RSF": (0.714795008912656, 0.5385214007782101),
        "DeepSurv": (0.7456922162804516, 0.4964980544747082),
    }
    for model, (old, new) in expected.items():
        assert abs(payload["scores"][model]["old"] - old) < 1e-12
        assert abs(payload["scores"][model]["new"] - new) < 1e-12
    assert payload["deepsurv"]["new_best_epoch"] == 38
    assert payload["deepsurv"]["stopped_epoch"] == 68
    assert payload["coxph"]["warnings"]
    assert payload["old_artifact_hashes"]["rsf_baseline.joblib"] == (
        "3ab8da5972812ec3411430419ec40862306316504b270a77a9b1a9ce96121574"
    )


def test_phase7_reports_have_no_patient_barcodes_or_audit_claims() -> None:
    reports = "\n".join(
        (ROOT / name).read_text(encoding="utf-8")
        for name in ("SPLIT_REDESIGN.md", "MODEL_RETRAIN_COMPARISON.md")
    )
    assert not re.search(r"\bTCGA-[A-Z0-9]{2}-[A-Z0-9]{4}\b", reports)
    assert "Faz 7.0'da F/S/U/G/G_signed/R" in reports
    assert "patient_overlap_zero" not in reports
