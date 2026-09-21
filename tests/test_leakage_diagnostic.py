import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_leakage_diagnostic_counts_and_scores() -> None:
    payload = json.loads(
        (ROOT / "results/leakage_diagnostic/diagnostic.json").read_text(encoding="utf-8")
    )
    assert "patient-overlapping (dirty) split" in payload["verdict"]
    assert payload["counts"]["overlap_unique_patients"] == 23
    assert payload["counts"]["overlap_train_rows"] == 25
    assert payload["counts"]["overlap_test_rows"] == 23
    assert payload["counts"]["disjoint_test_rows"] == 146
    assert payload["counts"]["disjoint_test_unique_patients"] == 141
    assert payload["deepsurv_internal_validation"]["overlap_unique_patients"] == 12
    expected = {
        "CoxPH": (0.7064765300059418, 0.700507614213198),
        "RSF": (0.714795008912656, 0.7292724196277496),
        "DeepSurv": (0.7456922162804516, 0.7639593908629442),
    }
    for model, (original, disjoint) in expected.items():
        assert abs(payload["scores"][model]["reproduced_original_169"] - original) < 1e-12
        assert abs(payload["scores"][model]["disjoint_146"] - disjoint) < 1e-12
    assert payload["hash_cross_checks_passed"] is True
    assert payload["patient_or_sample_identifiers_exported"] is False


def test_leakage_report_has_no_tcga_patient_barcode() -> None:
    report = (ROOT / "LEAKAGE_DIAGNOSTIC.md").read_text(encoding="utf-8")
    assert "Modeller örnek-düzeyinde ayrık fakat hasta-düzeyinde çakışmalı (kirli) split" in report
    assert not re.search(r"\bTCGA-[A-Z0-9]{2}-[A-Z0-9]{4}\b", report)

