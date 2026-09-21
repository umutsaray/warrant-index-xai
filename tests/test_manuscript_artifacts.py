import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_manuscript_consistency_report_passes() -> None:
    payload = json.loads(
        (ROOT / "results/manuscript_consistency/consistency.json").read_text(encoding="utf-8")
    )
    assert payload["passed"] is True
    assert payload["mismatch_count"] == 0
    assert payload["check_count"] >= 70
    assert payload["scientific_computation_performed"] is False


def test_manuscript_scope_and_source_warning_are_explicit() -> None:
    draft = (ROOT / "MANUSCRIPT_METHODS_TO_AVAILABILITY.md").read_text(encoding="utf-8")
    ledger = (ROOT / "RESULTS_LEDGER.md").read_text(encoding="utf-8")
    assert "## Abstract" not in draft
    assert "## 1. Introduction" not in draft
    assert "section2-warrant-index-taslak.md` is absent" in ledger
    assert "[TODO: Claude/kullanıcı ile netleştirilecek — public repo + lisans kararı bekliyor]" in draft

