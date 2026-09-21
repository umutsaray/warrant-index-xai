import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "manuscript_assets"


def test_all_publication_assets_exist_and_pass_format_checks() -> None:
    metadata = json.loads((ASSETS / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["new_experiment_or_model_call"] is False
    assert metadata["privacy"]["patient_or_sample_identifiers_exported"] is False
    for number in range(1, 5):
        for suffix in ("png", "pdf", "eps"):
            assert (ASSETS / f"figure_{number}.{suffix}").stat().st_size > 0
        output = metadata["outputs"][f"figure_{number}"]
        assert output["png_dpi"] == [600.0, 600.0]
        assert output["pdf_contains_embedded_truetype_font"] is True
        assert output["physical_size_mm"][1] < 225
    assert abs(metadata["outputs"]["figure_3"]["physical_size_mm"][0] - 85) < 0.1
    for number in (1, 2, 4):
        assert abs(metadata["outputs"][f"figure_{number}"]["physical_size_mm"][0] - 170) < 0.1


def test_tables_are_booktabs_and_numerically_frozen() -> None:
    table1 = (ASSETS / "table_1.tex").read_text(encoding="utf-8")
    table2 = (ASSETS / "table_2.tex").read_text(encoding="utf-8")
    for table in (table1, table2):
        assert all(command in table for command in (r"\toprule", r"\midrule", r"\botrule"))
        assert not re.search(r"\d,\d{3}", table)
        assert "color" not in table.lower()
    assert "22/77 & 8/77  & 8/77  & 36.4" in table1
    assert "20/77 & 10/77 & 10/77 & 50.0" in table1
    assert "8/77  & 0/77  & 0/77  & 0.0" in table1
    assert "CoxPH ($N=50$) & 0 & 0 & 0 & 1 & 76 & 77" in table2
    assert "RSF ($N=50$) & 0 & 0 & 0 & 5 & 72 & 77" in table2
    assert "DeepSurv ($N=50$) & 0 & 0 & 0 & 0 & 77 & 77" in table2

