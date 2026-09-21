"""Checks for the report-only Phase 7.3 feature-level exports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_phase73_full_tables_match_frozen_phase72_counts() -> None:
    expected = {
        "COXPH": (9, 8, 8),
        "RSF": (17, 12, 12),
        "DEEPSURV": (0, 0, 0),
    }
    for model, counts in expected.items():
        table = pd.read_csv(ROOT / f"PHASE73_{model}_FULLTABLE.csv")
        assert len(table) == 27
        assert table["Feature ID"].is_unique
        observed = tuple((table[column] == "evet").sum() for column in (
            "Point adequate", "BH certified", "BY certified"
        ))
        assert observed == counts


def test_phase73_intersection_and_directionality_are_complete() -> None:
    cox = pd.read_csv(ROOT / "PHASE73_COXPH_FULLTABLE.csv")
    rsf = pd.read_csv(ROOT / "PHASE73_RSF_FULLTABLE.csv")
    cset = set(cox.loc[cox["BH certified"] == "evet", "Feature ID"])
    rset = set(rsf.loc[rsf["BH certified"] == "evet", "Feature ID"])
    assert (len(cset - rset), len(cset & rset), len(rset - cset), len(cset | rset)) == (2, 6, 6, 14)
    directions = pd.read_csv(ROOT / "PHASE73_DIRECTIONALITY.csv")
    cross = pd.read_csv(ROOT / "PHASE73_DIRECTIONALITY_INTERSECTION.csv")
    assert len(directions) == len(cset) + len(rset) == 20
    assert len(cross) == len(cset & rset) == 6
    assert not cross["Same sign on matched queries"].any()
    assert not any(
        token in column.lower()
        for table in (cox, rsf, directions, cross)
        for column in table.columns
        for token in ("patient_id", "sample_id", "sample_key")
    )
