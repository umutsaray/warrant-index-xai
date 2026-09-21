"""Reproducible experiment utilities (Section C.3); no data-driven tolerances."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import scipy
import sklearn

from warrantxai.aggregation import adequate, warrant_index
from warrantxai.components import COMPONENT_NAMES

# assumption: illustrative tolerances fixed BEFORE all evaluation seeds/grids.
# They are not optimized to make the scenarios pass, nor universal thresholds.
CALIBRATION = {"a": 0.5, "b": 0.5, "u": 3.0, "eta": 0.02, "tau": 0.8, "alpha_U": 0.05}


def markdown_table(frame: pd.DataFrame) -> str:
    """Render results without an extra tabulate dependency (Section C.3).

    Returns
    -------
    str
        Portable Markdown table, with NA preserved.
    """

    def cell(value: object) -> str:
        """Format a result cell (Section C.3).

        Returns
        -------
        str
            Human-readable cell.
        """
        if isinstance(value, (float, np.floating)):
            return "NA" if np.isnan(value) else f"{value:.4f}"
        return str(value).replace("|", "/").replace("\n", " ")

    rows = ["| " + " | ".join(map(str, frame.columns)) + " |", "| " + " | ".join(["---"] * len(frame.columns)) + " |"]
    rows.extend("| " + " | ".join(cell(v) for v in row) + " |" for row in frame.itertuples(index=False, name=None))
    return "\n".join(rows)


def score_rows(z: np.ndarray, **labels: object) -> list[dict]:
    """Create feature rows without altering any WI formula (Section C.3).

    Returns
    -------
    list of dict
        Components, WI and non-inferential point adequacy.
    """
    wi = warrant_index(z)
    decision = adequate(z, CALIBRATION["tau"])
    return [
        {
            **labels,
            "feature": j,
            **dict(zip(COMPONENT_NAMES, z[j], strict=True)),
            "WI": wi[j],
            "point_adequate": decision[j],
        }
        for j in range(len(z))
    ]


def save_run(output: Path, stem: str, frame: pd.DataFrame, report: str, started: float, args: dict) -> None:
    """Write CSV, Markdown and machine-readable provenance (Section C.3/E).

    Returns
    -------
    None
        Creates reproducible local artifacts, never uploads them.
    """
    output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output / f"{stem}.csv", index=False)
    (output / f"{stem}.md").write_text(report, encoding="utf-8")
    package_root = Path(__file__).resolve().parents[1]
    hashes = {
        str(p.relative_to(package_root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for directory in ("warrantxai", "tests", "experiments")
        for p in sorted((package_root / directory).rglob("*.py"))
    }
    metadata = {
        "arguments": args,
        "fixed_calibration": CALIBRATION,
        "seconds": perf_counter() - started,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "sklearn": sklearn.__version__,
        "source_sha256": hashes,
        "scope": "synthetic smoke/validation study, not proof of comparative superiority",
    }
    (output / f"{stem}_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
