#!/usr/bin/env python3
"""Phase 19 preregistration and Phase 18 reproduction gate."""
from __future__ import annotations

import hashlib
import json
import subprocess
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase19"
OUT.mkdir(parents=True, exist_ok=True)
FORMULA_FILES = ("components.py", "aggregation.py", "certification.py")


def sha256(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest().upper()


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def main() -> None:
    started = time.perf_counter()
    prereg_commit = git("log", "-1", "--format=%H", "--", "PHASE19_PREREG.md")
    prereg_status = git("status", "--short", "--", "PHASE19_PREREG.md")
    assert prereg_commit == "d839855acef9e87a62801bbfe60bb3c1ad834825"
    assert prereg_status == ""

    entries = pd.read_csv(ROOT / "results/phase18/PHASE18_ADIM5_KENTER_81.csv")
    profile = pd.read_csv(ROOT / "results/phase18/PHASE18_ADIM5_KPROFILE.csv")
    checks: list[dict[str, object]] = []

    expected_entries = {"CoxPH-v2": 8, "RSF": 2, "DeepSurv": 0}
    for model, expected in expected_entries.items():
        observed = int((~entries.loc[entries.model == model, "right_censored_128"]).sum())
        checks.append({
            "check": "entered_by_k_le_128",
            "model": model,
            "component": "all",
            "expected": expected,
            "observed": observed,
            "passed": observed == expected,
        })
        assert observed == expected

    expected_binding = {
        "CoxPH-v2": {"R": 19, "N": 8},
        "RSF": {"N": 17, "R": 10},
        "DeepSurv": {"R": 17, "N": 10},
    }
    cell = profile[(profile.k == 64.0) & (profile.kappa == 0.25) & (profile.rho == 1.0)]
    assert len(cell) == 81
    for model, expected_counts in expected_binding.items():
        counts = cell.loc[cell.model == model, "bottleneck6"].value_counts().to_dict()
        for component, expected in expected_counts.items():
            observed = int(counts.get(component, 0))
            checks.append({
                "check": "k64_binding_count",
                "model": model,
                "component": component,
                "expected": expected,
                "observed": observed,
                "passed": observed == expected,
            })
            assert observed == expected
        assert sum(counts.values()) == 27
        assert set(counts) == set(expected_counts)

    output = pd.DataFrame(checks)
    output.to_csv(OUT / "PHASE19_ADIM0_GATE.csv", index=False)
    metadata = {
        "preregistration_commit": prereg_commit,
        "phase19_results_absent_before_preregistration": True,
        "checks": len(output),
        "all_passed": bool(output.passed.all()),
        "formula_sha256": {
            name: sha256(ROOT / "warrantxai" / name) for name in FORMULA_FILES
        },
        "elapsed_seconds": time.perf_counter() - started,
    }
    (OUT / "PHASE19_ADIM0_METADATA.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(output.to_string(index=False))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
