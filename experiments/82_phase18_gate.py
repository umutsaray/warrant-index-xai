#!/usr/bin/env python3
"""Phase 18 gate: independently reproduce Phase 17 primary and synthetic outputs."""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase18"
OUT.mkdir(parents=True, exist_ok=True)
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]

from _shared import CALIBRATION
from warrantxai.certification import certify


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def reproduce_primary() -> pd.DataFrame:
    primary = load("77_phase17_v2_certification.py", "phase18_gate_primary")
    cert = load("43_phase72_certification.py", "phase18_gate_cert")
    shared = load("12_rsf_certification.py", "phase18_gate_shared")
    remedy = load("11_rsf_remedy_calibration.py", "phase18_gate_remedy")
    rng = np.random.default_rng(cert.SEED + 725_000)
    indices = rng.integers(0, 50, size=(cert.DRAWS, 50))
    original = {key: CALIBRATION[key] for key in ("a", "b", "u", "eta")}
    rows = []
    try:
        for model, source in (
            ("CoxPH-v2", "phase72_cox"),
            ("RSF", "phase72_rsf"),
            ("DeepSurv", "phase72_deepsurv"),
        ):
            scale = primary.SCALES[model]
            CALIBRATION.update(a=scale, b=scale, u=6 * scale, eta=scale)
            if model == "CoxPH-v2":
                points, pvalues = cert.cox_statistics(indices)
                with np.load(ROOT / "results/phase72_cox/attributions_n50.npz") as saved:
                    full = saved["features"].astype(str).tolist()
                    features = saved["active_features"].astype(str).tolist()
                    selected = np.array([full.index(feature) for feature in features])
                    sites = [saved["source_phi"][:, :, selected], saved["external_phi"]]
            else:
                with np.load(ROOT / "results" / source / "attributions_n50.npz") as saved:
                    phi = saved["candidate"]
                    reference = saved["reference"]
                    features = saved["features"].astype(str).tolist()
                points, pvalues = shared.component_statistics(phi, reference, indices)
                sites = [phi[0], phi[1]]
            n_point, n_p = primary.n_stats(sites, 0.25 * scale, indices, shared, remedy)
            values = np.column_stack([points, n_point])
            component_p = np.column_stack([pvalues, n_p])
            bh = certify(component_p, 0.05, "bh")
            by = certify(component_p, 0.05, "by")
            for j, feature in enumerate(features):
                rows.append(
                    {
                        "model": model,
                        "feature": feature,
                        **{name: values[j, q] for q, name in enumerate(primary.COMP)},
                        "WI6": values[j].min(),
                        "p_joint": bh.joint_p[j],
                        "certified_BH": bh.rejected[j],
                        "certified_BY": by.rejected[j],
                    }
                )
    finally:
        CALIBRATION.update(original)
    result = pd.DataFrame(rows)
    frozen = pd.read_csv(OUT.parent / "phase17/PHASE17_ADIM3_V2_CERT_PRIMARY.csv")
    frozen = frozen.set_index(["model", "feature"]).loc[
        pd.MultiIndex.from_frame(result[["model", "feature"]])
    ].reset_index()
    numeric = list(primary.COMP) + ["WI6", "p_joint"]
    np.testing.assert_allclose(result[numeric], frozen[numeric], atol=1e-12, rtol=0)
    assert result["certified_BH"].tolist() == frozen["certified_BH"].tolist()
    assert result["certified_BY"].tolist() == frozen["certified_BY"].tolist()
    return result


def reproduce_synthetic() -> pd.DataFrame:
    syn = load("80_phase17_synthetic_feature_level.py", "phase18_gate_syn")
    shared = load("12_rsf_certification.py", "phase18_gate_syn_shared")
    remedy = load("11_rsf_remedy_calibration.py", "phase18_gate_syn_remedy")
    rng = np.random.default_rng(20261108)
    n, p, t, draws = 700, 30, 50, 999
    classes = np.array(["A"] * 8 + ["B"] * 6 + ["C"] * 6 + ["D"] * 5 + ["E"] * 5)
    beta0 = np.r_[np.full(8, 0.05), np.zeros(6), np.full(6, 0.008), np.full(5, 0.05), np.full(5, 0.03)]
    beta1 = beta0.copy()
    beta1[20:25] = -0.05
    covariance = np.full((p, p), 0.3)
    np.fill_diagonal(covariance, 1)
    checkpoint_csv = OUT / "PHASE18_ADIM0_GATE_SYNTH_CHECKPOINT.csv"
    checkpoint_state = OUT / "PHASE18_ADIM0_GATE_SYNTH_RNG_STATE.json"
    rows = []
    start_repeat = 0
    if checkpoint_csv.exists() and checkpoint_state.exists():
        checkpoint = pd.read_csv(checkpoint_csv)
        rows = checkpoint.to_dict("records")
        start_repeat = int(checkpoint["repeat"].max()) + 1
        rng.bit_generator.state = json.loads(checkpoint_state.read_text(encoding="utf-8"))
        print(f"resuming synthetic gate at {start_repeat + 1}/20", flush=True)
    for repeat in range(start_repeat, 20):
        x = [rng.multivariate_normal(np.zeros(p), covariance, n) for _ in range(2)]
        y = [x[0] @ beta0 + rng.normal(0, 0.1, n), x[1] @ beta1 + rng.normal(0, 0.1, n)]
        phi = np.empty((2, t, n, p))
        for site in range(2):
            for replica in range(t):
                sampled = rng.integers(0, n, n)
                coefficient = syn.ridge(x[site][sampled], y[site][sampled])
                phi[site, replica] = coefficient * (x[site] - x[site].mean(0))
        indices = rng.integers(0, t, size=(draws, t))
        _, _, v1_bh, v1_by = syn.stats(phi, indices, shared, remedy, 0.5, 0.5, 3, 0.02)
        importance = np.abs(phi[0]).mean(1).mean(0)
        scale = float(np.median(importance[importance > 0]))
        _, _, v2_bh, v2_by = syn.stats(
            phi, indices, shared, remedy, scale, scale, 6 * scale, scale, True, scale
        )
        for feature in range(p):
            rows.append(
                {
                    "repeat": repeat,
                    "feature": feature,
                    "truth_class": classes[feature],
                    "v1_BH": v1_bh[feature],
                    "v1_BY": v1_by[feature],
                    "v2_BH": v2_bh[feature],
                    "v2_BY": v2_by[feature],
                }
            )
        pd.DataFrame(rows).to_csv(checkpoint_csv, index=False)
        checkpoint_state.write_text(
            json.dumps(rng.bit_generator.state, indent=2), encoding="utf-8"
        )
        print(f"synthetic gate {repeat + 1}/20", flush=True)
    result = pd.DataFrame(rows)
    frozen = pd.read_csv(OUT.parent / "phase17/PHASE17_ADIM5_SYNTH.csv")
    columns = ["repeat", "feature", "truth_class", "v1_BH", "v1_BY", "v2_BH", "v2_BY"]
    pd.testing.assert_frame_equal(result[columns], frozen[columns], check_dtype=False)
    return result


def main() -> None:
    started = time.perf_counter()
    primary = reproduce_primary()
    synthetic = reproduce_synthetic()
    summary = []
    for method in ("v1_BH", "v1_BY", "v2_BH", "v2_BY"):
        for truth_class, group in synthetic.groupby("truth_class"):
            summary.append(
                {
                    "gate": "phase17_synthetic_power",
                    "method": method,
                    "truth_class": truth_class,
                    "successes": int(group[method].sum()),
                    "total": len(group),
                    "rate": float(group[method].mean()),
                    "exact_match": True,
                }
            )
    summary.append(
        {
            "gate": "phase17_primary_81",
            "method": "all",
            "truth_class": "NA",
            "successes": len(primary),
            "total": 81,
            "rate": 1.0,
            "exact_match": True,
        }
    )
    pd.DataFrame(summary).to_csv(OUT / "PHASE18_ADIM0_GATE.csv", index=False)
    print("gate passed; elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
