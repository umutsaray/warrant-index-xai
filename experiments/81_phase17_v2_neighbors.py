#!/usr/bin/env python3
"""Phase 17 preregistered neighboring-setting v2 certifications."""
from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase17"
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]

from _shared import CALIBRATION
from warrantxai.certification import certify

SCALES = {
    "CoxPH-v2": 0.023934215222165295,
    "RSF": 0.01707069437342852,
    "DeepSurv": 0.012307440008774696,
}
COMPONENTS = ("F", "S", "U", "G_claim", "R", "N")
SETTINGS = (
    (0.5, 1.0, 0.25),
    (2.0, 1.0, 0.25),
    (1.0, 1.0, 0.10),
    (1.0, 1.0, 0.50),
)


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    started = time.perf_counter()
    primary = load("77_phase17_v2_certification.py", "phase17_primary_helpers")
    cert = load("43_phase72_certification.py", "phase17_neighbor_cert")
    shared = load("12_rsf_certification.py", "phase17_neighbor_shared")
    remedy = load("11_rsf_remedy_calibration.py", "phase17_neighbor_remedy")
    rng = np.random.default_rng(cert.SEED + 725_000)
    indices = rng.integers(0, 50, size=(cert.DRAWS, 50))
    original = {key: CALIBRATION[key] for key in ("a", "b", "u", "eta")}
    rows: list[dict[str, object]] = []
    try:
        for model, source in (
            ("CoxPH-v2", "phase72_cox"),
            ("RSF", "phase72_rsf"),
            ("DeepSurv", "phase72_deepsurv"),
        ):
            scale = SCALES[model]
            if model == "CoxPH-v2":
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
                sites = [phi[0], phi[1]]

            for k, rho, kappa in SETTINGS:
                CALIBRATION.update(a=k * scale, b=k * scale, u=6 * k * scale, eta=rho * scale)
                if model == "CoxPH-v2":
                    points, pvalues = cert.cox_statistics(indices)
                else:
                    points, pvalues = shared.component_statistics(phi, reference, indices)
                n_point, n_p = primary.n_stats(
                    sites, kappa * scale, indices, shared, remedy
                )
                component_points = np.column_stack([points, n_point])
                component_p = np.column_stack([pvalues, n_p])
                bh = certify(component_p, 0.05, "bh")
                by = certify(component_p, 0.05, "by")
                for j, feature in enumerate(features):
                    values = component_points[j]
                    rows.append(
                        {
                            "model": model,
                            "feature": feature,
                            "k": k,
                            "rho": rho,
                            "kappa": kappa,
                            "scale_s_m": scale,
                            **{name: values[q] for q, name in enumerate(COMPONENTS)},
                            "WI6": values.min(),
                            "bottleneck6": COMPONENTS[int(values.argmin())],
                            "point_adequate6": values.min() > 0.8,
                            "p_joint": bh.joint_p[j],
                            "q_BH": bh.adjusted_p[j],
                            "q_BY": by.adjusted_p[j],
                            "certified_BH": bh.rejected[j],
                            "certified_BY": by.rejected[j],
                        }
                    )
    finally:
        CALIBRATION.update(original)
    output = pd.DataFrame(rows)
    output.to_csv(OUT / "PHASE17_ADIM3_NEIGHBORS.csv", index=False)
    print(output.groupby(["model", "k", "kappa"])[["point_adequate6", "certified_BH", "certified_BY"]].sum())
    print("elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
