#!/usr/bin/env python3
"""Phase 18 real-data point k-profile and entry-k diagnostics."""
from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase18"
OUT.mkdir(parents=True, exist_ok=True)
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]

SCALES = {
    "CoxPH-v2": 0.023934215222165295,
    "RSF": 0.01707069437342852,
    "DeepSurv": 0.012307440008774696,
}
COARSE_K = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
KAPPAS = (0.1, 0.25, 0.5)
FINE_K = np.geomspace(0.25, 128.0, 60)
COMPONENTS = ("F", "S", "U", "G_claim", "R", "N")


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def n_scores(item: dict, m0: float) -> np.ndarray:
    if "external" in item:
        sites = [item["phi"][0][:, :, item["targets"]], item["external"]]
    else:
        sites = [item["phi"][0], item["phi"][1]]
    leaves = [
        np.minimum(1.0, np.abs(site).mean(axis=1) / m0).mean(axis=0)
        for site in sites
    ]
    return np.minimum(leaves[0], leaves[1])


def main() -> None:
    started = time.perf_counter()
    phase12 = load("60_phase12_scale_tie_query_audit.py", "phase18_kprofile_phase12")
    data = phase12.load_inputs()
    phase77 = pd.read_csv(ROOT / "PHASE77_FULL_FEATURE_TABLE_FOLD0.csv")
    phase77 = phase77.set_index(["Model", "Feature ID"])
    profile_rows: list[dict[str, object]] = []
    entry_rows: list[dict[str, object]] = []

    for model, item in data.items():
        scale = SCALES[model]
        n_by_kappa = {kappa: n_scores(item, kappa * scale) for kappa in KAPPAS}
        for k in COARSE_K:
            scores = phase12.scores(item, k * scale, k * scale, 6 * k * scale, scale)
            for kappa in KAPPAS:
                n_value = n_by_kappa[kappa]
                for j, feature in enumerate(item["features"]):
                    values = np.r_[
                        scores.loc[feature, ["F", "S", "U", "G_claim", "R"]].to_numpy(float),
                        n_value[j],
                    ]
                    profile_rows.append(
                        {
                            "model": model,
                            "feature": feature,
                            "k": k,
                            "rho": 1.0,
                            "kappa": kappa,
                            **{name: values[q] for q, name in enumerate(COMPONENTS)},
                            "WI6": values.min(),
                            "point_adequate6": values.min() >= 0.8,
                            "bottleneck6": COMPONENTS[int(values.argmin())],
                        }
                    )

        n_primary = n_by_kappa[0.25]
        first_entry = {feature: np.nan for feature in item["features"]}
        entry_wi = {feature: np.nan for feature in item["features"]}
        for k in FINE_K:
            scores = phase12.scores(item, k * scale, k * scale, 6 * k * scale, scale)
            for j, feature in enumerate(item["features"]):
                if np.isfinite(first_entry[feature]):
                    continue
                values = np.r_[
                    scores.loc[feature, ["F", "S", "U", "G_claim", "R"]].to_numpy(float),
                    n_primary[j],
                ]
                if values.min() >= 0.8:
                    first_entry[feature] = float(k)
                    entry_wi[feature] = float(values.min())

        for feature in item["features"]:
            frozen = phase77.loc[(model, feature)]
            observed = first_entry[feature]
            entry_rows.append(
                {
                    "model": model,
                    "feature": feature,
                    "k_enter": observed,
                    "k_enter_label": "girmedi" if not np.isfinite(observed) else f"{observed:.12g}",
                    "right_censored_128": not np.isfinite(observed),
                    "rank_helper_k": observed if np.isfinite(observed) else 256.0,
                    "WI6_at_entry": entry_wi[feature],
                    "v1_WI": float(frozen["WI"]),
                    "v1_BH": str(frozen["BH certified"]).strip().lower() in {"evet", "true"},
                    "M_j": float(frozen["M_j mean absolute attribution"]),
                    "warrant_label": frozen["Informative/Degenerate label"],
                    "is_AJCC_Stage_I": "Stage I" in feature and "Stage II" not in feature and "Stage III" not in feature and "Stage IV" not in feature,
                }
            )

    profile = pd.DataFrame(profile_rows)
    entries = pd.DataFrame(entry_rows)
    profile.to_csv(OUT / "PHASE18_ADIM5_KPROFILE.csv", index=False)
    entries.to_csv(OUT / "PHASE18_ADIM5_KENTER_81.csv", index=False)

    print("coarse point counts")
    print(profile.groupby(["model", "kappa", "k"])["point_adequate6"].sum().to_string())
    for model, group in entries.groupby("model"):
        rho_wi = spearmanr(group["rank_helper_k"], group["v1_WI"]).statistic
        rho_m = spearmanr(group["rank_helper_k"], group["M_j"]).statistic
        observed = group.loc[~group["right_censored_128"], "k_enter"]
        print(model, "entered", len(observed), "median", observed.median(), "IQR", observed.quantile([.25, .75]).tolist(), "rho_WI", rho_wi, "rho_M", rho_m)
    pivot = entries.pivot(index="feature", columns="model", values="rank_helper_k")
    for left, right in (("CoxPH-v2", "RSF"), ("CoxPH-v2", "DeepSurv"), ("RSF", "DeepSurv")):
        print("Kendall", left, right, kendalltau(pivot[left], pivot[right]).statistic)
    print("elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
