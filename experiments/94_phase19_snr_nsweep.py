"""Phase 19 step 2: coefficient SNR and sample-size sweep."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from phase19_common import BETA0, BETA1, CLASSES, simulate_light
from warrantxai.components import pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase19"
OUT.mkdir(parents=True, exist_ok=True)
SEEDS = tuple(range(19000001, 19000021))
NS = (700, 2800, 11200)
KS = (1.0, 2.0)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def rank_scores(importance: np.ndarray, eta: float) -> np.ndarray:
    sites, t, p = importance.shape
    orders = pairwise_orders(importance, eta).astype(float)
    out = np.empty((sites, p))
    pairs = t * (t - 1) / 2
    for site in range(sites):
        for j in range(p):
            order = orders[site, :, j, :]
            loss = 0.0
            for first in range(t):
                loss += np.abs(order[first + 1:] - order[first]).sum() / (2 * (p - 1))
            out[site, j] = 1.0 - loss / pairs
    return out.min(axis=0)


def a_components(sim: dict, k: float) -> pd.DataFrame:
    coef = np.asarray(sim["coefficients"])
    xs = [np.asarray(sim["x0"]), np.asarray(sim["x1"])]
    centered = [x - x.mean(axis=0) for x in xs]
    abs_mean = np.stack([np.abs(x).mean(axis=0) for x in centered])
    abs_sd = np.stack([np.abs(x).std(axis=0, ddof=1) for x in centered])
    importance = np.abs(coef) * abs_mean[:, None, :]
    magnitude = importance[0].mean(axis=0)
    scale = float(np.median(magnitude[magnitude > 0]))
    a = b = eta = k * scale
    u = 6.0 * k * scale

    half = norm.ppf(.975) * np.abs(coef) * abs_sd[:, None, :] / np.sqrt(len(xs[0]))
    lower = np.maximum(0.0, importance - half)
    upper = importance + half
    u_point = (1.0 - np.minimum(1.0, (upper - lower) / u)).mean(axis=1).min(axis=0)
    r_point = rank_scores(importance, eta)
    g_mag = 1.0 - np.minimum(1.0, np.abs(importance[1].mean(axis=0) - importance[0].mean(axis=0)) / b)
    signed_means = [coef[d].mean(axis=0) * centered[d] for d in range(2)]
    ordered = [np.sort(value, axis=0) for value in signed_means]
    g_signed = 1.0 - np.minimum(1.0, np.abs(ordered[1] - ordered[0]) / a).mean(axis=0)
    g_claim = np.minimum(g_mag, g_signed)
    n_point = np.minimum(1.0, importance / (.25 * scale)).mean(axis=1).min(axis=0)

    # S only for class A; compute the exact frozen ordered-pair mean without phi.
    s_a = np.empty(8)
    t = coef.shape[1]
    for j in range(8):
        by_site = []
        absx = [np.abs(centered[d][:, j]) for d in range(2)]
        for d in range(2):
            total = 0.0
            for first in range(t):
                dc = np.abs(coef[d, first, j] - coef[d, :, j])[:, None]
                total += (np.minimum(1.0, dc * absx[d][None, :] / a) ** 2).mean(axis=1).sum()
            by_site.append(1.0 - total / (t * (t - 1)))
        s_a[j] = min(by_site)
    values = np.column_stack([np.ones(8), s_a, u_point[:8], g_claim[:8], r_point[:8], n_point[:8]])
    return pd.DataFrame({
        "feature": np.arange(8), "truth_class": "A", "k": k, "n": len(xs[0]),
        "s_m": scale, "F": values[:, 0], "S": values[:, 1], "U": values[:, 2],
        "G_claim": values[:, 3], "R": values[:, 4], "N": values[:, 5],
        "G_claim_pass": values[:, 3] >= .8, "WI6": values.min(axis=1),
        "WI6_pass": values.min(axis=1) >= .8,
    })


def snr_rows(seed: int, sim: dict) -> pd.DataFrame:
    coef = np.asarray(sim["coefficients"])
    rows = []
    for site, beta in enumerate((BETA0, BETA1)):
        y = np.asarray(sim[f"y{site}"])
        x = np.asarray(sim[f"x{site}"])
        residual = y - x @ beta
        r2 = 1.0 - float(np.sum(residual**2) / np.sum((y - y.mean())**2))
        for j in range(30):
            mean = float(coef[site, :, j].mean())
            sd = float(coef[site, :, j].std(ddof=1))
            rows.append({"seed": seed, "site": site, "feature": j, "truth_class": CLASSES[j],
                         "beta_true": beta[j], "beta_hat_mean": mean, "beta_hat_sd": sd,
                         "beta_hat_over_se": mean / sd if sd else np.nan,
                         "source_external_difference_sd": float((coef[1, :, j] - coef[0, :, j]).std(ddof=1)),
                         "noise_sd": sim["noise_sd"], "outcome_variance": float(y.var(ddof=1)),
                         "r2_true_dgp": r2})
    return pd.DataFrame(rows)


def main() -> None:
    started = time.perf_counter()
    snr_cp = OUT / "PHASE19_ADIM2_SNR_CHECKPOINT.csv"
    snr = pd.read_csv(snr_cp) if snr_cp.exists() else pd.DataFrame()
    done_snr = set(snr.seed.unique()) if not snr.empty else set()
    for seed in SEEDS:
        if seed not in done_snr:
            sim = simulate_light(seed, n=700)
            snr = pd.concat([snr, snr_rows(seed, sim)], ignore_index=True)
            snr.to_csv(snr_cp, index=False)
            print(f"P19 step2 SNR seed {seed} complete", flush=True)
    assert len(snr) == 20 * 2 * 30
    snr.to_csv(OUT / "PHASE19_ADIM2_SNR.csv", index=False)

    sweep_cp = OUT / "PHASE19_ADIM2_NSWEEP_CHECKPOINT.csv"
    sweep = pd.read_csv(sweep_cp) if sweep_cp.exists() else pd.DataFrame()
    done = set(zip(sweep.seed, sweep.n)) if not sweep.empty else set()
    for n in NS:
        for seed in SEEDS:
            if (seed, n) in done:
                continue
            sim = simulate_light(seed, n=n)
            cells = []
            for k in KS:
                cell = a_components(sim, k)
                cell.insert(0, "seed", seed)
                cells.append(cell)
            sweep = pd.concat([sweep, *cells], ignore_index=True)
            sweep.to_csv(sweep_cp, index=False)
            print(f"P19 step2 n={n} seed={seed} complete", flush=True)
    assert len(sweep) == 20 * 3 * 2 * 8
    sweep.to_csv(OUT / "PHASE19_ADIM2_NSWEEP.csv", index=False)

    summary = sweep.groupby(["n", "k"], as_index=False).agg(
        median_G_claim=("G_claim", "median"), G_claim_pass_rate=("G_claim_pass", "mean"),
        median_WI6=("WI6", "median"), WI6_point_pass_rate=("WI6_pass", "mean"))
    summary.to_csv(OUT / "PHASE19_ADIM2_NSWEEP_SUMMARY.csv", index=False)
    k1 = summary[summary.k == 1.0].sort_values("n")
    rates = k1.G_claim_pass_rate.to_numpy()
    h2 = bool(np.all(np.diff(rates) >= 0) and rates[-1] >= .95)
    meta = {
        "rows": {"SNR": len(snr), "NSWEEP": len(sweep)},
        "expected_rows": {"SNR": 1200, "NSWEEP": 960},
        "H2": {"k1_G_claim_pass_rates": dict(zip(k1.n.astype(str), rates)), "supported": h2},
        "dgp_code": "y_d = X_d beta_d + Normal(0, 0.1); R2=1-SSE_true/SST",
        "formula_sha256": {name: sha(ROOT / "warrantxai" / name) for name in ["components.py", "aggregation.py", "certification.py"]},
        "elapsed_seconds": time.perf_counter() - started,
    }
    (OUT / "PHASE19_ADIM2_METADATA.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
