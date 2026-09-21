#!/usr/bin/env python3
"""Phase 7.4 tau and eta sensitivity from frozen Phase 7.2 attribution tensors."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from _shared import CALIBRATION, markdown_table
from warrantxai.certification import certify
from warrantxai.components import importance_ci, rank_reproducibility

TAUS = (0.70, 0.75, 0.80, 0.85, 0.90)
ETAS = (0.01, 0.02, 0.05, 0.10)
COMPONENTS = ("F", "S", "U", "G_claim", "R")
DRAWS, SEED = 4999, 20260918
EXPECTED_HASHES = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def mc_p(point: np.ndarray, sampled: np.ndarray, tau: float) -> np.ndarray:
    return (1 + (sampled - point >= point - tau).sum(axis=0)) / (len(sampled) + 1)


def studentized_p(
    point: np.ndarray, sampled: np.ndarray, se: np.ndarray, sampled_se: np.ndarray, tau: float
) -> np.ndarray:
    observed = (point - tau) / np.maximum(se, 1e-12)
    bootstrap_t = (sampled - point) / np.maximum(sampled_se, 1e-12)
    return (1 + (bootstrap_t >= observed).sum(axis=0)) / (len(sampled) + 1)


def attribution_statistics(
    phi: np.ndarray, reference: np.ndarray, indices: np.ndarray, tau: float, phase53, pilot
) -> tuple[np.ndarray, np.ndarray]:
    """Recursive-studentized statistics for RSF/DeepSurv at an arbitrary tau."""
    sites, t, _, p = phi.shape
    draws = len(indices)
    counts = np.zeros((draws, t), dtype=np.int16)
    np.add.at(counts, (np.arange(draws)[:, None], indices), 1)
    points, pvalues = np.empty((p, 5)), np.empty((p, 5))

    f_by_rep = 1 - (np.minimum(1.0, np.abs((phi - reference) / CALIBRATION["a"])) ** 2).mean(axis=2)
    lower, upper = importance_ci(phi, CALIBRATION["alpha_U"])
    u_by_rep = 1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"])
    for by_rep, column in ((f_by_rep, 0), (u_by_rep, 2)):
        point_leaf = by_rep.mean(axis=1)
        sampled_leaf, sampled_se = phase53.resampled_means(by_rep, counts)
        se = by_rep.std(axis=1, ddof=1) / np.sqrt(t)
        leaf_p = studentized_p(
            point_leaf.T.reshape(-1),
            np.moveaxis(sampled_leaf, 1, 2).reshape(draws, -1),
            se.T.reshape(-1),
            np.moveaxis(sampled_se, 1, 2).reshape(draws, -1),
            tau,
        ).reshape(p, sites)
        pvalues[:, column], points[:, column] = leaf_p.max(axis=1), point_leaf.min(axis=0)

    targets = np.arange(p)
    importance = np.abs(phi).mean(axis=2)
    s_loss = pilot.pair_losses(phi, targets, CALIBRATION["a"])
    r_loss = pilot.order_pair_losses(importance, targets, CALIBRATION["eta"])
    for loss, column in ((s_loss, 1), (r_loss, 4)):
        point_leaf = 1 - loss.sum(axis=(-2, -1)) / (t * (t - 1))
        sampled_leaf = phase53.quadratic_scores(loss, counts, t)
        leaf_p = np.column_stack([mc_p(point_leaf[site], sampled_leaf[:, site], tau) for site in range(sites)])
        pvalues[:, column], points[:, column] = leaf_p.max(axis=1), point_leaf.min(axis=0)

    difference = importance[1] - importance[0]
    gm_point = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / CALIBRATION["b"])
    gm_mean, gm_mean_se = phase53.resampled_means(difference[None, ...], counts)
    gm_sample = 1 - np.minimum(1.0, np.abs(gm_mean[:, 0]) / CALIBRATION["b"])
    gm_se = difference.std(axis=0, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    p_gm = studentized_p(gm_point, gm_sample, gm_se, gm_mean_se[:, 0] / CALIBRATION["b"], tau)

    full_weights = np.ones((1, t), dtype=np.int16)
    gs_point = phase53.signed_scores(phi, full_weights)[0]
    full_leave = np.ones((t, t), dtype=np.int16) - np.eye(t, dtype=np.int16)
    gs_jack = phase53.signed_scores(phi, full_leave)
    gs_centered = gs_jack - gs_jack.mean(axis=0, keepdims=True)
    gs_se = np.sqrt((t - 1) / t * np.sum(gs_centered**2, axis=0))
    gs_sample, gs_sample_se = phase53.signed_bootstrap_with_jackknife_se(phi, indices)
    p_gs = studentized_p(gs_point, gs_sample, gs_se, gs_sample_se, tau)
    pvalues[:, 3], points[:, 3] = np.maximum(p_gm, p_gs), np.minimum(gm_point, gs_point)
    return points, pvalues


def cox_statistics(indices: np.ndarray, tau: float, phase53, pilot, oldcox) -> tuple[np.ndarray, np.ndarray]:
    """Clean Cox recursive-studentized statistics at an arbitrary tau."""
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as saved:
        source, external = saved["source_phi"], saved["external_phi"]
        coefficients_all = saved["coefficients"]
        all_features = saved["features"].astype(str).tolist()
        active = saved["active_features"].astype(str).tolist()
    targets = np.array([all_features.index(name) for name in active])
    coefficients = coefficients_all[:, targets]
    selected_source = source[:, :, targets]
    t, p = coefficients.shape
    draws = len(indices)
    counts = np.zeros((draws, t), dtype=np.int16)
    np.add.at(counts, (np.arange(draws)[:, None], indices), 1)
    points, pvalues = np.empty((p, 5)), np.empty((p, 5))

    source4 = source[None, ...]
    lower, upper = importance_ci(source4, CALIBRATION["alpha_U"])
    u_by_rep = (1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"]))[:, :, targets]
    for by_rep, column in ((np.ones((1, t, p)), 0), (u_by_rep, 2)):
        point_leaf = by_rep.mean(axis=1)
        sampled, sampled_se = phase53.resampled_means(by_rep, counts)
        se = by_rep.std(axis=1, ddof=1) / np.sqrt(t)
        pvalues[:, column] = studentized_p(
            point_leaf.T.reshape(-1), np.moveaxis(sampled, 1, 2).reshape(draws, -1),
            se.T.reshape(-1), np.moveaxis(sampled_se, 1, 2).reshape(draws, -1), tau,
        )
        points[:, column] = point_leaf[0]

    s_loss = pilot.pair_losses(source4, targets, CALIBRATION["a"])
    r_loss = pilot.order_pair_losses(np.abs(source).mean(axis=1)[None, ...], targets, CALIBRATION["eta"])
    for loss, column in ((s_loss, 1), (r_loss, 4)):
        point_leaf = 1 - loss.sum(axis=(-2, -1)) / (t * (t - 1))
        sampled = phase53.quadratic_scores(loss, counts, t)
        pvalues[:, column], points[:, column] = mc_p(point_leaf[0], sampled[:, 0], tau), point_leaf[0]

    difference = np.abs(external).mean(axis=1) - np.abs(selected_source).mean(axis=1)
    gm_point = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / CALIBRATION["b"])
    gm_mean, gm_mean_se = phase53.resampled_means(difference[None, ...], counts)
    gm_sample = 1 - np.minimum(1.0, np.abs(gm_mean[:, 0]) / CALIBRATION["b"])
    gm_se = difference.std(axis=0, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    p_gm = studentized_p(gm_point, gm_sample, gm_se, gm_mean_se[:, 0] / CALIBRATION["b"], tau)

    sorted_delta = []
    for j in range(p):
        row = int(np.argmax(np.abs(coefficients[:, j])))
        coefficient = coefficients[row, j]
        if coefficient == 0:
            sorted_delta.append(np.zeros(max(selected_source.shape[1], external.shape[1])))
        else:
            z_source, z_external = selected_source[row, :, j] / coefficient, external[row, :, j] / coefficient
            grid_n = max(len(z_source), len(z_external))
            sorted_delta.append(np.sort(np.abs(
                pilot.empirical_rank_grid(z_source, grid_n) - pilot.empirical_rank_grid(z_external, grid_n)
            )))
    mean = coefficients.mean(axis=0)
    gs_point = oldcox.signed_score_from_coefficient(mean, sorted_delta)[0]
    delete_means = (t * mean - coefficients) / (t - 1)
    gs_se = oldcox.jackknife_se(oldcox.signed_score_from_coefficient(delete_means, sorted_delta))
    sampled_means = coefficients[indices].mean(axis=1)
    gs_sample = oldcox.signed_score_from_coefficient(sampled_means, sorted_delta)
    gs_sample_se = np.empty((draws, p))
    for j in range(p):
        leave_means = (t * sampled_means[:, j, None] - coefficients[indices, j]) / (t - 1)
        leave_scores = oldcox.signed_score_from_coefficient(
            leave_means.reshape(-1, 1), [sorted_delta[j]]
        ).reshape(draws, t)
        gs_sample_se[:, j] = oldcox.jackknife_se(leave_scores.T)
    p_gs = studentized_p(gs_point, gs_sample, gs_se, gs_sample_se, tau)
    points[:, 3], pvalues[:, 3] = np.minimum(gm_point, gs_point), np.maximum(p_gm, p_gs)
    return points, pvalues


def feasibility() -> dict:
    artifacts = {}
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as saved:
        artifacts["CoxPH"] = {key: list(saved[key].shape) for key in ("source_phi", "external_phi", "coefficients")}
    for model, path in [
        ("RSF", ROOT / "results/phase72_rsf/attributions_n50.npz"),
        ("DeepSurv", ROOT / "results/phase72_deepsurv/attributions_n50.npz"),
    ]:
        with np.load(path, allow_pickle=False) as saved:
            artifacts[model] = {key: list(saved[key].shape) for key in ("candidate", "reference")}
    if artifacts["CoxPH"]["source_phi"][:2] != [50, 153] or any(
        artifacts[model]["candidate"][:3] != [2, 50, 12] for model in ("RSF", "DeepSurv")
    ):
        raise RuntimeError(f"Stored Phase 7.2 attribution primitives are insufficient: {artifacts}")
    return {
        "feasible_without_patient_bootstrap_or_model_explainer": True,
        "stored_monte_carlo_4999_draw_matrix": False,
        "reason": "B=4999 inferential resamples are exactly regenerable from frozen T=50 attribution tensors and seed",
        "artifacts": artifacts,
    }


def tau_scan() -> pd.DataFrame:
    phase53 = load_module("12_rsf_certification.py", "phase74_phase53")
    pilot = load_module("08_certification_pilot.py", "phase74_pilot")
    oldcox = load_module("13_cox_recursive_certification.py", "phase74_oldcox")
    rng = np.random.default_rng(SEED + 725_000)
    indices = rng.integers(0, 50, size=(DRAWS, 50))
    stored = {
        "CoxPH": pd.read_csv(ROOT / "results/phase72_cox/feature_report.csv"),
        "RSF": pd.read_csv(ROOT / "results/phase72_rsf/feature_report.csv"),
        "DeepSurv": pd.read_csv(ROOT / "results/phase72_deepsurv/feature_report.csv"),
    }
    certified = {
        "CoxPH": pd.read_csv(ROOT / "results/phase72_coxph_certification/feature_certification.csv"),
        "RSF": pd.read_csv(ROOT / "results/phase72_rsf_certification/feature_certification.csv"),
        "DeepSurv": pd.read_csv(ROOT / "results/phase72_deepsurv_certification/feature_certification.csv"),
    }
    attribution = {}
    for model, path in [
        ("RSF", ROOT / "results/phase72_rsf/attributions_n50.npz"),
        ("DeepSurv", ROOT / "results/phase72_deepsurv/attributions_n50.npz"),
    ]:
        with np.load(path, allow_pickle=False) as saved:
            attribution[model] = (saved["candidate"], saved["reference"])
    rows = []
    for tau in TAUS:
        model_results = {}
        cpoint, cp = cox_statistics(indices, tau, phase53, pilot, oldcox)
        model_results["CoxPH"] = (cpoint, cp)
        for model in ("RSF", "DeepSurv"):
            model_results[model] = attribution_statistics(*attribution[model], indices, tau, phase53, pilot)
        row = {"tau": tau, "reference_primary": tau == 0.80}
        for model, (points, pvalues) in model_results.items():
            np.testing.assert_allclose(points, stored[model][list(COMPONENTS)].to_numpy(), atol=2e-12)
            wi = points.min(axis=1)
            bh, by = certify(pvalues, alpha=0.05, method="bh"), certify(pvalues, alpha=0.05, method="by")
            if tau == 0.80:
                np.testing.assert_allclose(
                    pvalues,
                    certified[model][[f"p_{component}" for component in COMPONENTS]].to_numpy(),
                    atol=0,
                    rtol=0,
                )
                np.testing.assert_array_equal(bh.rejected, certified[model].certified_BH_alpha_0_05)
                np.testing.assert_array_equal(by.rejected, certified[model].certified_BY_alpha_0_05)
            row[f"{model}_point"] = int((wi > tau).sum())
            row[f"{model}_BH"] = int(bh.rejected.sum())
            row[f"{model}_BY"] = int(by.rejected.sum())
        rows.append(row)
    result = pd.DataFrame(rows)
    reference = result.loc[result.tau == 0.8].iloc[0]
    expected = {"CoxPH": (9, 8, 8), "RSF": (17, 12, 12), "DeepSurv": (0, 0, 0)}
    for model, counts in expected.items():
        observed = tuple(int(reference[f"{model}_{key}"]) for key in ("point", "BH", "BY"))
        if observed != counts:
            raise RuntimeError(f"tau=0.8 did not reproduce Phase 7.2 for {model}: {observed} != {counts}")
    return result


def eta_scan() -> tuple[pd.DataFrame, pd.DataFrame]:
    reports = {
        "CoxPH": pd.read_csv(ROOT / "results/phase72_cox/feature_report.csv"),
        "RSF": pd.read_csv(ROOT / "results/phase72_rsf/feature_report.csv"),
        "DeepSurv": pd.read_csv(ROOT / "results/phase72_deepsurv/feature_report.csv"),
    }
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as saved:
        source = saved["source_phi"]
        all_features = saved["features"].astype(str).tolist()
        active = saved["active_features"].astype(str).tolist()
    targets = np.array([all_features.index(name) for name in active])
    importance = {"CoxPH": np.abs(source).mean(axis=1)[None, ...]}
    selected = {"CoxPH": targets}
    for model, path in [
        ("RSF", ROOT / "results/phase72_rsf/attributions_n50.npz"),
        ("DeepSurv", ROOT / "results/phase72_deepsurv/attributions_n50.npz"),
    ]:
        with np.load(path, allow_pickle=False) as saved:
            importance[model] = np.abs(saved["candidate"]).mean(axis=2)
        selected[model] = np.arange(27)

    rows, feature_rows = [], []
    for eta in ETAS:
        for model in ("CoxPH", "RSF", "DeepSurv"):
            r = rank_reproducibility(importance[model], eta)[selected[model]]
            base = reports[model].copy()
            scores = base[["F", "S", "U", "G_claim"]].copy()
            scores["R"] = r
            bottleneck = scores.idxmin(axis=1)
            wi = scores.min(axis=1)
            rows.append({
                "eta": eta, "reference_primary": eta == 0.02, "model": model,
                "R_mean": float(r.mean()), "R_std": float(r.std(ddof=1)),
                "R_min": float(r.min()), "R_q25": float(np.quantile(r, 0.25)),
                "R_median": float(np.median(r)), "R_q75": float(np.quantile(r, 0.75)),
                "R_max": float(r.max()), "R_bottleneck_count": int((bottleneck == "R").sum()),
                "point_adequate_tau_0_8": int((wi > 0.8).sum()),
            })
            for feature, value, bottleneck_value, wi_value in zip(base.feature, r, bottleneck, wi, strict=True):
                feature_rows.append({"eta": eta, "model": model, "feature": feature,
                                     "R": value, "bottleneck": bottleneck_value, "WI": wi_value})
    result, feature_result = pd.DataFrame(rows), pd.DataFrame(feature_rows)
    reference = result[result.eta == 0.02].set_index("model")
    expected = {"CoxPH": 24, "RSF": 24, "DeepSurv": 27}
    observed = reference.R_bottleneck_count.astype(int).to_dict()
    if observed != expected:
        raise RuntimeError(f"eta=0.02 did not reproduce Phase 7.2 bottlenecks: {observed} != {expected}")
    for model, report in reports.items():
        actual = feature_result[(feature_result.eta == 0.02) & (feature_result.model == model)].R
        np.testing.assert_allclose(actual, report.R, atol=2e-12)
    return result, feature_result


def main() -> None:
    total_started = perf_counter()
    hashes = {name: sha(ROOT / "warrantxai" / name) for name in EXPECTED_HASHES}
    if hashes != EXPECTED_HASHES:
        raise RuntimeError(f"Frozen formula hashes changed: {hashes}")
    stage = {}
    started = perf_counter()
    feasible = feasibility()
    stage["ADIM 0 — feasibility"] = perf_counter() - started

    started = perf_counter()
    tau = tau_scan()
    stage["ADIM 1 — tau scan"] = perf_counter() - started
    tau.to_csv(ROOT / "PHASE74_TAU_SENSITIVITY.csv", index=False)
    deep70 = tau.loc[tau.tau == 0.70, ["DeepSurv_point", "DeepSurv_BH", "DeepSurv_BY"]].iloc[0]
    tau_note = f"""# Faz 7.4 — τ duyarlılık taraması

Ana analiz **τ=0.80** olarak kalır. Aynı Faz 7.2 T=50 attribution tensörleri ve aynı B={DRAWS} inferential resample akışı kullanılarak yalnız tek-yönlü null sınırı değiştirildi; yeni patient bootstrap, model fit veya explainer çağrısı yapılmadı. Point kararı istenen katı `WI>τ` kuralıdır.

{markdown_table(tau)}

DeepSurv τ=0.70 satırında point/BH/BY sayıları **{int(deep70.DeepSurv_point)}/{int(deep70.DeepSurv_BH)}/{int(deep70.DeepSurv_BY)}**. Bu, τ=0.80'deki 0/0/0 sonucunun düşük eşikte nasıl değiştiğini doğrudan gösterir.

DeepSurv τ=0.75'te 6/1/1'e, τ=0.70'te 21/6/6'ya çıktığı için ana τ=0.80'deki 0/0/0 sonucu aşağı yönlü eşik değişimine duyarlıdır. CoxPH ve RSF sayıları da τ yükseldikçe genel olarak azalır; ana sonuçlar değiştirilmemiştir.
"""
    (ROOT / "PHASE74_TAU_SENSITIVITY.md").write_text(tau_note, encoding="utf-8")

    started = perf_counter()
    eta, eta_features = eta_scan()
    stage["ADIM 2 — eta scan"] = perf_counter() - started
    eta.to_csv(ROOT / "PHASE74_ETA_SENSITIVITY.csv", index=False)
    eta_features.to_csv(ROOT / "PHASE74_ETA_FEATURE_VALUES.csv", index=False)
    eta10 = eta[eta.eta == 0.10].set_index("model")
    eta_note = f"""# Faz 7.4 — η rank-bottleneck duyarlılığı

Ana analiz **η=0.02** olarak kalır. R, saklı T=50 candidate attribution tensörlerinden η∈{{0.01,0.02,0.05,0.10}} için yeniden değerlendirildi; diğer dört bileşen değiştirilmedi.

{markdown_table(eta)}

η=0.10'da R-darboğaz sayıları CoxPH/RSF/DeepSurv için **{int(eta10.loc['CoxPH','R_bottleneck_count'])}/{int(eta10.loc['RSF','R_bottleneck_count'])}/{int(eta10.loc['DeepSurv','R_bottleneck_count'])}**. Tam feature düzeyi R/WI/bottleneck değerleri `PHASE74_ETA_FEATURE_VALUES.csv` içindedir.

R darboğaz baskınlığı η=0.10'da RSF ve DeepSurv için sürerken (25/27 ve 25/27), CoxPH'de 4/27'ye düşer. Bu nedenle η-seçimine dayanıklılık iki nonlinear modelde güçlü, CoxPH'de ise geniş tolerans değerinde zayıftır.
"""
    (ROOT / "PHASE74_ETA_SENSITIVITY.md").write_text(eta_note, encoding="utf-8")

    total = perf_counter() - total_started
    timing = pd.DataFrame([{"Step": key, "Seconds": value} for key, value in stage.items()])
    feasibility_table = pd.DataFrame([
        {"Model": model, **{key: str(value) for key, value in values.items()}}
        for model, values in feasible["artifacts"].items()
    ]).fillna("-")
    summary = f"""# Faz 7.4 — Eşik ve rank-tolerance duyarlılık özeti

## ADIM 0 — Fizibilite

**Yeniden patient-cluster bootstrap gerekmedi.** B=4999 draw matrisi ayrıca saklanmamış olsa da aynı inferential resample akışı, aşağıdaki dondurulmuş T=50 attribution tensörlerinden ve seed'den birebir yeniden üretildi. Model veya explainer çağrısı yapılmadı.

{markdown_table(feasibility_table)}

## Supplementary tablo S1 — τ taraması

Referans satır `reference_primary=True` olan **τ=0.80** satırıdır.

{markdown_table(tau)}

## Supplementary tablo S2 — η taraması

Referans satırlar `reference_primary=True` olan **η=0.02** satırlarıdır.

{markdown_table(eta)}

η=0.10'da R-darboğaz sayıları CoxPH/RSF/DeepSurv için **4/25/25** oldu. R baskınlığı RSF ve DeepSurv'de geniş η değerinde sürerken CoxPH'de belirgin biçimde azaldı. DeepSurv τ=0.70'te 21/6/6, τ=0.75'te 6/1/1 ve referans τ=0.80'de 0/0/0 olduğundan DeepSurv sertifikasyon sonucu aşağı yönlü τ değişimine duyarlıdır.

## Süre

{markdown_table(timing)}

Toplam wall-clock **{total:.3f} saniye**. Ana Faz 7.2 sonuçları ve artefaktları değiştirilmedi; formül SHA256 değerleri doğrulandı.
"""
    (ROOT / "PHASE74_SUMMARY.md").write_text(summary, encoding="utf-8")
    metadata = {
        "phase": "7.4_tau_eta_sensitivity",
        "feasibility": feasible,
        "tau_values": TAUS,
        "eta_values": ETAS,
        "bootstrap_training_replicates_reused": 50,
        "inferential_resample_draws": DRAWS,
        "seed": SEED,
        "model_or_explainer_calls": False,
        "new_patient_cluster_bootstrap_fits": False,
        "formula_sha256": hashes,
        "stage_seconds": stage,
        "wall_clock_seconds": total,
    }
    (ROOT / "results/phase74_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(tau.to_string(index=False), flush=True)
    print(eta.to_string(index=False), flush=True)
    print(f"wall={total:.3f}s", flush=True)


if __name__ == "__main__":
    main()
