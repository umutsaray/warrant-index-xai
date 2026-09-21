#!/usr/bin/env python3
"""Phase 7.2 full-family BH/BY certification after all fresh gates pass."""

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
COMPONENTS = ("F", "S", "U", "G_claim", "R")
SEED, DRAWS = 20260918, 4999

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from _shared import CALIBRATION, markdown_table
from warrantxai.certification import certify
from warrantxai.components import importance_ci


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def cox_statistics(indices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact clean Cox leaves with recursive studentization and full 610-feature R universe."""
    phase53 = load_module("12_rsf_certification.py", "phase72_cox_phase53")
    oldcox = load_module("13_cox_recursive_certification.py", "phase72_cox_old")
    pilot = load_module("08_certification_pilot.py", "phase72_cox_pilot_cert")
    remedy = load_module("11_rsf_remedy_calibration.py", "phase72_cox_remedy_cert")
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as saved:
        source = saved["source_phi"]
        external = saved["external_phi"]
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
    points = np.empty((p, 5))
    pvalues = np.empty_like(points)

    f_by_rep = np.ones((1, t, p))
    source4 = source[None, ...]
    lower, upper = importance_ci(source4, CALIBRATION["alpha_U"])
    u_by_rep = (1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"]))[:, :, targets]
    for by_rep, column in ((f_by_rep, 0), (u_by_rep, 2)):
        point_leaf = by_rep.mean(axis=1)
        sampled, sampled_se = phase53.resampled_means(by_rep, counts)
        se = by_rep.std(axis=1, ddof=1) / np.sqrt(t)
        pvalues[:, column] = remedy.studentized_p(
            point_leaf.T.reshape(-1),
            np.moveaxis(sampled, 1, 2).reshape(draws, -1),
            se.T.reshape(-1),
            np.moveaxis(sampled_se, 1, 2).reshape(draws, -1),
        )
        points[:, column] = point_leaf[0]

    s_loss = pilot.pair_losses(source4, targets, CALIBRATION["a"])
    r_loss = pilot.order_pair_losses(np.abs(source).mean(axis=1)[None, ...], targets, CALIBRATION["eta"])
    for loss, column in ((s_loss, 1), (r_loss, 4)):
        point_leaf = 1 - loss.sum(axis=(-2, -1)) / (t * (t - 1))
        sampled = phase53.quadratic_scores(loss, counts, t)
        pvalues[:, column] = remedy.mc_p(point_leaf[0], sampled[:, 0])
        points[:, column] = point_leaf[0]

    difference = np.abs(external).mean(axis=1) - np.abs(selected_source).mean(axis=1)
    gm_point = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / CALIBRATION["b"])
    gm_mean, gm_mean_se = phase53.resampled_means(difference[None, ...], counts)
    gm_sample = 1 - np.minimum(1.0, np.abs(gm_mean[:, 0]) / CALIBRATION["b"])
    gm_se = difference.std(axis=0, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    p_gm = remedy.studentized_p(gm_point, gm_sample, gm_se, gm_mean_se[:, 0] / CALIBRATION["b"])

    sorted_delta = []
    for j in range(p):
        row = int(np.argmax(np.abs(coefficients[:, j])))
        coefficient = coefficients[row, j]
        if coefficient == 0:
            sorted_delta.append(np.zeros(max(selected_source.shape[1], external.shape[1])))
            continue
        z_source = selected_source[row, :, j] / coefficient
        z_external = external[row, :, j] / coefficient
        grid_n = max(len(z_source), len(z_external))
        source_grid = pilot.empirical_rank_grid(z_source, grid_n)
        external_grid = pilot.empirical_rank_grid(z_external, grid_n)
        sorted_delta.append(np.sort(np.abs(source_grid - external_grid)))
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
    p_gs = remedy.studentized_p(gs_point, gs_sample, gs_se, gs_sample_se)
    points[:, 3] = np.minimum(gm_point, gs_point)
    pvalues[:, 3] = np.maximum(p_gm, p_gs)
    return points, pvalues


def write_one(model: str, source_dir: str, points: np.ndarray, component_p: np.ndarray, scope: str) -> dict:
    stored = pd.read_csv(ROOT / "results" / source_dir / "feature_report.csv")
    if "block" not in stored:
        stored.insert(1, "block", stored.feature.str.split("_").str[0])
    np.testing.assert_allclose(points, stored[list(COMPONENTS)].to_numpy(), atol=2e-12)
    bh, by = certify(component_p, alpha=0.05, method="bh"), certify(component_p, alpha=0.05, method="by")
    result = stored[["feature", "block", *COMPONENTS, "WI", "bottleneck"]].copy()
    result.insert(0, "inference_scope", scope)
    for k, component in enumerate(COMPONENTS):
        result[f"p_{component}"] = component_p[:, k]
    result["p_joint"] = bh.joint_p
    result["q_BH"] = bh.adjusted_p
    result["q_BY"] = by.adjusted_p
    result["certified_BH_alpha_0_05"] = bh.rejected
    result["certified_BY_alpha_0_05"] = by.rejected
    result["point_WI_ge_0_8"] = result.WI >= CALIBRATION["tau"]
    out = ROOT / "results" / f"phase72_{model.lower()}_certification"
    out.mkdir(parents=True, exist_ok=True)
    result.to_csv(out / "feature_certification.csv", index=False)
    counts = {
        "family_m": len(result),
        "point": int(result.point_WI_ge_0_8.sum()),
        "BH": int(result.certified_BH_alpha_0_05.sum()),
        "BY": int(result.certified_BY_alpha_0_05.sum()),
    }
    metadata = {
        "phase": f"7.2_{model}_full_family_certification",
        "model": model,
        "inference_scope": scope,
        "method": "recursive_studentized",
        "arguments": {"bootstrap_draws": DRAWS, "seed": SEED, "alpha": 0.05},
        "counts": counts,
        "gate_source": "results/phase72_calibration/metadata.json",
        "formula_sha256": json.loads((ROOT / "results/phase72_panel/metadata.json").read_text())["formula_sha256"],
        "privacy": {"identifiers_exported": False},
        "model_or_explainer_calls": False,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    pd.DataFrame([{"criterion": key, "count": value, "inference_scope": scope} for key, value in counts.items()]).to_csv(out / "summary.csv", index=False)
    return {"result": result, "counts": counts, "out": out}


def append_certification(report_path: Path, payload: dict) -> None:
    certified = payload["result"].loc[payload["result"].certified_BH_alpha_0_05,
                                      ["feature", "WI", "p_joint", "q_BH", "q_BY"]]
    text = report_path.read_text(encoding="utf-8")
    text += f"""

## ADIM 5 — Tam-aile sertifikasyon

Taze null kapısı geçtikten sonra m={payload['counts']['family_m']} ailenin tamamında B={DRAWS} recursive_studentized uygulandı. Nokta/BH/BY sayıları **{payload['counts']['point']}/{payload['counts']['BH']}/{payload['counts']['BY']}**.

{markdown_table(certified) if len(certified) else 'BH ile sertifikalanan feature yok.'}
"""
    report_path.write_text(text, encoding="utf-8")


def main() -> None:
    started = perf_counter()
    gate = json.loads((ROOT / "results/phase72_calibration/metadata.json").read_text())
    if not gate["all_models_gate_passed"]:
        raise RuntimeError("At least one Phase 7.2 null gate failed; certification is forbidden")
    rng = np.random.default_rng(SEED + 725_000)
    indices = rng.integers(0, 50, size=(DRAWS, 50))

    cox_points, cox_p = cox_statistics(indices)
    cox = write_one(
        "coxph", "phase72_cox", cox_points, cox_p,
        "153 TCGA + 1466 METABRIC fixed queries; exact closed-form Cox attribution",
    )
    phase53 = load_module("12_rsf_certification.py", "phase72_shared_cert")
    other = {}
    for model, source, scope in [
        ("rsf", "phase72_rsf", "fixed 12+12 queries; stored 8/128 averages; permutation-path uncertainty excluded"),
        ("deepsurv", "phase72_deepsurv", "fixed 12+12 queries; nested DeepExplainer k16/k128 backgrounds"),
    ]:
        with np.load(ROOT / "results" / source / "attributions_n50.npz", allow_pickle=False) as saved:
            phi, reference = saved["candidate"], saved["reference"]
        points, pvalues = phase53.component_statistics(phi, reference, indices)
        other[model] = write_one(model, source, points, pvalues, scope)

    append_certification(ROOT / "PHASE72_COXPH.md", cox)
    append_certification(ROOT / "PHASE72_RSF.md", other["rsf"])
    append_certification(ROOT / "PHASE72_DEEPSURV.md", other["deepsurv"])
    formula = {
        "components.py": sha(ROOT / "warrantxai/components.py"),
        "aggregation.py": sha(ROOT / "warrantxai/aggregation.py"),
    }
    expected = json.loads((ROOT / "results/phase72_panel/metadata.json").read_text())["formula_sha256"]
    if formula != expected:
        raise RuntimeError("Frozen formula SHA256 changed")
    elapsed = perf_counter() - started
    counts = {"CoxPH": cox["counts"], "RSF": other["rsf"]["counts"], "DeepSurv": other["deepsurv"]["counts"]}
    (ROOT / "results/phase72_calibration/certification_summary.json").write_text(
        json.dumps({"counts": counts, "wall_clock_seconds": elapsed, "formula_sha256": formula}, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(counts, indent=2), flush=True)


if __name__ == "__main__":
    main()
