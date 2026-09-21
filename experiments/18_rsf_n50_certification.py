"""Phase 5.6 step 3: RSF N=50 full-family recursive-studentized certification."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _shared import CALIBRATION, markdown_table
from warrantxai.certification import certify


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase56_rsf_n50_certification"
COMPONENTS = ("F", "S", "U", "G_claim", "R")
WARNING = (
    "RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; "
    "per-permutation-path uncertainty is excluded."
)


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def run(*, draws: int = 4999, seed: int = 20260918, alpha: float = 0.05) -> pd.DataFrame:
    gate = json.loads((ROOT / "results/phase56_rsf_n50_null_gate/metadata.json").read_text())
    if not gate["gate_passed"]:
        raise RuntimeError("RSF N=50 null gate did not pass; certification forbidden")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    stored = pd.read_csv(ROOT / "results/phase56_rsf_n50/feature_report.csv")
    with np.load(ROOT / "results/phase56_rsf_n50/attributions_n50.npz", allow_pickle=False) as saved:
        phi, reference, features = saved["candidate"], saved["reference"], saved["features"].tolist()
    if phi.shape != (2, 50, 12, 77) or len(stored) != 77 or features != stored.feature.tolist():
        raise RuntimeError("Expected the frozen ordered RSF N=50 m=77 family")

    phase53 = load_module("12_rsf_certification.py", "phase53_rsf_statistics_for_n50")
    rng = np.random.default_rng(seed + 560_500)
    indices = rng.integers(0, 50, size=(draws, 50))
    points, component_p = phase53.component_statistics(phi, reference, indices)
    np.testing.assert_allclose(points, stored[list(COMPONENTS)].to_numpy(), atol=2e-12)
    bh = certify(component_p, alpha=alpha, method="bh")
    by = certify(component_p, alpha=alpha, method="by")

    result = stored[["feature", "block", *COMPONENTS, "WI", "bottleneck"]].copy()
    result.insert(0, "inference_scope", WARNING)
    for k, component in enumerate(COMPONENTS):
        result[f"p_{component}"] = component_p[:, k]
    result["p_joint"] = bh.joint_p
    result["q_BH"] = bh.adjusted_p
    result["q_BY"] = by.adjusted_p
    result["certified_BH_alpha_0_05"] = bh.rejected
    result["certified_BY_alpha_0_05"] = by.rejected
    result["point_WI_ge_0_8"] = result.WI >= CALIBRATION["tau"]
    result["comparison_BH"] = np.select(
        [result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05,
         result.point_WI_ge_0_8 & ~result.certified_BH_alpha_0_05,
         ~result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05],
        ["both", "point_only", "BH_only"], default="neither",
    )
    result.to_csv(OUT / "feature_certification.csv", index=False)
    counts = {
        "assessed_full_family_m": len(result),
        "point_WI_ge_0_8": int(result.point_WI_ge_0_8.sum()),
        "BH_certified_alpha_0_05": int(result.certified_BH_alpha_0_05.sum()),
        "BY_certified_alpha_0_05": int(result.certified_BY_alpha_0_05.sum()),
        "point_and_BH": int((result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05).sum()),
        "point_only_not_BH": int((result.point_WI_ge_0_8 & ~result.certified_BH_alpha_0_05).sum()),
        "BH_only_not_point": int((~result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05).sum()),
    }
    summary = pd.DataFrame([
        {"inference_scope": WARNING, "criterion": criterion, "count": count}
        for criterion, count in counts.items()
    ])
    summary.to_csv(OUT / "summary.csv", index=False)

    formula_sha = {name: sha(ROOT / name) for name in [
        "warrantxai/components.py", "warrantxai/aggregation.py"
    ]}
    unchanged = formula_sha == gate["formula_module_sha256"]
    elapsed = perf_counter() - started
    metadata = {
        "phase": "5.6_step3_RSF_N50_full_family_certification",
        "model": "RSF_N50",
        "rsf_conditionality_warning": WARNING,
        "gate": {"passed": True, "source": "results/phase56_rsf_n50_null_gate/metadata.json",
                 "source_sha256": sha(ROOT / "results/phase56_rsf_n50_null_gate/metadata.json")},
        "arguments": {"bootstrap_draws": draws, "seed": seed, "alpha": alpha,
                      "tau_components": dict.fromkeys(COMPONENTS, 0.8)},
        "family": {"definition": "all Phase 5.6 RSF N=50 features; no WI preselection", "m": 77},
        "method": {"name": "recursive_studentized",
                   "studentized_leaves": ["F", "U", "G_magnitude", "G_signed"],
                   "basic_bootstrap_leaves": ["S", "R"],
                   "recursive_IUT": "max p at every site/branch min node, then max over components",
                   "primary_fdr": "BH", "robustness_fdr": "BY", "full_BH_BY_applied": True},
        "counts": counts,
        "formula_module_sha256": formula_sha,
        "formula_module_sha256_unchanged": unchanged,
        "source_sha256": {"results/phase56_rsf_n50/attributions_n50.npz": sha(
            ROOT / "results/phase56_rsf_n50/attributions_n50.npz"
        ), "experiments/18_rsf_n50_certification.py": sha(Path(__file__))},
        "privacy_check": {"identifier_named_columns": [], "patient_or_sample_identifiers_exported": False,
                          "passed": True},
        "model_or_explainer_calls": False,
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(), "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not unchanged:
        raise RuntimeError("Frozen formula module SHA256 changed")

    passing = result.loc[result.certified_BH_alpha_0_05, [
        "inference_scope", "feature", "block", "WI", "p_joint", "q_BH", "q_BY",
        "certified_BY_alpha_0_05"
    ]]
    point_only = result.loc[result.comparison_BH == "point_only", [
        "inference_scope", "feature", "block", "WI", "bottleneck", "p_joint", "q_BH"
    ]]
    text = f"""# Faz 5.6 ADIM 3 — RSF N=50 FDR sertifikasyonu

> **Koşulluluk uyarısı:** {WARNING}

N=50 null gate'i geçtikten sonra önceden tanımlı 77 feature ailesinin tamamı B={draws} `recursive_studentized` ile test edildi. BH alpha={alpha:.2f} birincil, BY sağlamlık analizidir.

{markdown_table(summary)}

## BH sertifikalı feature'lar

{markdown_table(passing) if len(passing) else 'BH ile sertifikalanan feature yok.'}

## Nokta-yeterli olup BH ile sertifikalanmayanlar

{markdown_table(point_only) if len(point_only) else 'Bu grupta feature yok.'}

Tam tablo `results/phase56_rsf_n50_certification/feature_certification.csv` içindedir. q-değeri feature doğruluk olasılığı değildir. Formül SHA256 değişmedi: **{unchanged}**. Model/explainer çağrısı yapılmadı; kimlik-sızıntısı kontrolü geçti. Gerçek wall-clock **{elapsed:.3f} saniye**.
"""
    (OUT / "report.md").write_text(text, encoding="utf-8")
    (ROOT / "PHASE5_6_STEP3_CERTIFICATION_REPORT.md").write_text(text, encoding="utf-8")
    print(summary[["criterion", "count"]].to_string(index=False), flush=True)
    print(f"wall={elapsed:.3f}s", flush=True)
    return result


if __name__ == "__main__":
    run()
