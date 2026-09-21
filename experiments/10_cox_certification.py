"""Phase 5.2A: full-family CoxPH N=50 IUT plus BH/BY certification.

The component p-values use the approved Phase 5 null-centred bootstrap.
This script evaluates every Phase 3.2 assessable feature and never calls a
model fit or explainer.
"""

from __future__ import annotations

import argparse
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
OUT = ROOT / "results/phase52_cox_certification"
COMPONENTS = ("F", "S", "U", "G_claim", "R")


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_pilot_module():
    spec = importlib.util.spec_from_file_location("phase5_pilot", ROOT / "experiments/08_certification_pilot.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("Phase 5 pilot implementation cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(*, draws: int = 4999, seed: int = 20260917, alpha: float = 0.05) -> pd.DataFrame:
    if draws < 99:
        raise ValueError("At least 99 bootstrap draws required")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0,1)")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    stored = pd.read_csv(ROOT / "results/phase3_brca/feature_report.csv")
    assessed = stored.loc[np.isfinite(stored[list(COMPONENTS)]).all(axis=1)].copy().reset_index(drop=True)
    selection = assessed[["feature"]]
    if len(selection) != 77:
        raise RuntimeError(f"Expected the frozen m=77 CoxPH family, found {len(selection)}")

    pilot = load_pilot_module()
    rng = np.random.default_rng(seed + 5000)
    indices = rng.integers(0, 50, size=(draws, 50))
    point, sampled, query_inventory = pilot.cox_bootstrap(selection, indices)
    np.testing.assert_allclose(point, assessed[list(COMPONENTS)].to_numpy(), atol=2e-12)

    component_p = np.column_stack(
        [pilot.centered_bootstrap_pvalue(point[:, k], sampled[:, :, k], CALIBRATION["tau"]) for k in range(5)]
    )
    bh = certify(component_p, alpha=alpha, method="bh")
    by = certify(component_p, alpha=alpha, method="by")
    np.testing.assert_allclose(bh.joint_p, by.joint_p)

    result = assessed[["feature", "block", *COMPONENTS, "WI", "bottleneck"]].copy()
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

    summary = pd.DataFrame(
        [
            {"criterion": "assessed_full_family_m", "count": len(result)},
            {"criterion": "point_WI_ge_0_8", "count": int(result.point_WI_ge_0_8.sum())},
            {"criterion": "BH_certified_alpha_0_05", "count": int(result.certified_BH_alpha_0_05.sum())},
            {"criterion": "BY_certified_alpha_0_05", "count": int(result.certified_BY_alpha_0_05.sum())},
            {"criterion": "point_and_BH", "count": int((result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05).sum())},
            {"criterion": "point_only_not_BH", "count": int((result.point_WI_ge_0_8 & ~result.certified_BH_alpha_0_05).sum())},
            {"criterion": "BH_only_not_point", "count": int((~result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05).sum())},
        ]
    )
    summary.to_csv(OUT / "summary.csv", index=False)

    formula_paths = [ROOT / "warrantxai/components.py", ROOT / "warrantxai/aggregation.py"]
    formula_hashes = {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path) for path in formula_paths}
    pilot_meta = json.loads((ROOT / "results/phase5_certification_pilot/metadata.json").read_text())
    formula_unchanged = formula_hashes == pilot_meta["formula_module_sha256"]
    elapsed = perf_counter() - started
    metadata = {
        "phase": "5.2A_CoxPH_full_family_certification",
        "model": "CoxPH_N50",
        "arguments": {"bootstrap_draws": draws, "seed": seed, "alpha": alpha, "tau_components": dict.fromkeys(COMPONENTS, 0.8)},
        "family": {"definition": "all Phase 3.2 assessable CoxPH features; no WI preselection", "m": len(result)},
        "method": {
            "component_p": "null-centred nonparametric basic bootstrap",
            "resampling_unit": "paired patient-cluster training-replicate index t, with replacement",
            "p_value": "(1 + sum[(theta_star-theta_hat) >= (theta_hat-tau)]) / (B+1)",
            "joint_p": "max over F,S,U,G_claim,R (intersection-union test)",
            "primary_fdr": "BH",
            "robustness_fdr": "BY",
            "full_BH_BY_applied": True,
        },
        "query_rows": dict(zip(query_inventory.site, query_inventory.query_rows, strict=True)),
        "counts": dict(zip(summary.criterion, summary["count"], strict=True)),
        "formula_module_sha256": formula_hashes,
        "formula_module_sha256_unchanged_from_phase5_pilot": formula_unchanged,
        "privacy_check": {"identifier_named_columns": [], "patient_or_sample_identifiers_exported": False, "passed": True},
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not formula_unchanged:
        raise RuntimeError("Frozen formula module hash changed")
    if any(any(token in column.lower() for token in ("patient", "sample_id", "sample_key")) for column in result.columns):
        raise RuntimeError("Identifier-like output column found")

    passing = result.loc[result.certified_BH_alpha_0_05, ["feature", "block", "WI", "p_joint", "q_BH", "q_BY", "certified_BY_alpha_0_05"]]
    report = f"""# Faz 5.2A — CoxPH tam-aile istatistiksel sertifikasyonu

## Kapsam ve yöntem

Faz 3.2'de değerlendirilebilir bulunan **{len(result)} feature'ın tamamı** önceden tanımlı test ailesine alındı; `WI≥0.8` filtresiyle seçim yapılmadı. Her feature ve F/S/U/G_claim/R bileşeni için onaylanan N=50 replika-bazlı, null-merkezlenmiş bootstrap testi B={draws} ile çalıştırıldı. Bileşen p-değerleri `p_joint=max_k(p_jk)` intersection-union testiyle birleştirildi. BH q-değeri α={alpha:.2f} için birincil sonuç, keyfi feature bağımlılığına izin veren BY q-değeri sağlamlık kontrolüdür.

Çıkarım sabit 141 TCGA ve 1466 METABRIC sorgusuna, preprocessing'e ve exact Cox attribution tanımına koşulludur. Faz 5.1 kalibrasyonunda CoxPH N=50 için Holm-düzeltilmiş Tip-I ihlali görülmemiştir.

## Nokta tahmini ve sertifikasyon karşılaştırması

{markdown_table(summary)}

BH ile sertifikalanan feature'lar:

{markdown_table(passing) if len(passing) else '_Hiçbir feature BH ile sertifikalanmadı._'}

BH birincil sonuçtur; BY sütunu dependence-robust sağlamlık kontrolüdür. FDR-düzeltilmiş q-değeri bir feature'ın doğru olma olasılığı değildir. Tam feature-bileşen tablosu `results/phase52_cox_certification/feature_certification.csv` içindedir.

## Dürüstlük kontrolleri

Formül modülü SHA256 değerleri Faz 5 pilotuyla aynıdır: **{formula_unchanged}**. Çıktılarda hasta veya örnek kimliği yoktur. Model yeniden fit edilmedi ve explainer çağrılmadı. Gerçek wall-clock süre **{elapsed:.3f} saniye** idi.
"""
    (ROOT / "PHASE5_2_COX_CERTIFICATION_REPORT.md").write_text(report, encoding="utf-8")
    print(summary.to_string(index=False))
    print(f"Wrote CoxPH m={len(result)} BH/BY certification in {elapsed:.3f}s")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=4999)
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()
    run(draws=args.draws, seed=args.seed, alpha=args.alpha)
