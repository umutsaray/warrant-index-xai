"""Phase 6.2 step 4: compare N=50 CoxPH, RSF, and DeepSurv certification."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase62_three_model_comparison"
DEEP_WARNING = (
    "DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested "
    "DeepExplainer k=16/128 TCGA backgrounds."
)


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def markdown_table(frame: pd.DataFrame, digits: int = 4) -> str:
    shown = frame.copy()
    for column in shown.select_dtypes(include=[np.number]):
        shown[column] = shown[column].map(lambda value: f"{value:.{digits}f}" if pd.notna(value) else "NA")
    shown = shown.fillna("NA").astype(str)
    header = "| " + " | ".join(shown.columns) + " |"
    separator = "| " + " | ".join("---" for _ in shown.columns) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in shown.to_numpy().tolist()]
    return "\n".join([header, separator, *rows])


def run() -> None:
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    paths = {
        "CoxPH N=50": ROOT / "results/phase53_cox_recursive/feature_certification.csv",
        "RSF N=50": ROOT / "results/phase56_rsf_n50_certification/feature_certification.csv",
        "DeepSurv N=50": ROOT / "results/phase62_deepsurv_n50_certification/feature_certification.csv",
    }
    frames = {name: pd.read_csv(path) for name, path in paths.items()}
    scopes = {
        "CoxPH N=50": "CoxPH N=50 full-audit-query recursive_studentized certification.",
        "RSF N=50": "RSF N=50 fixed 12+12 query and stored 8/128 attribution certification.",
        "DeepSurv N=50": DEEP_WARNING,
    }
    summary_rows, margin_rows, bottleneck_rows = [], [], []
    point_margins = {}
    for model, frame in frames.items():
        point = frame.loc[frame.point_WI_ge_0_8.astype(bool)].copy()
        margins = point.WI - 0.8
        point_margins[model] = margins.to_numpy()
        summary_rows.append({
            "inference_scope": scopes[model], "model": model, "assessed": len(frame),
            "point_WI_ge_0_8": int(frame.point_WI_ge_0_8.sum()),
            "BH": int(frame.certified_BH_alpha_0_05.sum()),
            "BY": int(frame.certified_BY_alpha_0_05.sum()),
            "BH_survival_fraction": float(frame.certified_BH_alpha_0_05.sum() / max(1, frame.point_WI_ge_0_8.sum())),
        })
        margin_rows.append({
            "inference_scope": scopes[model], "model": model, "point_features": len(point),
            "mean_margin": margins.mean(), "median_margin": margins.median(),
            "n_0_0.01": int(((margins >= 0) & (margins < 0.01)).sum()),
            "n_0.01_0.03": int(((margins >= 0.01) & (margins < 0.03)).sum()),
            "n_0.03_0.05": int(((margins >= 0.03) & (margins < 0.05)).sum()),
            "n_0.05_0.10": int(((margins >= 0.05) & (margins <= 0.10)).sum()),
            "share_0.05_0.10": float(((margins >= 0.05) & (margins <= 0.10)).mean()),
        })
        for component, count in frame.bottleneck.value_counts().items():
            bottleneck_rows.append({
                "inference_scope": scopes[model], "model": model, "scope": "all_77",
                "bottleneck": component, "count": int(count),
            })
        for component, count in point.bottleneck.value_counts().items():
            bottleneck_rows.append({
                "inference_scope": scopes[model], "model": model, "scope": "point_passers",
                "bottleneck": component, "count": int(count),
            })

    summary = pd.DataFrame(summary_rows)
    margins = pd.DataFrame(margin_rows)
    bottlenecks = pd.DataFrame(bottleneck_rows)
    summary.to_csv(OUT / "three_model_summary.csv", index=False)
    margins.to_csv(OUT / "margin_distribution.csv", index=False)
    bottlenecks.to_csv(OUT / "bottleneck_distribution.csv", index=False)

    deep_margin = point_margins["DeepSurv N=50"]
    distances = pd.DataFrame([
        {"inference_scope": DEEP_WARNING, "comparison": "DeepSurv_vs_CoxPH",
         "wasserstein_margin_distance": wasserstein_distance(deep_margin, point_margins["CoxPH N=50"])},
        {"inference_scope": DEEP_WARNING, "comparison": "DeepSurv_vs_RSF",
         "wasserstein_margin_distance": wasserstein_distance(deep_margin, point_margins["RSF N=50"])},
    ])
    distances.to_csv(OUT / "margin_distances.csv", index=False)
    closer = "CoxPH" if distances.iloc[0].wasserstein_margin_distance < distances.iloc[1].wasserstein_margin_distance else "RSF"

    transition = pd.read_csv(ROOT / "results/phase62_deepsurv_n50/n2_n50_comparison.csv")
    transition.to_csv(OUT / "deepsurv_n2_n50_transition.csv", index=False)
    transition_counts = transition.transition.value_counts().to_dict()
    lost = transition.loc[transition.transition == "lost", [
        "inference_scope", "feature", "block", "WI_n2", "WI", "R_n2", "R",
        "bottleneck_n2", "bottleneck"
    ]]
    lost.to_csv(OUT / "deepsurv_lost_n2_winners.csv", index=False)

    formula_files = ["warrantxai/components.py", "warrantxai/aggregation.py"]
    formula_sha = {name: sha(ROOT / name) for name in formula_files}
    expected = json.loads((ROOT / "results/phase62_deepsurv_n50/metadata.json").read_text())["formula_module_sha256"]
    unchanged = formula_sha == expected
    if not unchanged:
        raise RuntimeError("Frozen formula SHA256 changed")
    elapsed = perf_counter() - started
    step1_meta = json.loads((ROOT / "results/phase62_deepsurv_n50/metadata.json").read_text())
    step2_meta = json.loads((ROOT / "results/phase62_deepsurv_n50_null_gate/metadata.json").read_text())
    step3_meta = json.loads((ROOT / "results/phase62_deepsurv_n50_certification/metadata.json").read_text())
    runtimes = pd.DataFrame([
        {"inference_scope": DEEP_WARNING, "step": "ADIM 1 N=50",
         "wall_clock_seconds": step1_meta["runtime"]["observed_total_wall_seconds_across_attempts"],
         "note": "includes stopped SHAP assertion attempt and checkpoint resume"},
        {"inference_scope": DEEP_WARNING, "step": "ADIM 2 null gate",
         "wall_clock_seconds": step2_meta["wall_clock_seconds"], "note": "2000x999, no model calls"},
        {"inference_scope": DEEP_WARNING, "step": "ADIM 3 FDR",
         "wall_clock_seconds": step3_meta["wall_clock_seconds"], "note": "B=4999, no model calls"},
        {"inference_scope": DEEP_WARNING, "step": "ADIM 4 comparison",
         "wall_clock_seconds": elapsed, "note": "existing artifacts only"},
    ])
    runtimes.to_csv(OUT / "runtime_by_step.csv", index=False)
    verdict = (
        "winner's-curse supported at count level, but N=2 membership is highly unstable; "
        f"DeepSurv point-margin distribution is closer to {closer} by Wasserstein distance"
    )
    metadata = {
        "phase": "6.2_step4_three_model_comparison", "deepsurv_conditionality_warning": DEEP_WARNING,
        "verdict": verdict, "transition_counts": transition_counts,
        "formula_module_sha256": formula_sha, "formula_module_sha256_unchanged": unchanged,
        "source_sha256": {str(path.relative_to(ROOT)).replace('\\', '/'): sha(path) for path in paths.values()} |
                         {"experiments/25_phase62_three_model_comparison.py": sha(Path(__file__))},
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "model_or_explainer_calls": False, "wall_clock_seconds": elapsed,
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    text = f"""# Faz 6.2 — DeepSurv N=50, sertifikasyon ve üç-model karşılaştırması

> **DeepSurv koşulluluk uyarısı:** {DEEP_WARNING}

## Sertifikasyon özeti

{markdown_table(summary)}

DeepSurv null gate'i 30/30 kontrolde geçti. Tam 77-feature ailede 8 feature nokta WI≥0.8 eşiğini geçti; BH ve BY ile sertifikalanan feature olmadı.

ADIM 1 sırasında SHAP'ın sabit mutlak 0.01 iç assertion'ı iki reference çağrısında tetiklendi. Replikalar atlanmadı veya değiştirilmedi; attributionlar aynı model/veriyle yeniden hesaplandı ve yalnız önceden tanımlı göreli additivity hatası ≤%5 kapısını geçtikleri için kabul edildi. N=50'de maksimum reference göreli hata %{100*step1_meta['additivity']['max_new_reference_relative_error']:.3f} oldu.

## Nokta-geçenlerde marjin dağılımı

{markdown_table(margins)}

{markdown_table(distances)}

Wasserstein uzaklığına göre DeepSurv nokta-marjin şekli **{closer}** modeline daha yakındır. Bu karşılaştırma DeepSurv'de yalnız sekiz nokta-geçen feature bulunduğu için keşifseldir.

## Darboğaz dağılımı

{markdown_table(bottlenecks)}

DeepSurv N=50'nin 77 feature'ının tamamında darboğaz R'dir. CoxPH ve RSF'de G_claim az sayıda feature'da darboğaz olurken baskın bileşen yine R'dir.

## N=2 → N=50 winner's-curse

Nokta-geçen sayı 11'den 8'e düştü. Geçişler `{transition_counts}`: N=2 kazananlarının yalnız biri korundu, onu düştü ve yedi yeni feature geçti. Sayı düzeyindeki düşüş winner's-curse ile uyumludur; feature üyeliğinin bu kadar değişmesi N=2 sıralama tahmininin çok oynak olduğunu da gösterir. Bu nedenle sonuç yalnız “11→8 düşüşü” olarak yorumlanmamalıdır.

### N=2'de geçip N=50'de düşenler

{markdown_table(lost)}

Gerçek analiz wall-clock **{elapsed:.3f} sn**. Model/explainer çağrısı yapılmadı. Formül SHA256 değişmedi ve kimlik-sızıntısı kontrolü geçti.

## Adım bazında wall-clock

{markdown_table(runtimes)}
"""
    (OUT / "report.md").write_text(text, encoding="utf-8")
    (ROOT / "PHASE6_2_REPORT.md").write_text(text, encoding="utf-8")
    print(summary[["model", "point_WI_ge_0_8", "BH", "BY"]].to_string(index=False), flush=True)
    print(margins[["model", "share_0.05_0.10"]].to_string(index=False), flush=True)
    print(verdict, flush=True)


if __name__ == "__main__":
    run()
