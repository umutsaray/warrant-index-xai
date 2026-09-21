"""Create publication figures/tables from existing frozen artifacts only.

No model, explainer, bootstrap, or inferential computation is performed.  The
only derived quantities are counts of already stored bottleneck labels and
WI-0.8 margins for already stored point-passing features.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript_assets"
COMPONENTS = ["F", "S", "U", "G_claim", "R"]
MODELS = ["CoxPH N=50", "RSF N=50", "DeepSurv N=50"]
COLORS = ["#0072B2", "#E69F00", "#009E73"]  # Okabe-Ito subset
HATCHES = ["///", "\\\\", "..."]
MM_TO_INCH = 1 / 25.4


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "legend.fontsize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.linewidth": 0.75,
            "lines.linewidth": 1.0,
            "patch.linewidth": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
            "savefig.bbox": None,
            "savefig.pad_inches": 0,
        }
    )


def save_figure(fig: plt.Figure, number: int) -> dict[str, object]:
    paths = {suffix: OUT / f"figure_{number}.{suffix}" for suffix in ("png", "pdf", "eps")}
    fig.savefig(paths["png"], dpi=600, facecolor="white")
    fig.savefig(paths["pdf"], facecolor="white")
    fig.savefig(paths["eps"], facecolor="white")
    plt.close(fig)
    with Image.open(paths["png"]) as image:
        dpi = tuple(round(float(v), 1) for v in image.info.get("dpi", (0, 0)))
        pixels = image.size
    physical_mm = tuple(round(pixels[i] / dpi[i] * 25.4, 2) for i in range(2))
    pdf_bytes = paths["pdf"].read_bytes()
    return {
        "files": [str(path.relative_to(ROOT)) for path in paths.values()],
        "png_pixels": pixels,
        "png_dpi": dpi,
        "physical_size_mm": physical_mm,
        "pdf_contains_embedded_truetype_font": b"/FontFile2" in pdf_bytes,
    }


def add_box(ax, x, y, w, h, text, *, face="#F7F7F7", edge="#333333", fontsize=8) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.015",
            facecolor=face, edgecolor=edge, linewidth=0.9,
        )
    )
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize)


def arrow(ax, x1, x2, y=0.5) -> None:
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>", mutation_scale=10, lw=0.9, color="#333333"))


def figure_1() -> tuple[plt.Figure, str]:
    fig, ax = plt.subplots(figsize=(170 * MM_TO_INCH, 62 * MM_TO_INCH))
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.03, top=0.97)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    add_box(ax, 0.005, 0.18, 0.22, 0.64, "Component scores\n\n$F_j$  Faithfulness\n$S_j$  Stability\n$U_j$  Uncertainty\n$G^{claim}_j$  Generalizability\n$R_j$  Rank reproducibility", face="#E8F1F8", fontsize=7.0)
    add_box(ax, 0.265, 0.31, 0.17, 0.38, "$WI_j=\\min_k z_{jk}$\n\n$z_j=(F,S,U,$\n$G^{claim},R)$", face="#FFF4DE", fontsize=7.4)
    add_box(ax, 0.475, 0.20, 0.205, 0.60, "Point decision\n$\\mathbf{1}\\{z_{jk}>\\tau_k\\;\\forall k\\}$\n$\\tau_k=0.8$\n\nComponent tests\none-sided $p_{jk}$\n$p_j^{joint}=\\max_k p_{jk}$", face="#E8F6EF", fontsize=7.2)
    add_box(ax, 0.72, 0.29, 0.13, 0.42, "Full family\n$m=77$\nfeatures\n\nBH primary\nBY robustness", face="#F3EAF7", fontsize=7.2)
    add_box(ax, 0.875, 0.36, 0.11, 0.28, "Certified\nfeature set", face="#FDEDEC", fontsize=7.2)
    arrow(ax, 0.225, 0.265)
    arrow(ax, 0.435, 0.475)
    arrow(ax, 0.68, 0.72)
    arrow(ax, 0.85, 0.875)
    return fig, "Conceptual pipeline; no empirical values beyond frozen tau=0.8 and family size m=77."


def bar_style(rects, hatch: str) -> None:
    for rect in rects:
        rect.set_hatch(hatch)
        rect.set_edgecolor("#222222")
        rect.set_linewidth(0.75)


def figure_2(bottlenecks: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(170 * MM_TO_INCH, 78 * MM_TO_INCH))
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.16, top=0.93)
    x = np.arange(len(COMPONENTS))
    width = 0.23
    for i, model in enumerate(MODELS):
        vals = [int(bottlenecks.loc[model, component]) for component in COMPONENTS]
        bars = ax.bar(x + (i - 1) * width, vals, width, label=model.replace(" N=50", ""), color=COLORS[i])
        bar_style(bars, HATCHES[i])
    ax.set_xticks(x, ["F", "S", "U", r"$G^{claim}$", "R"])
    ax.set_ylabel("Features for which component is the bottleneck")
    ax.set_ylim(0, 82)
    ax.set_yticks(np.arange(0, 78, 10))
    ax.grid(axis="y", color="#D0D0D0", linewidth=0.5, linestyle=":")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=3, loc="upper left")
    return fig


def figure_3(summary: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(85 * MM_TO_INCH, 75 * MM_TO_INCH))
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.97)
    categories = ["Point WI≥0.8", "BH", "BY"]
    x = np.arange(len(MODELS))
    width = 0.23
    category_colors = ["#0072B2", "#E69F00", "#009E73"]
    category_hatches = ["///", "\\\\", "..."]
    for i, category in enumerate(categories):
        column = {"Point WI≥0.8": "point_WI_ge_0_8", "BH": "BH", "BY": "BY"}[category]
        vals = [int(summary.loc[model, column]) for model in MODELS]
        bars = ax.bar(x + (i - 1) * width, vals, width, label=category, color=category_colors[i])
        bar_style(bars, category_hatches[i])
        ax.bar_label(bars, labels=[str(v) for v in vals], padding=2, fontsize=6)
    ax.set_xticks(x, ["CoxPH", "RSF", "DeepSurv"])
    ax.set_ylabel("Feature count (of 77)")
    ax.set_ylim(0, 26)
    ax.grid(axis="y", color="#D0D0D0", linewidth=0.5, linestyle=":")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=1, loc="upper right")
    return fig


def figure_4(margins: dict[str, np.ndarray]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(170 * MM_TO_INCH, 78 * MM_TO_INCH))
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.16, top=0.97)
    positions = np.arange(1, 4)
    box = ax.boxplot(
        [margins[m] for m in MODELS], positions=positions, widths=0.42, patch_artist=True,
        showfliers=False, medianprops={"color": "#111111", "linewidth": 1.2},
        whiskerprops={"linewidth": 0.8}, capprops={"linewidth": 0.8}, boxprops={"linewidth": 0.8},
    )
    for patch, color, hatch in zip(box["boxes"], COLORS, HATCHES):
        patch.set_facecolor(color)
        patch.set_hatch(hatch)
        patch.set_edgecolor("#222222")
    for pos, model, color in zip(positions, MODELS, COLORS):
        values = np.sort(margins[model])
        offsets = np.linspace(-0.14, 0.14, len(values)) if len(values) > 1 else np.array([0.0])
        ax.scatter(pos + offsets, values, s=18, facecolor=color, edgecolor="#222222", linewidth=0.45, zorder=3)
    ax.axhline(0, color="#333333", linewidth=0.8, linestyle="--")
    ax.set_xticks(positions, ["CoxPH", "RSF", "DeepSurv"])
    ax.set_ylabel("Point-passing WI margin ($WI_j-0.8$)")
    ax.grid(axis="y", color="#D0D0D0", linewidth=0.5, linestyle=":")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def write_tables(summary: pd.DataFrame, bottlenecks: pd.DataFrame) -> None:
    table1 = r"""\begin{table}[t]
\caption{Point-estimate and multiplicity-adjusted Warrant Index results across the three model classes.}
\label{tab:three-model-comparison}
\centering
\begin{tabular}{lrrrr}
\toprule
Model & Point $WI\geq0.8$ & BH certified & BY certified & BH survival (\%) \\
\midrule
CoxPH ($N=50$)    & 22/77 & 8/77  & 8/77  & 36.4 \\
RSF ($N=50$)      & 20/77 & 10/77 & 10/77 & 50.0 \\
DeepSurv ($N=50$) & 8/77  & 0/77  & 0/77  & 0.0 \\
\botrule
\end{tabular}
\end{table}
"""
    (OUT / "table_1.tex").write_text(table1, encoding="utf-8")
    lines = [
        r"\begin{table}[t]",
        r"\caption{Distribution of the minimum-score (bottleneck) component across all 77 assessed features.}",
        r"\label{tab:bottleneck-distribution}",
        r"\centering",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r"Model & $F$ & $S$ & $U$ & $G^{\mathrm{claim}}$ & $R$ & Total \\",
        r"\midrule",
    ]
    for model in MODELS:
        label = model.replace(" N=50", r" ($N=50$)")
        vals = [int(bottlenecks.loc[model, c]) for c in COMPONENTS]
        lines.append(f"{label} & " + " & ".join(map(str, vals)) + r" & 77 \\")
    lines.extend([r"\botrule", r"\end{tabular}", r"\end{table}", ""])
    (OUT / "table_2.tex").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    configure()
    OUT.mkdir(exist_ok=True)
    summary_path = ROOT / "results/phase62_three_model_comparison/three_model_summary.csv"
    bottleneck_path = ROOT / "results/phase62_three_model_comparison/bottleneck_distribution.csv"
    summary = pd.read_csv(summary_path).set_index("model")
    stored_bottlenecks = pd.read_csv(bottleneck_path)
    stored_bottlenecks = stored_bottlenecks[stored_bottlenecks.scope == "all_77"]
    pivot = stored_bottlenecks.pivot(index="model", columns="bottleneck", values="count").fillna(0)
    bottlenecks = pivot.reindex(index=MODELS, columns=COMPONENTS, fill_value=0).astype(int)

    report_paths = {
        "CoxPH N=50": ROOT / "results/phase3_brca/feature_report.csv",
        "RSF N=50": ROOT / "results/phase56_rsf_n50/feature_report.csv",
        "DeepSurv N=50": ROOT / "results/phase62_deepsurv_n50/feature_report.csv",
    }
    certification_paths = {
        "CoxPH N=50": ROOT / "results/phase53_cox_recursive/feature_certification.csv",
        "RSF N=50": ROOT / "results/phase56_rsf_n50_certification/feature_certification.csv",
        "DeepSurv N=50": ROOT / "results/phase62_deepsurv_n50_certification/feature_certification.csv",
    }
    derived_counts = {}
    margins = {}
    for model, path in report_paths.items():
        frame = pd.read_csv(path)
        frame = frame[frame.WI.notna()]
        assert len(frame) == 77
        component_values = frame[COMPONENTS].to_numpy(float)
        recomputed = np.array(COMPONENTS)[np.argmin(component_values, axis=1)]
        assert np.array_equal(recomputed, frame.bottleneck.to_numpy())
        derived_counts[model] = pd.Series(recomputed).value_counts().reindex(COMPONENTS, fill_value=0)
    derived = pd.DataFrame(derived_counts).T.astype(int)
    pd.testing.assert_frame_equal(
        derived.loc[MODELS, COMPONENTS], bottlenecks.loc[MODELS, COMPONENTS], check_names=False
    )
    assert (bottlenecks.sum(axis=1) == 77).all()
    for model, path in certification_paths.items():
        frame = pd.read_csv(path)
        selected = frame[frame.point_WI_ge_0_8.astype(str).str.lower() == "true"]
        values = selected.WI.to_numpy(float) - 0.8
        assert len(values) == int(summary.loc[model, "point_WI_ge_0_8"])
        assert np.all(values >= -1e-12)
        margins[model] = values

    write_tables(summary, bottlenecks)
    generated = {}
    fig, note = figure_1()
    generated["figure_1"] = save_figure(fig, 1) | {"note": note}
    generated["figure_2"] = save_figure(figure_2(bottlenecks), 2)
    generated["figure_3"] = save_figure(figure_3(summary), 3)
    generated["figure_4"] = save_figure(figure_4(margins), 4)

    for name, details in generated.items():
        assert details["png_dpi"][0] >= 599 and details["png_dpi"][1] >= 599, (name, details)
        assert details["pdf_contains_embedded_truetype_font"], (name, "PDF font was not embedded")
    assert abs(generated["figure_3"]["physical_size_mm"][0] - 85) < 0.1
    for name in ("figure_1", "figure_2", "figure_4"):
        assert abs(generated[name]["physical_size_mm"][0] - 170) < 0.1
    for path in OUT.glob("figure_*.eps"):
        assert "/FontType 42 def" in path.read_text(encoding="latin-1", errors="ignore")
    tex = (OUT / "table_1.tex").read_text(encoding="utf-8") + (OUT / "table_2.tex").read_text(encoding="utf-8")
    assert not re.search(r"\d,\d{3}", tex)

    metadata = {
        "new_experiment_or_model_call": False,
        "threshold": 0.8,
        "family_size": 77,
        "table_1": summary.loc[MODELS, ["point_WI_ge_0_8", "BH", "BY", "BH_survival_fraction"]].to_dict("index"),
        "table_2": bottlenecks.loc[MODELS, COMPONENTS].to_dict("index"),
        "figure_4_point_counts": {model: len(values) for model, values in margins.items()},
        "outputs": generated,
        "privacy": {"patient_or_sample_identifiers_exported": False, "feature_labels_exported_in_figures": False},
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    provenance = f"""# Figure/Table Provenance Log

No new experiment, model fit, explainer call, bootstrap, threshold, or inferential calculation was performed. All visual elements use frozen aggregate or feature-level artifacts. Patient and sample identifiers are absent from every output; Figure 4 points are intentionally unlabeled.

## Table 1 — Three-model comparison

- Output: `manuscript_assets/table_1.tex`.
- Source: `results/phase62_three_model_comparison/three_model_summary.csv:2-4`.
- Values: CoxPH 22/77 point, 8/77 BH, 8/77 BY, 36.4%; RSF 20/77, 10/77, 10/77, 50.0%; DeepSurv 8/77, 0/77, 0/77, 0.0%.
- Cross-check: `RESULTS_LEDGER.md`, final three-model comparison.

## Table 2 — Bottleneck distribution

- Output: `manuscript_assets/table_2.tex`.
- Primary source: `results/phase62_three_model_comparison/bottleneck_distribution.csv`, rows with `scope=all_77`.
- Independent check: stored `bottleneck` labels in `results/phase3_brca/feature_report.csv`, `results/phase56_rsf_n50/feature_report.csv`, and `results/phase62_deepsurv_n50/feature_report.csv` exactly equal `argmin(F,S,U,G_claim,R)` for all 77 assessed features in each model.
- Values: CoxPH F=0, S=0, U=0, G_claim=1, R=76; RSF F=0, S=0, U=0, G_claim=5, R=72; DeepSurv F=0, S=0, U=0, G_claim=0, R=77. Every row sums to 77.
- Cross-check: `RESULTS_LEDGER.md`, final three-model comparison and bottleneck entries.

## Figure 1 — Warrant Index audit pipeline

- Outputs: `figure_1.png`, `figure_1.pdf`, `figure_1.eps`.
- Method source: verbatim frozen equations in `RESULTS_LEDGER.md` and Sections 2.2–2.4 of `MANUSCRIPT_METHODS_TO_AVAILABILITY.md`.
- Frozen elements shown: F/S/U/G_claim/R, $WI_j=min_k z_jk$, strict component threshold $tau_k=0.8$, $p_j^joint=max_k p_jk$, full family m=77, BH primary, and BY robustness.

## Figure 2 — Bottleneck distribution

- Outputs: `figure_2.png`, `figure_2.pdf`, `figure_2.eps`.
- Data: exactly Table 2 above; no interpolation or simulation.

## Figure 3 — Point pass and certification counts

- Outputs: `figure_3.png`, `figure_3.pdf`, `figure_3.eps`.
- Data: exactly Table 1 count columns from `three_model_summary.csv`; no interpolation or simulation.

## Figure 4 — Point-passing WI margins

- Outputs: `figure_4.png`, `figure_4.pdf`, `figure_4.eps`.
- Data availability decision: produced. Per-feature WI values and point-pass flags exist for all three models in the certification CSV files listed below.
- Sources: `results/phase53_cox_recursive/feature_certification.csv`, `results/phase56_rsf_n50_certification/feature_certification.csv`, and `results/phase62_deepsurv_n50_certification/feature_certification.csv`.
- Transformation: select stored `point_WI_ge_0_8=True` rows and plot the exact stored `WI-0.8` value. Counts are CoxPH={len(margins['CoxPH N=50'])}, RSF={len(margins['RSF N=50'])}, DeepSurv={len(margins['DeepSurv N=50'])}. Deterministic horizontal offsets prevent overplotting and encode no data.

## Output-format checks

- PNG files were saved at 600 dpi; recorded pixel dimensions and DPI are in `manuscript_assets/metadata.json`.
- PDF files contain embedded TrueType font streams (`/FontFile2`) under Matplotlib `pdf.fonttype=42`; EPS uses `ps.fonttype=42`.
- Figure canvases are exactly 170 mm wide (Figures 1, 2, and 4) or 85 mm wide (Figure 3), within pixel rounding at 600 dpi; all heights are below 225 mm.
- Lines and patch borders are at least 0.5 pt; data-bearing borders are at least 0.75 pt. Colorblind-safe colors are reinforced with distinct hatches.
- Tables use `booktabs`, no color/shading, and no thousands-separator commas.

## Skipped outputs

None. Table 2 was supported by both stored bottleneck labels and component scores. Figure 4 was supported by stored per-feature margins for all three models.
"""
    (ROOT / "FIGURES_TABLES_PROVENANCE.md").write_text(provenance, encoding="utf-8")
    print(json.dumps({"figures": 4, "tables": 2, "output": str(OUT), "table_2": True, "figure_4": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
