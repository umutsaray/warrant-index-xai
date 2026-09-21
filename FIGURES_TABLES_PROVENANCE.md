# Figure/Table Provenance Log

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
- Transformation: select stored `point_WI_ge_0_8=True` rows and plot the exact stored `WI-0.8` value. Counts are CoxPH=22, RSF=20, DeepSurv=8. Deterministic horizontal offsets prevent overplotting and encode no data.

## Output-format checks

- PNG files were saved at 600 dpi; recorded pixel dimensions and DPI are in `manuscript_assets/metadata.json`.
- PDF files contain embedded TrueType font streams (`/FontFile2`) under Matplotlib `pdf.fonttype=42`; EPS uses `ps.fonttype=42`.
- Figure canvases are exactly 170 mm wide (Figures 1, 2, and 4) or 85 mm wide (Figure 3), within pixel rounding at 600 dpi; all heights are below 225 mm.
- Lines and patch borders are at least 0.5 pt; data-bearing borders are at least 0.75 pt. Colorblind-safe colors are reinforced with distinct hatches.
- Tables use `booktabs`, no color/shading, and no thousands-separator commas.

## Skipped outputs

None. Table 2 was supported by both stored bottleneck labels and component scores. Figure 4 was supported by stored per-feature margins for all three models.
