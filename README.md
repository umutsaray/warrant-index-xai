# Warrant Index for Statistically Certified Survival-XAI Audits

This repository implements the Warrant Index (WI), a non-compensatory audit of feature-level explanation credibility, and the experiments used to evaluate it across CoxPH, random survival forests, and DeepSurv in breast-cancer survival analysis.

## Method

For feature (j), the primary score is

\[
WI_j = \min(F_j, S_j, U_j, G_j^{claim}, R_j),
\qquad G_j^{claim}=\min(G_j,G_j^{signed}).
\]

The five components measure faithfulness, stability, uncertainty control, cross-cohort generalizability, and rank reproducibility. Point adequacy does not imply statistical certification. The certification layer tests every component against the frozen threshold, combines component p-values with an intersection-union maximum, and controls the full feature family with Benjamini-Hochberg (primary) and Benjamini-Yekutieli (robustness) adjustments.

The frozen experimental settings are `a=b=0.5`, `u=3`, `eta=0.02`, and `tau=0.8`. These are experimental settings rather than clinically validated universal cutoffs.

## Repository contents

- `warrantxai/`: WI components, aggregation, certification, synthetic data generation, and public API.
- `experiments/`: synthetic experiments and the Phase 3-8 survival-XAI analysis scripts.
- `tests/`: unit, regression, privacy, artifact, and formula-integrity tests.
- Root-level `PHASE*.md` and aggregate CSV files: publication-facing reports and supplementary summaries that contain no patient identifiers.
- `manuscript_assets/`: publication figures and aggregate LaTeX tables.

Patient-level data, fitted models, checkpoints, raw cBioPortal responses, attribution tensors, and local result caches are deliberately excluded from version control.

## Installation

Python 3.10 or newer is required. The exact package versions used for the final local analysis are recorded in `requirements.txt`.

```bash
python -m venv .venv
# Windows
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python -m pip install -e .
# Linux/macOS
# .venv/bin/python -m pip install -r requirements.txt
# .venv/bin/python -m pip install -e .
```

Run the test suite:

```bash
python -m pytest -q
```

## High-level reproduction workflow

1. Run the synthetic validation scripts `experiments/01_*` through `05_*`.
2. Obtain TCGA-BRCA data from the NCI Genomic Data Commons and prepare the patient-grouped split with the companion preprocessing workflow.
3. Obtain METABRIC data through cBioPortal under its applicable terms, then run the Phase 3 panel and mapping scripts.
4. Run the model-specific audit scripts in numerical phase order. The Phase 7 scripts rebuild the clean patient-grouped reference audit; Phase 7.2 performs the N=50 audits, calibration gates, and BH/BY certification.
5. Run the Phase 7.3-8 scripts to reproduce the publication tables and robustness diagnostics from stored local artifacts.

The real-data scripts intentionally fail when required local data or model artifacts are absent. They never substitute fabricated data.

## Data availability

- **TCGA-BRCA:** [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/)
- **METABRIC:** [cBioPortal](https://www.cbioportal.org/), subject to cBioPortal and source-study terms

No raw or patient-level TCGA/METABRIC data are redistributed here. The tracked supplementary CSV files contain feature-level aggregate results only.

## Important interpretation limits

A high WI is a warrant for the audited explanation claim under the specified data, query set, estimator, tolerances, and replication design. It is not evidence of causality or clinical utility. A feature with identically zero attributions can be stable and reproducible; the Phase 7.6-7.7 reports therefore distinguish informative from degenerate certified warrants without changing WI or BH/BY decisions.

RSF and DeepSurv audit results are conditional on the fixed 12+12 query sets and their frozen explainer budgets. CoxPH faithfulness measures numerical self-consistency of the exact closed-form attribution rather than agreement between independent estimators.

## License

MIT License. See `LICENSE`.
