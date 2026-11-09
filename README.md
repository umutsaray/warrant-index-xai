# Warrant Index for Survival-XAI Audits: Statistical Certification (BH/BY, Tolerance-Conditional; Not Biological or Clinical Validation)

This repository implements the Warrant Index (WI), a non-compensatory audit of feature-level explanation credibility, and the experiments used to evaluate it across CoxPH, random survival forests, and DeepSurv in breast-cancer survival analysis.

## Method

For feature (j), the primary score is

\[
WI_j = \min(F_j, S_j, U_j, G_j^{claim}, R_j),
\qquad G_j^{claim}=\min(G_j,G_j^{signed}).
\]

The five components measure faithfulness, stability, uncertainty control, cross-cohort generalizability, and rank reproducibility. Point adequacy does not imply statistical certification. The certification layer tests every component against the frozen threshold, combines component p-values with an intersection-union maximum, and controls the full feature family with Benjamini-Hochberg (primary) and Benjamini-Yekutieli (robustness) adjustments.

The frozen experimental settings are `a=b=0.5`, `u=3`, `eta=0.02`, and `tau=0.8`. These are experimental settings rather than clinically validated universal cutoffs.

## Analysis tiers

- **v1 official analysis (pre-specified/frozen):** the five-component index uses the absolute tolerances `a=b=0.5`, `u=3`, `eta=0.02`, and `tau=0.8`. This is the official analysis reported in the main results.
- **v2 exploratory analysis (post hoc, then preregistered before v2 computation):** Phase 17 replaces the absolute tolerances with model-scale-normalized tolerances and adds the sixth exercised-signal component `N`. The frozen v2 definitions and grid are in [`PHASE17_PREREG.md`](PHASE17_PREREG.md); results are in [`PHASE17_V2.md`](PHASE17_V2.md).
- **Phase 21 v3 diagnostic replication (post hoc, preregistered before v3 computation):** the leakage-controlled screening and matched-family replication is documented under [`phase21_reviewer1/`](phase21_reviewer1/).
- **Tolerance profile and synthetic diagnostics (post hoc, then preregistered before computation):** Phase 18 profiles `k`, evaluates the data-generating-process scale, and tests synthetic discrimination. The plan is in [`PHASE18_PREREG.md`](PHASE18_PREREG.md) and the results are in [`PHASE18_KAPANIS.md`](PHASE18_KAPANIS.md). No standalone Phase 19 preregistration or report is present in this repository snapshot.

The frozen v1 analysis is the official analysis tier. The v2, tolerance-profile, and Phase 21 v3 analyses are post hoc exploratory or diagnostic replication layers and do **not** replace or modify the frozen v1 results.

## Key findings and limitations

- The v1 analysis produced 20 BH-certified model-feature pairs: CoxPH 8, RSF 12, and DeepSurv 0 ([Phase 17](PHASE17_V2.md)).
- The primary v2 setting produced 0 point-adequate and 0 BH/BY-certified pairs in every model ([Phase 17](PHASE17_V2.md)).
- In the synthetic analysis, the v1 absolute tolerance was approximately 19.6 times the simulated attribution scale and the v1 result was saturated across classes ([Phase 18](PHASE18_KAPANIS.md)).
- The synthetic mechanism/sign-shift class D was not separated from the reliable-signal class by the evaluated profile; its primary AUC against the zero-effect class was 0.952 ([Phase 18](PHASE18_KAPANIS.md)).
- Across `k <= 128`, the real-data profile admitted 8/27 CoxPH, 2/27 RSF, and 0/27 DeepSurv features ([Phase 18](PHASE18_KAPANIS.md)).
- These findings are conditional on the audited cohorts, query sets, explanation estimators, replication design, and experimental thresholds; they are not evidence of causality or clinical utility.

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

## Reproducing

1. Run the synthetic validation scripts `experiments/01_*` through `05_*`.
2. Obtain TCGA-BRCA data from the NCI Genomic Data Commons and prepare the patient-grouped split with the companion preprocessing workflow.
3. Obtain METABRIC data through cBioPortal under its applicable terms, then run the Phase 3 panel and mapping scripts.
4. Run the model-specific audit scripts in numerical phase order. The Phase 7 scripts rebuild the clean patient-grouped reference audit; Phase 7.2 performs the N=50 audits, calibration gates, and BH/BY certification.
5. Run the Phase 7.3-8 scripts to reproduce the publication tables and robustness diagnostics from stored local artifacts.
6. For the exploratory v2 layer, run `experiments/75_phase17_v1_gate.py` through `experiments/81_phase17_v2_neighbors.py` in numerical order. Measured reference times include 8.9 s for the primary certification, 83.6 s for the neighboring settings, and 2,029.6 s for the 20-repeat synthetic run.
7. For the tolerance-profile layer, run `experiments/82_phase18_gate.py` through `experiments/89_phase18_kstar_eval_boot.py` in numerical order. The successful Phase 18 main steps took approximately 2.44 h in the recorded environment; individual measured times are reported in [`PHASE18_KAPANIS.md`](PHASE18_KAPANIS.md).

The real-data scripts intentionally fail when required local data or model artifacts are absent. They never substitute fabricated data.

The repository requires Python 3.10 or newer. The current regression command is:

```bash
python -m pytest -q
```

The recorded complete suite contains 211 tests. Phase 17–18 real-data gates require the locally stored, non-versioned TCGA/METABRIC-derived matrices, fitted models, and attribution tensors. Synthetic-only Phase 17–18 steps do not download patient-level data. There are no Phase 19 scripts in this repository snapshot.

## Data access

- **TCGA-BRCA:** [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/)
- **METABRIC:** [cBioPortal METABRIC study](https://www.cbioportal.org/study/summary?id=brca_metabric), subject to cBioPortal and source-study terms

The Phase 20 repository scan found no raw or patient-level TCGA/METABRIC data in the tracked Git tree or its reachable history. Local patient-level caches, matrices, fitted models, checkpoints, and attribution tensors are excluded by `.gitignore`. The tracked supplementary CSV files contain feature-level aggregate or synthetic results only. See [`PHASE20_ADIM1_TARAMA.md`](PHASE20_ADIM1_TARAMA.md).

## Important interpretation limits

A high WI is a warrant for the audited explanation claim under the specified data, query set, estimator, tolerances, and replication design. It is not evidence of causality or clinical utility. A feature with identically zero attributions can be stable and reproducible; the Phase 7.6-7.7 reports therefore distinguish informative from degenerate certified warrants without changing WI or BH/BY decisions.

RSF and DeepSurv audit results are conditional on the fixed 12+12 query sets and their frozen explainer budgets. CoxPH faithfulness measures numerical self-consistency of the exact closed-form attribution rather than agreement between independent estimators.

## Licensing

- Source code under `warrantxai/`, `experiments/`, and `tests/` is licensed under the [MIT License](LICENSE).
- Repository-produced derived result tables, supplementary materials under `results/**` and `manuscript_assets/**`, and `PHASE*` reports are licensed under [Creative Commons Attribution 4.0 International](LICENSE-DATA).
- Source TCGA and METABRIC data are **not** included in this repository and are not relicensed here. They remain subject to the terms of the [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/) and the [cBioPortal METABRIC study](https://www.cbioportal.org/study/summary?id=brca_metabric) and its source study.
