# Warrant Index manuscript results ledger

This ledger is the numerical source of truth for `MANUSCRIPT_METHODS_TO_AVAILABILITY.md`. Values were re-read from stored artifacts; they were not copied from conversational memory.

## Source-integrity note

The requested file `section2-warrant-index-taslak.md` is absent from the workspace, the parent project, and the user-home filename search. The preserved Phase 2 kickoff attachment explicitly labels its equations as fixed Section 2 content. Accordingly, the verbatim formula block below is sourced from `C:/Users/UMUT/.codex/attachments/8129a6e2-5b8d-4427-88da-7b5e0036d95d/pasted-text.txt:19-46`; the signed extension is sourced from `PHASE2_REPORT_TR.md:72-74`. This missing-file condition must be resolved before submission by comparing the block against the canonical manuscript source.

## Verbatim frozen equations

```text
WI_j = min_k z_{jk}                              (Eq. 2.3.1)
Adequate_j = 1{z_jk > tau_k  for all k}           (Eq. 2.3.2)
WI_claim = min_{j in J*} WI_j                     (Eq. 2.3.3, J* önceden kayıtlı feature seti)

F_j = min_d (1/T) sum_t F_hat_{dtj}
  F_hat_{dtj} = 1 - (1/n) sum_i min(1, (phi_bar_{dtj}(x_i) - psi_tilde_j(f_dt,x_i))^2 / a_j^2)

S_j = min_d S_{dj}
  S_{dj} = 1 - (2/(T(T-1)n)) sum_{t<t'} sum_i min(1, (phi_bar_{dtj}(x_i)-phi_bar_{dt'j}(x_i))^2/a_j^2)

U_j = min_d U_{dj}
  U_{dj} = (1/T) sum_t [1 - min(1, (H_{dtj}-L_{dtj})/u_j)]

G_j = min_{d>=1} [1 - min(1, |I_hat_{dj} - I_hat_{0j}| / b_j)]

G_signed_j = min_d>=1 [1 - mean_i min(1, |phi_hat_dj(x_i)-phi_hat_0j(x_i)|/a_j)]
G_claim_j = min(G_j, G_signed_j)

R_j = min_d R_{dj}
  R_{dj} = 1 - (2/(T(T-1)(p-1))) sum_{t<t'} sum_{l!=j} |o_{dt,jl}-o_{dt',jl}|/2

Joint p-value: p_joint_j = max_k p_jk  (intersection-union test)
```

## Data and model inventory

| Ledger key | Verified value | Source |
|---|---:|---|
| `tcga_train_rows` | 673 | `results/phase3_brca/metadata.json:25` |
| `tcga_test_rows` | 169 | `results/phase3_brca/metadata.json:27` |
| `selected_features` | 616 | `results/phase3_brca/metadata.json:29` |
| `cox_model_features` | 614 | `results/phase3_brca/metadata.json:30` |
| `tcga_audit_patients` | 141 | `results/phase3_brca/metadata.json:33` |
| `metabric_complete_queries` | 1466 | `results/phase3_brca/metadata.json:165` |
| `common_features` | 77 = 72 RNA + 5 clinical | `PHASE3_REPORT_TR.md:17-27`; `PHASE4_DESIGN_NOTE.md:19` |
| `rsf_deepsurv_query_rows` | 12 TCGA + 12 METABRIC | `results/phase4_rsf/metadata.json:19-22`; `results/phase62_deepsurv_n50/metadata.json:3-10` |
| `cox_test_cindex` | 0.7064765300 | sibling `baseline_results.json`, `models.coxph.test_c_index` |
| `rsf_test_cindex` | 0.7147950089 | sibling `baseline_results.json`, `models.random_survival_forest.test_c_index` |
| `deepsurv_test_cindex` | 0.7456922163 | `results/phase6_deepsurv_pilot/metadata.json:16` |
| `deepsurv_best_epoch` | 26 | `results/phase6_deepsurv_pilot/metadata.json:17` |
| `deepsurv_dropout` | 0.3 | `results/phase6_deepsurv_pilot/metadata.json:19` |
| `fixed_tolerances` | a=b=0.5, u=3, eta=0.02, tau=0.8 | `results/phase3_brca/metadata.json:14-20` |

DeepSurv architecture: block encoders RNA 200–128–32, miRNA 200–64–16, methylation 200–128–32, clinical 16–32–8; head 88–64–16–1. Optimizer: AdamW, learning rate 0.001, weight decay 0.0001, maximum 300 epochs, patience 30 (`results/phase6_deepsurv_pilot/metadata.json:21-22`).

### CoxPH F computation provenance

CoxPH candidate and reference attributions are identical exact closed-form arrays. `experiments/06_brca_real.py:394-398` constructs `WarrantIndexEvaluator` without an `attributor`; `warrantxai/api.py:129-149` first sets the reference to `(anchors-background)*coefficients` and, when `attributor is None`, assigns that same array to the candidate. Thus F=1 is an exact-computation self-consistency check, not agreement with an independent numerical explainer. `results/phase3_brca/metadata.json:50` stores only `oracle_max_abs_error=2.842170943040401e-14`; it does not contain a candidate/reference-method field. The same interpretation is stated explicitly in `PHASE3_REPORT_TR.md:23` ("F=1 oracle'ın kendisiyle kontrolüdür; harici explainer doğrulaması değildir").

## CoxPH progression and final results

| Ledger key | Verified value | Source |
|---|---:|---|
| `cox_n6_point` | 27/77 | `results/phase3_brca/metadata.json:65` |
| `cox_n50_point` | 22/77 | `results/phase3_brca/metadata.json:51-52`; `results/phase53_cox_recursive/summary.csv:2-3` |
| `cox_n6_to_n50_retained` | 15/27 | `results/phase3_brca/metadata.json:53-54` |
| `cox_n6_R_median_IQR` | 0.0303834 | `results/phase3_brca/metadata.json:75` |
| `cox_n50_R_median_IQR` | 0.0193074 | `results/phase3_brca/metadata.json:76` |
| `cox_n50_bottlenecks` | R=76, G_claim=1 | `results/phase3_brca/metadata.json:70-73` |
| `cox_recursive_BH` | 8/77 | `results/phase53_cox_recursive/summary.csv:5` |
| `cox_recursive_BY` | 8/77 | `results/phase53_cox_recursive/summary.csv:6` |

The original basic-root and recursive-studentized methods produced the same eight CoxPH BH/BY certificates (`results/phase53_cox_recursive/summary.csv:4-9`).

## RSF progression, design decision, and final results

| Ledger key | Verified value | Source |
|---|---:|---|
| `rsf_fixed_bg_n2_point` | 43/77 | `results/phase4_rsf/metadata.json:35` |
| `rsf_sampled_bg_n2_point` | 35/77 | `results/phase4_rsf/metadata.json:36` |
| `rsf_n20_point` | 23/77 | `results/phase4_rsf/metadata.json:37` |
| `rsf_n20_bottlenecks` | R=72, G_claim=5 | `results/phase4_rsf/metadata.json:41-44` |
| `fixed_vs_sampled_max_delta_Gclaim` | 0.250446 | `PHASE4_DESIGN_NOTE.md:66` |
| `reference_128_to_256_max_delta_F` | 0.001367 | `PHASE4_DESIGN_NOTE.md:81` |
| `approved_rsf_background` | sampled k=16 unique TCGA patients; 8/128 paths | `PHASE4_DESIGN_NOTE.md:85`; `results/phase4_rsf/metadata.json:6-10` |
| `rsf_n50_point` | 20/77 | `results/phase56_rsf_n50/metadata.json:36` |
| `rsf_n20_to_n50_transition` | retained=16, lost=7, new=4 | `results/phase56_rsf_n50/metadata.json:37-39` |
| `rsf_n50_bottlenecks` | R=72, G_claim=5 | `results/phase56_rsf_n50/metadata.json:40-43` |
| `rsf_n50_BH` | 10/77 | `results/phase56_rsf_n50_certification/summary.csv:4` |
| `rsf_n50_BY` | 10/77 | `results/phase56_rsf_n50_certification/summary.csv:5` |

## Calibration and recursive-studentized inference

The first direct-root calibration used 2,000 outer simulations and 999 inner bootstrap resamples for two null configurations (`single_active` and `tie`), five components, and alpha in {0.01, 0.05, 0.10}. CoxPH had no Holm-corrected violation, whereas RSF N=20 failed in `single_active`: F at alpha=0.01 (0.0205), U at alpha=0.01 (0.0215), and G_claim at alpha=0.01/0.05 (0.0345/0.0730), from `PHASE5_1_NULL_CALIBRATION_REPORT.md:70-80,103-112`.

Three remedies were tested on the same battery. `recursive_basic` and `m_out_of_n_root` failed; `recursive_studentized` passed (`results/phase52_rsf_remedy/candidate_gate_summary.csv:1-4`). The accepted method studentizes the F, U, G_magnitude, and G_signed leaves, uses basic bootstrap tests for S and R leaves, and combines branch/site p-values by `max` at every `min` node; the five component p-values are again combined by the IUT maximum (`PHASE5_3_REPORT.md:7`; `PHASE5_3_COX_METHOD_COMPARISON_REPORT.md:53`).

Fresh N=50 gates both passed with 0 of 30 Holm-flagged controls: RSF (`results/phase56_rsf_n50_null_gate/gate_summary.csv:2`) and DeepSurv (`results/phase62_deepsurv_n50_null_gate/gate_summary.csv:2`). Each gate retained the 2,000×999 design.

## Winner's-curse and rejected rank-universe explanation

| Ledger key | Verified value | Source |
|---|---:|---|
| `rsf_selected23_mean_delta_R` | -0.0156899 | `results/phase56_comparison/R_delta_summary.csv:3` |
| `rsf_selected23_wilcoxon_p` | 0.00349458 | `results/phase56_comparison/R_delta_summary.csv:3` |
| `rsf_all77_mean_delta_R` | -0.00057505 | `results/phase56_comparison/R_delta_summary.csv:2` |
| `rsf_all77_wilcoxon_p` | 0.743291 | `results/phase56_comparison/R_delta_summary.csv:2` |
| `rank_universe_mean_delta_same22` | +0.00283 | `PHASE5_5_RANK_UNIVERSE_DIAGNOSTIC.md:7` |
| `rank_universe_W1_same22` | 0.05857→0.06076 | `PHASE5_5_RANK_UNIVERSE_DIAGNOSTIC.md:9`; `results/phase55_rank_universe_diagnostic/distance_to_rsf.csv:2` |
| `rank_universe_W1_all77` | 0.02174→0.02381 | `PHASE5_5_RANK_UNIVERSE_DIAGNOSTIC.md:9`; `results/phase55_rank_universe_diagnostic/distance_to_rsf.csv:3` |

The 614-versus-77 ranking-universe hypothesis was therefore not supported; restricting CoxPH to 77 competitors did not move its margin distribution toward RSF.

## DeepSurv progression and final results

| Ledger key | Verified value | Source |
|---|---:|---|
| `deepsurv_n2_point` | 11/77 | `results/phase61_deepsurv_explainer_revision/metadata.json:18` |
| `deepsurv_n2_bottlenecks` | R=75, G_claim=2 | `results/phase61_deepsurv_explainer_revision/metadata.json:19-22` |
| `deepsurv_design` | DeepExplainer candidate k=16/reference k=128 | `results/phase61_deepsurv_explainer_revision/metadata.json:9-11` |
| `deepsurv_n50_point` | 8/77 | `results/phase62_deepsurv_n50/metadata.json:21` |
| `deepsurv_n50_BH_BY` | 0/77, 0/77 | `results/phase62_deepsurv_n50_certification/summary.csv:4-5` |
| `deepsurv_n50_bottleneck` | R=77 | `results/phase62_deepsurv_n50/metadata.json:25-27` |
| `deepsurv_transition` | neither=59, lost=10, new=7, retained=1 | `results/phase62_three_model_comparison/metadata.json`, `transition_counts` |
| `deepsurv_margin_mean` | 0.00930518 | `results/phase62_three_model_comparison/margin_distribution.csv:4` |
| `deepsurv_margin_median` | 0.00946563 | `results/phase62_three_model_comparison/margin_distribution.csv:4` |
| `deepsurv_margin_share_005_010` | 0% | `results/phase62_three_model_comparison/margin_distribution.csv:4` |
| `deepsurv_null_gate` | 0/30 violations; pass | `results/phase62_deepsurv_n50_null_gate/gate_summary.csv:2` |

Two N=50 reference calls tripped SHAP's fixed absolute 0.01 internal assertion. No replicate was replaced or skipped; both recomputations passed the preregistered relative-error ≤5% gate, with maximum reference relative error 0.8075% (`results/phase62_deepsurv_n50/metadata.json:37-41`).

## Final three-model comparison

| Model | Point WI≥0.8 | BH | BY | BH survival |
|---|---:|---:|---:|---:|
| CoxPH N=50 | 22/77 | 8/77 | 8/77 | 36.4% |
| RSF N=50 | 20/77 | 10/77 | 10/77 | 50.0% |
| DeepSurv N=50 | 8/77 | 0/77 | 0/77 | 0% |

Source: `results/phase62_three_model_comparison/three_model_summary.csv:2-4`. Full-family bottlenecks were CoxPH R=76/G_claim=1, RSF R=72/G_claim=5, and DeepSurv R=77 (`results/phase62_three_model_comparison/bottleneck_distribution.csv:2-10`). The point-passing 0.05–0.10 margin shares were CoxPH 4.55%, RSF 45.0%, and DeepSurv 0% (`results/phase62_three_model_comparison/margin_distribution.csv:2-4`).

## Rounded presentation values and design constants used in the draft

This table makes the draft's rounded display literals traceable to their higher-precision sources. Section numbers, equation indices, and the algebraic constants 0, 1, and 2 inside the verbatim formulas are structural rather than empirical results.

| Draft literal(s) | Source value / meaning | Source |
|---|---|---|
| 0.7065; 0.7148 | CoxPH and RSF test C-indices | sibling `baseline_results.json`, `models.*.test_c_index` |
| 0.7457; epoch 26; seed 42 | DeepSurv test C-index, selected epoch, random state | `results/phase6_deepsurv_pilot/metadata.json`, `adim0` |
| 0.1; 200; 0.25 | CoxPH penalizer, retry steps, retry step size | `results/phase3_brca/metadata.json`, `arguments` |
| 2.84×10^-14 | CoxPH oracle reconstruction error (full value 2.842170943040401e-14) | `results/phase3_brca/metadata.json`, `stats.oracle_max_abs_error` |
| 500; 10; 5; 42 | RSF trees, split/leaf minima, random state | sibling `baseline_results.json`, `models.random_survival_forest` and stored model design |
| 539 | inactive nonlinear-model inputs (616−77) | `results/phase56_rsf_n50/metadata.json`, `design.inactive_background_integrated_features` |
| 256; 2,048 | first GradientExplainer budget passing the relative additivity gate | `results/phase61_deepsurv_explainer_revision/gradient_budget_curve.csv`; `PHASE6_1_EXPLAINER_REVISION.md` |
| 0.8075% | maximum DeepSurv N=50 reference relative error (full value 0.0080748124) | `results/phase62_deepsurv_n50/metadata.json`, `additivity` |
| 13; 11 | CoxPH remediated and numerical-retry replicate counts | `results/phase3_brca/metadata.json`, `stats` |
| 0.03038; 0.01931 | rounded N=6 and N=50 R median IQR widths | `results/phase3_brca/metadata.json`, `stats.R_*_median_IQR_width` |
| −0.01569; 0.00349 | rounded selected-23 mean ΔR and Wilcoxon p | `results/phase56_comparison/R_delta_summary.csv:3` |
| −0.000575; 39; 38; 0.74329 | rounded all-77 mean ΔR, decrease/increase counts, and Wilcoxon p | `results/phase56_comparison/R_delta_summary.csv:2` |
| 0.00931; 0.00947 | rounded DeepSurv point-passer mean and median margins | `results/phase62_three_model_comparison/margin_distribution.csv:4` |
| −1, 0, 1 | values of the η-tolerant pairwise-order indicator | frozen R definition and implementation |
| 20260917 | common bootstrap seed | N=50 model metadata files, `arguments.seed` |
| 4,999 | final full-family certification bootstrap draws | certification metadata and Phase 5 reports |
| 0.01; 0.03; 0.05; 0.10 | displayed calibration alpha levels and margin-bin boundaries | null-gate and three-model comparison artifacts |

## Ledger exclusions

- No patient or sample identifiers are reproduced.
- No external-survival performance claim is made for METABRIC; its outcome labels were not used for fitting or tuning.
- Abstract, Introduction, Related Work, literature citations, public repository URL, and final license statement remain outside this deliverable.
