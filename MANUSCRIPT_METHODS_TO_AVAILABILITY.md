# Manuscript draft: Methods through Availability

> Scope: first factual draft for later editorial integration. Abstract, Introduction, Related Work, and external-literature citations are intentionally omitted. All empirical values are tied to `RESULTS_LEDGER.md`.

## 2. Methods

### 2.1 Audit object and notation

We evaluated whether a feature-level attribution claim was simultaneously faithful, stable across training perturbations, precise, generalizable across cohorts, and reproducible in its relative rank. Site (d=0) denotes the source cohort and (d\geq1) an external cohort; (t=1,\ldots,T) indexes patient-cluster bootstrap training replicates; (i=1,\ldots,n) indexes audit queries; and (j=1,\ldots,p) indexes features. The signed candidate attribution is \(\bar\phi_{dtj}(x_i)\), its independent exact or higher-budget reference is \(\widetilde\psi_j(f_{dt},x_i)\), and \(\widehat I_{dj}\) is the corresponding cohort-level mean absolute importance. Tolerances were fixed before evaluation at (a=b=0.5), (u=3), (\eta=0.02), and (\tau=0.8).

### 2.2 Component scores

The following block reproduces the frozen Section 2.2 equations verbatim from the preserved specification. The canonical `section2-warrant-index-taslak.md` was unavailable in this workspace; this provenance issue is recorded in the results ledger and requires a final source comparison before submission.

```text
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
```

Here, (F) measures agreement with an independent reference explainer; (S) measures attribution stability across unordered pairs of training replicates; (U) penalizes wide conditional attribution-importance intervals; (G) compares absolute importance across cohorts; and (G^{signed}) detects signed mechanism shifts that magnitude-only importance can conceal. The primary cross-cohort component was (G^{claim}=\min(G,G^{signed})). Features without common external support were marked not assessable rather than assigned an extremal score. For (R), (o_{dt,j\ell}\in\{-1,0,1\}) is the (\eta)-tolerant pairwise order of features (j) and \(\ell\); persistent ties count as reproducible but do not establish rank identifiability.

### 2.3 Warrant Index and weakest-link uniqueness

The component vector was (z_j=(F_j,S_j,U_j,G_j^{claim},R_j)\in[0,1]^5). The frozen aggregation and decision rules were:

```text
WI_j = min_k z_{jk}                              (Eq. 2.3.1)
Adequate_j = 1{z_jk > tau_k  for all k}           (Eq. 2.3.2)
WI_claim = min_{j in J*} WI_j                     (Eq. 2.3.3, J* önceden kayıtlı feature seti)
```

The minimum is bounded, monotone, idempotent, continuous, and fully non-compensatory. The Section 2.3 uniqueness result follows from the conjunction of monotonicity, idempotency, and the full weakest-link requirement (A(z)\leq z_k) for every criterion (k). Let (m=\min_k z_k). The weakest-link requirement gives (A(z)\leq m). Because (m\mathbf 1\leq z), monotonicity and idempotency give (m=A(m\mathbf 1)\leq A(z)). Hence (A(z)=m). Thus, under these requirements, the minimum is the unique aggregator; geometric or additive alternatives relax full non-compensation and were not primary analyses.

### 2.4 Statistical certification

Point adequacy is not statistical certification. For each feature and component, we tested the one-sided null that the population component score was at most 0.8. Because adequacy requires every component, component-level p-values were combined by the intersection-union rule reproduced from the frozen specification:

```text
Joint p-value: p_joint_j = max_k p_jk  (intersection-union test)
```

The preregistered family contained all 77 assessable features, without filtering on observed WI. Benjamini-Hochberg correction at (\alpha=0.05) was the primary feature-level false-discovery-rate analysis; Benjamini-Yekutieli correction was the arbitrary-dependence robustness analysis.

The final `recursive_studentized` procedure decomposed nonsmooth minimum scores into their leaves. F, U, (G), and (G^{signed}) leaves used one-sided studentized bootstrap tests; S and R leaves used basic bootstrap tests. At each site or branch minimum, leaf p-values were combined by their maximum, and the five component p-values were again combined by their maximum. Certification used 4,999 bootstrap draws.

We calibrated this inferential layer before full-family testing. The calibration battery contained 2,000 independent outer simulations, 999 inner resamples, `single_active` and `tie` boundary configurations, five components, and nominal alpha levels 0.01, 0.05, and 0.10. Thirty checks per model/method were evaluated with Holm-corrected one-sided exact binomial tests. A direct-root method was anti-conservative for RSF N=20 in selected F, U, and (G^{claim}) cells. Of the remedies tested, `recursive_studentized` passed, whereas `recursive_basic` and `m_out_of_n_root` did not. Fresh N=50 calibrations then produced 0/30 Holm-flagged violations for both RSF and DeepSurv. CoxPH also showed no Holm-corrected violation, and its basic-root and recursive-studentized implementations produced the same eight certificates.

## 3. Experimental design and data

### 3.1 Cohorts, split, and common audit panel

The source data were TCGA-BRCA, using the existing 673-row training and 169-row test split with 616 selected inputs: 200 RNA, 200 miRNA, 200 methylation, and 16 clinical variables. The stored CoxPH model used 614 inputs after excluding two separation columns. Twenty-three patients occurred in both original train and test files; these overlaps were removed from the audit, and one sample per remaining test patient was retained, yielding 141 disjoint TCGA audit patients.

METABRIC was used as the external explanation cohort. The shared, fully observed panel contained 77 features: 72 RNA features, age at diagnosis, and four stage indicators. A total of 1,466 METABRIC patients had complete query vectors for these 77 variables. Cohort-specific standardization was used to align the prepared feature spaces. METABRIC outcome labels were not used for fitting, tuning, or an external-survival-performance claim.

CoxPH used the full 141 TCGA and 1,466 METABRIC query distributions for cross-cohort attribution comparisons. The computationally heavier RSF and DeepSurv audits used deterministic subsets of 12 TCGA and 12 METABRIC queries. Their results are therefore conditional on these fixed queries.

### 3.2 Model classes

The linear model was a penalized Cox proportional-hazards model implemented with `lifelines.CoxPHFitter` and penalizer 0.1. Its held-out test C-index was 0.7065. For this model, log-partial-hazard Shapley values were available in closed form as (\beta_j(x_j-\text{background}_j)); the maximum checked reconstruction error was (2.84\times10^{-14}). The CoxPH audit used this same exact closed-form tensor as both candidate and reference. Its F component therefore measures numerical self-consistency of the closed-form computation rather than agreement between independent attribution estimators; RSF instead used independent permutation budgets, while DeepSurv used separately budgeted nested backgrounds.

The tree model was `sksurv.ensemble.RandomSurvivalForest` with 500 trees, `min_samples_split=10`, `min_samples_leaf=5`, `max_features="sqrt"`, and random state 42. Its test C-index was 0.7148. Because the installed SHAP version did not support this estimator through TreeExplainer, we used model-agnostic permutation Shapley estimates. Candidate estimates used eight permutation paths and independent references used 128 paths.

The neural model was a PyTorch DeepSurv-type Cox network. Separate encoders mapped RNA 200–128–32, miRNA 200–64–16, methylation 200–128–32, and clinical 16–32–8; the concatenated 88-dimensional representation passed through a 64–16–1 head. The network used batch normalization, ReLU activations, dropout 0.3, AdamW with learning rate 0.001 and weight decay 0.0001, at most 300 epochs, patience 30, and seed 42. The selected checkpoint was epoch 26 and its test C-index was 0.7457.

### 3.3 Nonlinear-model explanation games and backgrounds

The external cohort supplied only 77 of 616 model inputs. We therefore defined a 77-feature explanation subgame for RSF and DeepSurv and integrated the remaining 539 inputs over sampled TCGA training backgrounds. A fixed-mean RSF background was rejected during design validation: across eight probe features it changed (G^{claim}) by as much as 0.250446 relative to sampled background. Sixteen deterministic, unique TCGA training patients were consequently used for candidate background integration. Increasing the RSF reference budget from 128 to 256 paths changed F by at most 0.001367, supporting the retained 8/128 path design.

For DeepSurv, DeepExplainer was technically supported and near-additive. The approved design used 16 nested TCGA background patients for the candidate and 128 for the reference, making background size the computation-budget axis. GradientExplainer required 256/2,048 samples before satisfying the prespecified relative additivity criterion, whereas DeepExplainer k=16/128 achieved effectively exact reconstruction at lower cost. Two of the N=50 reference calls exceeded SHAP's internal fixed absolute 0.01 assertion. The same attributions were recomputed without that internal assertion and retained only because they passed the prespecified relative additivity error limit of 5%; the maximum observed reference-relative error was 0.8075%. No bootstrap replicate was skipped or replaced.

### 3.4 Bootstrap design

All final model audits used 50 patient-cluster bootstrap training replicates with seed 20260917. Samples belonging to one patient were resampled together. The main fitted artifacts were not overwritten; bootstrap models were trained separately with the original model class and hyperparameters. CoxPH convergence remediation followed a fixed sequence: remove the separation column identified by `lifelines`, then, if needed, use `max_steps=200` and `step_size=0.25` while retaining the same sample and penalizer. Thirteen CoxPH replicates required remediation and 11 required the numerical retry; no final replicate retained a warning.

## 4. Results

### 4.1 CoxPH

At six bootstrap replicates, 27 of 77 features had point WI above 0.8. Expressed on the common denominator, 27/77 features passed at N=6, compared with 22/77 at N=50, and 15 of the original 27 remained above threshold. The median bootstrap-of-bootstrap IQR width for R decreased from 0.03038 to 0.01931. In the final N=50 audit, R was the point-estimate bottleneck for 76 features and (G^{claim}) for one. Full-family recursive-studentized inference certified eight features by BH and the same eight by BY.

### 4.2 Random survival forest

The initial N=2 comparison yielded 43/77 under fixed-mean background and 35/77 under sampled background. With the approved sampled design, 23/77 features passed at N=20 and 20/77 at N=50; 16 of the N=20 winners were retained, seven were lost, and four new features crossed the point threshold. The final bottleneck distribution was R=72 and (G^{claim}=5). Full-family inference certified 10/77 features by BH and 10/77 by BY.

The N=20 winner set showed a selected-subset maturation effect. R decreased in 19 of 23 selected features, with mean change -0.01569 and paired Wilcoxon (p=0.00349). Across all 77 features, however, the mean change was only -0.000575, with 39 decreases and 38 increases and Wilcoxon (p=0.74329). Thus, the evidence supports winner's curse among initially selected features rather than a general collapse of R.

### 4.3 DeepSurv

Under the approved DeepExplainer design, 11/77 features passed at N=2. At N=50, only 8/77 passed and R was the bottleneck for all 77 features. Membership was highly unstable: of the 77 features, 59 passed at neither depth, 10 were lost, seven were new at N=50, and only one was retained. The fresh DeepSurv N=50 null gate passed all 30 Holm-corrected controls.

Despite passing the calibration gate, none of the eight point-passing DeepSurv features survived BH correction, and none survived BY correction. Their mean WI margin above 0.8 was 0.00931 and the median was 0.00947; four margins were below 0.01, four were between 0.01 and 0.03, and none lay between 0.05 and 0.10. This zero-certification result is reported without post hoc threshold or method changes.

### 4.4 Three-model comparison

| Model | Point WI≥0.8 | BH certified | BY certified | BH survival among point passers |
|---|---:|---:|---:|---:|
| CoxPH N=50 | 22/77 | 8/77 | 8/77 | 36.4% |
| RSF N=50 | 20/77 | 10/77 | 10/77 | 50.0% |
| DeepSurv N=50 | 8/77 | 0/77 | 0/77 | 0% |

R dominated all three final audits: it was the bottleneck for 76/77 CoxPH features, 72/77 RSF features, and 77/77 DeepSurv features. Among point passers, the proportion with a WI margin between 0.05 and 0.10 was 4.55% for CoxPH, 45.0% for RSF, and 0% for DeepSurv. DeepSurv's point-pass margin distribution was closer to CoxPH than to RSF by Wasserstein distance, but this comparison is exploratory because DeepSurv contributed only eight point passers.

The apparently favorable RSF BH-survival fraction was not explained by weaker inference. Earlier matched-margin analysis found smaller, rather than larger, CoxPH p-values in the low-margin band. We also directly tested whether the ranking universe caused the difference: recomputing CoxPH R within 77 rather than 614 competitors changed R by only +0.00283 on average among the same 22 point passers. The Wasserstein distance to RSF increased from 0.05857 to 0.06076 for those 22 features and from 0.02174 to 0.02381 across all 77 matched features. The ranking-universe hypothesis was therefore rejected.

## 5. Discussion

The central empirical result is the recurrence of rank reproducibility as the weakest component across three structurally different survival models. This pattern was not confined to a linear Cox representation, an ensemble tree ranking, or a neural risk function. R accounted for 76/77, 72/77, and 77/77 bottlenecks in CoxPH, RSF, and DeepSurv, respectively. High fidelity, stability, uncertainty, or cross-cohort scores could not compensate for this weakness because the Warrant Index uses the minimum.

The DeepSurv result quantifies the separation between predictive discrimination and explanation warrant. DeepSurv had the highest stored test C-index (0.7457 versus 0.7148 for RSF and 0.7065 for CoxPH), yet it produced the fewest point-passing explanations and no BH- or BY-certified feature. This does not imply that neural survival models are intrinsically less explainable. It shows that predictive performance alone did not establish the joint explanation claim under this audit design, tolerance set, query sample, and explainer-background protocol.

Increasing bootstrap depth exposed selection instability in every model family. CoxPH point passers declined from 27 to 22, RSF from 23 to 20, and DeepSurv from 11 to eight. RSF supplied the clearest formal winner's-curse evidence: the selected N=20 subset declined significantly, whereas the full family did not. DeepSurv showed an even more unstable selected set, with only one of 11 N=2 winners retained, but N=2 contains only one unordered replicate pair and should be interpreted as a design pilot rather than a precise baseline.

The difference between CoxPH and RSF FDR survival remains incompletely explained. It was not attributable to weaker CoxPH p-values at matched low margins, and the 614-versus-77 R ranking universe did not account for it. The remaining difference may reflect model-specific attribution geometry, query conditioning, nonlinear background integration, or other dependencies not isolated here. We leave this as an open empirical question rather than selecting a post hoc explanation.

The calibration sequence also illustrates why inference was treated as a separate methodological layer. A plausible direct bootstrap was anti-conservative for several RSF boundary cases. Recursive leaf testing with studentization and IUT maximum combination resolved the tested violations and was freshly revalidated for N=50 RSF and DeepSurv. Passing these finite simulation gates supports use under the tested nuisance structures; it is not a universal proof of calibration.

## 6. Limitations

First, external generalizability was assessable only in the 77-feature TCGA–METABRIC intersection. The nonlinear-model results describe a 77-feature subgame with 539 other inputs background-integrated, not full 616-feature Shapley games. miRNA and methylation explanations could not receive a complete cross-cohort WI in this study.

Second, the empirical analysis covers one cancer type. Cohort-specific standardization, platform differences, mixed TCGA sample types, and clinical-variable harmonization limit biological interpretation. METABRIC was an explanation-distribution cohort; the study did not establish transported survival performance.

Third, RSF and DeepSurv inference is conditional on fixed sets of 12 TCGA and 12 METABRIC queries. These samples are much smaller than the CoxPH audit distributions of 141 and 1,466 patients, and query-sampling uncertainty is not included in the final p-values.

Fourth, the final DeepSurv protocol used deterministic DeepExplainer with nested k=16/128 backgrounds. GradientExplainer was examined during design calibration, but the final N=50 and certification results were not replicated with GradientExplainer, KernelSHAP, or other explanation families. The findings therefore do not establish explainer-method invariance.

The interpretation of F also differs across model classes. CoxPH compared the exact closed-form attribution tensor with itself, whereas the nonlinear audits compared estimates produced at distinct computation budgets. CoxPH F is consequently a numerical implementation check and is not directly comparable to the estimator-agreement test performed for RSF and DeepSurv.

Fifth, 50 training bootstrap replicates were the computational upper bound used here. Although increasing T reduced some rank-estimation uncertainty, finite-replicate variability remains, and the same bootstrap data were used to estimate the score and its inferential distribution.

Sixth, tolerances (a=b=0.5), (u=3), (\eta=0.02), and (\tau=0.8) were frozen experimental settings rather than clinically validated cutoffs. Results may change under independently justified tolerances.

Finally, the canonical Section 2 manuscript source was unavailable while this draft was generated. The equations were copied verbatim from the preserved fixed-specification prompt and signed-extension report, but they must be checked against the canonical source before submission.

## 7. Conclusion

The Warrant Index provided a non-compensatory audit across linear, tree-ensemble, and neural survival models. Rank reproducibility was the dominant limitation in every class, and increasing bootstrap depth reduced the number or stability of initially adequate features. Full-family inference substantially narrowed point-estimate claims: CoxPH retained eight BH/BY-certified features, RSF retained 10, and DeepSurv retained none. The DeepSurv result shows that comparatively strong predictive discrimination can coexist with an absence of statistically certified feature-level explanation warrants under the evaluated conditions. These findings support reporting component profiles, replication depth, calibration evidence, and multiplicity-adjusted decisions alongside any aggregate explanation-credibility score.

## 8. Availability

The implementation is provided as the `warrantxai` Python package, with experiment scripts and machine-readable result artifacts for the analyses reported here. Local installation uses the package configuration in `pyproject.toml`; optional real-data dependencies are defined by the project.

**[TODO: Claude/kullanıcı ile netleştirilecek — public repo + lisans kararı bekliyor]**
