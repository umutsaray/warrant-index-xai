# Faz 6.2 ADIM 3 — DeepSurv N=50 FDR sertifikasyonu

> **Koşulluluk uyarısı:** DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds.

N=50 null gate'i geçtikten sonra önceden tanımlı 77 feature ailesinin tamamı B=4999 `recursive_studentized` ile test edildi. BH alpha=0.05 birincil, BY sağlamlık analizidir.

| inference_scope | criterion | count |
| --- | --- | --- |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | assessed_full_family_m | 77 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | point_WI_ge_0_8 | 8 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | BH_certified_alpha_0_05 | 0 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | BY_certified_alpha_0_05 | 0 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | point_and_BH | 0 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | point_only_not_BH | 8 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | BH_only_not_point | 0 |

## BH sertifikalı feature'lar

BH ile sertifikalanan feature yok.

## Nokta-yeterli olup BH ile sertifikalanmayanlar

| inference_scope | feature | block | WI | bottleneck | p_joint | q_BH |
| --- | --- | --- | --- | --- | --- | --- |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000174899.12 | rna | 0.8072 | R | 0.4242 | 1.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000168917.9 | rna | 0.8144 | R | 0.1892 | 1.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000133800.9 | rna | 0.8158 | R | 0.1362 | 1.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000077616.11 | rna | 0.8077 | R | 0.4080 | 1.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000187531.14 | rna | 0.8138 | R | 0.3410 | 1.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000182173.13 | rna | 0.8113 | R | 0.3156 | 1.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000205186.4 | rna | 0.8009 | R | 0.5946 | 1.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000250334.6 | rna | 0.8034 | R | 0.5306 | 1.0000 |

Tam tablo `results/phase62_deepsurv_n50_certification/feature_certification.csv` içindedir. q-değeri feature doğruluk olasılığı değildir. Formül SHA256 değişmedi: **True**. Model/explainer çağrısı yapılmadı; kimlik-sızıntısı kontrolü geçti. Gerçek wall-clock **9.691 saniye**.
