# Faz 5.6 ADIM 3 — RSF N=50 FDR sertifikasyonu

> **Koşulluluk uyarısı:** RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

N=50 null gate'i geçtikten sonra önceden tanımlı 77 feature ailesinin tamamı B=4999 `recursive_studentized` ile test edildi. BH alpha=0.05 birincil, BY sağlamlık analizidir.

| inference_scope | criterion | count |
| --- | --- | --- |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | assessed_full_family_m | 77 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | point_WI_ge_0_8 | 20 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | BH_certified_alpha_0_05 | 10 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | BY_certified_alpha_0_05 | 10 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | point_and_BH | 10 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | point_only_not_BH | 10 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | BH_only_not_point | 0 |

## BH sertifikalı feature'lar

| inference_scope | feature | block | WI | p_joint | q_BH | q_BY | certified_BY_alpha_0_05 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000169800.14 | rna | 0.8549 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000109158.11 | rna | 0.8622 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000004846.16 | rna | 0.8670 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000116663.11 | rna | 0.8524 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000187758.8 | rna | 0.8430 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000204246.3 | rna | 0.8549 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | clin | 0.8571 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | clin | 0.8657 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | clin | 0.8635 | 0.0002 | 0.0015 | 0.0076 | True |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | clin | 0.8549 | 0.0002 | 0.0015 | 0.0076 | True |

## Nokta-yeterli olup BH ile sertifikalanmayanlar

| inference_scope | feature | block | WI | bottleneck | p_joint | q_BH |
| --- | --- | --- | --- | --- | --- | --- |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000006116.4 | rna | 0.8265 | R | 0.1064 | 0.6827 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000182591.6 | rna | 0.8478 | R | 0.0472 | 0.3304 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000108639.8 | rna | 0.8101 | R | 0.3390 | 1.0000 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000134531.10 | rna | 0.8096 | R | 0.3674 | 1.0000 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000205186.4 | rna | 0.8034 | R | 0.5030 | 1.0000 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000163145.13 | rna | 0.8176 | R | 0.2978 | 1.0000 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000121871.4 | rna | 0.8284 | R | 0.1274 | 0.7546 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000169764.16 | rna | 0.8088 | R | 0.4226 | 1.0000 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000163545.11 | rna | 0.8236 | R | 0.1632 | 0.8378 |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000187486.5 | rna | 0.8230 | R | 0.1590 | 0.8378 |

Tam tablo `results/phase56_rsf_n50_certification/feature_certification.csv` içindedir. q-değeri feature doğruluk olasılığı değildir. Formül SHA256 değişmedi: **True**. Model/explainer çağrısı yapılmadı; kimlik-sızıntısı kontrolü geçti. Gerçek wall-clock **10.495 saniye**.
