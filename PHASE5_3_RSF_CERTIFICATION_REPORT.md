# Faz 5.3A — RSF tam-aile recursive_studentized sertifikasyonu

> **Koşulluluk uyarısı:** RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

## Kapsam ve yöntem — RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

Faz 4.2'deki **77 feature'ın tamamı** aileye alındı; WI filtresi uygulanmadı. Gate'i geçen `recursive_studentized` yöntemi B=4999 ile kullanıldı: F/U/G_magnitude/G_signed yapraklarında studentized bootstrap, S/R yapraklarında basic bootstrap ve formül ağacındaki her `min` için recursive `max(p)` IUT. Bileşenler de `p_joint=max_k(p_jk)` ile birleştirildi. BH α=0.05 birincil, BY bağımlılığa dayanıklı sağlamlık analizidir.

## Nokta tahmini karşılaştırması — RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

| inference_scope | criterion | count |
| --- | --- | --- |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | assessed_full_family_m | 77 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | point_WI_ge_0_8 | 23 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | BH_certified_alpha_0_05 | 15 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | BY_certified_alpha_0_05 | 12 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | point_and_BH | 15 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | point_only_not_BH | 8 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | BH_only_not_point | 0 |

## BH sertifikalı feature'lar — RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

| inference_scope | feature | block | WI | p_joint | q_BH | q_BY | certified_BY_alpha_0_05 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000006116.4 | rna | 0.8631 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000169800.14 | rna | 0.8565 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000109158.11 | rna | 0.8631 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000004846.16 | rna | 0.8696 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000118407.15 | rna | 0.8506 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000116663.11 | rna | 0.8613 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000134531.10 | rna | 0.8545 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000187758.8 | rna | 0.8431 | 0.0088 | 0.0452 | 0.2226 | False |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000204246.3 | rna | 0.8565 | 0.0022 | 0.0130 | 0.0642 | False |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000163545.11 | rna | 0.8703 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000250334.6 | rna | 0.8652 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | clin | 0.8583 | 0.0060 | 0.0330 | 0.1626 | False |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | clin | 0.8648 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | clin | 0.8593 | 0.0002 | 0.0013 | 0.0063 | True |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | clin | 0.8565 | 0.0002 | 0.0013 | 0.0063 | True |

## Nokta-yeterli olup BH ile sertifikalanmayanlar — RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

| inference_scope | feature | block | WI | bottleneck | p_joint | q_BH |
| --- | --- | --- | --- | --- | --- | --- |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000178522.14 | rna | 0.8002 | R | 0.6212 | 1.0000 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000182591.6 | rna | 0.8574 | R | 0.1256 | 0.6044 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000198099.9 | rna | 0.8173 | R | 0.4398 | 1.0000 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000133800.9 | rna | 0.8259 | R | 0.3434 | 1.0000 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000187531.14 | rna | 0.8095 | R | 0.5046 | 1.0000 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000246430.7 | rna | 0.8022 | R | 0.6006 | 1.0000 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000121871.4 | rna | 0.8024 | R | 0.5990 | 1.0000 |
| RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna_ENSG00000187486.5 | rna | 0.8036 | R | 0.5574 | 1.0000 |

Bu gruptaki 8 feature'ın tamamında nokta-tahmin bottleneck'i `R` bileşenidir.

Tam feature-bileşen tablosu `results/phase53_rsf_certification/feature_certification.csv` içindedir. q-değeri feature doğruluk olasılığı değildir.

## Dürüstlük kontrolleri — RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

Formül SHA256 değerleri Faz 5.2 ile aynıdır: **True**. Kimlik-sızıntısı kontrolü geçti. Model yeniden fit edilmedi, explainer çağrılmadı. Gerçek wall-clock **3.937 saniye**.
