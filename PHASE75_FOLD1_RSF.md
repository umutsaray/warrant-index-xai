# Faz 7.5 fold 1 — RSF N=50 temiz audit

Hasta-gruplu fold-1 referans RSF modeli ve temiz 80-feature TCGA–METABRIC paneli üzerinde 50 patient-cluster bootstrap replika tamamlandı. **10/80** feature nokta tahminde WI≥0.8 geçti.

| block | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- |
| clin | 5 | 1.0000 | 0.9976 | 0.9964 | 0.9957 | 0.8297 | 0.8297 |
| rna | 75 | 0.9997 | 0.9437 | 0.9568 | 0.9027 | 0.7012 | 0.6829 |

Bottleneck dağılımı: `{'R': 76, 'G_claim': 4}`.

Bu sonuçlar sabit deterministik 12+12 query kümesine, saklanmamış permütasyon yollarına, sampled inactive background k=16 ve 8/128 permütasyon bütçesine koşulludur. Formül SHA256 değişmedi; hasta kimliği dışa aktarılmadı. Birikimli ölçülen hesap süresi **1883.2 sn**, bu invocation wall-clock **1884.6 sn**.


## ADIM 5 — Tam-aile sertifikasyon

Taze null kapısı geçtikten sonra m=80 ailenin tamamında B=4999 recursive_studentized uygulandı. Nokta/BH/BY sayıları **10/7/7**.

| feature | WI | p_joint | q_BH | q_BY |
| --- | --- | --- | --- | --- |
| rna_ENSG00000178522.14 | 0.8281 | 0.0002 | 0.0023 | 0.0113 |
| rna_ENSG00000204246.3 | 0.8581 | 0.0002 | 0.0023 | 0.0113 |
| rna_ENSG00000004846.16 | 0.8268 | 0.0002 | 0.0023 | 0.0113 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | 0.8576 | 0.0002 | 0.0023 | 0.0113 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | 0.8516 | 0.0002 | 0.0023 | 0.0113 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 0.8521 | 0.0002 | 0.0023 | 0.0113 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | 0.8581 | 0.0002 | 0.0023 | 0.0113 |
