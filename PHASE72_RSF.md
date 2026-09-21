# Faz 7.2 — RSF N=50 temiz audit

Hasta-gruplu fold-0 referans RSF modeli ve temiz 27-feature TCGA–METABRIC paneli üzerinde 50 patient-cluster bootstrap replika tamamlandı. **17/27** feature nokta tahminde WI≥0.8 geçti.

| block | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- |
| clin | 5 | 1.0000 | 0.9846 | 0.9883 | 0.9776 | 0.8901 | 0.8901 |
| rna | 22 | 1.0000 | 0.9714 | 0.9727 | 0.9410 | 0.8471 | 0.8304 |

Bottleneck dağılımı: `{'R': 24, 'G_claim': 3}`.

Bu sonuçlar sabit deterministik 12+12 query kümesine, saklanmamış permütasyon yollarına, sampled inactive background k=16 ve 8/128 permütasyon bütçesine koşulludur. Formül SHA256 değişmedi; hasta kimliği dışa aktarılmadı. Birikimli ölçülen hesap süresi **883.9 sn**, bu invocation wall-clock **885.2 sn**.


## ADIM 5 — Tam-aile sertifikasyon

Taze null kapısı geçtikten sonra m=27 ailenin tamamında B=4999 recursive_studentized uygulandı. Nokta/BH/BY sayıları **17/12/12**.

| feature | WI | p_joint | q_BH | q_BY |
| --- | --- | --- | --- | --- |
| rna_ENSG00000244362.3 | 0.9272 | 0.0002 | 0.0005 | 0.0019 |
| rna_ENSG00000205330.4 | 0.9038 | 0.0012 | 0.0027 | 0.0105 |
| rna_ENSG00000185926.1 | 0.9272 | 0.0002 | 0.0005 | 0.0019 |
| rna_ENSG00000169800.14 | 0.9272 | 0.0002 | 0.0005 | 0.0019 |
| rna_ENSG00000204246.3 | 0.9278 | 0.0002 | 0.0005 | 0.0019 |
| rna_ENSG00000186965.5 | 0.9272 | 0.0002 | 0.0005 | 0.0019 |
| rna_ENSG00000205126.2 | 0.9246 | 0.0002 | 0.0005 | 0.0019 |
| rna_ENSG00000254834.4 | 0.9272 | 0.0002 | 0.0005 | 0.0019 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | 0.9265 | 0.0002 | 0.0005 | 0.0019 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | 0.9240 | 0.0002 | 0.0005 | 0.0019 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 0.9116 | 0.0002 | 0.0005 | 0.0019 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | 0.9272 | 0.0002 | 0.0005 | 0.0019 |
