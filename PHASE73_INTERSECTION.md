# Faz 7.3 — CoxPH-v2 ↔ RSF sertifikalı özellik kesişimi

## BH kesişimi

Venn sayıları: **CoxPH-only=2, both=6, RSF-only=6**; birleşim=14.

CoxPH-only: Age at index, CA1  
RSF-only: AJCC Stage I, AJCC Stage II, AJCC Stage IV, OR4C46, ACCSL, OR6C1

| Feature | Feature ID | CoxPH WI | CoxPH p_joint | CoxPH q | RSF WI | RSF p_joint | RSF q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AJCC Stage III | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 0.8509 | 0.0002 | 0.0009 | 0.9116 | 0.0002 | 0.0005 |
| RBMY1F | rna_ENSG00000169800.14 | 0.8292 | 0.0002 | 0.0009 | 0.9272 | 0.0002 | 0.0005 |
| KRTAP19-2 | rna_ENSG00000186965.5 | 0.8291 | 0.0002 | 0.0009 | 0.9272 | 0.0002 | 0.0005 |
| OR13C3 | rna_ENSG00000204246.3 | 0.8205 | 0.0006 | 0.0023 | 0.9278 | 0.0002 | 0.0005 |
| KRTAP19-7 | rna_ENSG00000244362.3 | 0.8299 | 0.0002 | 0.0009 | 0.9272 | 0.0002 | 0.0005 |
| OR5M10 | rna_ENSG00000254834.4 | 0.8297 | 0.0002 | 0.0009 | 0.9272 | 0.0002 | 0.0005 |

## BY kesişimi

Venn sayıları: **CoxPH-only=2, both=6, RSF-only=6**; birleşim=14.

CoxPH-only: Age at index, CA1  
RSF-only: AJCC Stage I, AJCC Stage II, AJCC Stage IV, OR4C46, ACCSL, OR6C1

| Feature | Feature ID | CoxPH WI | CoxPH p_joint | CoxPH q | RSF WI | RSF p_joint | RSF q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AJCC Stage III | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 0.8509 | 0.0002 | 0.0035 | 0.9116 | 0.0002 | 0.0019 |
| RBMY1F | rna_ENSG00000169800.14 | 0.8292 | 0.0002 | 0.0035 | 0.9272 | 0.0002 | 0.0019 |
| KRTAP19-2 | rna_ENSG00000186965.5 | 0.8291 | 0.0002 | 0.0035 | 0.9272 | 0.0002 | 0.0019 |
| OR13C3 | rna_ENSG00000204246.3 | 0.8205 | 0.0006 | 0.0090 | 0.9278 | 0.0002 | 0.0019 |
| KRTAP19-7 | rna_ENSG00000244362.3 | 0.8299 | 0.0002 | 0.0035 | 0.9272 | 0.0002 | 0.0019 |
| OR5M10 | rna_ENSG00000254834.4 | 0.8297 | 0.0002 | 0.0035 | 0.9272 | 0.0002 | 0.0019 |

BH ve BY sertifikalı kümeleri iki modelde de **aynıdır**.
