# Faz 7.7 — Tekil-değer cancellation kontrolü

Metrikler, hiçbir agregasyon öncesinde saklı `(replica, query)` hücrelerinden hesaplandı. RSF için her feature 50×12=600 hücre; CoxPH-v2 için 50×153=7650 hücredir. Tarihsel kolon adı “All 600 cells” olmakla birlikte CoxPH satırlarının gerçek paydası `Cell denominator=7650` olarak açıkça saklanmıştır.

| Fold | Model | Feature | Feature ID | Phase 7.6 label | All cells exactly zero | M_j mean absolute attribution | Non-zero cell count | Cell denominator | P_j non-zero frequency | Q_j max absolute attribution |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | CoxPH-v2 | KRTAP19-7 | rna_ENSG00000244362.3 | informative | False | 0.0031657696466263377 | 4743 | 7650 | 0.62 | 1.482915856371915 |
| 0 | CoxPH-v2 | RBMY1F | rna_ENSG00000169800.14 | informative | False | 0.0011567857749820404 | 3825 | 7650 | 0.5 | 0.0030487974452089355 |
| 0 | CoxPH-v2 | CA1 | rna_ENSG00000133742.14 | informative | False | 0.010199643124597153 | 7650 | 7650 | 1.0 | 1.0721829649566568 |
| 0 | CoxPH-v2 | OR13C3 | rna_ENSG00000204246.3 | informative | False | 0.004172711044193061 | 7191 | 7650 | 0.94 | 2.2197815114596002 |
| 0 | CoxPH-v2 | KRTAP19-2 | rna_ENSG00000186965.5 | informative | False | 0.000786128236826272 | 4743 | 7650 | 0.62 | 0.2538327957353823 |
| 0 | CoxPH-v2 | OR5M10 | rna_ENSG00000254834.4 | informative | False | 0.0032016118040726744 | 4896 | 7650 | 0.64 | 0.7486836258977378 |
| 0 | CoxPH-v2 | Age at index | clin_numeric__age_at_index.demographic | informative | False | 0.13579757388852318 | 7650 | 7650 | 1.0 | 0.5744230924466924 |
| 0 | CoxPH-v2 | AJCC Stage III | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | informative | False | 0.09505526685596867 | 7650 | 7650 | 1.0 | 0.3723541574435621 |
| 0 | RSF | KRTAP19-7 | rna_ENSG00000244362.3 | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 0 | RSF | OR6C1 | rna_ENSG00000205330.4 | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 0 | RSF | OR4C46 | rna_ENSG00000185926.1 | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 0 | RSF | RBMY1F | rna_ENSG00000169800.14 | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 0 | RSF | OR13C3 | rna_ENSG00000204246.3 | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 0 | RSF | KRTAP19-2 | rna_ENSG00000186965.5 | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 0 | RSF | ACCSL | rna_ENSG00000205126.2 | informative | False | 0.004556248492826611 | 574 | 600 | 0.9566666666666667 | 0.019237500000000463 |
| 0 | RSF | OR5M10 | rna_ENSG00000254834.4 | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 0 | RSF | AJCC Stage I | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | informative | False | 0.00023856374007937017 | 44 | 600 | 0.07333333333333333 | 0.012000000000000455 |
| 0 | RSF | AJCC Stage II | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | informative | False | 0.0024314139239539007 | 150 | 600 | 0.25 | 0.05540912316849822 |
| 0 | RSF | AJCC Stage III | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | informative | False | 0.006888961370342453 | 100 | 600 | 0.16666666666666666 | 0.500118625288111 |
| 0 | RSF | AJCC Stage IV | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 1 | RSF | OR13C3 | rna_ENSG00000204246.3 | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |
| 1 | RSF | AJCC Stage IV | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | degenerate | True | 0.0 | 0 | 600 | 0.0 | 0.0 |

Ana assertion geçti: fold-0'daki sekiz ve fold-1'deki iki ayrı degenerate RSF sonucu dahil **10/10** satırda bütün 600 tekil değer exact `0.0`; dolayısıyla M_j=P_j=Q_j=0. Karşılaştırma satırlarının tamamında informative özellikler için exact-zero assertion `False`, M_j>0 ve P_j>0'dır. Süre **0.096 sn**.
