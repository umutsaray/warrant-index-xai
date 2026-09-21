# Faz 7.3 — Sertifikalı özellik yönleri

Yön, saklı candidate attribution'ların TCGA tarafındaki aritmetik ortalamasının işaretidir. CoxPH için tüm 153 fold-0 audit hastası ve 50 replika kullanıldı. RSF artefaktı yalnız dondurulmuş 12 TCGA query için attribution sakladığından RSF yönü bu 12-query alt kümesine ve 50 replikaya koşulludur; 153 hastaya genişletmek yeni explainer çağrısı gerektireceği için yapılmadı. Modeller-arası işaret kontrolünde query farkını gidermek için CoxPH'nin saklı 153-hasta tensöründen RSF ile aynı deterministik 12 satır ayrıca seçildi. Bu seçim yeni attribution hesabı değildir. DeepSurv BH kümesi boş olduğundan yön satırı yoktur. G_signed tek başına yönü saklamaz; bu nedenle ham attribution tensörleri okundu. Tam sıfır ortalama, pozitif veya negatif yön yerine `sıfır` olarak tutuldu.

## BH-sertifikalı özellikler

| Model | Feature | Feature ID | Mean signed attribution | Direction | Averaging scope |
| --- | --- | --- | --- | --- | --- |
| CoxPH-v2 | AJCC Stage III | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | -0.0051066449 | negatif | 153 TCGA audit patients × 50 replicas |
| CoxPH-v2 | Age at index | clin_numeric__age_at_index.demographic | 0.0027330933 | pozitif | 153 TCGA audit patients × 50 replicas |
| CoxPH-v2 | CA1 | rna_ENSG00000133742.14 | -2.3586326e-05 | negatif | 153 TCGA audit patients × 50 replicas |
| CoxPH-v2 | RBMY1F | rna_ENSG00000169800.14 | -0.0011567858 | negatif | 153 TCGA audit patients × 50 replicas |
| CoxPH-v2 | KRTAP19-2 | rna_ENSG00000186965.5 | -0.00025621341 | negatif | 153 TCGA audit patients × 50 replicas |
| CoxPH-v2 | OR13C3 | rna_ENSG00000204246.3 | 5.3203893e-05 | pozitif | 153 TCGA audit patients × 50 replicas |
| CoxPH-v2 | KRTAP19-7 | rna_ENSG00000244362.3 | -0.0015624577 | negatif | 153 TCGA audit patients × 50 replicas |
| CoxPH-v2 | OR5M10 | rna_ENSG00000254834.4 | -0.0015379201 | negatif | 153 TCGA audit patients × 50 replicas |
| RSF | AJCC Stage I | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | 0.00011641195 | pozitif | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | AJCC Stage II | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | 0.0023411669 | pozitif | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | AJCC Stage III | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 0.0068889614 | pozitif | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | AJCC Stage IV | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | 0 | sıfır | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | RBMY1F | rna_ENSG00000169800.14 | 0 | sıfır | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | OR4C46 | rna_ENSG00000185926.1 | 0 | sıfır | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | KRTAP19-2 | rna_ENSG00000186965.5 | 0 | sıfır | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | OR13C3 | rna_ENSG00000204246.3 | 0 | sıfır | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | ACCSL | rna_ENSG00000205126.2 | 0.00084130624 | pozitif | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | OR6C1 | rna_ENSG00000205330.4 | 0 | sıfır | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | KRTAP19-7 | rna_ENSG00000244362.3 | 0 | sıfır | stored deterministic 12 TCGA queries × 50 replicas |
| RSF | OR5M10 | rna_ENSG00000254834.4 | 0 | sıfır | stored deterministic 12 TCGA queries × 50 replicas |

## CoxPH–RSF BH kesişiminde işaret tutarlılığı

| Feature | Feature ID | CoxPH full-153 mean | CoxPH full direction | CoxPH matched-12 mean | CoxPH matched direction | RSF mean | RSF direction | Same sign on matched queries |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AJCC Stage III | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | -0.0051066449 | negatif | -0.022718597 | negatif | 0.0068889614 | pozitif | False |
| RBMY1F | rna_ENSG00000169800.14 | -0.0011567858 | negatif | -0.0011567858 | negatif | 0 | sıfır | False |
| KRTAP19-2 | rna_ENSG00000186965.5 | -0.00025621341 | negatif | 0.00026846678 | pozitif | 0 | sıfır | False |
| OR13C3 | rna_ENSG00000204246.3 | 5.3203893e-05 | pozitif | -6.6653038e-05 | negatif | 0 | sıfır | False |
| KRTAP19-7 | rna_ENSG00000244362.3 | -0.0015624577 | negatif | 0.00081768907 | pozitif | 0 | sıfır | False |
| OR5M10 | rna_ENSG00000254834.4 | -0.0015379201 | negatif | 0.00085417728 | pozitif | 0 | sıfır | False |

Aynı 12 query üzerinde aynı işaret sayısı **0/6**. Bu yalnız sayısal yön tutarlılığıdır; biyolojik yorum yapılmamıştır.
