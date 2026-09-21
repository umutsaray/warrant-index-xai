# Faz 5.2A — CoxPH tam-aile istatistiksel sertifikasyonu

## Kapsam ve yöntem

Faz 3.2'de değerlendirilebilir bulunan **77 feature'ın tamamı** önceden tanımlı test ailesine alındı; `WI≥0.8` filtresiyle seçim yapılmadı. Her feature ve F/S/U/G_claim/R bileşeni için onaylanan N=50 replika-bazlı, null-merkezlenmiş bootstrap testi B=4999 ile çalıştırıldı. Bileşen p-değerleri `p_joint=max_k(p_jk)` intersection-union testiyle birleştirildi. BH q-değeri α=0.05 için birincil sonuç, keyfi feature bağımlılığına izin veren BY q-değeri sağlamlık kontrolüdür.

Çıkarım sabit 141 TCGA ve 1466 METABRIC sorgusuna, preprocessing'e ve exact Cox attribution tanımına koşulludur. Faz 5.1 kalibrasyonunda CoxPH N=50 için Holm-düzeltilmiş Tip-I ihlali görülmemiştir.

## Nokta tahmini ve sertifikasyon karşılaştırması

| criterion | count |
| --- | --- |
| assessed_full_family_m | 77 |
| point_WI_ge_0_8 | 22 |
| BH_certified_alpha_0_05 | 8 |
| BY_certified_alpha_0_05 | 8 |
| point_and_BH | 8 |
| point_only_not_BH | 14 |
| BH_only_not_point | 0 |

BH ile sertifikalanan feature'lar:

| feature | block | WI | p_joint | q_BH | q_BY | certified_BY_alpha_0_05 |
| --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000133742.14 | rna | 0.8412 | 0.0002 | 0.0019 | 0.0095 | True |
| rna_ENSG00000169800.14 | rna | 0.8331 | 0.0002 | 0.0019 | 0.0095 | True |
| rna_ENSG00000133800.9 | rna | 0.8367 | 0.0002 | 0.0019 | 0.0095 | True |
| rna_ENSG00000109158.11 | rna | 0.8335 | 0.0002 | 0.0019 | 0.0095 | True |
| rna_ENSG00000187758.8 | rna | 0.8371 | 0.0002 | 0.0019 | 0.0095 | True |
| rna_ENSG00000204246.3 | rna | 0.8371 | 0.0002 | 0.0019 | 0.0095 | True |
| clin_numeric__age_at_index.demographic | clin | 0.9503 | 0.0002 | 0.0019 | 0.0095 | True |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | clin | 0.8670 | 0.0002 | 0.0019 | 0.0095 | True |

BH birincil sonuçtur; BY sütunu dependence-robust sağlamlık kontrolüdür. FDR-düzeltilmiş q-değeri bir feature'ın doğru olma olasılığı değildir. Tam feature-bileşen tablosu `results/phase52_cox_certification/feature_certification.csv` içindedir.

## Dürüstlük kontrolleri

Formül modülü SHA256 değerleri Faz 5 pilotuyla aynıdır: **True**. Çıktılarda hasta veya örnek kimliği yoktur. Model yeniden fit edilmedi ve explainer çağrılmadı. Gerçek wall-clock süre **176.233 saniye** idi.
