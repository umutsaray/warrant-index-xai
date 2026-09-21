# Faz 6.2 — DeepSurv N=50, sertifikasyon ve üç-model karşılaştırması

> **DeepSurv koşulluluk uyarısı:** DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds.

## Sertifikasyon özeti

| inference_scope | model | assessed | point_WI_ge_0_8 | BH | BY | BH_survival_fraction |
| --- | --- | --- | --- | --- | --- | --- |
| CoxPH N=50 full-audit-query recursive_studentized certification. | CoxPH N=50 | 77.0000 | 22.0000 | 8.0000 | 8.0000 | 0.3636 |
| RSF N=50 fixed 12+12 query and stored 8/128 attribution certification. | RSF N=50 | 77.0000 | 20.0000 | 10.0000 | 10.0000 | 0.5000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | DeepSurv N=50 | 77.0000 | 8.0000 | 0.0000 | 0.0000 | 0.0000 |

DeepSurv null gate'i 30/30 kontrolde geçti. Tam 77-feature ailede 8 feature nokta WI≥0.8 eşiğini geçti; BH ve BY ile sertifikalanan feature olmadı.

ADIM 1 sırasında SHAP'ın sabit mutlak 0.01 iç assertion'ı iki reference çağrısında tetiklendi. Replikalar atlanmadı veya değiştirilmedi; attributionlar aynı model/veriyle yeniden hesaplandı ve yalnız önceden tanımlı göreli additivity hatası ≤%5 kapısını geçtikleri için kabul edildi. N=50'de maksimum reference göreli hata %0.807 oldu.

## Nokta-geçenlerde marjin dağılımı

| inference_scope | model | point_features | mean_margin | median_margin | n_0_0.01 | n_0.01_0.03 | n_0.03_0.05 | n_0.05_0.10 | share_0.05_0.10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH N=50 full-audit-query recursive_studentized certification. | CoxPH N=50 | 22.0000 | 0.0293 | 0.0225 | 5.0000 | 8.0000 | 7.0000 | 1.0000 | 0.0455 |
| RSF N=50 fixed 12+12 query and stored 8/128 attribution certification. | RSF N=50 | 20.0000 | 0.0387 | 0.0454 | 3.0000 | 6.0000 | 2.0000 | 9.0000 | 0.4500 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | DeepSurv N=50 | 8.0000 | 0.0093 | 0.0095 | 4.0000 | 4.0000 | 0.0000 | 0.0000 | 0.0000 |

| inference_scope | comparison | wasserstein_margin_distance |
| --- | --- | --- |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | DeepSurv_vs_CoxPH | 0.0200 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | DeepSurv_vs_RSF | 0.0294 |

Wasserstein uzaklığına göre DeepSurv nokta-marjin şekli **CoxPH** modeline daha yakındır. Bu karşılaştırma DeepSurv'de yalnız sekiz nokta-geçen feature bulunduğu için keşifseldir.

## Darboğaz dağılımı

| inference_scope | model | scope | bottleneck | count |
| --- | --- | --- | --- | --- |
| CoxPH N=50 full-audit-query recursive_studentized certification. | CoxPH N=50 | all_77 | R | 76.0000 |
| CoxPH N=50 full-audit-query recursive_studentized certification. | CoxPH N=50 | all_77 | G_claim | 1.0000 |
| CoxPH N=50 full-audit-query recursive_studentized certification. | CoxPH N=50 | point_passers | R | 21.0000 |
| CoxPH N=50 full-audit-query recursive_studentized certification. | CoxPH N=50 | point_passers | G_claim | 1.0000 |
| RSF N=50 fixed 12+12 query and stored 8/128 attribution certification. | RSF N=50 | all_77 | R | 72.0000 |
| RSF N=50 fixed 12+12 query and stored 8/128 attribution certification. | RSF N=50 | all_77 | G_claim | 5.0000 |
| RSF N=50 fixed 12+12 query and stored 8/128 attribution certification. | RSF N=50 | point_passers | R | 20.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | DeepSurv N=50 | all_77 | R | 77.0000 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | DeepSurv N=50 | point_passers | R | 8.0000 |

DeepSurv N=50'nin 77 feature'ının tamamında darboğaz R'dir. CoxPH ve RSF'de G_claim az sayıda feature'da darboğaz olurken baskın bileşen yine R'dir.

## N=2 → N=50 winner's-curse

Nokta-geçen sayı 11'den 8'e düştü. Geçişler `{'neither': 59, 'lost': 10, 'new': 7, 'retained': 1}`: N=2 kazananlarının yalnız biri korundu, onu düştü ve yedi yeni feature geçti. Sayı düzeyindeki düşüş winner's-curse ile uyumludur; feature üyeliğinin bu kadar değişmesi N=2 sıralama tahmininin çok oynak olduğunu da gösterir. Bu nedenle sonuç yalnız “11→8 düşüşü” olarak yorumlanmamalıdır.

### N=2'de geçip N=50'de düşenler

| inference_scope | feature | block | WI_n2 | WI | R_n2 | R | bottleneck_n2 | bottleneck |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000198670.12 | rna | 0.8092 | 0.7215 | 0.8092 | 0.7215 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000090402.8 | rna | 0.8158 | 0.6858 | 0.8158 | 0.6858 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000270069.1 | rna | 0.8224 | 0.7327 | 0.8224 | 0.7327 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000169432.18 | rna | 0.8026 | 0.7198 | 0.8026 | 0.7198 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000157617.17 | rna | 0.8092 | 0.6941 | 0.8092 | 0.6941 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000004846.16 | rna | 0.8026 | 0.7706 | 0.8026 | 0.7706 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000108639.8 | rna | 0.8224 | 0.7202 | 0.8224 | 0.7202 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000134531.10 | rna | 0.8158 | 0.7596 | 0.8158 | 0.7596 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000092010.15 | rna | 0.8026 | 0.7962 | 0.8026 | 0.7962 | R | R |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna_ENSG00000124440.16 | rna | 0.8224 | 0.7326 | 0.8224 | 0.7326 | R | R |

Gerçek analiz wall-clock **0.014 sn**. Model/explainer çağrısı yapılmadı. Formül SHA256 değişmedi ve kimlik-sızıntısı kontrolü geçti.

## Adım bazında wall-clock

| inference_scope | step | wall_clock_seconds | note |
| --- | --- | --- | --- |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | ADIM 1 N=50 | 216.6175 | includes stopped SHAP assertion attempt and checkpoint resume |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | ADIM 2 null gate | 23.1611 | 2000x999, no model calls |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | ADIM 3 FDR | 9.6914 | B=4999, no model calls |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | ADIM 4 comparison | 0.0137 | existing artifacts only |
