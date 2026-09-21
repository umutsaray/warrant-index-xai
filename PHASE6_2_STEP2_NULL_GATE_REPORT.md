# Faz 6.2 ADIM 2 — DeepSurv N=50 null-kalibrasyon kapısı

> **Koşulluluk uyarısı:** DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds.

`recursive_studentized`, N=50 attribution varyans/bağımlılık yapısından türetilen nuisance ile aynı 2000×999, `single_active` + `tie`, üç alpha ve 30 kontrol içinde Holm düzeltmeli bataryadan geçirildi.

| model | method | controls | holm_flagged_controls | gate_passed | inference_scope |
| --- | --- | --- | --- | --- | --- |
| DeepSurv_N50 | recursive_studentized | 30 | 0 | True | DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. |

## Alpha=0.05

| inference_scope | null_scenario | component | empirical_type1 | exact_95ci_lower | exact_95ci_upper | holm_p_within_30 | anti_conservative_gate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | single_active | F | 0.0490 | 0.0400 | 0.0594 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | single_active | S | 0.0250 | 0.0186 | 0.0328 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | single_active | U | 0.0505 | 0.0413 | 0.0610 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | single_active | G_claim | 0.0435 | 0.0350 | 0.0534 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | single_active | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | tie | F | 0.0055 | 0.0027 | 0.0098 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | tie | S | 0.0010 | 0.0001 | 0.0036 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | tie | U | 0.0030 | 0.0011 | 0.0065 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | tie | G_claim | 0.0015 | 0.0003 | 0.0044 | 1.0000 | False |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | tie | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |

## Holm-düzeltilmiş ihlaller

İhlal yok.

Gate sonucu: **GEÇTİ**. Gerçek wall-clock 23.161 saniye. Formül SHA256 değişmedi; model/explainer çağrısı yapılmadı ve kimlik-sızıntısı kontrolü geçti.
