# Faz 7.5 fold 1 — Taze null-kalibrasyon kapıları

Her model kendi temiz N=50 replika varyans/bağımlılık yapısıyla yeniden kalibre edildi. Batarya model başına 2000×999, `single_active`+`tie`, üç alpha ve **30 kontrol içinde Holm** düzeltmesidir. Model veya explainer çağrılmadı.

| model | controls | holm_flagged_controls | gate_passed | inference_scope |
| --- | --- | --- | --- | --- |
| RSF_N50 | 30 | 0 | True | fixed 12+12 queries and stored 8/128 attribution averages; permutation-path uncertainty excluded |
| DeepSurv_N50 | 30 | 0 | True | fixed 12+12 queries and deterministic nested DeepExplainer k16/k128 backgrounds |

## Alpha=0.05

| model | null_scenario | component | empirical_type1 | exact_95ci_lower | exact_95ci_upper | holm_p_within_model_30 | anti_conservative_gate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RSF_N50 | single_active | F | 0.0385 | 0.0305 | 0.0479 | 1.0000 | False |
| RSF_N50 | single_active | S | 0.0220 | 0.0160 | 0.0294 | 1.0000 | False |
| RSF_N50 | single_active | U | 0.0420 | 0.0336 | 0.0517 | 1.0000 | False |
| RSF_N50 | single_active | G_claim | 0.0465 | 0.0377 | 0.0567 | 1.0000 | False |
| RSF_N50 | single_active | R | 0.0005 | 0.0000 | 0.0028 | 1.0000 | False |
| RSF_N50 | tie | F | 0.0005 | 0.0000 | 0.0028 | 1.0000 | False |
| RSF_N50 | tie | S | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| RSF_N50 | tie | U | 0.0035 | 0.0014 | 0.0072 | 1.0000 | False |
| RSF_N50 | tie | G_claim | 0.0030 | 0.0011 | 0.0065 | 1.0000 | False |
| RSF_N50 | tie | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| DeepSurv_N50 | single_active | F | 0.0410 | 0.0327 | 0.0506 | 1.0000 | False |
| DeepSurv_N50 | single_active | S | 0.0255 | 0.0190 | 0.0334 | 1.0000 | False |
| DeepSurv_N50 | single_active | U | 0.0485 | 0.0395 | 0.0588 | 1.0000 | False |
| DeepSurv_N50 | single_active | G_claim | 0.0445 | 0.0359 | 0.0545 | 1.0000 | False |
| DeepSurv_N50 | single_active | R | 0.0005 | 0.0000 | 0.0028 | 1.0000 | False |
| DeepSurv_N50 | tie | F | 0.0035 | 0.0014 | 0.0072 | 1.0000 | False |
| DeepSurv_N50 | tie | S | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| DeepSurv_N50 | tie | U | 0.0020 | 0.0005 | 0.0051 | 1.0000 | False |
| DeepSurv_N50 | tie | G_claim | 0.0025 | 0.0008 | 0.0058 | 1.0000 | False |
| DeepSurv_N50 | tie | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |

## Holm ihlalleri

İhlal yok.

Mevcut iki-model ortak gate sonucu: **GEÇTİ**. CoxPH, deterministik bootstrap replika 15'in dondurulmuş retry protokolünden sonra yakınsamaması nedeniyle bu kapıda yoktur ve CoxPH sertifikasyonu üretilmez. Gate kalırsa RSF/DeepSurv BH/BY de üretilmez. Formül SHA256 değişmedi; kimlik-sızıntısı kontrolü geçti. Wall-clock **48.701 saniye**.
