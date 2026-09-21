# Faz 7.2 — Taze null-kalibrasyon kapıları

Her model kendi temiz N=50 replika varyans/bağımlılık yapısıyla yeniden kalibre edildi. Batarya model başına 2000×999, `single_active`+`tie`, üç alpha ve **30 kontrol içinde Holm** düzeltmesidir. Model veya explainer çağrılmadı.

| model | controls | holm_flagged_controls | gate_passed | inference_scope |
| --- | --- | --- | --- | --- |
| CoxPH_N50 | 30 | 0 | True | 153 TCGA + 1466 METABRIC fixed queries; exact closed-form Cox attribution |
| RSF_N50 | 30 | 0 | True | fixed 12+12 queries and stored 8/128 attribution averages; permutation-path uncertainty excluded |
| DeepSurv_N50 | 30 | 0 | True | fixed 12+12 queries and deterministic nested DeepExplainer k16/k128 backgrounds |

## Alpha=0.05

| model | null_scenario | component | empirical_type1 | exact_95ci_lower | exact_95ci_upper | holm_p_within_model_30 | anti_conservative_gate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH_N50 | single_active | F | 0.0425 | 0.0341 | 0.0523 | 1.0000 | False |
| CoxPH_N50 | single_active | S | 0.0275 | 0.0208 | 0.0356 | 1.0000 | False |
| CoxPH_N50 | single_active | U | 0.0440 | 0.0354 | 0.0539 | 1.0000 | False |
| CoxPH_N50 | single_active | G_claim | 0.0420 | 0.0336 | 0.0517 | 1.0000 | False |
| CoxPH_N50 | single_active | R | 0.0005 | 0.0000 | 0.0028 | 1.0000 | False |
| CoxPH_N50 | tie | F | 0.0030 | 0.0011 | 0.0065 | 1.0000 | False |
| CoxPH_N50 | tie | S | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| CoxPH_N50 | tie | U | 0.0020 | 0.0005 | 0.0051 | 1.0000 | False |
| CoxPH_N50 | tie | G_claim | 0.0045 | 0.0021 | 0.0085 | 1.0000 | False |
| CoxPH_N50 | tie | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| RSF_N50 | single_active | F | 0.0480 | 0.0391 | 0.0583 | 1.0000 | False |
| RSF_N50 | single_active | S | 0.0295 | 0.0225 | 0.0379 | 1.0000 | False |
| RSF_N50 | single_active | U | 0.0475 | 0.0386 | 0.0578 | 1.0000 | False |
| RSF_N50 | single_active | G_claim | 0.0530 | 0.0436 | 0.0637 | 1.0000 | False |
| RSF_N50 | single_active | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| RSF_N50 | tie | F | 0.0015 | 0.0003 | 0.0044 | 1.0000 | False |
| RSF_N50 | tie | S | 0.0010 | 0.0001 | 0.0036 | 1.0000 | False |
| RSF_N50 | tie | U | 0.0025 | 0.0008 | 0.0058 | 1.0000 | False |
| RSF_N50 | tie | G_claim | 0.0015 | 0.0003 | 0.0044 | 1.0000 | False |
| RSF_N50 | tie | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| DeepSurv_N50 | single_active | F | 0.0530 | 0.0436 | 0.0637 | 1.0000 | False |
| DeepSurv_N50 | single_active | S | 0.0250 | 0.0186 | 0.0328 | 1.0000 | False |
| DeepSurv_N50 | single_active | U | 0.0365 | 0.0287 | 0.0457 | 1.0000 | False |
| DeepSurv_N50 | single_active | G_claim | 0.0525 | 0.0431 | 0.0632 | 1.0000 | False |
| DeepSurv_N50 | single_active | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| DeepSurv_N50 | tie | F | 0.0030 | 0.0011 | 0.0065 | 1.0000 | False |
| DeepSurv_N50 | tie | S | 0.0010 | 0.0001 | 0.0036 | 1.0000 | False |
| DeepSurv_N50 | tie | U | 0.0035 | 0.0014 | 0.0072 | 1.0000 | False |
| DeepSurv_N50 | tie | G_claim | 0.0035 | 0.0014 | 0.0072 | 1.0000 | False |
| DeepSurv_N50 | tie | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |

## Holm ihlalleri

İhlal yok.

Üç-model ortak gate sonucu: **GEÇTİ**. Gate kalırsa BH/BY üretilmez. Formül SHA256 değişmedi; kimlik-sızıntısı kontrolü geçti. Wall-clock **74.599 saniye**.
