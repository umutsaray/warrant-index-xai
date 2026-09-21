# Faz 5.6 ADIM 2 — RSF N=50 null-kalibrasyon kapısı

> **Koşulluluk uyarısı:** RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

`recursive_studentized`, N=50 attribution varyans/bağımlılık yapısından türetilen nuisance ile aynı 2000×999, `single_active` + `tie`, üç alpha ve 30 kontrol içinde Holm düzeltmeli bataryadan geçirildi.

| model | method | controls | holm_flagged_controls | gate_passed | inference_scope |
| --- | --- | --- | --- | --- | --- |
| RSF_N50 | recursive_studentized | 30 | 0 | True | RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. |

## Alpha=0.05

| inference_scope | null_scenario | component | empirical_type1 | exact_95ci_lower | exact_95ci_upper | holm_p_within_30 | anti_conservative_gate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | single_active | F | 0.0465 | 0.0377 | 0.0567 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | single_active | S | 0.0220 | 0.0160 | 0.0294 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | single_active | U | 0.0495 | 0.0404 | 0.0599 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | single_active | G_claim | 0.0435 | 0.0350 | 0.0534 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | single_active | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | tie | F | 0.0035 | 0.0014 | 0.0072 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | tie | S | 0.0015 | 0.0003 | 0.0044 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | tie | U | 0.0035 | 0.0014 | 0.0072 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | tie | G_claim | 0.0020 | 0.0005 | 0.0051 | 1.0000 | False |
| RSF N=50 inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | tie | R | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |

## Holm-düzeltilmiş ihlaller

İhlal yok.

Gate sonucu: **GEÇTİ**. Gerçek wall-clock 22.842 saniye. Formül SHA256 değişmedi; model/explainer çağrısı yapılmadı ve kimlik-sızıntısı kontrolü geçti.
