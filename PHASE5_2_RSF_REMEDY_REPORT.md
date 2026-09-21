# Faz 5.2B — RSF anti-conservative test düzeltmesi

> **Koşulluluk uyarısı:** RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

## Özet — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

Üç aday, Faz 5.1 ile aynı iki null senaryosu, üç alpha, 2000 dış simülasyon ve B=999 iç resample bataryasında sınandı. Her adayın 30 RSF kontrolü kendi içinde Holm ile düzeltildi. Frozen F/S/U/G/R/WI formülleri değiştirilmedi; yalnız test katmanı karşılaştırıldı. RSF ADIM 2 çalıştırılmadı.

| candidate | gate_failed | gate_passed |
| --- | --- | --- |
| m_out_of_n_root | True | False |
| recursive_basic | True | False |
| recursive_studentized | False | True |

## Alpha=0.05 — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

| candidate | null_scenario | component | empirical_type1 | exact_95ci_lower | exact_95ci_upper | anti_conservative_gate |
| --- | --- | --- | --- | --- | --- | --- |
| recursive_basic | single_active | F | 0.0615 | 0.0514 | 0.0729 | False |
| recursive_basic | single_active | S | 0.0290 | 0.0221 | 0.0373 | False |
| recursive_basic | single_active | U | 0.0660 | 0.0555 | 0.0778 | True |
| recursive_basic | single_active | G_claim | 0.0595 | 0.0495 | 0.0708 | False |
| recursive_basic | single_active | R | 0.0030 | 0.0011 | 0.0065 | False |
| recursive_basic | tie | F | 0.0040 | 0.0017 | 0.0079 | False |
| recursive_basic | tie | S | 0.0005 | 0.0000 | 0.0028 | False |
| recursive_basic | tie | U | 0.0060 | 0.0031 | 0.0105 | False |
| recursive_basic | tie | G_claim | 0.0020 | 0.0005 | 0.0051 | False |
| recursive_basic | tie | R | 0.0000 | 0.0000 | 0.0018 | False |
| recursive_studentized | single_active | F | 0.0365 | 0.0287 | 0.0457 | False |
| recursive_studentized | single_active | S | 0.0290 | 0.0221 | 0.0373 | False |
| recursive_studentized | single_active | U | 0.0350 | 0.0274 | 0.0440 | False |
| recursive_studentized | single_active | G_claim | 0.0290 | 0.0221 | 0.0373 | False |
| recursive_studentized | single_active | R | 0.0030 | 0.0011 | 0.0065 | False |
| recursive_studentized | tie | F | 0.0015 | 0.0003 | 0.0044 | False |
| recursive_studentized | tie | S | 0.0005 | 0.0000 | 0.0028 | False |
| recursive_studentized | tie | U | 0.0010 | 0.0001 | 0.0036 | False |
| recursive_studentized | tie | G_claim | 0.0005 | 0.0000 | 0.0028 | False |
| recursive_studentized | tie | R | 0.0000 | 0.0000 | 0.0018 | False |
| m_out_of_n_root | single_active | F | 0.0630 | 0.0527 | 0.0746 | False |
| m_out_of_n_root | single_active | S | 0.0480 | 0.0391 | 0.0583 | False |
| m_out_of_n_root | single_active | U | 0.0665 | 0.0560 | 0.0783 | True |
| m_out_of_n_root | single_active | G_claim | 0.0805 | 0.0689 | 0.0933 | True |
| m_out_of_n_root | single_active | R | 0.0315 | 0.0243 | 0.0401 | False |
| m_out_of_n_root | tie | F | 0.0420 | 0.0336 | 0.0517 | False |
| m_out_of_n_root | tie | S | 0.0070 | 0.0038 | 0.0117 | False |
| m_out_of_n_root | tie | U | 0.0375 | 0.0296 | 0.0468 | False |
| m_out_of_n_root | tie | G_claim | 0.0480 | 0.0391 | 0.0583 | False |
| m_out_of_n_root | tie | R | 0.0020 | 0.0005 | 0.0051 | False |

## Gate ihlalleri — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

| candidate | null_scenario | component | nominal_alpha | empirical_type1 | exact_95ci_lower | exact_95ci_upper | holm_p_within_candidate_30 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| recursive_basic | single_active | F | 0.0100 | 0.0195 | 0.0139 | 0.0266 | 0.0029 |
| recursive_basic | single_active | U | 0.0100 | 0.0200 | 0.0143 | 0.0271 | 0.0014 |
| recursive_basic | single_active | U | 0.0500 | 0.0660 | 0.0555 | 0.0778 | 0.0267 |
| m_out_of_n_root | single_active | F | 0.0100 | 0.0245 | 0.0182 | 0.0323 | 0.0000 |
| m_out_of_n_root | single_active | S | 0.0100 | 0.0230 | 0.0169 | 0.0306 | 0.0000 |
| m_out_of_n_root | single_active | U | 0.0100 | 0.0235 | 0.0173 | 0.0311 | 0.0000 |
| m_out_of_n_root | single_active | U | 0.0500 | 0.0665 | 0.0560 | 0.0783 | 0.0166 |
| m_out_of_n_root | single_active | G_claim | 0.0100 | 0.0405 | 0.0323 | 0.0501 | 0.0000 |
| m_out_of_n_root | single_active | G_claim | 0.0500 | 0.0805 | 0.0689 | 0.0933 | 0.0000 |
| m_out_of_n_root | single_active | G_claim | 0.1000 | 0.1255 | 0.1113 | 0.1408 | 0.0033 |

## Öneri — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

Gate'i geçen adaylar: recursive_studentized. Önerilen aday **recursive_studentized**; RSF gerçek-veri pilotu ve ADIM 2 için kullanıcı onayı gerektirir.

## Gerçek 5-feature remedy pilotu — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

| feature | selection_reason | F | S | U | G_claim | R | p_F_recursive_studentized | p_S_recursive_studentized | p_U_recursive_studentized | p_G_claim_recursive_studentized | p_R_recursive_studentized | p_joint_recursive_studentized |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000163545.11 | high_margin | 1.0000 | 0.9974 | 0.9949 | 0.9886 | 0.8703 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.0002 |
| rna_ENSG00000178522.14 | near_threshold | 1.0000 | 0.9918 | 0.9886 | 0.9805 | 0.8002 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.6154 | 0.6154 |
| rna_ENSG00000090402.8 | phase4.1_low_WI | 0.9977 | 0.5920 | 0.6834 | 0.0000 | 0.8260 | 0.0002 | 1.0000 | 0.9928 | 1.0000 | 0.2722 | 1.0000 |
| rna_ENSG00000231925.12 | seeded_random | 0.9997 | 0.9627 | 0.9599 | 0.8839 | 0.7576 | 0.0002 | 0.0002 | 0.0002 | 0.0004 | 0.9432 | 0.9432 |
| rna_ENSG00000187242.6 | seeded_random | 0.9998 | 0.8484 | 0.9307 | 0.7758 | 0.6899 | 0.0002 | 0.0752 | 0.0002 | 0.7570 | 1.0000 | 1.0000 |

Gerçek wall-clock: **19.07 saniye**. Formül SHA256 değişmedi: **True**. Kimlik-sızıntısı kontrolü geçti.
