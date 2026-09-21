# Faz 5.1 — null kalibrasyon / Tip-I hata kapısı

> **RSF koşulluluk uyarısı:** RSF nuisance calibration conditions on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

## Özet — RSF sabit 12+12 query ve saklı attribution ortalamalarına koşullu

Gerçek model veya explainer çağrısı yapılmadı. Her model ve null senaryosu için 2000 bağımsız dış simülasyon ve her simülasyonda B=999 null-merkezli iç bootstrap kullanıldı. Beş sentetik popülasyon hedefi analitik olarak tam 0.8'dir. `single_active` bir site/dal 0.8 ve diğeri 0.9 olan düzgün sınırı; `tie` iki site/dalın da 0.8 olduğu nonsmooth `min` sınırını test eder. Bileşen korelasyonu ve göreli gürültü ölçeği saklı CoxPH N=50 / RSF N=20 artefaktlarından alınmıştır.

Kapı kararı: **KALDI**. Anti-conservative etiketi, nominal oranın üzerindeki sapmanın 60 kontrol boyunca Holm-düzeltilmiş tek-yönlü exact binomial p<0.05 olmasıyla verildi.

## Alpha=0.05 özeti — RSF satırları sabit 12+12 query ve saklı attribution ortalamalarına koşullu

| model | null_scenario | component | nominal_alpha | empirical_type1 | exact_95ci_lower | exact_95ci_upper | holm_p_across_checks | anti_conservative_gate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH_N50 | single_active | F | 0.0500 | 0.0550 | 0.0454 | 0.0659 | 1.0000 | False |
| CoxPH_N50 | single_active | S | 0.0500 | 0.0295 | 0.0225 | 0.0379 | 1.0000 | False |
| CoxPH_N50 | single_active | U | 0.0500 | 0.0505 | 0.0413 | 0.0610 | 1.0000 | False |
| CoxPH_N50 | single_active | G_claim | 0.0500 | 0.0585 | 0.0486 | 0.0697 | 1.0000 | False |
| CoxPH_N50 | single_active | R | 0.0500 | 0.0005 | 0.0000 | 0.0028 | 1.0000 | False |
| CoxPH_N50 | tie | F | 0.0500 | 0.0280 | 0.0212 | 0.0362 | 1.0000 | False |
| CoxPH_N50 | tie | S | 0.0500 | 0.0125 | 0.0081 | 0.0184 | 1.0000 | False |
| CoxPH_N50 | tie | U | 0.0500 | 0.0360 | 0.0283 | 0.0451 | 1.0000 | False |
| CoxPH_N50 | tie | G_claim | 0.0500 | 0.0310 | 0.0238 | 0.0396 | 1.0000 | False |
| CoxPH_N50 | tie | R | 0.0500 | 0.0015 | 0.0003 | 0.0044 | 1.0000 | False |
| RSF_N20 | single_active | F | 0.0500 | 0.0570 | 0.0472 | 0.0681 | 1.0000 | False |
| RSF_N20 | single_active | S | 0.0500 | 0.0410 | 0.0327 | 0.0506 | 1.0000 | False |
| RSF_N20 | single_active | U | 0.0500 | 0.0615 | 0.0514 | 0.0729 | 0.6632 | False |
| RSF_N20 | single_active | G_claim | 0.0500 | 0.0730 | 0.0620 | 0.0853 | 0.0003 | True |
| RSF_N20 | single_active | R | 0.0500 | 0.0070 | 0.0038 | 0.0117 | 1.0000 | False |
| RSF_N20 | tie | F | 0.0500 | 0.0315 | 0.0243 | 0.0401 | 1.0000 | False |
| RSF_N20 | tie | S | 0.0500 | 0.0105 | 0.0065 | 0.0160 | 1.0000 | False |
| RSF_N20 | tie | U | 0.0500 | 0.0285 | 0.0217 | 0.0368 | 1.0000 | False |
| RSF_N20 | tie | G_claim | 0.0500 | 0.0265 | 0.0199 | 0.0345 | 1.0000 | False |
| RSF_N20 | tie | R | 0.0500 | 0.0020 | 0.0005 | 0.0051 | 1.0000 | False |

## Tip-I oranları — RSF satırları sabit 12+12 query ve saklı attribution ortalamalarına koşullu

| model | null_scenario | component | nominal_alpha | empirical_type1 | exact_95ci_lower | exact_95ci_upper | holm_p_across_checks | anti_conservative_gate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH_N50 | single_active | F | 0.0100 | 0.0115 | 0.0073 | 0.0172 | 1.0000 | False |
| CoxPH_N50 | single_active | F | 0.0500 | 0.0550 | 0.0454 | 0.0659 | 1.0000 | False |
| CoxPH_N50 | single_active | F | 0.1000 | 0.1080 | 0.0947 | 0.1224 | 1.0000 | False |
| CoxPH_N50 | single_active | S | 0.0100 | 0.0045 | 0.0021 | 0.0085 | 1.0000 | False |
| CoxPH_N50 | single_active | S | 0.0500 | 0.0295 | 0.0225 | 0.0379 | 1.0000 | False |
| CoxPH_N50 | single_active | S | 0.1000 | 0.0520 | 0.0427 | 0.0627 | 1.0000 | False |
| CoxPH_N50 | single_active | U | 0.0100 | 0.0105 | 0.0065 | 0.0160 | 1.0000 | False |
| CoxPH_N50 | single_active | U | 0.0500 | 0.0505 | 0.0413 | 0.0610 | 1.0000 | False |
| CoxPH_N50 | single_active | U | 0.1000 | 0.1005 | 0.0877 | 0.1145 | 1.0000 | False |
| CoxPH_N50 | single_active | G_claim | 0.0100 | 0.0130 | 0.0085 | 0.0190 | 1.0000 | False |
| CoxPH_N50 | single_active | G_claim | 0.0500 | 0.0585 | 0.0486 | 0.0697 | 1.0000 | False |
| CoxPH_N50 | single_active | G_claim | 0.1000 | 0.0965 | 0.0839 | 0.1103 | 1.0000 | False |
| CoxPH_N50 | single_active | R | 0.0100 | 0.0005 | 0.0000 | 0.0028 | 1.0000 | False |
| CoxPH_N50 | single_active | R | 0.0500 | 0.0005 | 0.0000 | 0.0028 | 1.0000 | False |
| CoxPH_N50 | single_active | R | 0.1000 | 0.0025 | 0.0008 | 0.0058 | 1.0000 | False |
| CoxPH_N50 | tie | F | 0.0100 | 0.0065 | 0.0035 | 0.0111 | 1.0000 | False |
| CoxPH_N50 | tie | F | 0.0500 | 0.0280 | 0.0212 | 0.0362 | 1.0000 | False |
| CoxPH_N50 | tie | F | 0.1000 | 0.0545 | 0.0450 | 0.0654 | 1.0000 | False |
| CoxPH_N50 | tie | S | 0.0100 | 0.0020 | 0.0005 | 0.0051 | 1.0000 | False |
| CoxPH_N50 | tie | S | 0.0500 | 0.0125 | 0.0081 | 0.0184 | 1.0000 | False |
| CoxPH_N50 | tie | S | 0.1000 | 0.0235 | 0.0173 | 0.0311 | 1.0000 | False |
| CoxPH_N50 | tie | U | 0.0100 | 0.0095 | 0.0057 | 0.0148 | 1.0000 | False |
| CoxPH_N50 | tie | U | 0.0500 | 0.0360 | 0.0283 | 0.0451 | 1.0000 | False |
| CoxPH_N50 | tie | U | 0.1000 | 0.0670 | 0.0564 | 0.0789 | 1.0000 | False |
| CoxPH_N50 | tie | G_claim | 0.0100 | 0.0075 | 0.0042 | 0.0123 | 1.0000 | False |
| CoxPH_N50 | tie | G_claim | 0.0500 | 0.0310 | 0.0238 | 0.0396 | 1.0000 | False |
| CoxPH_N50 | tie | G_claim | 0.1000 | 0.0580 | 0.0482 | 0.0692 | 1.0000 | False |
| CoxPH_N50 | tie | R | 0.0100 | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| CoxPH_N50 | tie | R | 0.0500 | 0.0015 | 0.0003 | 0.0044 | 1.0000 | False |
| CoxPH_N50 | tie | R | 0.1000 | 0.0065 | 0.0035 | 0.0111 | 1.0000 | False |
| RSF_N20 | single_active | F | 0.0100 | 0.0205 | 0.0148 | 0.0277 | 0.0013 | True |
| RSF_N20 | single_active | F | 0.0500 | 0.0570 | 0.0472 | 0.0681 | 1.0000 | False |
| RSF_N20 | single_active | F | 0.1000 | 0.0970 | 0.0844 | 0.1108 | 1.0000 | False |
| RSF_N20 | single_active | S | 0.0100 | 0.0175 | 0.0122 | 0.0243 | 0.0774 | False |
| RSF_N20 | single_active | S | 0.0500 | 0.0410 | 0.0327 | 0.0506 | 1.0000 | False |
| RSF_N20 | single_active | S | 0.1000 | 0.0670 | 0.0564 | 0.0789 | 1.0000 | False |
| RSF_N20 | single_active | U | 0.0100 | 0.0215 | 0.0156 | 0.0289 | 0.0003 | True |
| RSF_N20 | single_active | U | 0.0500 | 0.0615 | 0.0514 | 0.0729 | 0.6632 | False |
| RSF_N20 | single_active | U | 0.1000 | 0.1025 | 0.0895 | 0.1166 | 1.0000 | False |
| RSF_N20 | single_active | G_claim | 0.0100 | 0.0345 | 0.0269 | 0.0435 | 0.0000 | True |
| RSF_N20 | single_active | G_claim | 0.0500 | 0.0730 | 0.0620 | 0.0853 | 0.0003 | True |
| RSF_N20 | single_active | G_claim | 0.1000 | 0.1210 | 0.1070 | 0.1361 | 0.0719 | False |
| RSF_N20 | single_active | R | 0.0100 | 0.0060 | 0.0031 | 0.0105 | 1.0000 | False |
| RSF_N20 | single_active | R | 0.0500 | 0.0070 | 0.0038 | 0.0117 | 1.0000 | False |
| RSF_N20 | single_active | R | 0.1000 | 0.0075 | 0.0042 | 0.0123 | 1.0000 | False |
| RSF_N20 | tie | F | 0.0100 | 0.0065 | 0.0035 | 0.0111 | 1.0000 | False |
| RSF_N20 | tie | F | 0.0500 | 0.0315 | 0.0243 | 0.0401 | 1.0000 | False |
| RSF_N20 | tie | F | 0.1000 | 0.0630 | 0.0527 | 0.0746 | 1.0000 | False |
| RSF_N20 | tie | S | 0.0100 | 0.0030 | 0.0011 | 0.0065 | 1.0000 | False |
| RSF_N20 | tie | S | 0.0500 | 0.0105 | 0.0065 | 0.0160 | 1.0000 | False |
| RSF_N20 | tie | S | 0.1000 | 0.0175 | 0.0122 | 0.0243 | 1.0000 | False |
| RSF_N20 | tie | U | 0.0100 | 0.0075 | 0.0042 | 0.0123 | 1.0000 | False |
| RSF_N20 | tie | U | 0.0500 | 0.0285 | 0.0217 | 0.0368 | 1.0000 | False |
| RSF_N20 | tie | U | 0.1000 | 0.0625 | 0.0523 | 0.0740 | 1.0000 | False |
| RSF_N20 | tie | G_claim | 0.0100 | 0.0090 | 0.0053 | 0.0142 | 1.0000 | False |
| RSF_N20 | tie | G_claim | 0.0500 | 0.0265 | 0.0199 | 0.0345 | 1.0000 | False |
| RSF_N20 | tie | G_claim | 0.1000 | 0.0500 | 0.0409 | 0.0605 | 1.0000 | False |
| RSF_N20 | tie | R | 0.0100 | 0.0000 | 0.0000 | 0.0018 | 1.0000 | False |
| RSF_N20 | tie | R | 0.0500 | 0.0020 | 0.0005 | 0.0051 | 1.0000 | False |
| RSF_N20 | tie | R | 0.1000 | 0.0080 | 0.0046 | 0.0130 | 1.0000 | False |

## Yorum — RSF değerlendirmesi sabit 12+12 query ve saklı attribution ortalamalarına koşullu

Aşağıdaki hücreler anti-conservative kapıyı geçti; mevcut yöntemle ADIM 2'ye geçilmemelidir:

| model | null_scenario | component | nominal_alpha | empirical_type1 | exact_95ci_lower | exact_95ci_upper | holm_p_across_checks | anti_conservative_gate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RSF_N20 | single_active | F | 0.0100 | 0.0205 | 0.0148 | 0.0277 | 0.0013 | True |
| RSF_N20 | single_active | U | 0.0100 | 0.0215 | 0.0156 | 0.0289 | 0.0003 | True |
| RSF_N20 | single_active | G_claim | 0.0100 | 0.0345 | 0.0269 | 0.0435 | 0.0000 | True |
| RSF_N20 | single_active | G_claim | 0.0500 | 0.0730 | 0.0620 | 0.0853 | 0.0003 | True |

RSF N=20 `single_active` sınırında F ve U için alpha=0.01 ret oranları sırasıyla 0.0205 ve 0.0215'tir. G_claim daha belirgin biçimde alpha=0.01'de 0.0345, alpha=0.05'te 0.0730 ret vermiştir. Tie senaryoları anti-conservative değildir ve çoğunlukla muhafazakârdır; dolayısıyla yalnız tie null'u kullanmak sorunu gizlerdi. R her iki modelde belirgin muhafazakârdır; bu Tip-I ihlali değil, olası güç kaybıdır.

Olası düzeltmeler: F/U için studentized bootstrap-t veya bounded-score tabanlı muhafazakâr bir tek-yönlü bound; nonsmooth G_claim için null-constrained bootstrap, m-out-of-n/subsampling veya iki alt bileşen için ayrı geçerli testlerin intersection-union birleşimi. Hangi yaklaşımın seçileceği ayrı tasarım ve yeniden kalibrasyon gerektirir; kullanıcı onayı olmadan hiçbiri uygulanmadı.

Formül SHA256 değişmedi: **True**. Kimlik-sızıntısı kontrolü: **geçti**. Gerçek wall-clock: **15.36 saniye**. Tam BH/BY uygulanmadı.
