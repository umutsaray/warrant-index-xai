# Faz 7.4 — τ duyarlılık taraması

Ana analiz **τ=0.80** olarak kalır. Aynı Faz 7.2 T=50 attribution tensörleri ve aynı B=4999 inferential resample akışı kullanılarak yalnız tek-yönlü null sınırı değiştirildi; yeni patient bootstrap, model fit veya explainer çağrısı yapılmadı. Point kararı istenen katı `WI>τ` kuralıdır.

| tau | reference_primary | CoxPH_point | CoxPH_BH | CoxPH_BY | RSF_point | RSF_BH | RSF_BY | DeepSurv_point | DeepSurv_BH | DeepSurv_BY |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.7000 | False | 22 | 17 | 15 | 26 | 25 | 22 | 21 | 6 | 6 |
| 0.7500 | False | 16 | 9 | 9 | 22 | 16 | 12 | 6 | 1 | 1 |
| 0.8000 | True | 9 | 8 | 8 | 17 | 12 | 12 | 0 | 0 | 0 |
| 0.8500 | False | 2 | 1 | 1 | 13 | 12 | 12 | 0 | 0 | 0 |
| 0.9000 | False | 1 | 1 | 1 | 13 | 10 | 10 | 0 | 0 | 0 |

DeepSurv τ=0.70 satırında point/BH/BY sayıları **21/6/6**. Bu, τ=0.80'deki 0/0/0 sonucunun düşük eşikte nasıl değiştiğini doğrudan gösterir.

DeepSurv τ=0.75'te 6/1/1'e, τ=0.70'te 21/6/6'ya çıktığı için ana τ=0.80'deki 0/0/0 sonucu aşağı yönlü eşik değişimine duyarlıdır. CoxPH ve RSF sayıları da τ yükseldikçe genel olarak azalır; ana sonuçlar değiştirilmemiştir.
