# Faz 7.6 — Özet

Bu faz yalnız saklı attribution tensörlerini ve fitted fold-0 RSF ağaçlarını okudu. Yeni model fit, explainer, bootstrap, kalibrasyon veya sertifikasyon çalıştırılmadı; mevcut BH/BY kararları değiştirilmedi.

| Adım | Sonuç | Süre |
|---|---|---:|
| Fold-0 sınıflandırma | CoxPH-v2 8/0; RSF 4/8 informative/degenerate; assertion geçti | 0.083 sn |
| Fold-0 RSF split kullanımı | 8 degenerate feature'ın 6'sı kullanılmamış; OR6C1 10, AJCC Stage IV 25 ağaçta kullanılmış | 1.097 sn |
| Fold-1 RSF sınıflandırma | 5/2 informative/degenerate | 0.006 sn |

Toplam wall-clock **1.195 sn**. Formül SHA256 değerleri değişmedi: `components.py` `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa`, `aggregation.py` `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68`.
