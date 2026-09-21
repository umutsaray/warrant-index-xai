# Faz 7.5 — Özet

## Kapsam ve seçim

Ek audit için fold 1 seçildi: RSF C-index 0.789609 ve DeepSurv C-index 0.719193 ile her iki modelin beş-fold maksimumudur. Ortak medyan fold olmadığı ve ikinci tam tekrar isteğe bağlı olduğu için yalnız fold 1 çalıştırıldı. Fold 0 resmi audit referansı olarak kaldı.

## Adım sonuçları

| Adım | Sonuç | Ölçülen süre |
|---|---|---:|
| Model refit/doğrulama | CoxPH 0.705882; RSF 0.789609; DeepSurv checkpoint 0.719193 | 29.9 sn |
| Panel ve METABRIC | 80 feature; 153 TCGA audit hastası; 1466 METABRIC hastası; 12+12 query | 23.3 sn |
| DeepSurv N=50 | 7/80 point; R:80 | 407.3 sn |
| RSF N=50 | 10/80 point; R:76, G_claim:4 | 1884.6 sn |
| CoxPH N=50 | 15 replika tamamlandı; replika 15 iki kez aynı frozen retry sonrası başarısız | 922.2 sn başarılı replika + 406.7 sn son tanı |
| Null gate | RSF 0/30, DeepSurv 0/30 Holm ihlali; ikisi de geçti | 48.7 sn |
| B=4999 sertifikasyon | RSF 7/7 BH/BY; DeepSurv 0/0; CoxPH NA | 20.3 sn |

İlk artefakt ile son sertifikasyon arasındaki gerçek elapsed süre yaklaşık **35.8 dakika**dir; üç N=50 işi paralel çalıştığı için satır süreleri toplanmamalıdır.

## CoxPH beklenmeyen sapması

Deterministik patient-cluster bootstrap replika 15 (sıfır tabanlı indeks) ilk `max_steps=50` denemesinden ve dondurulmuş `max_steps=200, step_size=0.25` retry'sinden sonra aynı “Newton-Raphson failed to converge sufficiently” uyarısını verdi. Ön tarama 12 koşullu-sıfır-varyans sütununu dışladı; lifelines ayrıca isimli bir separation sütunu işaretlemedi. Replika atlanmadı, değiştirilmedi ve protokol genişletilmedi. Bu nedenle CoxPH N=50, taze gate ve BH/BY sonucu üretilmedi.

## Bütünlük

`components.py` SHA256 `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa`, `aggregation.py` SHA256 `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68` olarak kaldı. Audit çıktılarında hasta/örnek kimliği yoktur ve train–test hasta kesişimi sıfırdır.
