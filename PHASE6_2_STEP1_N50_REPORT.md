# Faz 6.2 ADIM 1 — DeepSurv N=50 tam koşu

> **Koşulluluk uyarısı:** DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds.

Faz 6.1'in DeepExplainer candidate k=16 / reference k=128 tasarımı aynı 12+12 query ve 77 feature ile N=50'ye çıkarıldı. İlk iki deterministik attribution korundu; aynı patient-cluster bootstrap seed akışındaki 48 model ayrıca fit edildi. Ana model değiştirilmedi.

## Block agregatları

| inference_scope | block | features | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | clin | 5 | 5 | 0.9997 | 0.9911 | 0.9936 | 0.9887 | 0.7476 | 0.7476 |
| DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds. | rna | 72 | 72 | 0.9981 | 0.9724 | 0.9873 | 0.9582 | 0.7398 | 0.7398 |

N=2'de 11/77, N=50'de **8/77** feature WI≥0.8 geçti. N=2 kazananlarının 1 tanesi korundu, 10 tanesi düştü; 7 yeni feature geçti. N=50 bottleneck dağılımı `{'R': 77}`.

Yeni 48 replikanın birikimli eğitim+açıklama süresi **224.314 sn**; ilk assertion ile duran invocation ve checkpoint-resume birlikte gerçek wall-clock **216.618 sn** sürdü (duran ilk deneme 177.045 sn, resume 39.573 sn). Yeni replikalarda maksimum göreli additivity hatası candidate/reference için 1.77e-06/0.00807. SHAP'ın sabit mutlak 0.01 iç assertion'ı 2 çağrıda tetiklendi; bu çağrılar açıkça kaydedildi ve yalnız önceden tanımlı göreli hata ≤%5 olduğunda kabul edildi. Formül SHA256 değişmedi ve kimlik-sızıntısı kontrolü geçti.

