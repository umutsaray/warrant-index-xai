# Faz 7.5 — Fold seçimi

Faz 7.1'de performans görülmeden seçilmiş resmi audit referansı fold 0 olarak kalır. Ek sağlamlık auditi için fold 1 seçildi; bu fold hem RSF hem DeepSurv için beş fold içindeki en yüksek dış-test C-index'ine sahiptir.

| Fold | CoxPH | RSF | DeepSurv | RSF sırası | DeepSurv sırası |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.667704 | 0.538521 | 0.496498 | 5 | 5 |
| 1 | 0.705882 | 0.789609 | 0.719193 | 1 | 1 |
| 2 | 0.681159 | 0.636232 | 0.650000 | 4 | 3 |
| 3 | 0.655113 | 0.700471 | 0.531878 | 2 | 4 |
| 4 | 0.650173 | 0.650748 | 0.662831 | 3 | 2 |

İkinci ek fold bu çalıştırmada seçilmedi. RSF'nin medyanı fold 4, DeepSurv'ün medyanı fold 2 olduğundan ortak bir “medyan performanslı fold” yoktur; bunlardan birini seçmek modele göre post-hoc tercih gerektirirdi. Kullanıcının “zaman/bütçe izin verirse” koşulu altında, tam Faz 7.2 tekrarının ölçülen maliyeti ve bu belirsizlik nedeniyle önceden tanımlı karşıt uç olan fold 1 ile tek ek-fold analizi yapılır.

