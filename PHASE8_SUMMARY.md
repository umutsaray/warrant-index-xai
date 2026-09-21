# Faz 8 — Özet

**Sorgu-varyansı hipotezi dört null-attribution fold-feature çiftinin tamamında reddedildi:** fiilen kullanılan 24 sorgunun varyansı sıfırdan büyüktür; split kullanımı ile tam-sıfır attribution arasındaki açıklama açık kalır.

| Adım | Sonuç | Süre |
|---|---|---:|
| Sorgu varyansı | Fold-0 Stage IV/OR6C1 ve fold-1 Stage IV/OR13C3: 0/4 varyans=0, 4/4 varyans>0 | 0.010 sn |
| CoxPH retry kontrolü | 28 retry vs 22 direct replika; 5/27 feature nominal `|fark|>1.96×SE` flag'i | 1.243 sn |
| GitHub hazırlık ve yayın | 210 dosya, ~3.1 MB; yasak tracked dosya yok; public `main` push doğrulandı | yaklaşık 4.5 dk |

CoxPH retry kontrolündeki beş flag yeni veya çoklu-test düzeltilmiş bir sertifikasyon sonucu değildir; yalnız tanımlayıcı sağlamlık sinyalidir. WI, component skor tanımları, eşikler ve mevcut BH/BY kararları değiştirilmedi.

Public depo: https://github.com/umutsaray/warrant-index-survival-xai

Formül SHA256 değerleri değişmedi: `components.py` `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa`; `aggregation.py` `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68`.
