# Faz 10 özeti

**ADIM 1 sorgu yeniden örneklemesi RSF için özellik düzeyinde duyarlı, DeepSurv için kararlı çıktı.** RSF'nin point-adequate sayısı orijinal 12+12 sorguda 17/27 iken 24+24'te seed'e göre 11–15/27, 48+48'te 12–18/27 oldu; 10/27 özelliğin point durumu ve 5/27 özelliğin darboğazı en az bir konfigürasyonda değişti. DeepSurv bütün boyut ve seed'lerde 0/27 point-adequate kaldı ve hiçbir darboğaz değişmedi. Faz 7.6/7.7'de TCGA/source 12-sorgu kapsamında sıfır sınıflanan sekiz RSF özelliğinin yalnız 4/8'i (`KRTAP19-7`, `KRTAP19-2`, `RBMY1F`, `OR5M10`) daha büyük TCGA sorgularında bütün seed'lerde sıfır kaldı. Dolayısıyla sekizli null-attribution örüntüsü sorgu örneklemesine duyarlıdır. Bu nokta analizi yeni BH/BY sertifikasyonu değildir ve resmi kararları değiştirmez.

**ADIM 2'de tam fold-1 harmonizasyonu mevcut dondurulmuş modellerle hesaplanamaz.** Fold-0'ın 27 özelliğinin yalnız 10'u fold-1'in dondurulmuş 616-feature model girdisinde bulunur; 17 özellik model girdisi değildir ve NA bırakılmıştır. Eksik girdilere sıfır attribution atanmamış, model yeniden fit edilmemiştir. Ortak 10 özellik için mevcut fold-1 RSF/DeepSurv değerleri yan yana verilmiştir; bunların BH/BY kararları özgün 80-feature ailesine aittir. Fold-1 CoxPH, deterministik replika 15 dondurulmuş retry sonrasında yakınsamadığı için ortak özelliklerde de NA'dır.

**ADIM 3'te retry replikalarını dışlamak iki sertifikasyon kararını değiştirdi.** Resmi tam-50 CoxPH-v2 sonucunda 9/27 point-adequate ve 8/27 BH=BY sertifikalı özellik varken, yalnız 22 doğrudan yakınsayan replikada point sayı 9/27, BH=BY sayı 6/27 oldu. `CA1` ve `OR13C3` sertifikadan düştü; yeni özellik eklenmedi. Retry grubunun katsayı L2 normu daha küçüktü (fark −0.122733, %95 GA [−0.190674, −0.054792], Welch p=0.001291); penalize edilmemiş log-partial-likelihood farkı istatistiksel olarak belirgin değildi (−18.3251, %95 GA [−41.1012, 4.4510], p=0.122420). Bu ek duyarlılık analizi resmi tam-50 sonucunu değiştirmez.

## Çıktılar

- `PHASE10_QUERY_RESAMPLING_SENSITIVITY.csv/.md`: 2 model × 27 özellik × 7 sorgu konfigürasyonu = 378 satır; F/S/U/G_signed/G_claim/R/WI, darboğaz, point kararı ve attribution büyüklükleri.
- `PHASE10_FOLD1_HARMONIZED_27FEATURE.csv/.md`: üç model × 27 özellik = 81 satır; 10 ortak özellik ve 17 açık NA.
- `PHASE10_COXPH_RETRY_BIAS_DEPTH.csv/.md`: 27 özellik için tam-50 ile direct-22 kararları ve iki grup tanı istatistikleri.

## Süre ve bütünlük

| Adım | Wall-clock |
| --- | ---: |
| ADIM 1 RSF | 8,136.0 sn |
| ADIM 1 DeepSurv | 2,276.8 sn |
| ADIM 1 toplam | 10,412.8 sn (2 sa 53 dk 32.8 sn) |
| ADIM 2 | 0.759 sn |
| ADIM 3 | 1.959 sn |

`components.py` SHA256 değeri `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa`, `aggregation.py` SHA256 değeri `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68` olarak değişmeden doğrulandı. ADIM 1 sorguları mevcut TCGA/METABRIC audit havuzlarından, mevcut sızıntı dışlama kurallarıyla çekildi. Makale dosyaları ve resmi Table 3–8 sonuçları değiştirilmedi.
