# Faz 5.3 — RSF sertifikasyonu ve CoxPH yöntem tutarlılığı

> **RSF koşulluluk uyarısı:** RSF sonuçları sabit 12 TCGA + 12 METABRIC sorgusuna ve saklanmış 8/128 attribution ortalamalarına koşulludur; tekil permütasyon-yolu belirsizliğini içermez.

## RSF ADIM 2 — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

Faz 4.2'nin önceden tanımlı 77 feature ailesinin tamamı, Faz 5.2'de null gate'i geçen `recursive_studentized` yöntemiyle test edildi. F/U/G_magnitude/G_signed yapraklarında studentized bootstrap, S/R yapraklarında basic bootstrap kullanıldı; her `min` düğümünde p-değerleri `max` ile recursive IUT olarak birleştirildi.

| Ölçüt — RSF sabit 12+12 query ve saklı attribution ortalamalarına koşullu | Sayı |
| --- | ---: |
| Değerlendirilen tam aile | 77 |
| Nokta tahmini WI≥0.8 | 23 |
| BH sertifikalı, α=0.05 | 15 |
| BY sertifikalı, α=0.05 | 12 |
| Nokta + BH | 15 |
| Yalnız nokta | 8 |
| Yalnız BH | 0 |

BH'ye giremeyen sekiz nokta-yeterli feature'ın tamamında nokta-tahmin darboğazı R'dir. Tam sonuçlar `results/phase53_rsf_certification/feature_certification.csv` dosyasındadır.

## CoxPH yöntem tutarlılığı

Aynı 77 feature, Faz 5.2A basic-root yöntemine ek olarak `recursive_studentized` ile çalıştırıldı.

| Yöntem | BH | BY |
| --- | ---: | ---: |
| Basic-root | 8 | 8 |
| Recursive-studentized | 8 | 8 |

İki yöntemin BH sertifika kümeleri birebir aynıdır. F/S/U/R p-değerleri değişmedi. Recursive yöntem yalnız üç feature'ın G_claim p-değerini büyüttü: ENSG00000121351.8 için 0.0436→0.1112, ENSG00000090402.8 için 0.0002→0.0384 ve ENSG00000178522.14 için 0.0002→0.0022. Bu farkların hiçbiri sertifika kümesini değiştirmedi.

## Yöntem önerisi

Makale için ortak birincil yöntem olarak `recursive_studentized` önerilir: RSF N=20'de 30 kontrollük null gate'i geçen yöntemdir ve CoxPH N=50'de kalibre edilmiş basic-root yöntemle aynı sertifika kümesini verir. CoxPH basic-root sonucu robustness analizi olarak korunmalıdır. Bu sonuç, RSF'nin küçük ve sabit query kümesi ile saklanmamış permütasyon yolları sınırlamasını kaldırmaz.

Frozen formüller değişmedi. Model veya explainer çağrısı yapılmadı. RSF tam-aile hesabı 3.937 saniye, CoxPH recursive tutarlılık hesabı 11.673 saniye sürdü; kimlik-sızıntısı kontrolleri geçti.
