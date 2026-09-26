# Faz 10 — CoxPH-v2 retry önyargısı derinleştirmesi

50 dondurulmuş replika yeniden fit edilmeden aynı 28 retry / 22 direct ayrımıyla incelendi. Objective, her replikaya ait saklı katsayı vektörü ve deterministik olarak yeniden kurulan hasta-cluster bootstrap örneği üzerinde hesaplanan penalize edilmemiş Cox log-partial-likelihood değeridir.

| metric | retry_n | direct_n | retry_mean | direct_mean | difference_retry_minus_direct | combined_SE | CI95_lower | CI95_upper | Welch_t | Welch_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| coefficient_L2 | 28 | 22 | 1.4321586670754785 | 1.554891636479284 | -0.12273296940380551 | 0.034663928015841255 | -0.19067426831485437 | -0.054791670492756656 | -3.5406538274519006 | 0.0012913850242353524 |
| log_partial_likelihood | 28 | 22 | -290.42422800387254 | -272.09910239528017 | -18.325125608592373 | 11.620452402281991 | -41.10121231706508 | 4.4509610998803275 | -1.5769717885504815 | 0.12242048697380295 |

## Yalnız 22 doğrudan-yakınsayan replika

`recursive_studentized`, B=4999, τ=0.8 ve aynı BH/BY prosedürü yalnız 22 direct replika üzerinde yeniden uygulandı. Tam-50 sonuçları değiştirilmedi.

| Analiz | Point | BH | BY |
| --- | ---: | ---: | ---: |
| Tam 50 | 9 | 8 | 8 |
| Direct 22 | 9 | 6 | 6 |

### Sertifikalı kümelerin kesin dökümü

- **Tam 50 replika:** 27 özellikten **8** tanesi BH ve BY ile sertifikalıdır: `KRTAP19-7, RBMY1F, CA1, OR13C3, KRTAP19-2, OR5M10, Age at index, AJCC Stage III`.
- **Yalnız direct 22 replika:** 27 özellikten **6** tanesi BH ve BY ile sertifikalıdır: `KRTAP19-7, RBMY1F, KRTAP19-2, OR5M10, Age at index, AJCC Stage III`.
- **Sertifikadan düşenler:** `CA1, OR13C3`.
- **Direct-22 analizinde yeni eklenenler:** `yok`.

Dolayısıyla CA1 ve OR13C3'ün **ikisi de** tam-50 kümesinde sertifikalıyken direct-22 analizinde sertifikadan düşmüştür; biri düşüp diğeri eklenmemiştir. BH ve BY kümeleri her iki senaryoda da kendi içinde aynıdır.

Point-adequate sonuç tam-50 için **9/27**, direct-22 için **9/27**'dir. BH/BY kararı değişen özellik sayısı **2/27**.

Süre **1.959 sn**. Yeni bootstrap/model fit üretilmedi; resmi sertifikasyon korunmuştur.
