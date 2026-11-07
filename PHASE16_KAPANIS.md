# Faz 16 — Faz 15 sonrası üç kapanış tanısı

Bu çalışma yalnız tanısal/keşifseldir. Yeni model eğitilmedi, panel yeniden seçilmedi ve mevcut resmî `results/**` artefaktları değiştirilmedi. Yeni ham çıktılar yalnız `results/phase16/` altındadır.

## ADIM 1 — Tam RNA havuzunda alternatif sıralama

Fold-0 train matrisinde seçim öncesi RNA havuzu 60.660 sütundur; 57.782'si fit edilebilir resmî Wald değeri taşır. İlk deneme 256'lık blok kullanınca beş p-değeri `rtol=1e-6` sınırını çok az aştığı için kabul edilmedi. Resmî kodun 512'lik blok boyutuyla yeniden koşuldu; kayıtlı 60.660 Wald p-değeri `rtol=1e-6` içinde ve ilk-200 kümesi/sırası birebir eşleşti.

| Kriter | Wald ilk-200 ile kesişim | Jaccard | Aktif 22 RNA ilk-200 | İlk-200'de >%90 sıfır | Havuzda p<.05 |
|---|---:|---:|---:|---:|---:|
| Wald (resmî) | 200 | 1.0000 | 22 | 166 | 7.784 |
| LRT | 3 | 0.0076 | 1 | 14 | 7.284 |
| Skor | 45 | 0.1268 | 6 | 116 | 8.123 |
| ±5 winsorized Wald | 15 | 0.0390 | 3 | 61 | 7.663 |

| Gen | Wald p/sıra | LRT p/sıra/BH-q | Skor p/sıra | ±5 Wald p/sıra |
|---|---:|---:|---:|---:|
| RBMY1F | 1.36e-55 / 61 | .02399 / 4.559 / .3192 | 7.64e-7 / 106 | .07824 / 10.034 |
| OR6C1 | 1.84e-78 / 35 | .03470 / 5.795 / .3633 | 7.95e-4 / 1.083 | .15068 / 14.855 |
| OR4C46 | 4.81e-63 / 51 | .03113 / 5.399 / .3497 | 2.04e-4 / 653 | .00776 / 2.655 |
| OR13C3 | 1.15e-43 / 85 | .02850 / 5.079 / .3404 | 3.26e-5 / 341 | .02968 / 5.631 |
| OR5M10 | 4.89e-11 / 197 | .08564 / 10.249 / .5066 | .00518 / 2.589 | .53443 / 34.837 |
| KRTAP19-2 | 2.82e-17 / 174 | .04740 / 7.047 / .4080 | 9.10e-5 / 470 | .35530 / 26.059 |
| KRTAP19-7 | 1.05e-118 / 15 | .07294 / 9.241 / .4787 | .00248 / 1.826 | .40061 / 28.298 |
| KRTAP11-1 | 1.07e-37 / 102 | .04599 / 6.907 / .4038 | .00346 / 2.148 | .51591 / 33.946 |

Sekiz odak genin 8/8'i LRT ve winsorized Wald ilk-200'ünün, 7/8'i skor ilk-200'ünün dışındadır; hiçbir odak LRT BH q<.05 değildir. Bütün 60.660 özellik ve aktif 22 RNA'nın sıraları ham CSV'lerdedir.

**Sonuç:** Önceden yazılı hipotez desteklenmiştir: >%90 sıfır genler ilk-200'de Wald 166'dan LRT/skor/winsorized kriterlerinde 14/116/61'e düşmüş; sekiz odak genin çoğu değil, LRT ve winsorized kriterlerinde tamamı ilk 200 dışına çıkmıştır. Bu tanı paneli değiştirmez.

## ADIM 2 — METABRIC platform kayması ve konum/ölçek ayrıştırması

Faz 15'in doğru `metabric_phase72` cache eşlemesi kullanıldı; yeni indirme yapılmadı. Cox'un 1.466 tam sorgusu ile özellik başına tüm mevcut METABRIC değerleri ayrı scope satırları olarak raporlandı. Resmî ve yalnız-TCGA biçimleri Faz 15 ile `atol=1e-12` eşleşti.

| Biçim | Point-adequate | Darboğaz | Medyan G_claim | Medyan WI |
|---|---:|---|---:|---:|
| Resmî konum+ölçek | 9/27 | R:24, G:3 | .9583 | .7586 |
| Yalnız TCGA ölçeği | 3/27 | G:23, R:4 | .0000 | .0000 |
| Yalnız konum düzeltmesi | 7/27 | R:15, G:12 | .8980 | .7431 |
| Yalnız ölçek düzeltmesi | 3/27 | G:24, R:3 | .0000 | .0000 |

Tam sorgu kümesinde medyan platform tanısı:

| Tür | Medyan METABRIC−TCGA kayması (TCGA SD) | Medyan SD oranı | TCGA train aralığı dışı oran |
|---|---:|---:|---:|
| RNA (22) | 147.64 | 4.396 | 1.000 |
| Yaş (1) | 0.222 | 0.999 | 0.0048 |
| AJCC (4) | -0.022 | 0.858 | 0.000 |

Resmî eksi yalnız-TCGA G_claim düşüşü RNA'da min/medyan/maks = .7515/.9574/.9913; yaşta .0433/.0433/.0433; AJCC'de -.0696/-.0113/-.0021'dir. Düşüşün `|kayma|` ile Spearman korelasyonu tüm 27 özellikte ρ=.4280 (p=.0260), `|log SD oranı|` ile ρ=.3553 (p=.0689); yalnız 22 RNA'da sırasıyla ρ=-.0401 ve -.0559'dur. RNA içindeki düşük korelasyon, tüm RNA kaymalarının aşırı ve G'nin çoğunda tabana doymasıyla birliktedir; bu nedensel bir kanıt değildir.

**Sonuç:** Hipotez kısmen/güçlü yapısal biçimde desteklenmiştir: düşüş RNA'da yoğunlaşır ve yaş/evrede küçüktür; konum düzeltmesi tek başına point sayısını 3→7 geri getirirken yalnız ölçek düzeltmesi 3'te kalır. Ancak RNA içindeki düşüş büyüklüğü `|kayma|` veya SD oranıyla monoton ilişkili değildir.

## ADIM 3 — R null ve sınır-yakın alternatifler

Faz 15 ile aynı küçültülmüş tasarım (`outer=1000`, `inner=499`, `seed=20261108`) kullanıldı; bu resmî 2000×999 kalibrasyon kapısı değildir. R merkezleri önceden .80 null ve .81/.82/.83 yakın alternatifler olarak sabitlendi; doymuş .85/.90 kontrol satırları korundu.

| Senaryo / merkez | α | Resmî payda güç/ret | Farklı-indeks güç/ret |
|---|---:|---:|---:|
| single_active / .80 null | .01/.05/.10 | .000/.000/.002 | .000/.000/.006 |
| tie / .80 null | .01/.05/.10 | .000/.000/.000 | .000/.000/.000 |
| single_active / .81 | .01/.05/.10 | .000/.002/.031 | .000/.009/.162 |
| tie / .81 | .01/.05/.10 | .000/.000/.000 | .000/.000/.025 |
| single_active / .82 | .01/.05/.10 | .001/.467/1.000 | .021/.808/1.000 |
| tie / .82 | .01/.05/.10 | .000/.210/.999 | .000/.661/1.000 |
| single_active / .83 | .01/.05/.10 | .390/1.000/1.000 | .644/1.000/1.000 |
| tie / .83 | .01/.05/.10 | .168/1.000/1.000 | .428/1.000/1.000 |
| iki senaryo / .85 ve .90 | tüm α | 1.000 | 1.000 |

Exact %95 binom GA'ları ham CSV'dedir.

**Sonuç:** Farklı-indeks düzeltmesi R null'unda nominali aşmamıştır ve test her iki paydada belirgin muhafazakârdır. Yakın alternatifte fark en büyük merkez .82, α=.05'te görülür: güç single_active için .467→.808, tie için .210→.661 yükselir.

## ADIM 4 — İsteğe bağlı yaprak değerleri

Yapılmadı. Faz 15'te ayrıştıran ağaç ve terminal splitler yeniden üretilebilse de 8-yollu attribution oyununda hangi background/permutation koalisyon durumunun ilgili yaprak çiftini oluşturduğu ara yol düzeyinde saklanmamıştır. “Sorgu yolu” ve “baz yolu” için tek bir terminal çıktı seçmek bu nedenle benzersiz değildir; tahminle doldurulmadı. Faz 15'in 128-yollu sıfır sonucu korunur.

## Komutlar ve bütünlük

```text
..\.venv\Scripts\python.exe experiments/72_phase16_platform_shift.py
..\.venv\Scripts\python.exe experiments/73_phase16_full_pool_ranking.py
..\.venv\Scripts\python.exe experiments/74_phase16_r_pair_sim.py
```

Başlangıç ve bitiş formül SHA256 değerleri aynıdır:

| Dosya | SHA256 |
|---|---|
| `components.py` | `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa` |
| `aggregation.py` | `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68` |
| `certification.py` | `76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487` |

Bu analiz dondurulmuş ayarları, Tablo 3–9'u, BH/BY kararlarını, saklı atıf tensörlerini, öznitelik panelini ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.
