# Faz 18 ön-kayıt — tolerans profili, sentetik ölçek ve ayırt etme

**Kayıt zamanı:** hiçbir Faz 18 sonucu hesaplanmadan önce.  
**Etiket:** post hoc/keşifsel; v1 ve Faz 17 v2 birincil analizleri değişmez.

Bu belge Faz 18 ızgaralarını, tohumlarını, seçim/doğrulama ayrımını, karar kuralını ve hipotezlerini sonuçlar görülmeden sabitler. Sonradan yalnız açıkça `EK` başlığı altında genişletilebilir; mevcut tanımlar değiştirilemez.

## Donmuş analizler ve bütünlük kapıları

- V1: `tau=0.8`, `eta=0.02`, `a=b=0.5`, `u=3`.
- V2 birincil: `k=1`, `rho=1`, `kappa=0.25`; `a=b=k*s_m`, `u=6*k*s_m`, `eta=rho*s_m`, `m0=kappa*s_m`.
- Tablo 3–9, resmî BH/BY kararları, saklı attribution tensörleri, özellik paneli ve CoxPH-v2 retry protokolü değişmez.
- `components.py`, `aggregation.py`, `certification.py` ve Faz 17 modülleri değiştirilemez.
- Faz 17 v2 birincil 81 nokta satırı `atol=1e-12`; sentetik 20 tekrar A–E v1/v2 BH güçleri birebir yeniden üretilemezse Faz 18 durur.

## Ayrık tohum kümeleri

Faz 17 sentetik tohumu `20261108` ve onun tek-RNG zinciridir. Aşağıdaki kümeler hem birbirinden hem Faz 17'den ayrıdır.

`CAL_SEEDS` (ilk 20'si sınırlı bootstrap hücrelerinde, 60'ının tamamı nokta/AUC analizinde kullanılır):

```text
18000001, 18000002, 18000003, 18000004, 18000005, 18000006, 18000007, 18000008, 18000009, 18000010,
18000011, 18000012, 18000013, 18000014, 18000015, 18000016, 18000017, 18000018, 18000019, 18000020,
18000021, 18000022, 18000023, 18000024, 18000025, 18000026, 18000027, 18000028, 18000029, 18000030,
18000031, 18000032, 18000033, 18000034, 18000035, 18000036, 18000037, 18000038, 18000039, 18000040,
18000041, 18000042, 18000043, 18000044, 18000045, 18000046, 18000047, 18000048, 18000049, 18000050,
18000051, 18000052, 18000053, 18000054, 18000055, 18000056, 18000057, 18000058, 18000059, 18000060
```

`EVAL_SEEDS` (k* seçildikten sonra yalnız doğrulamada kullanılır):

```text
18010001, 18010002, 18010003, 18010004, 18010005, 18010006, 18010007, 18010008, 18010009, 18010010,
18010011, 18010012, 18010013, 18010014, 18010015, 18010016, 18010017, 18010018, 18010019, 18010020,
18010021, 18010022, 18010023, 18010024, 18010025, 18010026, 18010027, 18010028, 18010029, 18010030,
18010031, 18010032, 18010033, 18010034, 18010035, 18010036, 18010037, 18010038, 18010039, 18010040,
18010041, 18010042, 18010043, 18010044, 18010045, 18010046, 18010047, 18010048, 18010049, 18010050,
18010051, 18010052, 18010053, 18010054, 18010055, 18010056, 18010057, 18010058, 18010059, 18010060
```

## Sentetik DGP ve ızgaralar

Faz 17'nin `n=700`, `p=30`, eş-korelasyon `0.3`, iki kohort, T=50 ridge-OLS ve kapalı-form attribution DGP'si aynen kullanılır. Sınıflar A=8 güvenilir sinyal, B=6 sıfır etki, C=6 zayıf/gürültülü, D=5 mekanizma/işaret kayması, E=5 yakın-eşit sıralardır.

Nokta ızgarası:

- `k in {0.25, 0.5, 1, 2, 4, 8, 16}`
- `kappa in {0.1, 0.25, 0.5}`
- `rho=1`
- ayrıca N'siz WI5 ve v1 mutlak tolerans satırları.

Gerçek veri k-profili:

- kaba `k in {0.25, 0.5, 1, 2, 4, 8, 16, 32, 64}`;
- `kappa in {0.1, 0.25, 0.5}`, `rho=1`;
- giriş `k` değeri için 0.25–128 arasında doğal-log aralıklı 60 nokta.

Sınırlı bootstrap: ilk hücrenin süresinden toplam süre tahmin edilir. Tam plan 6 saati aşarsa yalnız `(k,kappa)=(1,.25),(2,.25),(4,.25)`, ilk 20 `CAL_SEEDS`, B=999 ile çalıştırılır; bu resmî B=4,999 kapısı değildir.

## Eşikten bağımsız ayırt etme

Her CAL tekrarı ve hücresinde WI6, WI5 ve v1 WI skorları için A/B/C/D/E sınıfını sıfır-etkili B'den ayırma ROC AUC'si ortalama sıralarla hesaplanır. %95 güven aralıkları tekrarların bootstrap'ı üzerinden hesaplanır; bootstrap tohumu `18180001`, 9,999 çekimdir. B-vs-B tanımlanmaz; raporda B referans sınıfı olarak gösterilir.

## Önceden yazılı işletim noktası kuralı

`k*`, `kappa=0.25`, `rho=1` altında CAL_SEEDS ile B sınıfının özellik-tekrar başına nokta yanlış-geçiş oranı `<=0.05` olan ızgaradaki **en küçük k**'dır. Hiçbiri sağlamazsa k* seçilmez. Kural sonuca göre değiştirilemez.

K* varsa 60 EVAL_SEEDS ile nokta gücü ve B yanlış-geçişi doğrulanır; bootstrap maliyeti izin verirse ilk 20 EVAL seed ile çalıştırılır. Gerçek veriye uygulama yalnız keşifseldir ve sentetik kalibrasyondan geldiği açıkça etiketlenir.

## Eşleşmiş ölçek

V1'in mutlak toleransları sentetik `s_m_sim` cinsinden ifade edilir. V2, v1 `a` toleransına karşılık gelen `k_v1=0.5/s_m_sim` noktasında; `rho_v1=0.02/s_m_sim` ve `u/a=6` korunarak ayrıca çalıştırılır. N dahil sonuç `kappa=0.25`; N'siz sonuç da ayrı raporlanır. Bu, v1 yüksek sertifikasyonunun ölçekten gelip gelmediğini tanısal olarak sınar.

## Önceden yazılı hipotezler

- **H1:** Sentetikte v1 tolerans seti `s_m_sim`'in en az 10 katıdır ve v2 formülü aynı mutlak noktada v1'in yüksek sertifikasyonunu yeniden üretir.
- **H2:** V2'nin düşük A gücünü baskın olarak tek bileşen (R veya N) belirler.
- **H3:** Küçük k (`<=2`) sıfır-etkili B'yi A/E'den ayırır; D AUC yaklaşık 0.5 kalır.
- **H4:** B yanlış-geçişini `<=.05` yapan bir k* vardır ve gerçek veride model başına en çok 5/27 özelliği point düzeyinde geçirir.
- **H5:** K-profil sıralamasının modeller arası tutarlılığı zayıftır (`Kendall tau < .3`).

## Karar ve raporlama kuralı

Birincil v1/v2 sonuçları ve k* kuralı değiştirilmeyecek; tüm hücreler ve olumsuz bulgular raporlanacaktır. Gerçek veri yeniden fit edilmeyecek. Yapılamayan alt adım `yapılamadı: <neden>` biçiminde kaydedilecektir.
