# Faz 19 ön-kayıt — G_signed, D kör noktası, A tutarlılığı ve k-profili tanıları

**Kayıt zamanı:** hiçbir Faz 19 sonucu hesaplanmadan önce.
**Etiket:** post hoc/tanısal; v1 ve v2 birincil analizleri değiştirmez.

## Donmuş kapsam ve bütünlük

- Frozen v1, v2 birincil (`k=1`, `rho=1`, `kappa=.25`), Tablo 3–9, BH/BY kararları, saklı tensörler, panel ve CoxPH-v2 retry protokolü değiştirilemez.
- `warrantxai/components.py`, `warrantxai/aggregation.py`, `warrantxai/certification.py` ve Faz 17–18 modülleri değiştirilemez.
- Yeni betikler `experiments/90_phase19_*.py` ve devamı; yerel ham çıktılar `results/phase19/` altında oluşturulur.
- Bu faz tanısaldır. Yeni bir resmî `G_signed` formülü önerilmez veya ana analize uygulanmaz. Paired ve göreli varyantlar yalnız açıkça etiketli tanılardır.
- Gerçek veride yeniden model fit veya explainer çağrısı yapılmaz.
- Yapılamayan alt adım `yapılamadı: <neden>` biçiminde raporlanır; tahminle hücre doldurulmaz.

Başlangıç SHA256:

| Dosya | SHA256 |
|---|---|
| `warrantxai/components.py` | `7887CF447958726DF0FC2FD17FFC5BC4D70B00DE92B0888F36D3C4F49057E7AA` |
| `warrantxai/aggregation.py` | `10179864A63104ECA5AEEC8145983F0D25B2539A191A5B2E2095EB4B3BF5FF68` |
| `warrantxai/certification.py` | `76AB49716CC3EAD6DEC71656BBD60C6D71C1D76CAB8A83B2F0662AF6D6A05487` |

## Tohumlar

`P19_SEEDS = 19000001..19000040` (uçlar dahil). Bu küme Faz 17 ve Faz 18 tohumlarından ayrıdır. ADIM 1, ADIM 2 ve P19 tekrarları istenen yerde ilk 20 tohumu (`19000001..19000020`) kullanır. Hiçbir tohum atlanmaz veya sonuçlara göre değiştirilmez.

## ADIM 0 — yeniden üretim kapısı

Faz 18 ADIM 5 gerçek-veri k-profili girdilerinden, `atol=1e-12` ile şu iki ayrı sağlama yeniden üretilir:

1. `kappa=.25`, `rho=1`, ince `k` profili altında en az bir kez giren özellik sayıları: CoxPH-v2 8/27, RSF 2/27, DeepSurv 0/27.
2. `k=64`, `kappa=.25`, `rho=1` bağlayıcı dağılımları: CoxPH R 19/N 8; RSF N 17/R 10; DeepSurv R 17/N 10.

Bu kapı geçmeden ADIM 1–4 başlamaz.

## ADIM 1 — G_signed mekaniği ve D

1. Faz 17/18 sentetik modülünde `G_signed` girdisinin iki kohortta nasıl hizalandığı kod referansıyla kaydedilir: aynı sorgu/anchor, `np.sort` veya başka yöntem.
2. V2 birincil tolerans altında P19 ilk 20 tekrarda A ve D için `G`, `G_signed`, `G_claim` medyan ve IQR'ları raporlanır.
3. Tanısal paired varyantta iki kohort aynı 700 sorgu noktasında değerlendirilir ve mevcut `generalizability_signed` eşleşmiş `(sites,n,p)` girdisiyle çağrılır. A ve D özetleri ve B'ye karşı AUC raporlanır.
4. İkinci DGP kontrolünde yalnız D özellikleri çarpık, standartlaştırılmış log-normal dağılıma alınır; diğer A–E parametreleri korunur. Sıralanmış `G_signed` özetleri ve ayırt etme sonucu raporlanır.

**H1:** Sıralanmış `G_signed` D'yi ayıramaz, paired varyant ayırır: `AUC_D(paired) >= .9` ve paired `G_signed_D` medyanı A medyanından en az `.3` düşüktür.

## ADIM 2 — A sınıfı SNR ve n-süpürme

1. P19 ilk 20 tekrarda A–E için kaynak ve dış kohort `beta_hat` ortalaması, SD'si, `beta_hat/SE` z oranı, kohortlar arası `beta_hat` farkının SD'si, DGP sonuç varyansı ve R-kare raporlanır.
2. `n in {700, 2800, 11200}`, `p=30`, `k in {1,2}`, `kappa=.25`, v2 birincil diğer ayarlar ve P19 ilk 20 tekrar kullanılır. A sınıfı için `G_claim` medyanı, tek-bileşen geçiş oranı ve WI6 point geçişi hesaplanır. B=999 bootstrap yapılmaz.

Beklenen n-süpürme hücre sayısı: `3 n x 2 k x 20 seed = 120` seed-hücre özeti; özellik-düzeyi ham çıktı tutulursa A'nın 8 özelliğiyle `960` satır.

**H2:** Fark `beta_hat` örnekleme gürültüsünden geliyorsa A'nın `G_claim` geçiş oranı n ile monoton artar ve `n=11200` için en az `.95` olur. Artmazsa gözlenen değerlerle `artmadı` diye raporlanır.

## ADIM 3 — gerçek veri rho süpürmesi ve RSF N

Önce Faz 18 ham CSV'lerinden şu değerler yeniden üretilir:

- `k=64`, `kappa=.25`, `rho=1`: CoxPH R 19/N 8; RSF N 17/R 10; DeepSurv R 17/N 10.
- DeepSurv R aralığı `.625–.783`; R>=.8 olan özellik yoktur.
- `k_enter*s_m/M_j`: Cox medyan 4.41 (minimum .76 yaş, maksimum 46.7); RSF 8.67 ve 43.0.

Yeni rho ızgarası: üç model, `k=64`, `kappa=.25`, `rho in {1,2,4,8,16}`, 27 özellik. Beklenen özellik-düzeyi satır sayısı `3 x 5 x 27 = 405`tir. Model-rho özet tablosu tutulursa ayrıca 15 satırdır.

RSF'de N bağlayıcı çiftler için kohort başına `I_hat_dj` ve eşik `kappa*s_m` raporlanır. Uyumlu 36-sorgu çıktısı bulunursa N duyarlılığı hesaplanır; bulunmazsa neden belirtilerek yapılmadı yazılır.

**H3:** `rho=16`, `k=64` altında DeepSurv'de en az bir özellik point-adequate olur.

**H4:** Cox giriş oranı medyanı en az 3'tür. Bu hipotez mevcut ham veriden beklenen 4.41 değeriyle yalnız yeniden üretim sağlamasına tabidir.

## ADIM 4 — mutlak tolerans, G ve M_j

P19 ilk 20 tekrarda ve Faz 18 `CAL_SEEDS` 60 tekrarda şu tanılar yeniden üretilir:

- B'ye karşı A–E AUC: `G_claim` için A .235, C .335, D .230, E .269; `N` için A/D/E 1.000 ve C .806; `M_j` için A/D/E 1.000 ve C .805; `WI6` için .958/.804/.969/.976.
- A–D AUC her bileşende yaklaşık .5; `G_claim` .504.

Tanısal göreli tolerans varyantı yalnız yeni Faz 19 işlevinde hesaplanır:

`G_rel_j = 1 - mean_i min(1, abs(delta_phi_ij) / (c * Mbar_j))`,

burada `Mbar_j` iki kohortun ortalama mutlak atıf büyüklüğüdür ve `c in {.25,.5,1}`. A–E için medyan ve B'ye karşı AUC raporlanır. Mevcut formül modülleri değiştirilmez.

P19 ham özellik-tekrar-c tablosu tutulursa beklenen satır sayısı `20 seed x 30 feature x 3 c = 1,800`dür. CAL yeniden üretim tablosu özellik-tekrar düzeyindeyse `60 x 30 = 1,800` satırdır.

**H5:** `G_rel` için A–E vs B AUC değerleri en az `.5` ve A'nın geçiş oranı mutlak G'den yüksektir.

## Öncelik ve kesinti kuralı

ADIM 0 kapısından sonra çalışma önceliği: **ADIM 3 → ADIM 4 → ADIM 1 → ADIM 2**.

Her tamamlanan adım ham CSV satır sayısı, beklenen ızgara/tohum sayısı, komut, süre, assertion ve kesintilerle raporlanır. Checkpoint kullanan uzun adımlar RNG durumunu veya deterministik tamamlanan-seed listesini saklar. Tamamlanmış adımlar yeniden hesaplanmaz. Kesinti olursa checkpoint'ten devam edilir; checkpoint yoksa yalnız kesilen adım aynı tohumlarla baştan çalıştırılır.

## Karar ve raporlama

H1–H5 desteklenmese dahi tüm sonuçlar raporlanır. Kötü görünen hücre veya tohum atlanmaz. Yeni formül/varyantlar resmî WI'ye aktarılmaz. Faz 19 kapanış raporu şu cümleyle biter:

> Bu analiz v1 ve v2 birincil ayarlarını, Tablo 3–9'u, BH/BY kararlarını, saklı atıf tensörlerini, öznitelik panelini ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.
