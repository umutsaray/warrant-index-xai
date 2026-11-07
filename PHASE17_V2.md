# Faz 17 — ölçeğe normalize v2, exercised-signal bileşeni ve sentetik doğrulama

## Statü ve kapsam

Bu paket **v2 post hoc/keşifsel** analizidir; v1 sonuçları görüldükten sonra tasarlanmıştır. Birincil ayar ve 45 hücreli ızgara sonuçlara bakılmadan [ön-kayıtta](PHASE17_PREREG.md) sabitlenmiştir. Gerçek modellere yeniden fit uygulanmadı, özellik paneli değiştirilmedi ve v1 sonuç dosyalarının üzerine yazılmadı.

## ADIM 0 — ön-kayıt ve v1 yeniden üretim kapısı

Ön-kayıt `37f5749` commit'inde, herhangi bir v2 hesabından önce kaydedildi. Ölçekler Faz 12 tanımıyla yeniden üretildi:

| Model | `s_m` |
|---|---:|
| CoxPH-v2 | 0.023934215222165295 |
| RSF | 0.01707069437342852 |
| DeepSurv | 0.012307440008774696 |

Komut: `python experiments/75_phase17_v1_gate.py`. Yeni kod yolu v1 ayarları (`a=b=0.5`, `u=3`, `eta=0.02`, N yok) ile çalıştırıldı. `PHASE17_ADIM0_V1_REPRODUCTION.csv` içindeki 81/81 model-özellik satırında beş bileşen, `p_joint`, BH ve BY kararları resmî Faz 7.2 dosyalarıyla `atol=1e-12` düzeyinde eşleşti. Kapı geçti.

## ADIM 1 — atıf oyunundaki skaler çıktılar

Ayrıntılı kod referansları `results/phase17/PHASE17_ADIM1_OUTPUT_SCALARS.md` içindedir.

| Model | Skaler çıktı | Birim/yorum |
|---|---|---|
| CoxPH-v2 | `predict_log_partial_hazard`; kapalı form `beta_j * (z_j - background_j)` | log kısmi hazard |
| RSF | `RandomSurvivalForest.predict` | eğitim olay zamanlarındaki tahmini kümülatif hazard toplamına eşdeğer boyutsuz risk skoru |
| DeepSurv | ağın tek skaler çıktısı | log-risk / log göreli hazard |

Dolayısıyla üç modelin atıf büyüklükleri aynı fiziksel ölçekte değildir; model-özgü tek `s_m` kullanımı bu farkı azaltır, modeller arası birim eşitliği iddia etmez.

## ADIM 2 — 45 hücreli v2 nokta ızgarası

Komut: `python experiments/76_phase17_point_grid.py`. `PHASE17_ADIM2_POINT_GRID.csv`, `3 model x 45 ayar x 27 özellik = 3645` ham satır içerir.

Birincil ayar (`k=1`, `rho=1`, `kappa=0.25`):

| Model | WI5 point | WI6 point | WI6==R | WI6 darboğazları |
|---|---:|---:|---:|---|
| CoxPH-v2 | 1 | 0 | 0/27 | G_claim 20, N 5, S 2 |
| RSF | 10 | 0 | 0/27 | N 15, G_claim 9, S 3 |
| DeepSurv | 0 | 0 | 0/27 | G_claim 20, S 5, F 1, N 1 |

Tüm ızgarada WI6 point sayısı CoxPH-v2 için 0–5, RSF için 0–3, DeepSurv için 0–1 aralığındadır. R, 45 hücrenin yalnız CoxPH-v2'de 3'ünde en sık darboğazdır; RSF ve DeepSurv'de hiçbir hücrede en sık darboğaz değildir. Böylece v1'deki 75/81 `WI==R` deseni model-ölçekli toleranslarda kaybolur; bu, v1 R-baskınlığının tolerans ölçeğine duyarlı olduğunu gösterir.

## ADIM 3 — v2 sertifikasyon

### Birincil ayar

Komut: `python experiments/77_phase17_v2_certification.py`. B=4,999, v1 tohum zinciri ve `recursive_studentized` kullanıldı. Birincil hesap 8.9 saniye sürdü.

| Model | v1 point | v1 BH/BY | v2 point | v2 BH/BY |
|---|---:|---:|---:|---:|
| CoxPH-v2 | 9 | 8/8 | 0 | 0/0 |
| RSF | 17 | 12/12 | 0 | 0/0 |
| DeepSurv | 0 | 0/0 | 0 | 0/0 |

V1'de sertifikalı 20 model-özellik çiftinin tamamı v2'de düştü; yeni sertifika yoktur. Buna RSF'nin sekiz exact-zero/unexercised warrant'ı ve AJCC Stage I dahildir. N, sıfır sinyali doğrudan engellediği için `Q_j>0` post hoc etiketi v2 tanımında kavramsal olarak gereksizdir; ancak birincil v2'de hiç sertifika bulunmadığından “informative v2 sertifikası” gözlemsel olarak doğrulanamamıştır. C-index sırası 0.6677/0.5385/0.4965 iken v2 sertifika sayıları 0/0/0'dır; aynı sıralamayı üretmez.

### Ön-kayıtlı komşu ayarlar

Komut: `python experiments/81_phase17_v2_neighbors.py`; 83.6 saniye. `PHASE17_ADIM3_NEIGHBORS.csv` 324 satır içerir.

| Ayar `(k,rho,kappa)` | Cox point/BH/BY | RSF point/BH/BY | DeepSurv point/BH/BY |
|---|---:|---:|---:|
| (0.5,1,0.25) | 0/0/0 | 0/0/0 | 0/0/0 |
| (2,1,0.25) | 1/0/0 | 0/0/0 | 0/0/0 |
| (1,1,0.10) | 0/0/0 | 0/0/0 | 0/0/0 |
| (1,1,0.50) | 0/0/0 | 0/0/0 | 0/0/0 |

Komşu ayarlar birincil 0/0/0 sertifikasyon sonucunu değiştirmedi.

## ADIM 4 — taban çizgileri

Komut: `python experiments/78_phase17_baselines.py`. B1, 50 replikanın ortalama mutlak atıf yüzdelik %95 GA'sıdır; B2 yalnız bir **rank belirsizliği yaklaşımıdır**, Goldwasser–Hooker yönteminin uygulaması değildir; B3 sertifikasız nokta indeksidir.

| Model | B1 önemli | B2 top-10 rank-kararlı | B3 v1 point | B3 v2 point | v1 BH | v2 BH |
|---|---:|---:|---:|---:|---:|---:|
| CoxPH-v2 | 19 | 2 | 9 | 0 | 8 | 0 |
| RSF | 13 | 4 | 17 | 0 | 12 | 0 |
| DeepSurv | 18 | 3 | 0 | 0 | 0 | 0 |

RSF'deki sekiz unexercised warrant'ın tamamı B1 ve B2 tarafından işaretlenmedi. AJCC Stage I de B1/B2 tarafından işaretlenmedi; Stage II/III B1 tarafından işaretlendi. B1 sıfır sinyali elemekte yararlı olsa da modelden modele 13–19 özellik işaretleyerek claim-level kararlılık/taşınabilirlik sertifikası sağlamaz. V2 ise bu veri ve birincil ayarda aşırı seçicidir: hiçbir özelliği point veya FDR düzeyinde geçirmez.

## ADIM 6 — küçültülmüş N null kalibrasyonu

Öncelik sırasına uygun olarak sentetik adımdan önce çalıştırıldı. Komut: `python experiments/79_phase17_n_calibration.py`. Bu, resmî kapı olmayan `outer=1000`, `inner=499` küçültülmüş bataryadır. Üç model, iki senaryo ve üç alpha için 18 kontrolün hiçbirinde Holm-düzeltilmiş anti-conservative flag yoktur.

| Model | single_active Tip-I (`.01/.05/.10`) | tie Tip-I (`.01/.05/.10`) |
|---|---|---|
| CoxPH-v2 | .009 / .043 / .093 | .000 / .003 / .012 |
| RSF | .005 / .049 / .098 | .000 / .004 / .011 |
| DeepSurv | .007 / .036 / .091 | .000 / .000 / .007 |

Bu adım yeni N yaprağını v2 amplitüdünde sınar; beş eski bileşenin yeni toleranslar altındaki bütün ortak dağılımını yeniden kalibre eden resmî 2000x999 kapı değildir. N için nominal aşım gözlenmedi, tie senaryosu belirgin biçimde muhafazakârdır.

## ADIM 5 — gerçek-değerli sentetik özellik doğrulaması

Komutlar:

```text
python experiments/80_phase17_synthetic_feature_level.py --repeats 20
python experiments/80_phase17_synthetic_feature_level.py --repeats 100  # maliyet kapısında durduruldu
```

Gauss doğrusal DGP'de `n=700`, `p=30`, eş-korelasyon 0.3, iki kohort, T=50 ridge-OLS bootstrap ve tam 700-sorgu kullanıldı. İç bootstrap B=999'dur ve resmî 4,999 kapısı değildir. 20 tekrar 2,029.6 saniye sürdü. 100 tekrarlı genişletme 10/100 noktasındaki gerçek hızla 4 saati aşacak biçimde projekte edildiği için ön-kayıtlı maliyet kuralıyla durduruldu; nihai dosya eksiksiz 20-tekrarlı, 600 satırlık çıktıdır.

Oranlar pooled özellik-tekrar düzeyindedir; parantezler Clopper–Pearson exact %95 GA'dır.

| Yöntem | A güç | B yanlış | C yanlış | D yanlış | E (belirsiz) |
|---|---:|---:|---:|---:|---:|
| v1 BH | 1.000 (.977–1) | 1.000 (.970–1) | 1.000 (.970–1) | 1.000 (.964–1) | 1.000 (.964–1) |
| v1 BY | 1.000 (.977–1) | 1.000 (.970–1) | 1.000 (.970–1) | 1.000 (.964–1) | 1.000 (.964–1) |
| v2 BH | .231 (.168–.304) | .000 (.000–.030) | .008 (.000–.046) | .180 (.110–.269) | .220 (.143–.314) |
| v2 BY | .094 (.053–.150) | .000 (.000–.030) | .000 (.000–.030) | .090 (.042–.164) | .120 (.064–.200) |
| B1 | 1.000 (.977–1) | 1.000 (.970–1) | 1.000 (.970–1) | 1.000 (.964–1) | 1.000 (.964–1) |
| B2 | .175 (.120–.243) | .000 (.000–.030) | .000 (.000–.030) | .110 (.056–.188) | .000 (.000–.036) |

Pooled özellik-düzeyi FDR, yanlış sınıflar B+C+D'nin tüm sertifikalar içindeki payı olarak v1 BH/BY'de 340/600=.567; v2 BH'de 19/78=.244 ve v2 BY'de 9/36=.250'dir. E belirsiz sınıfı FDR payına katılmamıştır.

Önceden yazılı hipotez yalnız kısmen desteklenmiştir: v1 B'yi daima sertifikalarken v2 B'yi hiç sertifikalamamıştır; fakat v2 A gücü v1'e göre kayda değer biçimde düşmüştür ve D yanlış sertifikaları sürmüştür. Simetrik Gauss sorgularda işaret terslemesi atıf dağılımını simetrik bıraktığından quantile-coupled signed ölçü D'yi güvenilir biçimde ayıramamıştır. Dolayısıyla bu sentetik sonuç WI v2 üstünlüğünü doğrulamaz; N sıfır-sinyal kusurunu giderirken belirgin güç kaybı ve mekanizma-kayması kör noktası yaratır/gösterir.

## ADIM 7 — deletion testi

**Yapılamadı:** saklı 27-feature tensörleri, 616 girdili modellerde ortak bir deletion yolu tanımlamaz. CoxPH kohort-içi sıfır-z bazını, RSF 16 hastalık sampled inactive-background oyununu, DeepSurv ise candidate/reference için 16/128 background bütçelerini kullanır. Bu üç farklı oyunu tek bir deletion AUC'sinde birleştirmek yeni ve ön-kayıtsız bir baseline seçimi gerektirirdi. F'nin tanımı değiştirilmedi.

## Bütünlük, çıktılar ve değişiklik kaydı

| Dosya | Başlangıç SHA256 | Bitiş SHA256 |
|---|---|---|
| `warrantxai/components.py` | `7887CF447958726DF0FC2FD17FFC5BC4D70B00DE92B0888F36D3C4F49057E7AA` | aynı |
| `warrantxai/aggregation.py` | `10179864A63104ECA5AEEC8145983F0D25B2539A191A5B2E2095EB4B3BF5FF68` | aynı |
| `warrantxai/certification.py` | `76AB49716CC3EAD6DEC71656BBD60C6D71C1D76CAB8A83B2F0662AF6D6A05487` | aynı |

Ham çıktılar `results/phase17/` altındadır. Faz commitleri ayrı ayrı push edilmiştir: ön-kayıt `37f5749`, v1 kapısı `aa5aba5`, skalerler `23076fb`, ızgara `6693f77`, birincil sertifikasyon `7a4d8f1`, taban çizgileri `f02a2b1`, N kalibrasyonu `0bc7bc0`, sentetik `3c46aef`, komşular `62b1d22`.

Regresyon doğrulaması: Faz 17'ye özgü 5/5 test geçti. Tam paket ilk çağrıda 206 test geçtikten sonra yalnız Windows kullanıcı Temp klasörü izin hatasında durdu; çalışma alanındaki temiz `--basetemp` ile yinelenen nihai çağrıda **207/207 test geçti**.

Bu analiz v1'in dondurulmuş ayarlarını, Tablo 3–9'u, resmi BH/BY kararlarını, saklı atıf tensörlerini, öznitelik panelini ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR; v2 sonuçları ayrı, post hoc etiketli, keşifsel bir analizdir.
