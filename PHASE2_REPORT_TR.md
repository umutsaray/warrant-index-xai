# Warrant Index — Faz 2 sonuç raporu

17 Eylül 2026. Araştırma sürümü: `warrantxai 0.1.0`.

İstenen paket, 11 parametreli sentetik senaryo, beş çalıştırılabilir deney ve sklearn-style değerlendirme arayüzü tamamlandı. F/S/U/G/R ve WI formülleri değiştirilmedi. Gerçek veri, SHAP/LIME bağımlılığı veya uzaktaki bir servise yayın eklenmedi.

## Doğrulama

- **93 test geçti.** Paket modüllerinde 670 statement ve 182 branch için %100 kapsam; deney betikleri bu kapsam oranına dahil değildir. Kapsam, bilimsel geçerliliğin kanıtı değildir.
- Hypothesis ile min aksiyomları; geometrik ortalama/WI karşı örneği; elle hesaplanan bileşenler; eksik destek; oracle; IUT/BH/BY; subsampling/cache; API testleri mevcut. Ruff kontrolü geçti.
- Wheel ve kaynak dağıtımı üretildi. Ayrı, sistem paketlerini kullanmayan Windows/Python 3.14.7 ortamında wheel kuruldu; `pip check` temiz ve README örneği çalışıyor. Örnekte WI_claim=0.950326. Diğer Python sürümleri bu çalışmada denenmedi.
- Çıktılar `results/` içinde CSV, Markdown ve argüman/sürüm/kaynak SHA256 içeren JSON olarak saklandı. MIT lisansı eklendi.

## Deney 01: temel kalibrasyon çalışıyor

Beş seed, n={250,1000,4000}, rho={0,0.5,0.9}, p=5. Aşağıdaki değerler feature/seed/rho ortalamalarıdır.

| Senaryo | WI, n=250 | WI, n=1000 | WI, n=4000 |
|---|---:|---:|---:|
| Oracle, invariant | 0.947929 | 0.973799 | 0.986910 |
| Stabil ama yanlış | 0.208608 | 0.204446 | 0.204586 |
| Replikalar arasında bozulan | 0.472405 | 0.469631 | 0.469724 |

Oracle senaryosu tüm değerlendirmelerde 0.8 eşiğini geçti; iki bozulmuş senaryo geçmedi. Stabil ama yanlış açıklamada S/G/R=1 kalırken F düşüyor. Exact coalition hesabı ile kapalı form arasındaki en büyük fark 1.776e-15. Detay: `results/01_calibration.md`.

## Deney 02: beklentiler ve önemli istisnalar

20 seed, üç site, altı eğitim replikası, 1000 anchor ve beş feature kullanıldı. Sabit toleranslar a=b=0.5, u=3, eta=0.02, tau=0.8; sonuçlara göre ayarlanmadı. Bunlar örnek toleranslardır, evrensel eşikler değildir.

| Senaryo | Ortalama WI | Belirgin sonuç |
|---|---:|---|
| Oracle invariant | 0.973833 | F/S/G/R=1 |
| Stabil ama yanlış | 0.206514 | F düşük, S/R yüksek |
| Noisy | 0.080635 | U=0.089252; S=1 |
| Training instability | 0.145372 | S=0.245898, R=0.497500 |
| Covariate shift | 0.973833 | Ortak Q altında G=1 |
| Mechanism shift: etki kaybı | 0.783781 | Değişen feature için G=0; diğerleri geçebilir |
| Spurious shortcut | 0.039061 | G=0.039061 |
| Near-tied ranks | 0.843571 | R=0.847333; tüm feature'lar reddedilmiyor |
| Redundant duplicates | 0.599404 | Bireysel paylar değişirken grup toplamı sabit |
| OOD | NA | Ortak destek yok; değerlendirme yapılmıyor |
| Zero model | 1.000000 | Tekrarlanabilir sıfır açıklama; önem veya fayda kanıtı değil |

Önceden belirlenen yönsel/yapısal kontroller her senaryoda 20/20 geçti. Bu, C.2'deki her yorumun doğrulandığı veya her sorunlu feature'ın yakalandığı anlamına gelmez. Near-tied senaryosunda feature kararlarının %46'sı geçiyor; ortalama R tek bir belirsiz karşılaştırmayı gizleyebilir. Kalıcı eşitlikler R=1 üretebilir.

**Kritik karşı örnek: gerçek mekanizmada salt işaret değişimi yakalanmıyor.** Dış sitenin gerçek katsayısı ters çevrildiğinde, açıklama kendi modeline sadık kaldığı için F=1; G mutlak önem üzerinden tanımlandığı için G=1; S/R de 1 kalıyor. Değişen feature WI≈0.950396 ile eşiği geçiyor. Dolayısıyla C.2'nin işaret değişiminin G'yi düşüreceği beklentisi mevcut formülle sağlanmıyor. Formül değiştirilmedi. Sonraki insan kararı: iddiayı mutlak önem taşınabilirliğiyle sınırlamak veya ayrı, açıkça tanımlanmış yön duyarlılığı ölçüsü araştırmak.

Detaylar: `results/02_scenarios.md`, `02_pattern_checks.csv`, `02_sign_flip.csv`.

## Baseline, maliyet ve certification

Deney 03'te önceden sabitlenen claim seti tüm beş feature'dır; aşağıdakiler feature ortalaması değil, claim minimumunun seed ortalamasıdır. WI oracle için 0.950258, stabil ama yanlış için 0.065542; ERI-min iki durumda da 0.747535. Bu sentetik dönüşümde WI'nin fidelity kontrolü ayrım sağlıyor. Ancak ERI eksenleri ayrı hesaplandı ve kullanılan protokol küçük bir lineer stress testidir; ERI-Bench replikasyonu veya WI'nin genel üstünlüğünün kanıtı değildir. Tam tablo `results/03_baselines.md` içinde.

Deney 04'ün son koşusunda p=500 için R pair subsampling medyan hızlanması **47.74 kat**, ortalama mutlak hata 0.040279; S hızlanması **4.07 kat**, hata 0.040403. S'de 0.8 eşiği kararlarının medyan %7'si değişti; R'de bu özel eşikte değişiklik yoktu. Beş tekrarın medyanları ayrı ayrı raporlanır; medyan hızlanma, medyan sürelerin oranıyla birebir aynı olmak zorunda değildir. Süreler tek makineye özgü, model eğitimi/explainer maliyeti hariçtir. Bellek ölçümü process RSS değil tracemalloc tahsisleridir. Site minimumu, unbiased site-mean tahmininden sonra bile örnekleme yanlılığı taşıyabilir.

Deney 05'te geçerli sentetik p-değerleri altında 200000 IUT tekrarında Type I oranı **0.049860**; %95 binom aralığı [0.048910, 0.050823]. 10000 ailede bağımsız BH FDR=0.025270, bağımsız BY=0.007117 ve ortak döndürülmüş-uniform bağımlılık altında BY=0.007273. Bunlar certification motorunun bilinen geçerli girdilerle kontrolüdür; WI puanlarından otomatik geçerli p-değeri üretildiğini göstermez. Point adequacy ile istatistiksel certification ayrıdır.

## Varsayımlar ve kapsam

assumption: oracle hedefi marginal lineer Shapley; logistic modelde probability değil logit ölçeği. Ortak anchor Q site karşılaştırmasını standartlaştırır. API'de normal CI, sabit modele koşullu bağımsız anchor varsayımı kullanır; noisy deneyinde ayrıca bilinen Monte Carlo belirsizliği hesaba katılır. Harici açıklayıcılar uygun aralıklarını doğrudan bileşen API'sine sağlamalıdır. OOD desteği bu sentetik deneylerde tasarımdan bilinir.

Arayüz `fit(models, data_by_site)`, `evaluate()` ve `feature_report()` sağlar; modeller önceden eğitilmiş olmalıdır. sklearn parametre/clone sözleşmesini destekler; standart supervised Pipeline transformer olduğu iddia edilmez. Bağımsız-block Hoeffding yardımcısı gerçek bağımsız uçtan uca çalışmalar gerektirir; aynı verinin bootstrap replikaları bunun yerine geçmez.

Çalıştırma komutları ve hızlı örnek README'dedir. Sonraki bilimsel karar işaret değişimi sınırıdır. PyPI/GitHub yayını, isim rezervasyonu ve gerçek veri seçimi yapılmadı.

## Faz 2.1 eki — signed/claim generalizability

17 Eylül 2026. Yukarıdaki Faz 2 sonuçları magnitude-only sürümün tarihsel kaydıdır; işaret değişimiyle ilgili eski sınır varsayılan claim modu için aşağıdaki sonuçla güncellenmiştir. Diğer formüller değiştirilmedi.

Verilen tanımlar birebir uygulandı:

`G_signed_j = min_d>=1 [1 - mean_i min(1, |phi_hat_dj(x_i)-phi_hat_0j(x_i)|/a_j)]`

`G_claim_j = min(G_j, G_signed_j)`

`component_scores` ve `WarrantIndexEvaluator` artık varsayılan `generalizability_mode="claim"` ile F/S/U/G_claim/R üretir. `"magnitude_only"` eski G'yi açıkça etiketlenen sensitivity olarak korur. Ortak destek yoksa her iki modda da NaN / not assessable döner. assumption: ensemble tensöründe phi_hat, aynı anchor üzerindeki eğitim replikalarının ortalamasıdır; bu ortalama mutlak fark ve kırpmadan önce alınır.

Deney 02 yeniden çalıştırıldı (20 seed; ek sign-flip kontrolü seed=20260916, n=1000). Değişen ilk feature için:

| Mod | Magnitude G | G_signed | G_claim | WI | WI > 0.8 |
|---|---:|---:|---:|---:|---|
| claim (varsayılan) | 1.000000 | 0.052177 | 0.052177 | 0.052177 | Hayır |
| magnitude_only sensitivity | 1.000000 | 0.052177 | 0.052177 | 0.950396 | Evet |

İkinci satırdaki G_claim tanısal değerdir; o modun WI hesabına magnitude G girer. F/S/R=1 ve U=0.950396 değişmedi. Beklenen sign-flip reddi gerçekleşti; formül değiştirmeyi gerektiren beklenmedik sonuç çıkmadı. 11 senaryonun yönsel kontrolleri yine 20/20 geçti. Training-instability senaryosunda G ortalaması 0.282204'ten 0.176421'e, WI 0.145372'den 0.128389'a düştü; daha kısıtlayıcı G_claim tanımının beklenen etkisidir.

Yeni çıktılar eski sonuçlarla karışmaması için `results/phase21/` altında saklandı; `02_sign_flip.csv` iki modu ve G'nin üç tanımını yan yana verir. Tekrar komutu: `python experiments/02_failure_scenarios.py --output results/phase21`.

Eski 93 testin tamamı, sign-flip testinin yeni varsayılanı açıklayan yorumla güncellenmesiyle geçti; üç yeni testle toplam **96 test**, paket kodunda **%100 satır/dal kapsamı** elde edildi. Elle hesaplanan signed örnek, anchor öncesi clipping, en kötü site, support eksikliği, replica-mean sırası, mode doğrulaması ve API üzerinden eski sensitivity davranışı kontrol edildi.

Sınır: G_signed toleransa göre küçük signed farkları veya replikalar ortalamasında birbirini götüren etkileri yine kabul edebilir; her işaret değişiminin koşulsuz reddini garanti etmez. Ortak anchor/background ve destek varsayımları geçerliliğini korur. Point adequacy istatistiksel certification değildir. Bu ek kaynak kodu ve Deney 02 kapsamındadır; önceki `dist/` arşivleri Faz 2 sürümüdür, Faz 2.1 kurulumu için güncel kaynak kullanılmalıdır.
