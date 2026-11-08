# Faz 19 — kapanış raporu

## Kapsam, ön-kayıt ve bütünlük

Bu faz tanısaldır. Ön-kayıt, herhangi bir Faz 19 sonucu üretilmeden önce `d839855` commit'i ile kaydedilip `main` dalına gönderildi. P19 tohumları `19000001..19000040` olarak ayrıldı; hesap gerektiren ADIM 1, 2 ve 4'te ön-kayıt uyarınca ilk 20 tohum kullanıldı. Gerçek veri yeniden fit edilmedi; v1/v2 sonuçları ve saklı tensörler salt okunur kullanıldı.

Başlangıç ve bitiş SHA256 değerleri aynıdır:

| Dosya | SHA256 |
|---|---|
| `warrantxai/components.py` | `7887CF447958726DF0FC2FD17FFC5BC4D70B00DE92B0888F36D3C4F49057E7AA` |
| `warrantxai/aggregation.py` | `10179864A63104ECA5AEEC8145983F0D25B2539A191A5B2E2095EB4B3BF5FF68` |
| `warrantxai/certification.py` | `76AB49716CC3EAD6DEC71656BBD60C6D71C1D76CAB8A83B2F0662AF6D6A05487` |

## ADIM 0 — kapı

Komut: `python experiments/90_phase19_gate.py`.

Dokuz kontrolün dokuzu geçti. Faz 18 k-profili girdilerinden CoxPH-v2=8, RSF=2 ve DeepSurv=0 point-adequate sayıları yeniden üretildi; seçilen iki zorunlu kontrol dahil tüm karşılaştırmaların mutlak hatası `1e-12` sınırının altındaydı. Çıktı 9/9 satırdır: `results/phase19/PHASE19_ADIM0_GATE.csv`.

## ADIM 1 — G_signed hizalama mekaniği

Komut: `python experiments/93_phase19_gsigned_alignment.py`.

Kod hizalaması: gerçek sertifikasyon yolundaki `experiments/12_rsf_certification.py:signed_scores`, önce replikaları ortalayıp sorgu eksenini `np.sort` ile sıralar. Faz 18 sentetik `point_components` aynı işlevi çağırır. Dolayısıyla iki yol da eşlenmiş anchor değil, quantile-coupled sıralı atıf dağılımlarını karşılaştırmaktadır.

Gaussian DGP özeti:

| Sınıf | G medyanı (IQR) | sorted G_signed medyanı (IQR) | paired G_signed medyanı (IQR) |
|---|---:|---:|---:|
| A | 0.855 (0.727–0.931) | 0.835 (0.725–0.889) | 0.834 (0.725–0.927) |
| D | 0.869 (0.765–0.941) | 0.841 (0.756–0.890) | 0.099 (0.091–0.106) |

Sorted Gaussian A-vs-D AUC 0.465 iken aynı 700 kaynak anchor'ı üzerindeki paired varyantta A-vs-D AUC 1.000 ve medyan ayrım 0.735 oldu. D sütunları standartlaştırılmış log-normal yapıldığında sorted varyant A-vs-D AUC 1.000'e çıktı (D medyan G_signed=0.190). Sonuç: H1 desteklendi; kör nokta `np.sort` ile simetrik özellik dağılımının etkileşimidir, frozen formül değiştirilmemiştir.

Sağlamalar: normal 600/600, skewed-D 600/600 ve AUC 40/40 satır. Çıktılar `PHASE19_ADIM1_GSIGNED_ALIGNMENT.csv`, `PHASE19_ADIM1_PAIRED_AUC.csv`, `PHASE19_ADIM1_SKEWED_D.csv`.

## ADIM 2 — A sınıfı SNR ve n-süpürmesi

Komut: `python experiments/94_phase19_snr_nsweep.py`.

DGP kodu `y_d = X_d beta_d + Normal(0, 0.1)`; raporlanan gerçek-DGP R² değeri `1-SSE_true/SST`'dir. Kaynak/dış site ortalamaları:

| Sınıf | beta_hat kaynak / dış | ort. bootstrap SD kaynak / dış | medyan z kaynak / dış | kaynak-dış fark SD |
|---|---:|---:|---:|---:|
| A | 0.05001 / 0.04999 | 0.00451 / 0.00459 | 11.20 / 11.07 | 0.00647 |
| B | 0.00061 / 0.00009 | 0.00454 / 0.00453 | 0.10 / 0.00 | 0.00645 |
| C | 0.00784 / 0.00776 | 0.00455 / 0.00461 | 1.64 / 1.78 | 0.00643 |
| D | 0.04970 / -0.04989 | 0.00459 / 0.00452 | 10.96 / -11.23 | 0.00645 |
| E | 0.03026 / 0.02968 | 0.00453 / 0.00452 | 6.78 / 6.54 | 0.00640 |

Sonuç varyansı kaynakta ortalama 0.2526, dış sitede 0.0706; gerçek-DGP R² sırasıyla 0.9596 ve 0.8585'tir.

| n | k | medyan G_claim | G_claim geçiş oranı | medyan WI6 | WI6 point geçiş oranı |
|---:|---:|---:|---:|---:|---:|
| 700 | 1 | 0.8351 | 0.6063 | 0.8351 | 0.6063 |
| 700 | 2 | 0.9169 | 0.9188 | 0.9162 | 0.9188 |
| 2800 | 1 | 0.9027 | 0.9375 | 0.9027 | 0.9375 |
| 2800 | 2 | 0.9514 | 1.0000 | 0.9514 | 1.0000 |
| 11200 | 1 | 0.9568 | 1.0000 | 0.9568 | 1.0000 |
| 11200 | 2 | 0.9784 | 1.0000 | 0.9784 | 1.0000 |

H2 desteklendi: k=1 G_claim geçişi n ile monoton arttı ve n=11200'de 1.00 oldu. Büyük-n hesabı 4-D tensörü belleğe almadan aynı cebirle yapıldı; n=700 tam-tensör referansına karşı F/S/U/G_claim/R/N/WI6 maksimum mutlak farkı `6.11e-15` idi. SNR 1200/1200, n-süpürmesi 960/960 satırdır.

## ADIM 3 — gerçek veri rho süpürmesi

Komut: `python experiments/91_phase19_real_rho.py`.

κ=.25, ρ=1 ve k=64 bağlayıcı dağılımları yeniden üretildi: CoxPH-v2 R=19/N=8; RSF N=17/R=10; DeepSurv R=17/N=10. DeepSurv R aralığı 0.624992–0.783438 olup 0/27 özellik 0.8'i geçti.

| Model | ρ=1 | ρ=2 | ρ=4 | ρ=8 | ρ=16 |
|---|---:|---:|---:|---:|---:|
| CoxPH-v2 point-adequate | 8 | 15 | 19 | 19 | 19 |
| RSF point-adequate | 2 | 3 | 6 | 7 | 10 |
| DeepSurv point-adequate | 0 | 1 | 9 | 13 | 13 |

ρ=16'da bağlayıcı dağılımları CoxPH-v2 N=20/G_claim=7; RSF N=22/G_claim=4/R=1; DeepSurv N=24/G_claim=2/S=1 oldu. H3 desteklendi: DeepSurv'de 13 özellik girdi. H4 yeniden doğrulandı: Cox giriş oranı medyanı 4.4167 (min 0.7654, maks 46.7562); RSF'de iki giren özellik 8.6541 ve 42.9331 (medyan 25.7936). RSF N-bağlayıcı 17 çift için Î_dj ve κs_m `PHASE19_ADIM3_RSF_N.csv` içindedir.

36-sorgu duyarlılığı yapılamadı: aynı dondurulmuş panel ve 8/128 bütçeyle uyumlu 36-sorguluk attribution tensörü mevcut değildir. Satır sağlamaları rho süpürmesi 405/405, RSF-N 17/17'dir.

## ADIM 4 — mutlak tolerans ve M_j

Komut: `python experiments/92_phase19_relative_g.py`.

Faz 18 özgün-20 değerleri yeniden üretildi: G_claim A/C/D/E-vs-B AUC'leri 0.235/0.335/0.230/0.269; N 1.000/0.806/1.000/1.000; M_j 1.000/0.805/1.000/1.000; WI6 0.958/0.804/0.969/0.976. G_claim A-vs-D AUC=0.504'tür. CAL60 ve bağımsız P19-first20 sonuçları aynı yönü korudu.

Göreli-tolerans tanısında A-vs-B AUC c=.25/.5/1 için 0.831/0.827/0.826 oldu; D için 0.844/0.843/0.841, E için 0.771/0.769/0.770. Ancak C her üçünde 0.453/0.453/0.450 ile 0.5'in altında kaldı. A geçiş oranı yalnız c=1'de mutlak G_claim oranını geçti (0.8063 vs 0.6063); c=.25 ve .5'te 0.0688 ve 0.4875 idi. Bu nedenle H5 bütünüyle desteklenmedi: göreli ölçekleme A/D/E ayrımını tersine çevirdi, fakat C'yi düzeltmedi ve sonuç c'ye duyarlıydı.

Satır sağlamaları: P19 600/600, CAL60 1800/1800, göreli G 1800/1800, AUC 120/120.

## Kesintiler ve yeniden başlatmalar

Başlangıç denetiminde önceki Faz 19 çalışmasına ait sonuç veya checkpoint bulunmadı; bütün adımlar ön-kayıt commit'inden sonra başladı. ADIM 4'te tüm seed'ler tamamlandıktan sonra A-vs-D satırlarının çok-indeks seçimine karışması rapor doğrulamasını durdurdu; checkpoint'ler kullanılarak hesap tekrarlanmadan yalnız filtre düzeltildi. ADIM 1'de yine tüm seed'ler tamamlandıktan sonra `median` sütun adının pandas yöntemiyle çakışması rapor üretimini durdurdu; checkpoint'ten yalnız raporlama yeniden çalıştırıldı. ADIM 2 kesintisiz tamamlandı. Hiçbir hücre atlanmadı veya başka tohumla değiştirilmedi.

## Hipotez kararı

| Hipotez | Karar | Sayısal gerekçe |
|---|---|---|
| H1 | Desteklendi | Paired A-vs-D AUC=1.000, medyan fark=0.735; sorted Gaussian AUC=0.465 |
| H2 | Desteklendi | k=1 A G_claim geçişi 0.606→0.938→1.000 |
| H3 | Desteklendi | DeepSurv ρ=16'da 13/27 point-adequate |
| H4 | Yeniden doğrulandı | Cox k_enter·s_m/M_j medyan=4.4167 |
| H5 | Desteklenmedi | C-vs-B AUC tüm c değerlerinde <0.5; A geçiş üstünlüğü yalnız c=1 |

Kod doğrulaması: altı Faz 19 betiği `py_compile` kontrolünü geçti; tam test paketi çalışma alanı içindeki açık `--basetemp` ile `211 passed` verdi. İlk test çağrısındaki tek kurulum hatası sistem geçici klasöründeki `WinError 5` erişim izniydi; kod/test başarısızlığı değildi.

Bu analiz v1 ve v2 birincil ayarlarını, Tablo 3–9'u, BH/BY kararlarını, saklı atıf tensörlerini, öznitelik panelini ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.
