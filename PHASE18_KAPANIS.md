# Faz 18 — tolerans profili, sentetik ölçek ve eşikten bağımsız ayırt etme

## Kapsam ve statü

Bu çalışma v1 ve Faz 17 v2 birincil analizlerini değiştirmeyen, post hoc/keşifsel bir tanıdır. Izgaralar, 60'ar CAL/EVAL tohumu, k* kuralı ve H1–H5 herhangi bir Faz 18 sonucu hesaplanmadan `099e186` commit'inde ön-kayıtlandı. Sağ-sansürlü sıra ve iç-bootstrap tohum ayrıntıları sonuçlarından önce `EK` olarak ayrıca commit edildi.

## ADIM 0 — yeniden üretim kapısı

Komut: `python experiments/82_phase18_gate.py`.

- Faz 17 birincil v2: 81/81 model-özellik satırında altı bileşen, WI6 ve `p_joint` `atol=1e-12`; BH/BY kararları tam eşleşti.
- Faz 17 sentetik: aynı `20261108` RNG zinciriyle 20 tekrar gerçekten yeniden çalıştırıldı; 600/600 özellik-tekrar v1/v2 BH/BY kararı birebir eşleşti.
- İlk süreç 5/20 sonrasında Windows dış süreç kesintisiyle durdu; bilimsel assertion hatası yoktu. Tekrar/RNG checkpoint'i eklenip aynı zincir baştan çalıştırıldı ve kapı geçti.
- Başarılı kapı duvar süresi: **3.437,7 s**.

Kapı özeti `results/phase18/PHASE18_ADIM0_GATE.csv` içindedir.

## ADIM 1 — DGP ölçeği ve bileşen davranışı

Komut: `python experiments/84_phase18_synth_breakdown.py`; süre **53,9 s**.

### DGP

Faz 17 kodu aynen kullanıldı: iki kohort, her kohortta `n=700`, `p=30`, eş-korelasyon `0.3`; ridge-OLS (`lambda=.1`) ve kapalı-form doğrusal attribution `beta_hat_j*(x_j-mean(x_j))`; T=50 bootstrap fit. Katsayılar:

| Sınıf | n özellik | Kaynak beta | Dış beta | Tanım |
|---|---:|---:|---:|---|
| A | 8 | .05 | .05 | güvenilir sinyal |
| B | 6 | 0 | 0 | sıfır etki / unexercised |
| C | 6 | .008 | .008 | zayıf/gürültülü |
| D | 5 | .05 | -.05 | işaret mekanizma kayması |
| E | 5 | .03 | .03 | yakın-eşit sıralar |

Sentetik `s_m_sim` medyanı **0.02553** (20 tekrar ortalaması 0.02534; min–maks 0.01993–0.02821). Sınıf medyan `M_j`: A .03995, B .00367, C .00622, D .03900, E .02454.

| Ölçek | v1 `a/s_m` | v1 `u/s_m` | v1 `eta/s_m` | v2 `a/s_m` |
|---|---:|---:|---:|---:|
| Sentetik, medyan | 19.59 | 117.53 | .784 | 1 |
| CoxPH-v2 gerçek | 20.89 | 125.34 | .836 | 1 |
| RSF gerçek | 29.29 | 175.74 | 1.172 | 1 |
| DeepSurv gerçek | 40.63 | 243.76 | 1.625 | 1 |

Dolayısıyla v1 toleransı sentetikte `s_m`'nin yaklaşık **19,6 katı**, gerçek veride yaklaşık **21–41 katıdır**.

### Bileşen dökümü

Aşağıdaki değerler medyan (IQR) ve tek-bileşen `>=.8` geçiş oranıdır:

| Sınıf | v2 G_claim | geçiş | v2 N | geçiş | baskın v2 darboğaz |
|---|---:|---:|---:|---:|---|
| A | .839 (.747–.891) | .612 | 1 (1–1) | 1.000 | G_claim 151/160 |
| B | .926 (.848–.961) | .892 | .474 (.432–.535) | .000 | N 119/120 |
| C | .857 (.785–.925) | .708 | .651 (.526–.805) | .275 | N 100/120 |
| D | .841 (.760–.882) | .640 | 1 (1–1) | 1.000 | G_claim 98/100 |
| E | .832 (.727–.912) | .570 | 1 (1–1) | 1.000 | G_claim 81/100 |

F, S, U ve R bütün sınıflarda tek başına %100 geçti. V1'de beş bileşenin tamamı bütün A–E satırlarında geçti ve her satırın darboğazı R idi. V2 A gücünü R veya N değil, ağırlıkla **G_claim** düşürdü; H2 reddedildi.

Ham değerler `PHASE18_ADIM1_DGP_SCALE.csv` ve `PHASE18_ADIM1_COMPONENT_BREAKDOWN.csv` içindedir.

## ADIM 2 — sentetik k x kappa profili

### Nokta düzeyi

Komut: `python experiments/85_phase18_synth_point_grid.py`; 60 ayrık CAL seed, **743,4 s**, 52.200 ham satır.

Birincil `kappa=.25` altında point-geçiş oranları:

| k | A | B | C | D | E |
|---:|---:|---:|---:|---:|---:|
| .25 | 0 | 0 | 0 | 0 | 0 |
| .5 | .169 | 0 | .075 | .137 | .200 |
| 1 | .660 | .006 | .214 | .637 | .653 |
| 2 | .944 | .006 | .244 | .920 | .943 |
| 4 | 1.000 | .006 | .244 | 1.000 | 1.000 |
| 8 | 1.000 | .006 | .244 | 1.000 | 1.000 |
| 16 | 1.000 | .006 | .244 | 1.000 | 1.000 |

N'siz WI5 k>=4'te bütün A–E sınıflarını %100 geçirir. Kappa=.1 B'yi k>=4'te .328 oranında geçirirken kappa=.5 B'yi tüm hücrelerde sıfırda tutar; N eşiği sonuçta belirleyicidir.

### Eşleşmiş mutlak ölçek

Her tekrarda `k_v1=.5/s_m_sim` ve `rho_v1=.02/s_m_sim` kullanıldığında v2-WI5, v1 bileşenlerini `atol=1e-12` ile aynen verdi ve A–E'nin tamamını %100 geçirdi. N eklendiğinde A/D/E yine %100, B 2/360, C 88/360 geçti. V1'in sentetik doygunluğu esas olarak mutlak tolerans ölçeği artefaktıdır; N yalnız sıfır/zayıf sinyali ayırır, D mekanizma kaymasını ayırmaz. H1 desteklendi.

### Sınırlı B=999 bootstrap

21 hücrelik tam planın ölçülen projeksiyonu altı saati aştığı için ön-kayıtlı dar plan kullanıldı: ilk 20 CAL seed, k=1/2/4, kappa=.25. Komut `python experiments/88_phase18_boot_cells.py`; **3.336,5 s**.

| k | A BH/BY | B BH/BY | C BH/BY | D BH/BY | E BH/BY |
|---:|---:|---:|---:|---:|---:|
| 1 | .325/.169 | 0/0 | .033/0 | .310/.180 | .310/.180 |
| 2 | .900/.863 | 0/0 | .083/.017 | .890/.850 | .870/.790 |
| 4 | .994/.994 | 0/0 | .083/.017 | 1/1 | .980/.970 |

Gevşek k A gücünü artırır, fakat D/E'yi de aynı hızda sertifikalar. Ham sonuç `PHASE18_ADIM2_BOOT_CELLS.csv` içindedir.

## ADIM 3 — eşikten bağımsız AUC

Komut: `python experiments/86_phase18_auc.py`; 9.999 repeat-cluster bootstrap, seed 18180001, **4,6 s**.

Birincil k=1, kappa=.25 için B'ye karşı özellik-tekrar AUC:

| Pozitif sınıf | AUC | %95 bootstrap GA |
|---|---:|---:|
| A | .968 | .957–.978 |
| C | .780 | .750–.809 |
| D | .952 | .934–.968 |
| E | .970 | .958–.981 |

Küçük/orta k'da A ve E, B'den ayrılır; ancak D de A ile neredeyse aynı ayrılır. H3'ün D için AUC≈.5 beklentisi reddedildi. Profil sıfır-etkiyi buluyor fakat mekanizma kaymasını “kötü” olarak ayırmıyor.

## ADIM 4 — önceden yazılı k* kuralı

CAL kuralı, kappa=.25 altında B yanlış-geçişi <=.05 olan en küçük grid değerini seçti: **k*=.25**. CAL B oranı k=.25/.5'te 0, k>=1'de .0056 idi; kural güç içermediği için en katı ilk hücreyi seçti.

Komut `python experiments/87_phase18_kstar.py`; 60 held-out EVAL seed, **83,8 s**. EVAL point-geçişi A–E sınıflarının tamamında 0'dı. Gerçek veriye keşifsel transferde CoxPH-v2/RSF/DeepSurv sonuçları 0/27, 0/27, 0/27; v1'in 20 sertifikasının 0'ı kaldı.

Maliyet izin verdiği için ilk 20 EVAL seed'de B=999 bootstrap da yapıldı (`89_phase18_kstar_eval_boot.py`; **965,5 s**): A–E'nin tamamında BH=BY=0.

H4 biçimsel olarak sağlandı (B yanlış-geçişi 0 ve gerçek veri <=5/27), fakat seçilen nokta bütün sinyali bastıran trivial bir çözümdür; savunulabilir/pratik bir işletim noktası değildir.

## ADIM 5 — gerçek veri k-profili ve giriş k'sı

Komut: `python experiments/83_phase18_kprofile.py`; frozen tensörler, yeni fit yok, **137,8 s**.

Kappa=.25 point-geçiş eğrisi:

| k | CoxPH-v2 | RSF | DeepSurv |
|---:|---:|---:|---:|
| .25 | 0 | 0 | 0 |
| .5 | 0 | 0 | 0 |
| 1 | 0 | 0 | 0 |
| 2 | 1 | 0 | 0 |
| 4 | 2 | 0 | 0 |
| 8 | 5 | 0 | 0 |
| 16 | 6 | 0 | 0 |
| 32 | 8 | 2 | 0 |
| 64 | 8 | 2 | 0 |

0.25–128 arasındaki 60 noktalı ince profilde CoxPH-v2 8/27 özellik için girdi (`k_enter` medyan 4.83, IQR 4.23–18.18); RSF 2/27 (medyan 30.75, IQR 29.94–31.57); DeepSurv 0/27. V1'in CoxPH sekiz sertifikasının üçü girdi, RSF on iki sertifikasının hiçbiri girmedi. Sekiz RSF unexercised warrant ve AJCC Stage I üç modelde de girmedi.

| Model | Spearman k_enter–v1 WI | Spearman k_enter–M_j |
|---|---:|---:|
| CoxPH-v2 | -.478 | -.278 |
| RSF | .127 | -.267 |
| DeepSurv | tanımsız | tanımsız |

CoxPH–RSF aynı-özellik Kendall tau **.178**'dir. DeepSurv'ün 27/27 çifti sağ-sansürlü olduğundan diğer iki Kendall tau tanımsızdır. H5 Cox–RSF için desteklendi, fakat üç-model iddiası DeepSurv sansürü nedeniyle tam sınanamadı. K-profili, hangi özelliklerin hangi tolerans gevşekliğinde girdiğini gösterir; sıralama Cox–RSF arasında zayıf ve DeepSurv için gözlenemezdir.

Ham dosyalar `PHASE18_ADIM5_KPROFILE.csv` ve `PHASE18_ADIM5_KENTER_81.csv` içindedir.

## ADIM 6 — deletion fizibilitesi

**Yapılamadı:** frozen CoxPH-v2 modeli ve 1.466 METABRIC özellik satırı mevcuttur, fakat bu satırlarla hizalı METABRIC süre/event hedefleri ve kimlik eşlemesi saklanmamıştır. C-index düşüşü outcome olmadan tanımlı değildir. Yeni dış veri indirme/eşleme kapsam dışı olduğundan `PHASE18_ADIM6_DELETION.csv` üretilmedi. RSF/DeepSurv için de yeni tahmin oyunu kurulmadı.

## Hipotez kararı

| Hipotez | Karar | Sayısal gerekçe |
|---|---|---|
| H1 | desteklendi | v1 a/s_m=19.59; matched WI5 A–E %100 |
| H2 | reddedildi | A darboğazı G_claim 151/160; N ve R değil |
| H3 | kısmen/reddedildi | A/E AUC yüksek, fakat D AUC=.952, ≈.5 değil |
| H4 | biçimsel destek, pratikte trivial | k*=.25; B=0 ama A–E de 0 |
| H5 | kısmi destek | Cox–RSF tau=.178; DeepSurv tau tanımsız |

## Bütünlük ve süre

| Dosya | Başlangıç/bitiş SHA256 |
|---|---|
| `components.py` | `7887CF447958726DF0FC2FD17FFC5BC4D70B00DE92B0888F36D3C4F49057E7AA` |
| `aggregation.py` | `10179864A63104ECA5AEEC8145983F0D25B2539A191A5B2E2095EB4B3BF5FF68` |
| `certification.py` | `76AB49716CC3EAD6DEC71656BBD60C6D71C1D76CAB8A83B2F0662AF6D6A05487` |

Faz 17 modülleri `75`–`81` için `git diff f8640d3` boştur; dosyalar değiştirilmedi. Başarılı ana adımların toplam ölçülen süresi yaklaşık **2,44 saat**tir (başarısız ilk kapı süreci hariç). Bütün ham Faz 18 sonuçları `results/phase18/` altındadır.

Regresyon doğrulaması: Faz 18'e özgü 4/4 test ve çalışma alanı `--basetemp` ile tam paket **211/211 test** geçti.

Bu analiz v1 ve v2 birincil ayarlarını, Tablo 3–9'u, BH/BY kararlarını, saklı atıf tensörlerini, öznitelik panelini ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.
