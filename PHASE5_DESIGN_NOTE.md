# Faz 5 — istatistiksel sertifikasyon tasarımı ve pilotu

> **RSF koşulluluk uyarısı:** Bütün RSF p-değerleri sabit 12+12 sorgu kümesine ve saklanmış 8/128 attribution ortalamalarına koşulludur; tekil permütasyon yollarının belirsizliğini içermez.

## ADIM 0 — artefakt envanteri

CoxPH N=50 için 50×614 katsayı dizisi, feature sırası ve replika fit kayıtları saklıdır. Attribution tensörü ayrıca yazılmamıştır; ancak frozen exact Cox oracle `beta_j*z_j` ve diskteki kimliksiz sorgu girdileriyle bütün bileşen girdileri model/explainer çağrısı olmadan yeniden kurulabilir. TCGA sorgusu 141, METABRIC sorgusu 1466 satırdır.

RSF N=20 için candidate ve bağımsız yüksek bütçeli reference attribution tensörlerinin ikisi de 2×20×12×77 olarak saklıdır. Bütün frozen bileşenler bunlardan yeniden hesaplanabilir. Tek tek 8/128 permütasyon yolları saklı değildir; bu nedenle explainer Monte Carlo varyansı ayrıca resample edilemez.

## ADIM 1 — önerilen pilot yöntem

Hedef, gözlenen sabit audit sorgu kümeleri, preprocessing, ana model ve RSF explainer ortalamalarına koşullu; patient-cluster bootstrap eğitim dağılımı üzerindeki bileşen değeridir. Site sayısı iki sabit, önceden belirlenmiş kohort bağlamıdır ve site resample edilmez. Her bootstrap çekiminde replika indisi `t` replacement ile N kez seçilir, aynı indis iki sitede ve beş bileşende birlikte kullanılır, sonra frozen F/S/U/G_claim/R fonksiyonelleri aynen yeniden hesaplanır. S ve R içindeki bütün çiftler yeni çekimde yeniden oluşur.

Her bileşen için ortak frozen eşik `tau_k=0.8` kullanılır. Sınır-null altında basic bootstrap merkezi `0.8 + (theta* - theta_hat)` olur; sağ-kuyruk p-değeri `(1 + #[(theta* - theta_hat) >= (theta_hat - 0.8)])/(B+1)` olarak okunur. B=4999, dolayısıyla en küçük raporlanabilir değer 0.0002'dir. Studentization seçilmedi: N=20 içinde iç içe SE bootstrap'ı zayıf, ayrıca min, clipping ve eta-toleranslı sıralama fonksiyonelleri düzgün değildir. Bu pilot yöntem asimptotik bir öneridir; gerçek-verideki uç-durum sanity check sonlu-örneklem Tip-I geçerliliğinin kanıtı değildir. Tam sertifikasyondan önce null altında kapsama/Tip-I simülasyonu ayrıca onay kriteri olmalıdır.

Yöntem seçimi, Efron ve Tibshirani'nin bootstrap hipotez testi çerçevesindeki null-merkezleme fikrine dayanır ([Chapter 16](https://djvu.online/file/DrRxMumhAaOZY)). Monte Carlo kuyruğunun sıfır raporlanmaması için +1 düzeltmesi kullanılır; bu tercih Phipson ve Smyth'in random-resampling p-değerleri uyarısıyla uyumludur ([2010 makalesi](https://pubmed.ncbi.nlm.nih.gov/21044043/)). İlerideki çoklu test aşamasında BH yalnız independence/PRDS varsayımıyla primary, BY keyfi bağımlılık robustness kontrolü olacaktır ([BH 1995](https://rss.onlinelibrary.wiley.com/doi/pdf/10.1111/j.2517-6161.1995.tb02031.x); [BY 2001](http://projecteuclid.org/download/pdf_1/euclid.aos/1013699998)). Bu kaynaklar pilotun bileşen fonksiyonelleri için özel sonlu-örneklem geçerlilik ispatı sağlamaz.

## 5+5 feature kalibrasyon pilotu — RSF sabit 12+12 query ve saklı attribution ortalamalarına koşullu

| model | feature | selection_reason | F | S | U | G_claim | R | p_F | p_S | p_U | p_G_claim | p_R | p_joint |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH_N50 | clin_numeric__age_at_index.demographic | high_margin | 1.0000 | 0.9889 | 0.9909 | 0.9503 | 0.9690 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.0002 |
| CoxPH_N50 | rna_ENSG00000139874.6 | near_threshold | 1.0000 | 0.9970 | 0.9987 | 0.9869 | 0.7989 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.6442 | 0.6442 |
| CoxPH_N50 | rna_ENSG00000090402.8 | phase4.1_low_WI | 1.0000 | 0.9937 | 0.9985 | 0.8790 | 0.7641 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.9180 | 0.9180 |
| CoxPH_N50 | rna_ENSG00000198670.12 | seeded_random | 1.0000 | 0.9973 | 0.9989 | 0.9173 | 0.7463 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.9996 | 0.9996 |
| CoxPH_N50 | rna_ENSG00000121871.4 | seeded_random | 1.0000 | 0.9921 | 0.9980 | 0.9829 | 0.7811 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.8540 | 0.8540 |
| RSF_N20 | rna_ENSG00000163545.11 | high_margin | 1.0000 | 0.9974 | 0.9949 | 0.9886 | 0.8703 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.0002 |
| RSF_N20 | rna_ENSG00000178522.14 | near_threshold | 1.0000 | 0.9918 | 0.9886 | 0.9805 | 0.8002 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.6278 | 0.6278 |
| RSF_N20 | rna_ENSG00000090402.8 | phase4.1_low_WI | 0.9977 | 0.5920 | 0.6834 | 0.0000 | 0.8260 | 0.0002 | 1.0000 | 0.9956 | 1.0000 | 0.2684 | 1.0000 |
| RSF_N20 | rna_ENSG00000231925.12 | seeded_random | 0.9997 | 0.9627 | 0.9599 | 0.8839 | 0.7576 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 0.9452 | 0.9452 |
| RSF_N20 | rna_ENSG00000187242.6 | seeded_random | 0.9998 | 0.8484 | 0.9307 | 0.7758 | 0.6899 | 0.0002 | 0.0660 | 0.0002 | 0.7672 | 1.0000 | 1.0000 |

Yüksek-marjin seçimlerde F ve p_F (RSF satırı sabit 12+12 query ve saklı attribution ortalamalarına koşullu):

| model | F | p_F |
| --- | --- | --- |
| CoxPH_N50 | 1.0000 | 0.0002 |
| RSF_N20 | 1.0000 | 0.0002 |

R eşiğine en yakın seçimlerde R ve p_R (RSF satırı sabit 12+12 query ve saklı attribution ortalamalarına koşullu):

| model | R | p_R |
| --- | --- | --- |
| CoxPH_N50 | 0.7989 | 0.6442 |
| RSF_N20 | 0.8002 | 0.6278 |

Sanity check beklenen yöndedir: iki yüksek-marjin feature için p_F Monte Carlo tabanı olan 0.0002'dir. Eşiğe en yakın R değerleri için p_R sırasıyla 0.6442 ve 0.6278'dir; sıfıra yakın değildir ve 0.5 çevresindeki belirsiz bölgededir. RSF düşük-WI kontrolünde G_claim=0 ve p_G_claim=1 çıkmıştır. B=4999 için p tahmininin en kötü durum binomial Monte Carlo SE'si 0.0071'tür.


Formül dosyaları Faz 4.2 hash'leriyle aynıdır: **True**. Çıktılarda hasta/örnek kimliği yoktur. Gerçek pilot wall-clock: **15.114 saniye**.

## Karar noktası

BH/BY uygulanmadı ve tam feature ailesine geçilmedi. Pilot, yalnız koşullu bootstrap davranışını gösterir. Tam ADIM 2 için yöntem onayı yanında null kalibrasyon simülasyonunun zorunlu tutulması önerilir; RSF sonucu ayrıca 12+12 sabit sorgu ve saklanmamış permütasyon-yolu belirsizliğine koşulludur.
