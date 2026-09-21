# Faz 7.1 — Beş-fold hasta-gruplu CV özeti

## Tasarım

Beş `StratifiedGroupKFold` fold'unun her birinde preprocessing ve feature selection yalnız ilgili train partisyonunda öğrenildi. CoxPH, fit uyarısı varsa aynı `penalizer=0.1` korunarak `max_steps=200, step_size=0.25` ile retry edildi; fold 0 ve 2 retry gerektirdi ve nihai beş fit uyarısızdır. RSF ayarları değişmedi. DeepSurv her fold'da dondurulmuş mimari, dropout 0.3, AdamW (`lr=0.001`, `weight_decay=0.0001`), 300 maksimum epoch, patience 30 ve seed 42 ile eğitildi. Fold 0 DeepSurv yeniden eğitilmedi; resmi referans checkpoint'i yeniden kullanıldı.

## Fold skorları

| Fold | Train/test | CoxPH | RSF | DeepSurv | DeepSurv val | Seçilen epoch |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 674/168 | 0.667704 | 0.538521 | 0.496498 | 0.852814 | 38 |
| 1 | 674/168 | 0.705882 | 0.789609 | 0.719193 | 0.823666 | 53 |
| 2 | 674/168 | 0.681159 | 0.636232 | 0.650000 | 0.885370 | 37 |
| 3 | 673/169 | 0.655113 | 0.700471 | 0.531878 | 0.855314 | 5 |
| 4 | 673/169 | 0.650173 | 0.650748 | 0.662831 | 0.963492 | 38 |

| Model | Beş skor | Ortalama | Minimum | Maksimum | Aralık genişliği |
|---|---|---:|---:|---:|---:|
| CoxPH | .6677, .7059, .6812, .6551, .6502 | **0.672006** | 0.650173 | 0.705882 | 0.055710 |
| RSF | .5385, .7896, .6362, .7005, .6507 | **0.663116** | 0.538521 | 0.789609 | 0.251088 |
| DeepSurv | .4965, .7192, .6500, .5319, .6628 | **0.612080** | 0.496498 | 0.719193 | 0.222695 |

“Ortalama±aralık” raporlaması açık biçimde şöyle okunmalıdır: CoxPH **0.6720 [0.6502, 0.7059]**, RSF **0.6631 [0.5385, 0.7896]**, DeepSurv **0.6121 [0.4965, 0.7192]**. Standart sapma yerine talep edilen min–maks aralığı verilmiştir.

## Yorum

CoxPH retry sonrası foldlar arasında görece kararlıdır. RSF ve DeepSurv geniş fold aralıkları gösterir; fold 0 her ikisinde de minimumdur. Bu durum fold 0'ın audit referansı olma kararını değiştirmez, çünkü referans performansa bakılmadan önceden sabitlenmiştir.

DeepSurv'ün iç validation C-index'i her fold'da dış testten yüksektir. Val–test farkları fold sırasıyla 0.3563, 0.1045, 0.2354, 0.3234 ve 0.3007'dir. Böylece aşırı-uyum/genelleme boşluğu beş fold'un tamamında görülür; aynı zamanda dış test skorunun 0.4965–0.7192 arasında değişmesi önemli fold oynaklığını gösterir. Mimari veya regularizasyon post-hoc değiştirilmedi.

Tam tablo [cv_fold_scores.csv](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase71_cv/cv_fold_scores.csv), DeepSurv ayrıntıları [deepsurv_fold_details.csv](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase71_cv/deepsurv_fold_details.csv), ayarlar ve hash kontrolleri [cv_summary.json](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase71_cv/cv_summary.json) içindedir.

## Süre ve kapsam

Faz 7.1'in bu prompt içindeki ek çalışma süresi 34.33 saniyedir: referans CoxPH-v2 üretimi/doğrulaması 16.17 saniye, CV tamamlama 18.16 saniye. Yeniden kullanılan Faz 7.0b preprocessing/selection/CoxPH/RSF fold taraması 1059.17 saniye, yeniden kullanılan fold-2 Cox retry 17.53 saniyedir. Baştan sona bu CV kanıt zincirinin ölçülen hesap maliyeti yaklaşık **1111.03 saniye (18.52 dakika)** olmuştur.

F/S/U/G/G_signed/R, N=50 bootstrap, null kalibrasyon, BH/BY ve METABRIC çalıştırılmadı. Fold-0 resmi referans artefaktlarının hashleri CV öncesi/sonrası değişmedi.
