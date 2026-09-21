# Faz 5.4 — BH hayatta-kalma farkı tanısı

Bu analiz yalnız mevcut Faz 5.3 sertifikasyon tablolarını okur. Yeni bootstrap, model fit'i, explainer çağrısı veya formül hesabı yapılmadı.

## Ana sonuç

Ters-sezgisel oran farkının ana açıklaması **güç üstünlüğü değil, marjin dağılımıdır**. RSF'nin nokta-geçen 23 feature'ından 16'sı 0.03'ten büyük marjine sahipken CoxPH'de bu sayı 9'dur. RSF'de 0.03'ün üzerindeki 16 feature'ın 15'i BH'den geçer; CoxPH'de 9 feature'ın 8'i geçer. Buna karşılık 0–0.03 bandında iki modelde de hiçbir feature sertifikalanmaz.

Aynı düşük-marjin bandında CoxPH testi daha zayıf görünmez: 0–0.03 için ortalama darboğaz p-değeri CoxPH'de 0.2644, RSF'de 0.5237'tür. Bu gözlem N=50'nin aynı marjin civarında daha fazla kesinlik sağlamasıyla uyumludur; model, estimand ve replika yapıları farklı olduğu için nedensel bir N etkisi olarak yorumlanamaz.

## Kontrol 1 — marjinler

| model | BH_status | n | margin_mean | margin_median | margin_q25 | margin_q75 | margin_min | margin_max | p_bottleneck_mean | p_joint_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH N=50 | point only | 14.0000 | 0.0148 | 0.0138 | 0.0066 | 0.0217 | 0.0036 | 0.0308 | 0.2582 | 0.2582 |
| CoxPH N=50 | BH certified | 8.0000 | 0.0545 | 0.0371 | 0.0359 | 0.0476 | 0.0331 | 0.1503 | 0.0002 | 0.0002 |
| RSF N=20 | point only | 8.0000 | 0.0148 | 0.0066 | 0.0023 | 0.0195 | 0.0002 | 0.0574 | 0.4583 | 0.4739 |
| RSF N=20 | BH certified | 15.0000 | 0.0595 | 0.0593 | 0.0565 | 0.0640 | 0.0431 | 0.0703 | 0.0008 | 0.0013 |

CoxPH'nin 14 kaçıranının ortalama marjini 0.0148, medyanı 0.0138; RSF'nin sekiz kaçıranında sırasıyla 0.0148 ve 0.0066'tür. Ortalamalar neredeyse aynıdır; RSF medyanı daha küçüktür ve 0.0574 marjinli tek aykırı kaçıran ortalamayı yükseltir. Dolayısıyla CoxPH kaçıranları sistematik olarak eşiğe daha yakın değildir. WI bir bileşen minimumu olduğundan, tekil nokta-darboğaz durumlarında WI marjini ile darboğaz bileşen marjini aynıdır; CSV iki sütunu da açıkça içerir.

## Kontrol 2 — darboğazlar

| model | bottleneck | n |
| --- | --- | --- |
| CoxPH N=50 | R | 14 |
| RSF N=20 | R | 8 |

CoxPH'nin kaçırdığı 14 feature'ın tamamının nokta-tahmin ve inferans darboğazı R'dir. RSF'nin sekiz kaçıranının da nokta-tahmin darboğazı R'dir; bunların yedisinde en büyük bileşen p-değeri yine R, `rna_ENSG00000182591.6` için ise U'dur. Bu feature'da p_R=0.0002 olmasına rağmen p_U=p_joint=0.1256 olduğu için recursive IUT sertifikayı engeller.

## Kontrol 3 — aynı marjin bandında p-değerleri

| model | margin_band | n | margin_mean | p_bottleneck_mean | p_bottleneck_median | p_joint_mean | p_joint_median | BH_certified |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH N=50 | 0.01-0.03 | 8.0000 | 0.0186 | 0.1400 | 0.1616 | 0.1400 | 0.1616 | 0.0000 |
| CoxPH N=50 | 0.03-0.05 | 7.0000 | 0.0356 | 0.0255 | 0.0002 | 0.0255 | 0.0002 | 6.0000 |
| CoxPH N=50 | <0.01 | 5.0000 | 0.0057 | 0.4634 | 0.4712 | 0.4634 | 0.4712 | 0.0000 |
| CoxPH N=50 | >=0.10 | 1.0000 | 0.1503 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 1.0000 |
| CoxPH N=50 | 0.05-0.10 | 1.0000 | 0.0670 | 0.0002 | 0.0002 | 0.0002 | 0.0002 | 1.0000 |
| RSF N=20 | <0.01 | 5.0000 | 0.0036 | 0.5766 | 0.5990 | 0.5766 | 0.5990 | 0.0000 |
| RSF N=20 | 0.05-0.10 | 15.0000 | 0.0605 | 0.0002 | 0.0002 | 0.0091 | 0.0002 | 14.0000 |
| RSF N=20 | 0.01-0.03 | 2.0000 | 0.0216 | 0.3916 | 0.3916 | 0.3916 | 0.3916 | 0.0000 |
| RSF N=20 | 0.03-0.05 | 1.0000 | 0.0431 | 0.0088 | 0.0088 | 0.0088 | 0.0088 | 1.0000 |

Önceden sorulan 0–0.03 bandının birleştirilmiş özeti:

| model | n | mean_margin | mean_p_bottleneck | median_p_bottleneck | mean_p_joint |
| --- | --- | --- | --- | --- | --- |
| CoxPH N=50 | 13.0000 | 0.0136 | 0.2644 | 0.2192 | 0.2644 |
| RSF N=20 | 7.0000 | 0.0087 | 0.5237 | 0.5574 | 0.5237 |

Marjin ile darboğaz p-değeri güçlü ters ilişkilidir: Spearman rho CoxPH için -0.954, RSF için -0.851. CoxPH'nin 0.01–0.03 bandındaki ortalama p'si RSF'den daha düşüktür. Bu nedenle mevcut veride “aynı marjin için CoxPH daha az güçlü” hipotezini destekleyen işaret yoktur. 0.03 üzerindeki bantlarda örnek sayıları dengesiz olduğundan model-etkisini ayrı bir katsayı gibi tahmin etmek güvenilir değildir.

![Marjin ve p-değeri ilişkisi](results/phase54_bh_survival_diagnostic/p_vs_margin.png)

## Kontrol 4 — R'nin T ve p ile ilişkisi

R'deki T(T−1) paydası T×(T−1)/2 replika çifti ortalamasını normalleştirir. Çiftler aynı replikaları paylaştığı için bağımsız değildir; bu nedenle etkili bilgi T² hızında büyümez. T=50 yine T=20'ye göre daha fazla replika bilgisi sağlar, fakat salt paydaya bakarak standart hatanın oranı çıkarılamaz.

`p−1` terimi de karşılaştırılan rakip feature'lar üzerindeki ortalamayı normalleştirir. CoxPH R sıralaması 614 model feature'ı, RSF R sıralaması 77 aktif feature'ı kullanır. Daha çok rakip, bağımsız olsalardı ortalamayı yumuşatabilirdi; burada bütün karşılaştırmalar aynı hedef feature'ı ve replika sıralamalarını paylaştığı için güçlü bağımlılık vardır. Ayrıca p değişince R'nin estimandı da değişir: CoxPH ve RSF aynı sıralama evrenini ölçmez. Yakın bağlar ve eta toleransı hangi rakiplerin oynaklığa katkı verdiğini belirler. Dolayısıyla p−1 farkı tek başına gözlenen sertifika oranını açıklamaz.

## Değerlendirme

En olası açıklama, nokta-geçen ailelerin farklı marjin bileşimidir. RSF'de sertifikalanan 15 feature, 0.0431–0.0703 aralığındaki yüksek-marjin kümesinde yer alır; CoxPH'de 14 feature düşük-marjin bölgede kalır. R her iki modelde de düşük-marjin başarısızlıklarının merkezindedir, fakat aynı 0–0.03 bandında RSF p-değerleri daha büyüktür. Bu nedenle sonuç, RSF N=20 prosedürünün CoxPH N=50'den daha güçlü olduğu şeklinde yorumlanmamalıdır. RSF'nin dar 12+12 query kümesine ve saklı 8/128 attribution ortalamalarına koşulluluk uyarısı aynen geçerlidir.

Analiz 0.358 saniye sürdü. Formül SHA256 değerleri sabit kaldı ve kaynak sertifikasyon CSV'lerinin SHA256 değerleri metadata'ya yazıldı. Kimlik-sızıntısı kontrolü geçti.
