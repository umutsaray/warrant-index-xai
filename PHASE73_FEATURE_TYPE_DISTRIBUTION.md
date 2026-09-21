# Faz 7.3 — Özellik-tipi dağılımı

## Tam aile

| Subset | Total | RNA | clinical (age) | clinical (AJCC stage indicator) | Clinical total | RNA % | Clinical % |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 27-feature family | 27 | 22 | 1 | 4 | 5 | 81.4815 | 18.5185 |

## Model alt kümeleri

| Model | Subset | Total | RNA | clinical (age) | clinical (AJCC stage indicator) | Clinical total | RNA % | Clinical % | Clinical share minus family pp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH-v2 | point adequate | 9 | 6 | 1 | 2 | 3 | 66.6667 | 33.3333 | 14.8148 |
| CoxPH-v2 | BH certified | 8 | 6 | 1 | 1 | 2 | 75.0000 | 25.0000 | 6.4815 |
| RSF | point adequate | 17 | 13 | 0 | 4 | 4 | 76.4706 | 23.5294 | 5.0109 |
| RSF | BH certified | 12 | 8 | 0 | 4 | 4 | 66.6667 | 33.3333 | 14.8148 |
| DeepSurv | point adequate | 0 | 0 | 0 | 0 | 0 | NA | NA | NA |
| DeepSurv | BH certified | 0 | 0 | 0 | 0 | 0 | NA | NA | NA |

`Clinical share minus family pp`, alt kümedeki klinik pay ile tam ailenin %18.52 klinik payı arasındaki yüzde-puan farkıdır. Pozitif değer klinik özelliklerin sayısal olarak daha yüksek temsilini, negatif değer daha düşük temsilini gösterir; biyolojik yorum yapılmamıştır. Boş DeepSurv alt kümelerinde oran tanımsızdır.

Sayısal karşılaştırmada klinik pay CoxPH point/BH alt kümelerinde sırasıyla %33.33/%25.00 (+14.81/+6.48 yüzde puan), RSF point/BH alt kümelerinde %23.53/%33.33 (+5.01/+14.81 yüzde puan) oldu. Bu dört dolu alt küme tam aileye göre klinik özellikler yönünde kaymıştır; örneklem sayıları küçük olduğundan bu ifade yalnız betimseldir.
