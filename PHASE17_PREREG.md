# Faz 17 v2 ön-kayıt

**Kayıt zamanı:** v2 sonuçları hesaplanmadan önce.  
**Etiket:** v2 post hoc/keşifsel; v1 sonuçları görüldükten sonra tasarlandı.

Bu belge birincil tanımları, duyarlılık ızgarasını ve karar kuralını sonuçlar görülmeden sabitler. Sonradan değiştirilmeyecek; gerekirse yalnız açıkça `EK` olarak işaretlenen eklemeler yapılacaktır.

## Ölçek

`s_m`, fold-0 TCGA kohortunda, sıfır olmayan `M_j` özniteliklerinde `M_j`'nin medyanıdır (Faz 12 tanımı). Saklı Faz 12 girdilerinden yeniden doğrulanan değerler:

| Model | s_m |
|---|---:|
| CoxPH-v2 | 0.023934215222165295 |
| RSF | 0.01707069437342852 |
| DeepSurv | 0.012307440008774696 |

Bu değerler hedef 0.0239/0.0171/0.0123 ile ±1e-4 içinde eşleşir.

## v2 birincil tolerans seti

- `a_m = b_m = k · s_m`
- `u_m = 6k · s_m`; v1'deki `u/a=6` oranı korunur.
- `η_m = ρ · s_m`
- **Birincil:** `k=1`, `ρ=1`.

Gerekçe: tolerans, “tipik bir atıf birimi” büyüklüğünde bir sapmanın tam kayıp sayılması anlamına gelir. Ölçek modele özgü tek bir sayıdır; öznitelik-başına normalizasyon değildir, çünkü sıfır atıflı özniteliklerde böyle bir ölçek tanımsızdır.

## Altıncı bileşen N — exercised signal

Her kohort `d`, replika `t` ve özellik `j` için:

- `Î_dtj = mean_i |φ̄_dtj(x_i)|`, sorgular üzerindeki ortalama mutlak atıf.
- `N_dtj = min(1, Î_dtj / m0)`.
- `N_dj = (1/T) Σ_t N_dtj`.
- `N_j = min_d N_dj`; kaynak dahil bütün kohortlarda sinyal bulunmalıdır.
- `m0 = κ · s_m`.
- **Birincil:** `κ=0.25`.

N, F/U ile aynı yapıda sınırlı bir replika ortalamasıdır; çıkarım U ile aynı `recursive_studentized` pivotunu kullanır. Null hipotezi `H0: N_j ≤ τ`; `τ=0.8` korunur.

## Altı bileşenli indeks ve sertifikasyon

- `WI6 = min(F, S, U, G^claim, R, N)`.
- Sertifikasyon altı bileşenli intersection-union testidir.
- Aile model başına 27 özelliktir.
- BH ve BY düzeltmeleri v1 ile aynıdır.

## Duyarlılık ızgarası

Tüm hücreler raporlanacak, hiçbiri sonuçlara göre seçilmeyecektir:

- `k ∈ {0.25, 0.5, 1, 2, 4}`
- `ρ ∈ {0.5, 1, 2}`
- `κ ∈ {0.1, 0.25, 0.5}`

Toplam `5×3×3=45` hücre/model vardır.

## Karar kuralı

v2 sonuçları, birincil ayar ve bütün duyarlılık ızgarası dahil, sonuca bakılmaksızın tam raporlanacaktır. Birincil ayar sonuçlara göre değiştirilmeyecektir. v1 resmî analiz olarak korunur; v2 ayrı, post hoc ve keşifsel bir analizdir.
