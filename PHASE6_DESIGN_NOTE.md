# Faz 6 — DeepSurv keşif, açıklanabilirlik ve N=2 tasarım pilotu

> **Koşulluluk uyarısı:** DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result.

## ADIM 0 — hazır model keşfi

Hazır `explainable_dl_model.pt` bulundu ve artifact içinden doğrulandı. Model aynı 673 train / 169 test ve 616-feature panelini kullanır. Dört bloklu mimari: RNA 200→128→32, miRNA 200→64→16, methylation 200→128→32, clinical 16→32→8; birleşik head 88→64→16→1. Encoder katmanları Linear–BatchNorm–ReLU–Dropout, dropout 0.3'tür. AdamW lr=0.001, weight decay=0.0001, en çok 300 epoch, patience=30 ve seed=42 kullanılmıştır. En iyi epoch 26, validation C-index 0.872659; test C-index artifacttan yeniden üretilerek **0.745692** bulundu. Ana model yeniden eğitilmedi.

Karşılaştırılabilirlik sınırlaması: ana DeepSurv modeli 673 satırlık train içinde %15 satır-stratified validation kullanır. Mixed cohortta aynı hastanın birden çok örneği bu iç split'i aşabilir; bu, mevcut model artifactının kayıtlı metodolojik sınırlamasıdır.

## ADIM 1 — gerçek explainer testi

| explainer | supported | shape | finite | max_abs_additivity_error | elapsed_seconds |
| --- | --- | --- | --- | --- | --- |
| DeepExplainer | True | [3, 616] | True | 0.0000 | 0.0168 |
| GradientExplainer | True | [3, 616] | True | 0.0248 | 0.0398 |

Her iki yöntem gerçek checkpoint ve gerçek 616-feature tensörlerinde çalıştı. DeepExplainer maksimum additivity hatası 1.77e-07 ile near-exact destek kontrolüdür. Pilot için GradientExplainer seçildi; çünkü `nsamples` parametresi bağımsız candidate/reference 8/128 bütçesini uygular. DeepExplainer deterministik olduğundan bu F bütçe karşılaştırmasını dejenere eder. GradientExplainer sonucu yaklaşık ve stokastiktir; N=2 candidate/ref maksimum additivity hataları sırasıyla 0.4392/0.0777'tür.

## ADIM 2 — sampled-background sağlık kontrolü

Sabit-ortalama background yeniden denenmedi. Diğer 539 feature, RSF ile aynı deterministik 16 benzersiz TCGA train hastası üzerinden prediction düzeyinde entegre edildi; aktif 77-feature baseline'ı TCGA train ortalamasıdır. Sekiz çeşitli feature için nested sampled k=8→16 ve reference 128→256 kontrolü:

| inference_scope | feature | block | WI_k16_ref128 | F_k16_ref128 | G_claim_k16 | F_k8_ref128 | G_claim_k8 | delta_F_k16_minus_k8 | delta_G_k16_minus_k8 | F_k16_ref256 | delta_F_ref256_minus_128 | candidate_attr_mean_abs_delta_k16_vs_k8 | reference_attr_mean_abs_delta_k16_vs_k8 | reference_attr_mean_abs_delta_256_vs_128 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000124343.14 | rna | 0.4145 | 0.9988 | 0.9318 | 0.9983 | 0.9280 | 0.0005 | 0.0038 | 0.9987 | -0.0001 | 0.0033 | 0.0030 | 0.0006 |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000100911.16 | rna | 0.4211 | 0.9998 | 0.9615 | 0.9998 | 0.9593 | 0.0001 | 0.0021 | 0.9999 | 0.0000 | 0.0039 | 0.0037 | 0.0005 |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000233452.7 | rna | 0.4342 | 0.9998 | 0.9555 | 0.9998 | 0.9570 | -0.0000 | -0.0016 | 0.9998 | 0.0000 | 0.0040 | 0.0038 | 0.0007 |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000250334.6 | rna | 0.8224 | 1.0000 | 0.9867 | 1.0000 | 0.9885 | 0.0000 | -0.0018 | 1.0000 | -0.0000 | 0.0032 | 0.0033 | 0.0004 |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000122025.15 | rna | 0.8289 | 0.9998 | 0.8481 | 0.9997 | 0.8481 | 0.0001 | 0.0000 | 0.9998 | 0.0000 | 0.0047 | 0.0047 | 0.0007 |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000134531.10 | rna | 0.8355 | 1.0000 | 0.9857 | 1.0000 | 0.9849 | 0.0000 | 0.0008 | 1.0000 | -0.0000 | 0.0034 | 0.0034 | 0.0004 |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000104327.7 | rna | 0.7632 | 0.9998 | 0.9653 | 0.9999 | 0.9660 | -0.0000 | -0.0007 | 0.9998 | -0.0000 | 0.0017 | 0.0016 | 0.0004 |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000116663.11 | rna | 0.5724 | 0.9999 | 0.9865 | 0.9999 | 0.9856 | 0.0000 | 0.0009 | 0.9999 | 0.0000 | 0.0019 | 0.0017 | 0.0004 |

Maksimum |ΔF| k8→k16 **0.000512**, maksimum |ΔG_claim| **0.003799**. Reference 128→256 maksimum |ΔF| **0.000098**. Bunlar pilot sağlık ölçüleridir; eşik kalibrasyonu değildir.

## ADIM 3 — N=2 pilot

> **Koşulluluk uyarısı:** DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result.

| inference_scope | block | features | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | clin | 5 | 5 | 0.9999 | 0.9922 | 0.9903 | 0.9729 | 0.6316 | 0.6316 |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna | 72 | 72 | 0.9994 | 0.9738 | 0.9825 | 0.9418 | 0.6832 | 0.6807 |

77/77 feature değerlendirildi; 15 feature WI≥0.8 geçti. Bottleneck dağılımı `{'R': 75, 'G_claim': 2}`.

### Örnek feature'lar

| inference_scope | feature | block | F | S | U | G_claim | R | WI | bottleneck |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000124343.14 | rna | 0.9988 | 0.9071 | 0.9658 | 0.9318 | 0.4145 | 0.4145 | R |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000100911.16 | rna | 0.9998 | 0.9400 | 0.9846 | 0.9615 | 0.4211 | 0.4211 | R |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000233452.7 | rna | 0.9998 | 0.8917 | 0.9814 | 0.9555 | 0.4342 | 0.4342 | R |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000072657.9 | rna | 1.0000 | 0.9766 | 0.9897 | 0.9546 | 0.4408 | 0.4408 | R |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000134531.10 | rna | 1.0000 | 0.9868 | 0.9916 | 0.9857 | 0.8355 | 0.8355 | R |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000122025.15 | rna | 0.9998 | 0.9147 | 0.9644 | 0.8481 | 0.8289 | 0.8289 | R |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000090402.8 | rna | 0.9999 | 0.9885 | 0.9719 | 0.8586 | 0.8224 | 0.8224 | R |
| DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 attribution draws; N=2 is a design pilot, not a certification result. | rna_ENSG00000108639.8 | rna | 1.0000 | 0.9956 | 0.9938 | 0.9795 | 0.8224 | 0.8224 | R |

Onaylı ana N=2 eğitim+candidate+reference süresi 2.635 sn. Doğrusal projeksiyon N=20 için 26.3 sn (0.44 dk), N=50 için 65.9 sn (1.10 dk). Bütün sensitivity kontrolleri ve raporlama dahil gerçek invocation 6.478 sn sürdü.

## Karar ve sınırlar

Sampled k=16 background sağlık kontrolü ve reference=128 bütçesi onaylanabilir: k8→k16 maksimum |ΔG_claim| 0.0038, reference 128→256 maksimum |ΔF| 9.8e-5'tir. Buna karşılık GradientExplainer candidate=8 maksimum additivity hatası 0.439 risk-skoru birimidir; bu nedenle tasarım kapısı **REVİZYON GEREKTİRİYOR**. N≥20 mevcut 8/128 ayarıyla başlatılmamalıdır. Sonraki ucuz karar turunda candidate bütçesi 32/64'e yükseltilerek additivity ve F yakınsaklığı ölçülebilir veya near-exact DeepExplainer'ın F katmanındaki rolü yeniden tanımlanabilir. Bu iki seçenek frozen WI bileşen formüllerini değiştirmez, fakat explainer protokolü için kullanıcı kararı gerektirir.

N≥20 koşusu, null kalibrasyon ve FDR başlatılmadı. Query kümesi sabit 12+12'dir; GradientExplainer örnekleme belirsizliği daha sonraki inferansta ayrıca ele alınmalıdır.

Formül SHA256 değişmedi: **True**. Kimlik-sızıntısı kontrolü geçti. Runtime: `C:\moenv\Scripts\python.exe`.
