# Faz 3 — TCGA-BRCA / METABRIC CoxPH yöntem doğrulaması

Gerçek veri üzerinde teknik yöntem doğrulaması; klinik/biyolojik keşif iddiası yok. 50 hasta-kümesi bootstrap replikası; sabit penalizer=0.1.

## Manifest ve veri kaynakları

Train 673/106 event, test 169/26 event eşleşiyor. **SAPMA: 623 değil 616 feature** (200 RNA, 200 miRNA, 200 methylation, 16 clinical). Kayıtlı model iki separation sütunu çıkarılarak 614 feature ile eğitilmiş; model/scaler/features içeren joblib yüklendi. Ana model yeniden eğitilmedi; bootstrap modelleri aynı penalizer/max_steps ile eğitildi. Detay ve tam liste kökte MANIFEST.md.

**Hasta sızıntısı:** orijinal splitte 23 ortak hasta var. Audit'te bunlar çıkarıldı, kalan hastalarda dosya sırasındaki ilk örnek seçildi: 141 hasta/16 event. Train 628 benzersiz hasta; cluster bootstrap aynı hastanın örneklerini birlikte örnekliyor. Önceden seçilmiş panel ve scaler sabit: S tüm seçim/imputation hattının belirsizliğini ölçmüyor.

[cBioPortal REST API](https://www.cbioportal.org/api/swagger-ui/index.html): brca_metabric, brca_metabric_mrna (Illumina HT-12 v3); clinical + seçilmiş RNA paneli için hedefli sorguların ham JSON yanıtları data/external/metabric altında. Tam transkriptom indirilmedi. 1980 expression örneği, 72 eşleşen gen. OS_MONTHS/OS_STATUS mevcut: 1979 kullanılabilir etiket, 1143 event. Etiketler yalnız envanter içindir; METABRIC ile model eğitimi/eşik seçimi yapılmadı. CNA/mutation/methylation profilleri indirilmedi; profil listesi yalnız metadata.

Ensembl sürüm eki kaldırıldı; MyGene.info gene-symbol anotasyonu ile cBioPortal HUGO symbol birebir eşlendi; alias veya belirsiz çoklu eşleşme kullanılmadı. AGE_AT_DIAGNOSIS yaşa, TUMOR_STAGE 1–4 TCGA Stage I–IV sütunlarına eşlendi. assumption: bu genel stage sınıfları karşılaştırılabilir; AJCC sürümü/ölçüm ayrıntılarının eşdeğerliği doğrulanmadı. Grade stage yerine kullanılmadı. ER/PR/HER2 ve grade TCGA seçilmiş panelinde yok. Cinsiyet bu görevin izin verdiği klinik eşleme listesi dışında tutuldu; dış kohortta ayrıca sabit. Irk/etnisite eşlenmedi.

## Harmonizasyon, oracle ve G'nin sınırı

assumption: ortak 72 RNA + yaş + dört 0/1 evre feature'ının her biri her kohortta ayrı mean/std (ddof=0) ile z-score. TCGA train'de kayıtlı scaler öğrenilmiş; temiz test aynı scaler ile dönüştürüldü. METABRIC'te 77 feature'ın tamamı gözlenen tam hasta vektörleri üzerinde kendi ortalama/standart sapmasıyla standardize edildi. RNA'nın hazırlanmış ölçeği korundu; ham RNA-seq/microarray değerleri karşılaştırılmadı. z-score platform/batch etkilerini tamamen giderdiğini kanıtlamaz.

**assumption: background=0**, TCGA tarafında TCGA train ortalamasını, METABRIC tarafında METABRIC tam-query kohort ortalamasını temsil eder. Bunlar iki ayrı z-score uzayındaki iki ayrı referanstır; ortak, kohort-bağımsız biyolojik referans değildir.

**Kohort kapsamı:** kök veri primary-tumor-only değildir. Train: 609 primary solid tumor, 60 solid tissue normal, 4 metastatic; test: 153/15/1. [GDC örnek türü kodları](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) 01/11/06 ile doğrulandı. İstenen model/split korunmuştur; başka primary-only dizine sessiz geçiş yapılmadı. Bu karışım biyolojik genelleme için ek sınırlamadır.

Oracle, [lifelines CoxPH](https://lifelines.readthedocs.io/en/latest/fitters/regression/CoxPHFitter.html) log-partial-hazard için exact marginal Shapley: beta_j*(z_j-background_j), background=0. Toplam attribution + base ile predict_log_partial_hazard farkı en fazla 2.84e-14. Exact ifadesi survival probability veya conditional/causal SHAP için kullanılmaz. F=1 oracle'ın kendisiyle kontrolüdür; harici explainer doğrulaması değildir.

Faz 3'te kullanılan aynı-girdi sonucu silinmedi: `G_claim_range_overlap_diagnostic` ve `WI_phase3_range_overlap` kolonlarında korunur. Bu tanıda METABRIC aralığına düşen TCGA x_i iki siteye birden verildiği için G≈1 matematiksel olarak trivialdir; WI'nin varsayılan girdisi değildir.

Faz 3.1 G/G_signed, site0'da 141 temiz TCGA audit hastasının; site1'de ortak 77 feature'ı eksiksiz olan gerçek METABRIC hasta vektörlerinin attribution dağılımlarını kullanır. Faz 3.2'de aynı query-set ve beta_j*(z_j-0) oracle mantığı 50 bootstrap katsayı dizisine uygulanır. G, iki gerçek sorgu kümesindeki mean_t mean_i |phi| değerlerinden hesaplanır. G_signed için hasta eşleştirmesi yapılmaz: her kohortun gözlenen signed attribution order-statistic'leri max(n_TCGA,n_METABRIC) büyüklüğündeki ortak ampirik rank gridine tekrarlanır; interpolasyon veya sentetik hasta yoktur. Böylece dondurulmuş Faz 2.1 anchor ortalaması ve clipping formülü değişmeden iki farklı boyuttaki dağılıma uygulanır. Bu deterministik kuantil coupling bir metodoloji varsayımıdır; başka dağılımsal coupling sonuçları değiştirebilir.

F/S/U/R TCGA temiz audit kümesinden gelir. U sabit modele koşullu normal mean-|phi| aralığıdır ve coverage garantisi yoktur. Eksik omik modalitelerle tam METABRIC hazard/prediction üretilmedi. Bu koşu dış survival performansını veya mekanizma invariance'ını sınamaz.

Toleranslar Faz 2'den aynen: a=b=0.5,u=3,eta=0.02,tau=0.8. Gerçek veri için kalibre edilmiş klinik eşikler değildir. Point adequacy istatistiksel certification değildir; bootstrap replikalarından IUT p-değeri üretilmedi.

Bootstrap eğitiminde sabit/separation sütunları her replikada train verisinden dışlandı. İlk koşunun altıncı replikasında meth_cg04831594 için düşük koşullu varyans ve yakınsamama uyarısı çıktı. Aynı örnekleme/penalizer=0.1 korunarak lifelines'ın işaretlediği separation sütunu çıkarıldı. 50 adımda kalan sayısal yakınsamama için max_steps=200 ve step_size=0.25 ile aynı optimizasyon hedefi tekrar çözüldü; başarısız replikalar başarılı olanlarla değiştirilmedi. N=50 koşusunda 13 replika birden çok fit denemesi, 11 replika sayısal retry gerektirdi. Separation uyarıları: `{'meth_cg04831594': 9}`. Nihai hiçbir replikada warning bırakılmadı. Fit denemeleri ve dışlanan feature'lar bootstrap_summary.json içinde; dışlanan sütun katsayısı o replikada sıfırdır.

## Faz 3.1 düzeltmesi — önce/sonra

METABRIC 77-feature query matrisi: 1466 tam gerçek hasta; 514 eksik vektör çıkarıldı. Çıktılarda kimlik yoktur.

### Faz 3/Faz 3.1 block özeti

| method | block | features | assessed | mean_WI | F | S | U | G_claim | R |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| phase3_range_overlap_diagnostic | clin | 16 | 5 | 0.8563 | 1.0000 | 0.9923 | 0.9957 | 1.0000 | 0.8094 |
| phase3_range_overlap_diagnostic | meth | 200 | 0 | NA | 1.0000 | 0.9944 | 0.9981 | NA | 0.7387 |
| phase3_range_overlap_diagnostic | mirna | 200 | 0 | NA | 1.0000 | 0.9934 | 0.9968 | NA | 0.7590 |
| phase3_range_overlap_diagnostic | rna | 200 | 72 | 0.7593 | 1.0000 | 0.9952 | 0.9972 | 1.0000 | 0.7780 |
| phase3.1_actual_METABRIC_queries_n6 | clin | 16 | 5 | 0.8529 | 1.0000 | 0.9923 | 0.9957 | 0.9368 | 0.8094 |
| phase3.1_actual_METABRIC_queries_n6 | meth | 200 | 0 | NA | 1.0000 | 0.9944 | 0.9981 | NA | 0.7387 |
| phase3.1_actual_METABRIC_queries_n6 | mirna | 200 | 0 | NA | 1.0000 | 0.9934 | 0.9968 | NA | 0.7590 |
| phase3.1_actual_METABRIC_queries_n6 | rna | 200 | 72 | 0.7593 | 1.0000 | 0.9952 | 0.9972 | 0.9707 | 0.7780 |

Yeni G_claim trivial değildir: N=6 değerlendirilen feature'larda minimum=0.833009, ortalama=0.968535, maksimum=0.997988. Faz 3 range-overlap tanısına göre WI 2 feature'da değişti.

## Faz 3.2 — replika sayısı artırma

Yalnız patient-cluster bootstrap replikası 6'dan 50'ye çıkarıldı. Seed=20260917 ve deterministik RNG sırası korundu; ilk altı checkpoint yeniden kullanıldı, 7–50 aynı akışın devamıdır. CoxPH ana model, query setleri, toleranslar ve tüm formüller değişmedi.

### N=6 ve N=50 R dağılımı — önceki 27 geçer feature

| feature | R_n6 | R_n6_resample_std | R_n6_q025 | R_n6_q25 | R_n6_q75 | R_n6_q975 | R_n50 | R_n50_resample_std | R_n50_q025 | R_n50_q25 | R_n50_q75 | R_n50_q975 | IQR_width_n6 | IQR_width_n50 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000133742.14 | 0.8347 | 0.0220 | 0.8334 | 0.8439 | 0.8763 | 0.9113 | 0.8412 | 0.0044 | 0.8346 | 0.8416 | 0.8473 | 0.8519 | 0.0324 | 0.0058 |
| rna_ENSG00000169800.14 | 0.8358 | 0.0176 | 0.8404 | 0.8497 | 0.8738 | 0.9051 | 0.8331 | 0.0016 | 0.8337 | 0.8354 | 0.8374 | 0.8399 | 0.0241 | 0.0020 |
| rna_ENSG00000090402.8 | 0.8505 | 0.0165 | 0.8556 | 0.8649 | 0.8844 | 0.9143 | 0.7641 | 0.0293 | 0.7112 | 0.7486 | 0.7905 | 0.8267 | 0.0195 | 0.0419 |
| rna_ENSG00000139874.6 | 0.8505 | 0.0186 | 0.8505 | 0.8606 | 0.8840 | 0.9203 | 0.7989 | 0.0143 | 0.7753 | 0.7935 | 0.8128 | 0.8321 | 0.0233 | 0.0193 |
| rna_ENSG00000198099.9 | 0.8605 | 0.0157 | 0.8605 | 0.8707 | 0.8928 | 0.9174 | 0.7421 | 0.0178 | 0.7152 | 0.7348 | 0.7596 | 0.7841 | 0.0220 | 0.0247 |
| rna_ENSG00000270069.1 | 0.8266 | 0.0216 | 0.8246 | 0.8379 | 0.8695 | 0.9023 | 0.7762 | 0.0169 | 0.7472 | 0.7702 | 0.7921 | 0.8120 | 0.0316 | 0.0220 |
| rna_ENSG00000174899.12 | 0.8014 | 0.0209 | 0.8063 | 0.8199 | 0.8454 | 0.8838 | 0.7832 | 0.0176 | 0.7552 | 0.7769 | 0.8002 | 0.8241 | 0.0255 | 0.0233 |
| rna_ENSG00000168917.9 | 0.8134 | 0.0274 | 0.8073 | 0.8230 | 0.8656 | 0.9016 | 0.7522 | 0.0177 | 0.7229 | 0.7452 | 0.7693 | 0.7903 | 0.0426 | 0.0241 |
| rna_ENSG00000270547.6 | 0.8007 | 0.0335 | 0.7873 | 0.8071 | 0.8599 | 0.9039 | 0.8200 | 0.0175 | 0.7865 | 0.8118 | 0.8355 | 0.8564 | 0.0529 | 0.0236 |
| rna_ENSG00000133800.9 | 0.8334 | 0.0224 | 0.8288 | 0.8420 | 0.8754 | 0.9085 | 0.8367 | 0.0087 | 0.8222 | 0.8354 | 0.8476 | 0.8548 | 0.0334 | 0.0122 |
| rna_ENSG00000077616.11 | 0.8467 | 0.0162 | 0.8517 | 0.8610 | 0.8803 | 0.9117 | 0.8118 | 0.0129 | 0.7902 | 0.8069 | 0.8244 | 0.8382 | 0.0193 | 0.0175 |
| rna_ENSG00000226674.11 | 0.8579 | 0.0155 | 0.8605 | 0.8695 | 0.8883 | 0.9211 | 0.7657 | 0.0220 | 0.7257 | 0.7532 | 0.7832 | 0.8091 | 0.0188 | 0.0300 |
| rna_ENSG00000109158.11 | 0.8350 | 0.0183 | 0.8367 | 0.8485 | 0.8734 | 0.9051 | 0.8335 | 0.0044 | 0.8279 | 0.8340 | 0.8401 | 0.8444 | 0.0249 | 0.0062 |
| rna_ENSG00000023892.11 | 0.8171 | 0.0216 | 0.8171 | 0.8284 | 0.8613 | 0.8939 | 0.8036 | 0.0109 | 0.7822 | 0.8007 | 0.8150 | 0.8267 | 0.0329 | 0.0143 |
| rna_ENSG00000118407.15 | 0.8001 | 0.0441 | 0.7698 | 0.8035 | 0.8794 | 0.9139 | 0.8050 | 0.0124 | 0.7810 | 0.7999 | 0.8169 | 0.8307 | 0.0760 | 0.0169 |
| rna_ENSG00000108639.8 | 0.8290 | 0.0230 | 0.8264 | 0.8382 | 0.8717 | 0.9093 | 0.7781 | 0.0172 | 0.7518 | 0.7696 | 0.7955 | 0.8142 | 0.0335 | 0.0258 |
| rna_ENSG00000076043.11 | 0.8249 | 0.0325 | 0.8096 | 0.8278 | 0.8745 | 0.9313 | 0.7228 | 0.0192 | 0.6965 | 0.7136 | 0.7407 | 0.7691 | 0.0467 | 0.0271 |
| rna_ENSG00000182173.13 | 0.8566 | 0.0167 | 0.8569 | 0.8681 | 0.8897 | 0.9174 | 0.8222 | 0.0128 | 0.8024 | 0.8165 | 0.8356 | 0.8497 | 0.0216 | 0.0190 |
| rna_ENSG00000240108.1 | 0.8081 | 0.0287 | 0.8020 | 0.8177 | 0.8573 | 0.9041 | 0.7945 | 0.0170 | 0.7612 | 0.7878 | 0.8102 | 0.8317 | 0.0396 | 0.0224 |
| rna_ENSG00000187758.8 | 0.8524 | 0.0158 | 0.8576 | 0.8648 | 0.8844 | 0.9161 | 0.8371 | 0.0084 | 0.8237 | 0.8349 | 0.8471 | 0.8559 | 0.0195 | 0.0122 |
| rna_ENSG00000104327.7 | 0.8103 | 0.0264 | 0.8040 | 0.8212 | 0.8629 | 0.8992 | 0.7821 | 0.0100 | 0.7684 | 0.7796 | 0.7933 | 0.8078 | 0.0417 | 0.0137 |
| rna_ENSG00000204246.3 | 0.8324 | 0.0177 | 0.8355 | 0.8475 | 0.8700 | 0.9010 | 0.8371 | 0.0025 | 0.8354 | 0.8386 | 0.8420 | 0.8454 | 0.0225 | 0.0034 |
| rna_ENSG00000205186.4 | 0.8426 | 0.0168 | 0.8475 | 0.8560 | 0.8781 | 0.9095 | 0.8228 | 0.0112 | 0.8024 | 0.8199 | 0.8347 | 0.8466 | 0.0221 | 0.0148 |
| clin_numeric__age_at_index.demographic | 0.9674 | 0.0141 | 0.9529 | 0.9601 | 0.9880 | 0.9961 | 0.9690 | 0.0049 | 0.9605 | 0.9664 | 0.9734 | 0.9786 | 0.0279 | 0.0069 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | 0.8424 | 0.0433 | 0.8100 | 0.8349 | 0.9130 | 0.9448 | 0.8670 | 0.0156 | 0.8374 | 0.8596 | 0.8811 | 0.9007 | 0.0781 | 0.0215 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | 0.9114 | 0.0189 | 0.8982 | 0.9110 | 0.9414 | 0.9646 | 0.7496 | 0.0336 | 0.6929 | 0.7321 | 0.7764 | 0.8242 | 0.0304 | 0.0443 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 0.8694 | 0.0500 | 0.8260 | 0.8463 | 0.9343 | 0.9834 | 0.8308 | 0.0283 | 0.7785 | 0.8154 | 0.8564 | 0.8880 | 0.0880 | 0.0410 |

500 replikalı bootstrap-of-bootstrap özetinde medyan IQR genişliği N=6 için 0.030383, N=50 için 0.019307; oran=0.6355. Böylece koşullu R tahmin bandı daraldı; %36.5'lik daralma yetersiz replika sayısından gelen tahmin gürültüsü hipotezini destekler. Bununla birlikte önceki 27 feature'ın nokta R aralığı 0.8001–0.9674'den 0.7228–0.9690'a genişledi ve 12 feature eşiği kaybetti. Dolayısıyla artan N oynaklığı azaltmış, fakat rank heterojenliğini ortadan kaldırmamıştır; bu koşu kalan heterojenliği tek başına biyolojik sinyal olarak tanımlayamaz. Bu bantlar gözlenen eğitim replikalarını yeniden örnekler; hasta düzeyi veya biyolojik güven aralığı değildir.

### N=6 → N=50 önceki 27 feature

| feature | block | G_claim_n6 | R_n6 | WI_n6 | G_claim | R | WI | bottleneck | point_adequate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000133742.14 | rna | 0.9579 | 0.8347 | 0.8347 | 0.9757 | 0.8412 | 0.8412 | R | 1.0000 |
| rna_ENSG00000169800.14 | rna | 0.9850 | 0.8358 | 0.8358 | 0.9801 | 0.8331 | 0.8331 | R | 1.0000 |
| rna_ENSG00000090402.8 | rna | 0.9150 | 0.8505 | 0.8505 | 0.8790 | 0.7641 | 0.7641 | R | 0.0000 |
| rna_ENSG00000139874.6 | rna | 0.9904 | 0.8505 | 0.8505 | 0.9869 | 0.7989 | 0.7989 | R | 0.0000 |
| rna_ENSG00000198099.9 | rna | 0.9859 | 0.8605 | 0.8605 | 0.9748 | 0.7421 | 0.7421 | R | 0.0000 |
| rna_ENSG00000270069.1 | rna | 0.9968 | 0.8266 | 0.8266 | 0.9954 | 0.7762 | 0.7762 | R | 0.0000 |
| rna_ENSG00000174899.12 | rna | 0.9920 | 0.8014 | 0.8014 | 0.9897 | 0.7832 | 0.7832 | R | 0.0000 |
| rna_ENSG00000168917.9 | rna | 0.9976 | 0.8134 | 0.8134 | 0.9875 | 0.7522 | 0.7522 | R | 0.0000 |
| rna_ENSG00000270547.6 | rna | 0.9753 | 0.8007 | 0.8007 | 0.9715 | 0.8200 | 0.8200 | R | 1.0000 |
| rna_ENSG00000133800.9 | rna | 0.9960 | 0.8334 | 0.8334 | 0.9967 | 0.8367 | 0.8367 | R | 1.0000 |
| rna_ENSG00000077616.11 | rna | 0.9892 | 0.8467 | 0.8467 | 0.9919 | 0.8118 | 0.8118 | R | 1.0000 |
| rna_ENSG00000226674.11 | rna | 0.9913 | 0.8579 | 0.8579 | 0.9815 | 0.7657 | 0.7657 | R | 0.0000 |
| rna_ENSG00000109158.11 | rna | 0.9568 | 0.8350 | 0.8350 | 0.9605 | 0.8335 | 0.8335 | R | 1.0000 |
| rna_ENSG00000023892.11 | rna | 0.9975 | 0.8171 | 0.8171 | 0.9974 | 0.8036 | 0.8036 | R | 1.0000 |
| rna_ENSG00000118407.15 | rna | 0.9933 | 0.8001 | 0.8001 | 0.9950 | 0.8050 | 0.8050 | R | 1.0000 |
| rna_ENSG00000108639.8 | rna | 0.9970 | 0.8290 | 0.8290 | 0.9930 | 0.7781 | 0.7781 | R | 0.0000 |
| rna_ENSG00000076043.11 | rna | 0.9706 | 0.8249 | 0.8249 | 0.9769 | 0.7228 | 0.7228 | R | 0.0000 |
| rna_ENSG00000182173.13 | rna | 0.9967 | 0.8566 | 0.8566 | 0.9980 | 0.8222 | 0.8222 | R | 1.0000 |
| rna_ENSG00000240108.1 | rna | 0.8910 | 0.8081 | 0.8081 | 0.8959 | 0.7945 | 0.7945 | R | 0.0000 |
| rna_ENSG00000187758.8 | rna | 0.9825 | 0.8524 | 0.8524 | 0.9814 | 0.8371 | 0.8371 | R | 1.0000 |
| rna_ENSG00000104327.7 | rna | 0.9833 | 0.8103 | 0.8103 | 0.9791 | 0.7821 | 0.7821 | R | 0.0000 |
| rna_ENSG00000204246.3 | rna | 0.9572 | 0.8324 | 0.8324 | 0.9390 | 0.8371 | 0.8371 | R | 1.0000 |
| rna_ENSG00000205186.4 | rna | 0.9914 | 0.8426 | 0.8426 | 0.9709 | 0.8228 | 0.8228 | R | 1.0000 |
| clin_numeric__age_at_index.demographic | clin | 0.9513 | 0.9674 | 0.9513 | 0.9503 | 0.9690 | 0.9503 | G_claim | 1.0000 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | clin | 0.8888 | 0.8424 | 0.8424 | 0.8817 | 0.8670 | 0.8670 | R | 1.0000 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | clin | 0.9974 | 0.9114 | 0.9114 | 0.9977 | 0.7496 | 0.7496 | R | 0.0000 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | clin | 0.8682 | 0.8694 | 0.8682 | 0.8744 | 0.8308 | 0.8308 | R | 1.0000 |

N=6'da geçen 27 feature'ın 15 tanesi N=50'de hâlâ geçiyor; 12 tanesi artık geçmiyor.

### N=50 block bazında güncel agregat

| method | block | features | assessed | mean_WI | F | S | U | G_claim | R |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| phase3.2_actual_METABRIC_queries_n50 | clin | 16 | 5 | 0.8285 | 1.0000 | 0.9927 | 0.9953 | 0.9353 | 0.8037 |
| phase3.2_actual_METABRIC_queries_n50 | meth | 200 | 0 | NA | 1.0000 | 0.9944 | 0.9981 | NA | 0.7408 |
| phase3.2_actual_METABRIC_queries_n50 | mirna | 200 | 0 | NA | 1.0000 | 0.9936 | 0.9967 | NA | 0.7626 |
| phase3.2_actual_METABRIC_queries_n50 | rna | 200 | 72 | 0.7619 | 1.0000 | 0.9956 | 0.9974 | 0.9714 | 0.7825 |

### Çalışma süresi

N=6 için tarihsel tam timer saklanmamıştı: checkpoint ilk/son yazımı alt sınırı 199.8 sn, ardışık checkpoint medyan aralığı eklenmiş dürüst tahmin 233.7 sn (3.9 dk). N=50 koşusu bir kesintiden sonra checkpoint'lerden sürdürüldüğü için tek-process tam timer yoktur; N=6 tahmini ile 7–50. replikaların kaydedilmiş fit sürelerinin toplamı 2690.6 sn (44.8 dk) karşılaştırılabilir training wall-clock tahminidir. Aynı medyan replika maliyetiyle N=100 projeksiyonu 3943.8 sn (65.7 dk). Bunlar tek makineye özgüdür; veri okuma/raporlama/test toplamı değildir.

Tabloda F/S/U/R, modelde bulunan tüm blok feature'ları üzerinden; WI yalnız dış destek açısından değerlendirilebilen altküme üzerinden ortalanır. Yüksek WI yüksek feature önemi veya klinik yarar anlamına gelmez.

77 feature için WI hesaplandı, 22 feature tüm eşikleri geçti. Eksik G skorları ne sıfır ne bir ile dolduruldu.

Güncel bottleneck dağılımı: `{'R': 76, 'G_claim': 1}`. Faz 2.1 formula modüllerinin SHA256 değerleri değişmedi. 101 test geçti; paket satır/dal kapsamı %100, deney betikleri bu kapsam yüzdesine dahil değildir. Bu kontroller yöntemin klinik geçerliliğini kanıtlamaz.

## Dış kohortta değerlendirilebilen feature'lar

| feature | F | S | U | G_claim | R | WI | bottleneck | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000121351.8 | 1.0000 | 0.9903 | 0.9878 | 0.8629 | 0.6919 | 0.6919 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000178522.14 | 1.0000 | 0.9952 | 0.9989 | 0.8876 | 0.7309 | 0.7309 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000006116.4 | 1.0000 | 0.9769 | 0.9940 | 0.9573 | 0.7229 | 0.7229 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000080910.13 | 1.0000 | 0.9972 | 0.9986 | 0.8937 | 0.8146 | 0.8146 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000133742.14 | 1.0000 | 0.9988 | 0.9992 | 0.9757 | 0.8412 | 0.8412 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000113889.14 | 1.0000 | 0.9883 | 0.9948 | 0.8844 | 0.6862 | 0.6862 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000182591.6 | 1.0000 | 0.9890 | 0.9967 | 0.9383 | 0.7725 | 0.7725 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000169800.14 | 1.0000 | 1.0000 | 1.0000 | 0.9801 | 0.8331 | 0.8331 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000198670.12 | 1.0000 | 0.9973 | 0.9989 | 0.9173 | 0.7463 | 0.7463 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000090402.8 | 1.0000 | 0.9937 | 0.9985 | 0.8790 | 0.7641 | 0.7641 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000167654.18 | 1.0000 | 0.9943 | 0.9970 | 0.9848 | 0.8189 | 0.8189 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000233452.7 | 1.0000 | 0.9871 | 0.9947 | 0.9910 | 0.7274 | 0.7274 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000139874.6 | 1.0000 | 0.9970 | 0.9987 | 0.9869 | 0.7989 | 0.7989 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000198099.9 | 1.0000 | 0.9941 | 0.9978 | 0.9748 | 0.7421 | 0.7421 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000169760.17 | 1.0000 | 0.9949 | 0.9969 | 0.9687 | 0.7033 | 0.7033 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000231861.2 | 1.0000 | 0.9930 | 0.9966 | 0.9535 | 0.7190 | 0.7190 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000270069.1 | 1.0000 | 0.9952 | 0.9987 | 0.9954 | 0.7762 | 0.7762 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000174899.12 | 1.0000 | 0.9970 | 0.9988 | 0.9897 | 0.7832 | 0.7832 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000100911.16 | 1.0000 | 0.9961 | 0.9981 | 0.9824 | 0.7171 | 0.7171 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000113389.16 | 1.0000 | 0.9902 | 0.9975 | 0.9835 | 0.6871 | 0.6871 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000169432.18 | 1.0000 | 0.9979 | 0.9980 | 0.9676 | 0.7369 | 0.7369 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000157617.17 | 1.0000 | 0.9954 | 0.9971 | 0.9807 | 0.7020 | 0.7020 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000158122.12 | 1.0000 | 0.9934 | 0.9979 | 0.9917 | 0.7475 | 0.7475 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000168917.9 | 1.0000 | 0.9970 | 0.9985 | 0.9875 | 0.7522 | 0.7522 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000270547.6 | 1.0000 | 0.9950 | 0.9942 | 0.9715 | 0.8200 | 0.8200 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000133800.9 | 1.0000 | 0.9983 | 0.9993 | 0.9967 | 0.8367 | 0.8367 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000077616.11 | 1.0000 | 0.9975 | 0.9990 | 0.9919 | 0.8118 | 0.8118 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000226674.11 | 1.0000 | 0.9893 | 0.9965 | 0.9815 | 0.7657 | 0.7657 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000064042.18 | 1.0000 | 0.9963 | 0.9965 | 0.9838 | 0.7232 | 0.7232 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000138795.10 | 1.0000 | 0.9920 | 0.9981 | 0.9944 | 0.7585 | 0.7585 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000100908.14 | 1.0000 | 0.9975 | 0.9984 | 0.9852 | 0.7715 | 0.7715 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000266964.6 | 1.0000 | 0.9972 | 0.9984 | 0.9711 | 0.7401 | 0.7401 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000109158.11 | 1.0000 | 0.9990 | 0.9994 | 0.9605 | 0.8335 | 0.8335 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000060709.15 | 1.0000 | 0.9969 | 0.9970 | 0.9722 | 0.7748 | 0.7748 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000004846.16 | 1.0000 | 0.9982 | 0.9989 | 0.9748 | 0.8056 | 0.8056 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000023892.11 | 1.0000 | 0.9970 | 0.9990 | 0.9974 | 0.8036 | 0.8036 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000214575.9 | 1.0000 | 0.9886 | 0.9972 | 0.9906 | 0.6857 | 0.6857 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000231925.12 | 1.0000 | 0.9960 | 0.9984 | 0.9883 | 0.7427 | 0.7427 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000183747.12 | 1.0000 | 0.9956 | 0.9982 | 0.9649 | 0.8049 | 0.8049 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000172901.20 | 1.0000 | 0.9989 | 0.9991 | 0.9802 | 0.8253 | 0.8253 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000118407.15 | 1.0000 | 0.9970 | 0.9988 | 0.9950 | 0.8050 | 0.8050 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000187531.14 | 1.0000 | 0.9966 | 0.9984 | 0.9873 | 0.7416 | 0.7416 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000108639.8 | 1.0000 | 0.9954 | 0.9985 | 0.9930 | 0.7781 | 0.7781 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000116663.11 | 1.0000 | 0.9916 | 0.9983 | 0.9975 | 0.7291 | 0.7291 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000124343.14 | 1.0000 | 0.9925 | 0.9983 | 0.9925 | 0.7307 | 0.7307 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000076043.11 | 1.0000 | 0.9916 | 0.9962 | 0.9769 | 0.7228 | 0.7228 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000134531.10 | 1.0000 | 0.9942 | 0.9984 | 0.9952 | 0.7756 | 0.7756 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000182173.13 | 1.0000 | 0.9974 | 0.9991 | 0.9980 | 0.8222 | 0.8222 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000172725.14 | 1.0000 | 0.9968 | 0.9982 | 0.9870 | 0.7482 | 0.7482 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000240108.1 | 1.0000 | 0.9968 | 0.9973 | 0.8959 | 0.7945 | 0.7945 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000187758.8 | 1.0000 | 0.9981 | 0.9991 | 0.9814 | 0.8371 | 0.8371 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000092010.15 | 1.0000 | 0.9966 | 0.9988 | 0.9979 | 0.8129 | 0.8129 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000104327.7 | 1.0000 | 0.9916 | 0.9965 | 0.9791 | 0.7821 | 0.7821 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000204246.3 | 1.0000 | 0.9982 | 0.9984 | 0.9390 | 0.8371 | 0.8371 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000064961.19 | 1.0000 | 0.9972 | 0.9989 | 0.9893 | 0.7888 | 0.7888 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000124440.16 | 1.0000 | 0.9985 | 0.9971 | 0.9534 | 0.7378 | 0.7378 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000205186.4 | 1.0000 | 0.9988 | 0.9992 | 0.9709 | 0.8228 | 0.8228 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000163145.13 | 1.0000 | 0.9960 | 0.9972 | 0.9865 | 0.7018 | 0.7018 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000181788.4 | 1.0000 | 0.9933 | 0.9980 | 0.9925 | 0.6940 | 0.6940 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000246430.7 | 1.0000 | 0.9959 | 0.9979 | 0.9654 | 0.7465 | 0.7465 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000072657.9 | 1.0000 | 0.9949 | 0.9981 | 0.9775 | 0.7951 | 0.7951 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000121871.4 | 1.0000 | 0.9921 | 0.9980 | 0.9829 | 0.7811 | 0.7811 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000169764.16 | 1.0000 | 0.9964 | 0.9974 | 0.9776 | 0.7275 | 0.7275 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000140332.17 | 1.0000 | 0.9966 | 0.9976 | 0.9839 | 0.7095 | 0.7095 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000231672.8 | 1.0000 | 0.9873 | 0.9968 | 0.9917 | 0.7284 | 0.7284 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000163545.11 | 1.0000 | 0.9957 | 0.9982 | 0.9879 | 0.7281 | 0.7281 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000187242.6 | 1.0000 | 0.9936 | 0.9970 | 0.9488 | 0.7782 | 0.7782 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000250334.6 | 1.0000 | 0.9943 | 0.9975 | 0.9856 | 0.7237 | 0.7237 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000187486.5 | 1.0000 | 0.9975 | 0.9979 | 0.9766 | 0.7177 | 0.7177 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000095203.14 | 1.0000 | 0.9918 | 0.9967 | 0.9752 | 0.6895 | 0.6895 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000149557.14 | 1.0000 | 0.9937 | 0.9985 | 0.9961 | 0.7442 | 0.7442 | R | assessed (different real cohort query distributions) |
| rna_ENSG00000122025.15 | 1.0000 | 0.9943 | 0.9952 | 0.9631 | 0.8096 | 0.8096 | R | assessed (different real cohort query distributions) |
| clin_numeric__age_at_index.demographic | 1.0000 | 0.9889 | 0.9909 | 0.9503 | 0.9690 | 0.9503 | G_claim | assessed (different real cohort query distributions) |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | 1.0000 | 0.9940 | 0.9928 | 0.8817 | 0.8670 | 0.8670 | R | assessed (different real cohort query distributions) |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | 1.0000 | 0.9846 | 0.9990 | 0.9977 | 0.7496 | 0.7496 | R | assessed (different real cohort query distributions) |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 1.0000 | 0.9881 | 0.9950 | 0.8744 | 0.8308 | 0.8308 | R | assessed (different real cohort query distributions) |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | 1.0000 | 0.9946 | 0.9840 | 0.9725 | 0.7447 | 0.7447 | R | assessed (different real cohort query distributions) |

## Not assessable — dış kohort/modalite yok veya modelden dışlanmış

miRNA ve methylation için G/G_signed/WI NA; modelde bulunan sütunların TCGA F/S/U/R değerleri aşağıda korunur. Modelden separation nedeniyle dışlanan feature'ların tüm bileşenleri NA'dır.

| feature | F | S | U | G_claim | R | WI | bottleneck | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000229086.5 | 1.0000 | 0.9993 | 0.9997 | NA | 0.8237 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000252076.1 | 1.0000 | 0.9997 | 0.9995 | NA | 0.8372 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000201919.1 | 1.0000 | 0.9980 | 0.9989 | NA | 0.8339 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000234863.1 | 1.0000 | 0.9951 | 0.9953 | NA | 0.8105 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000236485.1 | 1.0000 | 0.9981 | 0.9992 | NA | 0.8118 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000213958.2 | 1.0000 | 0.9871 | 0.9947 | NA | 0.7340 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000253940.1 | 1.0000 | 0.9997 | 1.0000 | NA | 0.8257 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000207788.1 | 1.0000 | 0.9999 | 1.0000 | NA | 0.8344 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000285569.2 | 1.0000 | 0.9940 | 0.9984 | NA | 0.7831 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000230161.1 | 1.0000 | 0.9933 | 0.9948 | NA | 0.8063 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000240320.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8326 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000286371.1 | 1.0000 | 0.9963 | 0.9976 | NA | 0.8294 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000224250.1 | 1.0000 | 0.9896 | 0.9965 | NA | 0.7379 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000201523.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8334 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000278687.1 | 1.0000 | 0.9999 | 1.0000 | NA | 0.8400 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000199220.1 | 1.0000 | 0.9954 | 0.9961 | NA | 0.8275 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000285684.1 | 1.0000 | 0.9922 | 0.9953 | NA | 0.7954 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000249994.1 | 1.0000 | 0.9952 | 0.9986 | NA | 0.7624 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000252978.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8336 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000227143.1 | 1.0000 | 0.9941 | 0.9964 | NA | 0.7943 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000278672.1 | 1.0000 | 0.9905 | 0.9919 | NA | 0.7920 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000266448.1 | 1.0000 | 0.9974 | 0.9978 | NA | 0.8186 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000224935.2 | 1.0000 | 0.9956 | 0.9981 | NA | 0.8136 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000215368.5 | 1.0000 | 0.9996 | 1.0000 | NA | 0.8320 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000259694.1 | 1.0000 | 0.9890 | 0.9852 | NA | 0.7138 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000286720.1 | 1.0000 | 0.9989 | 0.9996 | NA | 0.8064 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000206891.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8335 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000212532.1 | 1.0000 | 0.9977 | 0.9974 | NA | 0.7938 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000251264.1 | 1.0000 | 0.9967 | 0.9983 | NA | 0.7799 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000237242.1 | 1.0000 | 0.9968 | 0.9987 | NA | 0.8257 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000261340.1 | 1.0000 | 0.9995 | 0.9996 | NA | 0.8387 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000286103.1 | 1.0000 | 0.9961 | 0.9985 | NA | 0.8181 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000232395.1 | 1.0000 | 0.9950 | 0.9991 | NA | 0.7219 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000260620.1 | 1.0000 | 0.9999 | 1.0000 | NA | 0.8381 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000238077.1 | 1.0000 | 0.9964 | 0.9989 | NA | 0.7889 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000260639.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8333 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000253776.1 | 1.0000 | 0.9989 | 0.9994 | NA | 0.8309 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000270123.4 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8365 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000236656.1 | 1.0000 | 0.9894 | 0.9493 | NA | 0.7051 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000229834.1 | 1.0000 | 0.9958 | 0.9953 | NA | 0.8480 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000258803.1 | 1.0000 | 0.9980 | 0.9990 | NA | 0.8087 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000200434.1 | 1.0000 | 0.9985 | 0.9991 | NA | 0.8161 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000203876.10 | 1.0000 | 0.9956 | 0.9987 | NA | 0.7663 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000253488.1 | 1.0000 | 0.9975 | 0.9990 | NA | 0.8232 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000260413.1 | 1.0000 | 0.9999 | 1.0000 | NA | 0.8405 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000274244.1 | 1.0000 | 0.9948 | 0.9981 | NA | 0.8042 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000243995.3 | 1.0000 | 0.9928 | 0.9971 | NA | 0.8130 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000232728.1 | 1.0000 | 0.9961 | 0.9989 | NA | 0.7819 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000201752.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8334 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000285990.1 | 1.0000 | 0.9961 | 0.9928 | NA | 0.7367 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000241777.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8334 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000250923.1 | 1.0000 | 0.9998 | 1.0000 | NA | 0.8303 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000213607.4 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8330 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000255618.2 | 1.0000 | 0.9976 | 0.9984 | NA | 0.8373 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000254751.3 | 1.0000 | 0.9996 | 1.0000 | NA | 0.8493 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000256484.1 | 1.0000 | 0.9915 | 0.9931 | NA | 0.7538 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000249479.1 | 1.0000 | 0.9965 | 0.9972 | NA | 0.8077 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000226695.1 | 1.0000 | 0.9968 | 0.9975 | NA | 0.8462 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000224922.1 | 1.0000 | 0.9998 | 1.0000 | NA | 0.8390 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000252891.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8332 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000250391.1 | 1.0000 | 0.9953 | 0.9973 | NA | 0.8001 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000258115.1 | 1.0000 | 0.9979 | 0.9988 | NA | 0.8281 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000285831.1 | 1.0000 | 0.9965 | 0.9972 | NA | 0.7155 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000286252.1 | 1.0000 | 0.9983 | 0.9982 | NA | 0.8425 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000251338.1 | 1.0000 | 0.9997 | 0.9995 | NA | 0.8368 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000276118.1 | 1.0000 | 0.9995 | 1.0000 | NA | 0.8465 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000207097.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8380 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000212383.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8332 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000279376.1 | 1.0000 | 0.9990 | 0.9996 | NA | 0.8047 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000267125.2 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8329 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000226374.3 | 1.0000 | 0.9982 | 0.9993 | NA | 0.8189 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000201412.1 | 1.0000 | 0.9977 | 0.9985 | NA | 0.8287 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000285598.1 | 1.0000 | 0.9974 | 0.9987 | NA | 0.8431 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000267603.1 | 1.0000 | 0.9976 | 0.9985 | NA | 0.8040 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000248328.1 | 1.0000 | 0.9951 | 0.9977 | NA | 0.8088 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000206871.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8330 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000261103.3 | 1.0000 | 0.9891 | 0.9975 | NA | 0.7297 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000272855.1 | 1.0000 | 0.9946 | 0.9953 | NA | 0.8121 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000225399.4 | 1.0000 | 0.9970 | 0.9981 | NA | 0.7153 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000231426.7 | 1.0000 | 0.9904 | 0.9919 | NA | 0.6965 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000253346.1 | 1.0000 | 0.9889 | 0.9899 | NA | 0.7228 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000239438.1 | 1.0000 | 0.9999 | 1.0000 | NA | 0.8393 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000248869.7 | 1.0000 | 0.9966 | 0.9988 | NA | 0.8081 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000237560.7 | 1.0000 | 0.9965 | 0.9969 | NA | 0.7328 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000254403.1 | 1.0000 | 0.9938 | 0.9973 | NA | 0.7710 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000176134.5 | 1.0000 | 0.9949 | 0.9943 | NA | 0.7428 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000232821.1 | 1.0000 | 0.9984 | 0.9991 | NA | 0.8473 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000233884.2 | 1.0000 | 0.9978 | 0.9984 | NA | 0.8227 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000252969.1 | 1.0000 | 0.9979 | 0.9954 | NA | 0.8131 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000215474.8 | 1.0000 | 0.9984 | 0.9994 | NA | 0.8111 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000229237.2 | 1.0000 | 0.9932 | 0.9972 | NA | 0.6841 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000201809.1 | 1.0000 | 0.9975 | 0.9970 | NA | 0.8449 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000255162.1 | 1.0000 | 0.9874 | 0.9877 | NA | 0.7133 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000286892.1 | 1.0000 | 0.9988 | 0.9991 | NA | 0.8211 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000147223.6 | 1.0000 | 0.9918 | 0.9957 | NA | 0.6954 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000234723.1 | 1.0000 | 0.9966 | 0.9977 | NA | 0.8454 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000287616.1 | 1.0000 | 0.9969 | 0.9990 | NA | 0.7895 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000279204.1 | 1.0000 | 0.9988 | 0.9990 | NA | 0.8287 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000253177.2 | 1.0000 | 0.9979 | 0.9988 | NA | 0.7976 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000261472.1 | 1.0000 | 0.9975 | 0.9991 | NA | 0.8061 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000253105.6 | 1.0000 | 0.9959 | 0.9987 | NA | 0.7886 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000226622.6 | 1.0000 | 0.9989 | 0.9992 | NA | 0.8072 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000225131.2 | 1.0000 | 0.9938 | 0.9964 | NA | 0.6948 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000237870.6 | 1.0000 | 0.9959 | 0.9981 | NA | 0.7427 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000230147.2 | 1.0000 | 0.9933 | 0.9965 | NA | 0.7205 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000267583.5 | 1.0000 | 0.9927 | 0.9980 | NA | 0.6976 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000226592.1 | 1.0000 | 0.9967 | 0.9993 | NA | 0.7757 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000225056.1 | 1.0000 | 0.9954 | 0.9973 | NA | 0.7395 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000163632.14 | 1.0000 | 0.9833 | 0.9975 | NA | 0.7034 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000233648.1 | 1.0000 | 0.9903 | 0.9930 | NA | 0.7053 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000271776.1 | 1.0000 | 0.9958 | 0.9986 | NA | 0.7759 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000224365.1 | 1.0000 | 0.9869 | 0.9893 | NA | 0.7110 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000229498.1 | 1.0000 | 0.9960 | 0.9988 | NA | 0.7910 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000236466.2 | 1.0000 | 0.9857 | 0.9975 | NA | 0.6858 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000286750.1 | 1.0000 | 0.9971 | 0.9989 | NA | 0.8162 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000269124.1 | 1.0000 | 0.9974 | 0.9989 | NA | 0.8043 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000224609.8 | 1.0000 | 0.9932 | 0.9973 | NA | 0.7268 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000287179.1 | 1.0000 | 0.9949 | 0.9972 | NA | 0.8353 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000237202.1 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8338 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000227748.2 | 1.0000 | 0.9862 | 0.9867 | NA | 0.7095 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000229950.1 | 1.0000 | 0.9934 | 0.9939 | NA | 0.8539 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000249102.3 | 1.0000 | 0.9931 | 0.9956 | NA | 0.7082 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000250993.1 | 1.0000 | 0.9981 | 0.9987 | NA | 0.8491 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000254135.3 | 1.0000 | 0.9970 | 0.9989 | NA | 0.8106 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000279499.1 | 1.0000 | 0.9952 | 0.9986 | NA | 0.7768 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000280304.1 | 1.0000 | 0.9935 | 0.9985 | NA | 0.7568 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000287785.1 | 1.0000 | 0.9907 | 0.9921 | NA | 0.7408 | NA | not assessable | not assessable (no external cohort) |
| rna_ENSG00000237417.1 | 1.0000 | 0.9906 | 0.9979 | NA | 0.7314 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4533 | 1.0000 | 0.9916 | 0.9941 | NA | 0.7098 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-548ac | 1.0000 | 1.0000 | 1.0000 | NA | 0.8343 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4760 | 1.0000 | 0.9995 | 0.9990 | NA | 0.8406 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4493 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8323 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-640 | 1.0000 | 0.9997 | 1.0000 | NA | 0.8370 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4788 | 1.0000 | 0.9968 | 0.9982 | NA | 0.7134 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4776-1 | 1.0000 | 0.9944 | 0.9951 | NA | 0.7997 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-22 | 1.0000 | 0.9923 | 0.9969 | NA | 0.6860 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6715a | 1.0000 | 0.9915 | 0.9975 | NA | 0.6966 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-452 | 1.0000 | 0.9971 | 0.9991 | NA | 0.8090 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-891a | 1.0000 | 0.9921 | 0.9977 | NA | 0.7117 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-129-2 | 1.0000 | 0.9935 | 0.9967 | NA | 0.6980 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-190b | 1.0000 | 0.9955 | 0.9981 | NA | 0.6936 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-874 | 1.0000 | 0.9913 | 0.9969 | NA | 0.6899 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-224 | 1.0000 | 0.9945 | 0.9986 | NA | 0.7553 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-219a-2 | 1.0000 | 0.9911 | 0.9967 | NA | 0.6841 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-342 | 1.0000 | 0.9948 | 0.9984 | NA | 0.7822 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6883 | 1.0000 | 0.9896 | 0.9962 | NA | 0.7047 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-133a-2 | 1.0000 | 0.9982 | 0.9991 | NA | 0.8213 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-23c | 1.0000 | 0.9903 | 0.9941 | NA | 0.6975 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-618 | 1.0000 | 0.9943 | 0.9954 | NA | 0.7532 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1468 | 1.0000 | 0.9950 | 0.9983 | NA | 0.7564 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-490 | 1.0000 | 0.9934 | 0.9970 | NA | 0.7874 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-378a | 1.0000 | 0.9932 | 0.9985 | NA | 0.7560 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-486-1 | 1.0000 | 0.9973 | 0.9989 | NA | 0.8284 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-141 | 1.0000 | 0.9977 | 0.9979 | NA | 0.7858 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-486-2 | 1.0000 | 0.9971 | 0.9989 | NA | 0.8026 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6742 | 1.0000 | 0.9941 | 0.9960 | NA | 0.7074 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-192 | 1.0000 | 0.9970 | 0.9977 | NA | 0.7145 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-891b | 1.0000 | 0.9965 | 0.9945 | NA | 0.7525 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-221 | 1.0000 | 0.9921 | 0.9979 | NA | 0.7368 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-193a | 1.0000 | 0.9924 | 0.9976 | NA | 0.7044 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-133a-1 | 1.0000 | 0.9977 | 0.9991 | NA | 0.8159 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-200c | 1.0000 | 0.9975 | 0.9979 | NA | 0.7778 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5094 | 1.0000 | 0.9878 | 0.9909 | NA | 0.7627 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-378c | 1.0000 | 0.9948 | 0.9959 | NA | 0.7681 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4435-2 | 1.0000 | 0.9898 | 0.9953 | NA | 0.7009 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4724 | 1.0000 | 0.9958 | 0.9984 | NA | 0.7537 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-466 | 1.0000 | 0.9895 | 0.9916 | NA | 0.7597 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-129-1 | 1.0000 | 0.9950 | 0.9984 | NA | 0.7655 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-33b | 1.0000 | 0.9929 | 0.9975 | NA | 0.6916 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3199-2 | 1.0000 | 0.9891 | 0.9981 | NA | 0.7134 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4458 | 1.0000 | 0.9954 | 0.9976 | NA | 0.8305 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4751 | 1.0000 | 0.9964 | 0.9952 | NA | 0.7569 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-16-1 | 1.0000 | 0.9961 | 0.9980 | NA | 0.7240 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-135a-2 | 1.0000 | 0.9937 | 0.9959 | NA | 0.7117 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5591 | 1.0000 | 0.9995 | 1.0000 | NA | 0.8306 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-375 | 1.0000 | 0.9904 | 0.9972 | NA | 0.6706 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-96 | 1.0000 | 0.9975 | 0.9981 | NA | 0.7357 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-365b | 1.0000 | 0.9963 | 0.9986 | NA | 0.7990 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-664b | 1.0000 | 0.9922 | 0.9923 | NA | 0.9354 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-365a | 1.0000 | 0.9962 | 0.9986 | NA | 0.8002 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4738 | 1.0000 | 0.9935 | 0.9947 | NA | 0.8162 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-139 | 1.0000 | 0.9973 | 0.9973 | NA | 0.7194 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-378d-2 | 1.0000 | 0.9924 | 0.9982 | NA | 0.7380 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-16-2 | 1.0000 | 0.9967 | 0.9982 | NA | 0.7345 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-2355 | 1.0000 | 0.9932 | 0.9978 | NA | 0.7224 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6798 | 1.0000 | 0.9934 | 0.9958 | NA | 0.7139 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4678 | 1.0000 | 0.9986 | 0.9993 | NA | 0.8205 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-549a | 1.0000 | 0.9909 | 0.9938 | NA | 0.8735 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-222 | 1.0000 | 0.9941 | 0.9984 | NA | 0.7588 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4469 | 1.0000 | 0.9884 | 0.9953 | NA | 0.6873 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-182 | 1.0000 | 0.9964 | 0.9959 | NA | 0.7486 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-424 | 1.0000 | 0.9956 | 0.9967 | NA | 0.7213 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-636 | 1.0000 | 0.9933 | 0.9928 | NA | 0.9139 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-548at | 1.0000 | 0.9945 | 0.9954 | NA | 0.7523 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3170 | 1.0000 | 0.9940 | 0.9957 | NA | 0.8207 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-605 | 1.0000 | 0.9912 | 0.9968 | NA | 0.6897 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-556 | 1.0000 | 0.9950 | 0.9986 | NA | 0.7851 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4732 | 1.0000 | 0.9974 | 0.9989 | NA | 0.8209 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-892a | 1.0000 | 0.9822 | 0.9837 | NA | 0.7670 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4528 | 1.0000 | 0.9915 | 0.9943 | NA | 0.8216 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-let-7b | 1.0000 | 0.9943 | 0.9965 | NA | 0.6976 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4771-2 | 1.0000 | 0.9912 | 0.9929 | NA | 0.7056 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6874 | 1.0000 | 0.9938 | 0.9987 | NA | 0.7521 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-431 | 1.0000 | 0.9945 | 0.9982 | NA | 0.7403 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-191 | 1.0000 | 0.9975 | 0.9989 | NA | 0.8071 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1293 | 1.0000 | 0.9956 | 0.9979 | NA | 0.7414 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-548o | 1.0000 | 0.9876 | 0.9944 | NA | 0.8166 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3683 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8355 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5683 | 1.0000 | 0.9918 | 0.9981 | NA | 0.7142 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4684 | 1.0000 | 0.9917 | 0.9947 | NA | 0.6957 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6780a | 1.0000 | 0.9909 | 0.9914 | NA | 0.7317 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-21 | 1.0000 | 0.9965 | 0.9985 | NA | 0.7653 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5682 | 1.0000 | 0.9995 | 0.9988 | NA | 0.8548 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-585 | 1.0000 | 0.9902 | 0.9974 | NA | 0.6931 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5694 | 1.0000 | 0.9949 | 0.9979 | NA | 0.7431 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3688-1 | 1.0000 | 0.9868 | 0.9962 | NA | 0.7727 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4655 | 1.0000 | 0.9902 | 0.9918 | NA | 0.7661 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6796 | 1.0000 | 0.9824 | 0.9912 | NA | 0.7507 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-28 | 1.0000 | 0.9942 | 0.9983 | NA | 0.7559 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4453 | 1.0000 | 0.9904 | 0.9935 | NA | 0.7857 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-24-2 | 1.0000 | 0.9954 | 0.9987 | NA | 0.7868 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5191 | 1.0000 | 0.9962 | 0.9974 | NA | 0.7784 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-147b | 1.0000 | 0.9943 | 0.9985 | NA | 0.7646 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3907 | 1.0000 | 0.9968 | 0.9975 | NA | 0.8416 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-144 | 1.0000 | 0.9960 | 0.9988 | NA | 0.7948 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-24-1 | 1.0000 | 0.9954 | 0.9987 | NA | 0.7955 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4803 | 1.0000 | 0.9925 | 0.9957 | NA | 0.8159 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5691 | 1.0000 | 0.9945 | 0.9982 | NA | 0.7267 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-548f-1 | 1.0000 | 0.9975 | 0.9988 | NA | 0.8226 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-200a | 1.0000 | 0.9983 | 0.9963 | NA | 0.7402 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1245b | 1.0000 | 0.9856 | 0.9916 | NA | 0.7738 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3688-2 | 1.0000 | 0.9873 | 0.9891 | NA | 0.7614 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4532 | 1.0000 | 0.9863 | 0.9944 | NA | 0.7737 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5588 | 1.0000 | 0.9966 | 0.9939 | NA | 0.8262 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1225 | 1.0000 | 0.9890 | 0.9937 | NA | 0.7107 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-335 | 1.0000 | 0.9948 | 0.9983 | NA | 0.7916 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4658 | 1.0000 | 0.9980 | 0.9968 | NA | 0.7597 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3661 | 1.0000 | 0.9808 | 0.9865 | NA | 0.6937 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3622a | 1.0000 | 0.9840 | 0.9951 | NA | 0.7073 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-146a | 1.0000 | 0.9899 | 0.9930 | NA | 0.8650 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5100 | 1.0000 | 0.9997 | 1.0000 | NA | 0.8406 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-662 | 1.0000 | 0.9951 | 0.9966 | NA | 0.8439 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4796 | 1.0000 | 0.9941 | 0.9950 | NA | 0.7404 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-135a-1 | 1.0000 | 0.9920 | 0.9977 | NA | 0.7339 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-574 | 1.0000 | 0.9954 | 0.9982 | NA | 0.7312 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3162 | 1.0000 | 0.9911 | 0.9960 | NA | 0.8220 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1258 | 1.0000 | 0.9946 | 0.9986 | NA | 0.7632 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-145 | 1.0000 | 0.9978 | 0.9964 | NA | 0.7468 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6851 | 1.0000 | 0.9910 | 0.9940 | NA | 0.8052 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-659 | 1.0000 | 0.9946 | 0.9983 | NA | 0.7556 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4670 | 1.0000 | 0.9957 | 0.9986 | NA | 0.7716 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5705 | 1.0000 | 0.9958 | 0.9969 | NA | 0.8433 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-8073 | NA | NA | NA | NA | NA | NA | not assessable | not assessable (excluded from Cox model: separation) |
| mirna_hsa-mir-6716 | 1.0000 | 0.9922 | 0.9980 | NA | 0.7280 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4714 | 1.0000 | 0.9919 | 0.9976 | NA | 0.7094 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5192 | 1.0000 | 0.9997 | 1.0000 | NA | 0.8540 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-217 | 1.0000 | 0.9880 | 0.9977 | NA | 0.7098 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-20b | 1.0000 | 0.9960 | 0.9980 | NA | 0.7096 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6499 | 1.0000 | 0.9881 | 0.9971 | NA | 0.7233 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-30a | 1.0000 | 0.9929 | 0.9961 | NA | 0.7147 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-142 | 1.0000 | 0.9967 | 0.9987 | NA | 0.7727 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-429 | 1.0000 | 0.9974 | 0.9972 | NA | 0.7183 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6717 | 1.0000 | 0.9790 | 0.9935 | NA | 0.6734 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1257 | 1.0000 | 0.9942 | 0.9970 | NA | 0.8091 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6746 | 1.0000 | 0.9925 | 0.9974 | NA | 0.7166 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-106a | 1.0000 | 0.9943 | 0.9980 | NA | 0.7224 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6759 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8360 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-548f-4 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8360 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6892 | 1.0000 | 0.9918 | 0.9964 | NA | 0.6862 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-133b | 1.0000 | 0.9978 | 0.9992 | NA | 0.8265 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4720 | 1.0000 | 0.9924 | 0.9906 | NA | 0.7296 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-200b | 1.0000 | 0.9972 | 0.9967 | NA | 0.7258 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3139 | 1.0000 | 0.9931 | 0.9914 | NA | 0.7507 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4520-2 | 1.0000 | 0.9982 | 0.9979 | NA | 0.8452 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-let-7i | 1.0000 | 0.9886 | 0.9979 | NA | 0.7083 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-148a | 1.0000 | 0.9932 | 0.9976 | NA | 0.6977 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-125b-1 | 1.0000 | 0.9962 | 0.9971 | NA | 0.7104 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-541 | 1.0000 | 0.9940 | 0.9987 | NA | 0.7424 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6715b | 1.0000 | 0.9928 | 0.9977 | NA | 0.7624 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-320b-2 | 1.0000 | 0.9865 | 0.9977 | NA | 0.6980 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6830 | 1.0000 | 0.9843 | 0.9914 | NA | 0.7168 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-922 | 1.0000 | 0.9966 | 0.9961 | NA | 0.8448 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-600 | 1.0000 | 0.9918 | 0.9970 | NA | 0.7074 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-125b-2 | 1.0000 | 0.9961 | 0.9971 | NA | 0.7103 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4668 | 1.0000 | 0.9962 | 0.9989 | NA | 0.7821 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1915 | 1.0000 | 0.9953 | 0.9979 | NA | 0.8193 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1-2 | 1.0000 | 0.9957 | 0.9988 | NA | 0.7819 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4450 | 1.0000 | 0.9958 | 0.9960 | NA | 0.8395 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5189 | 1.0000 | 0.9930 | 0.9954 | NA | 0.7378 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3619 | 1.0000 | 0.9942 | 0.9967 | NA | 0.7181 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1246 | 1.0000 | 0.9982 | 0.9986 | NA | 0.8362 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4716 | 1.0000 | 0.9937 | 0.9930 | NA | 0.7801 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-143 | 1.0000 | 0.9944 | 0.9978 | NA | 0.7043 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4502 | 1.0000 | 0.9962 | 0.9974 | NA | 0.8483 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-2114 | 1.0000 | 0.9931 | 0.9951 | NA | 0.7167 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4510 | 1.0000 | 0.9928 | 0.9952 | NA | 0.7806 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4641 | 1.0000 | 0.9956 | 0.9981 | NA | 0.8358 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6087 | 1.0000 | 0.9850 | 0.9967 | NA | 0.6774 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5693 | 1.0000 | 0.9873 | 0.9953 | NA | 0.7904 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-548ai | 1.0000 | 0.9949 | 0.9914 | NA | 0.7888 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-194-2 | 1.0000 | 0.9935 | 0.9983 | NA | 0.7425 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-326 | 1.0000 | 0.9930 | 0.9986 | NA | 0.7490 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6885 | 1.0000 | 0.9893 | 0.9973 | NA | 0.6912 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4435-1 | 1.0000 | 0.9941 | 0.9974 | NA | 0.8079 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6756 | 1.0000 | 0.9960 | 0.9986 | NA | 0.8036 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-140 | 1.0000 | 0.9943 | 0.9982 | NA | 0.7503 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-584 | 1.0000 | 0.9964 | 0.9967 | NA | 0.7110 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-586 | 1.0000 | 0.9952 | 0.9952 | NA | 0.8394 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6511b-2 | 1.0000 | 0.9948 | 0.9977 | NA | 0.6958 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3115 | 1.0000 | 0.9915 | 0.9967 | NA | 0.7665 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6825 | 1.0000 | 0.9896 | 0.9964 | NA | 0.7237 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5186 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8355 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-183 | 1.0000 | 0.9970 | 0.9955 | NA | 0.7816 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-5187 | 1.0000 | 0.9808 | 0.9968 | NA | 0.6904 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4452 | 1.0000 | 0.9977 | 0.9972 | NA | 0.8478 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-451a | 1.0000 | 0.9975 | 0.9980 | NA | 0.7350 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-550a-2 | 1.0000 | 0.9947 | 0.9946 | NA | 0.8804 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6852 | 1.0000 | 0.9927 | 0.9987 | NA | 0.7299 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-32 | 1.0000 | 0.9972 | 0.9989 | NA | 0.8256 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-592 | 1.0000 | 0.9932 | 0.9984 | NA | 0.7375 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4444-2 | 1.0000 | 0.9875 | 0.9971 | NA | 0.6827 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-6800 | 1.0000 | 0.9778 | 0.9906 | NA | 0.6853 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-31 | 1.0000 | 0.9925 | 0.9942 | NA | 0.8441 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4681 | 1.0000 | 1.0000 | 1.0000 | NA | 0.8329 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-181d | 1.0000 | 0.9875 | 0.9960 | NA | 0.6980 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-4770 | 1.0000 | 0.9977 | 0.9989 | NA | 0.8178 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-3936 | 1.0000 | 0.9856 | 0.9946 | NA | 0.7010 | NA | not assessable | not assessable (no external cohort) |
| mirna_hsa-mir-1294 | 1.0000 | 0.9904 | 0.9963 | NA | 0.7458 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23972735 | 1.0000 | 0.9944 | 0.9979 | NA | 0.6970 | NA | not assessable | not assessable (no external cohort) |
| meth_cg09856092 | 1.0000 | 0.9974 | 0.9931 | NA | 0.8715 | NA | not assessable | not assessable (no external cohort) |
| meth_cg09229620 | 1.0000 | 0.9936 | 0.9977 | NA | 0.6932 | NA | not assessable | not assessable (no external cohort) |
| meth_cg03854238 | 1.0000 | 0.9959 | 0.9973 | NA | 0.7097 | NA | not assessable | not assessable (no external cohort) |
| meth_cg01375382 | 1.0000 | 0.9967 | 0.9980 | NA | 0.8354 | NA | not assessable | not assessable (no external cohort) |
| meth_cg06738031 | 1.0000 | 0.9935 | 0.9974 | NA | 0.7033 | NA | not assessable | not assessable (no external cohort) |
| meth_cg13420364 | 1.0000 | 0.9906 | 0.9959 | NA | 0.8033 | NA | not assessable | not assessable (no external cohort) |
| meth_cg02467858 | 1.0000 | 0.9962 | 0.9988 | NA | 0.7418 | NA | not assessable | not assessable (no external cohort) |
| meth_cg12134633 | 1.0000 | 0.9970 | 0.9987 | NA | 0.7990 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18618429 | 1.0000 | 0.9954 | 0.9989 | NA | 0.7483 | NA | not assessable | not assessable (no external cohort) |
| meth_cg21616785 | 1.0000 | 0.9923 | 0.9983 | NA | 0.6855 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07205203 | 1.0000 | 0.9983 | 0.9992 | NA | 0.8239 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22077197 | 1.0000 | 0.9964 | 0.9983 | NA | 0.7136 | NA | not assessable | not assessable (no external cohort) |
| meth_cg09737197 | 1.0000 | 0.9957 | 0.9989 | NA | 0.7412 | NA | not assessable | not assessable (no external cohort) |
| meth_cg03532263 | 1.0000 | 0.9960 | 0.9987 | NA | 0.7186 | NA | not assessable | not assessable (no external cohort) |
| meth_cg15037823 | 1.0000 | 0.9905 | 0.9975 | NA | 0.6988 | NA | not assessable | not assessable (no external cohort) |
| meth_cg14669379 | 1.0000 | 0.9983 | 0.9993 | NA | 0.8213 | NA | not assessable | not assessable (no external cohort) |
| meth_cg01444716 | 1.0000 | 0.9967 | 0.9989 | NA | 0.7816 | NA | not assessable | not assessable (no external cohort) |
| meth_cg10665891 | 1.0000 | 0.9911 | 0.9966 | NA | 0.7087 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00778807 | 1.0000 | 0.9963 | 0.9990 | NA | 0.7551 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04723534 | 1.0000 | 0.9914 | 0.9982 | NA | 0.7060 | NA | not assessable | not assessable (no external cohort) |
| meth_cg26389380 | 1.0000 | 0.9952 | 0.9983 | NA | 0.7017 | NA | not assessable | not assessable (no external cohort) |
| meth_cg16962463 | 1.0000 | 0.9967 | 0.9990 | NA | 0.7788 | NA | not assessable | not assessable (no external cohort) |
| meth_cg16672562 | 1.0000 | 0.9974 | 0.9981 | NA | 0.7302 | NA | not assessable | not assessable (no external cohort) |
| meth_cg09136878 | 1.0000 | 0.9954 | 0.9987 | NA | 0.7501 | NA | not assessable | not assessable (no external cohort) |
| meth_cg01973456 | 1.0000 | 0.9967 | 0.9988 | NA | 0.7750 | NA | not assessable | not assessable (no external cohort) |
| meth_cg15429127 | 1.0000 | 0.9973 | 0.9928 | NA | 0.8689 | NA | not assessable | not assessable (no external cohort) |
| meth_cg26607031 | 1.0000 | 0.9973 | 0.9987 | NA | 0.7452 | NA | not assessable | not assessable (no external cohort) |
| meth_cg21729992 | 1.0000 | 0.9925 | 0.9975 | NA | 0.6946 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07904452 | 1.0000 | 0.9977 | 0.9993 | NA | 0.8191 | NA | not assessable | not assessable (no external cohort) |
| meth_cg19418318 | 1.0000 | 0.9941 | 0.9987 | NA | 0.7440 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07109453 | 1.0000 | 0.9961 | 0.9990 | NA | 0.7826 | NA | not assessable | not assessable (no external cohort) |
| meth_cg17648080 | 1.0000 | 0.9823 | 0.9972 | NA | 0.6760 | NA | not assessable | not assessable (no external cohort) |
| meth_cg27598107 | 1.0000 | 0.9919 | 0.9982 | NA | 0.7043 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18519308 | 1.0000 | 0.9968 | 0.9977 | NA | 0.7072 | NA | not assessable | not assessable (no external cohort) |
| meth_cg19880942 | 1.0000 | 0.9932 | 0.9981 | NA | 0.6828 | NA | not assessable | not assessable (no external cohort) |
| meth_cg14981189 | 1.0000 | 0.9954 | 0.9985 | NA | 0.7455 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23497383 | 1.0000 | 0.9962 | 0.9989 | NA | 0.7994 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00887153 | 1.0000 | 0.9926 | 0.9973 | NA | 0.6926 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00382999 | 1.0000 | 0.9948 | 0.9985 | NA | 0.7341 | NA | not assessable | not assessable (no external cohort) |
| meth_cg10296205 | 1.0000 | 0.9938 | 0.9948 | NA | 0.8861 | NA | not assessable | not assessable (no external cohort) |
| meth_cg27063372 | 1.0000 | 0.9908 | 0.9965 | NA | 0.7154 | NA | not assessable | not assessable (no external cohort) |
| meth_cg27494647 | 1.0000 | 0.9977 | 0.9992 | NA | 0.8306 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23824801 | 1.0000 | 0.9907 | 0.9971 | NA | 0.7248 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04178787 | 1.0000 | 0.9976 | 0.9992 | NA | 0.8027 | NA | not assessable | not assessable (no external cohort) |
| meth_cg14769786 | 1.0000 | 0.9921 | 0.9981 | NA | 0.7088 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22891070 | 1.0000 | 0.9969 | 0.9980 | NA | 0.7211 | NA | not assessable | not assessable (no external cohort) |
| meth_cg26010218 | 1.0000 | 0.9949 | 0.9985 | NA | 0.6887 | NA | not assessable | not assessable (no external cohort) |
| meth_cg06600725 | 1.0000 | 0.9949 | 0.9985 | NA | 0.7514 | NA | not assessable | not assessable (no external cohort) |
| meth_cg08047907 | 1.0000 | 0.9959 | 0.9986 | NA | 0.7636 | NA | not assessable | not assessable (no external cohort) |
| meth_cg08436419 | 1.0000 | 0.9941 | 0.9984 | NA | 0.7091 | NA | not assessable | not assessable (no external cohort) |
| meth_cg25181170 | 1.0000 | 0.9975 | 0.9991 | NA | 0.7731 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18142262 | 1.0000 | 0.9955 | 0.9987 | NA | 0.7431 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24504442 | 1.0000 | 0.9900 | 0.9965 | NA | 0.7174 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00892742 | 1.0000 | 0.9954 | 0.9985 | NA | 0.7338 | NA | not assessable | not assessable (no external cohort) |
| meth_cg27408285 | 1.0000 | 0.9899 | 0.9972 | NA | 0.7018 | NA | not assessable | not assessable (no external cohort) |
| meth_cg11157127 | 1.0000 | 0.9975 | 0.9988 | NA | 0.7502 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22675660 | 1.0000 | 0.9958 | 0.9988 | NA | 0.7765 | NA | not assessable | not assessable (no external cohort) |
| meth_cg19412295 | 1.0000 | 0.9940 | 0.9987 | NA | 0.7211 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24181266 | 1.0000 | 0.9933 | 0.9985 | NA | 0.7434 | NA | not assessable | not assessable (no external cohort) |
| meth_cg08923594 | 1.0000 | 0.9903 | 0.9962 | NA | 0.7120 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04831594 | 1.0000 | 0.9898 | 0.9940 | NA | 0.8460 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24492202 | 1.0000 | 0.9924 | 0.9987 | NA | 0.7379 | NA | not assessable | not assessable (no external cohort) |
| meth_cg17239057 | 1.0000 | 0.9927 | 0.9982 | NA | 0.7366 | NA | not assessable | not assessable (no external cohort) |
| meth_cg16610086 | 1.0000 | 0.9951 | 0.9990 | NA | 0.7681 | NA | not assessable | not assessable (no external cohort) |
| meth_cg08511651 | 1.0000 | 0.9972 | 0.9991 | NA | 0.7857 | NA | not assessable | not assessable (no external cohort) |
| meth_cg10932018 | 1.0000 | 0.9903 | 0.9981 | NA | 0.7001 | NA | not assessable | not assessable (no external cohort) |
| meth_cg03359666 | 1.0000 | 0.9917 | 0.9949 | NA | 0.8428 | NA | not assessable | not assessable (no external cohort) |
| meth_cg17165836 | 1.0000 | 0.9931 | 0.9982 | NA | 0.6943 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04682791 | 1.0000 | 0.9954 | 0.9989 | NA | 0.7809 | NA | not assessable | not assessable (no external cohort) |
| meth_cg19126300 | 1.0000 | 0.9942 | 0.9983 | NA | 0.6881 | NA | not assessable | not assessable (no external cohort) |
| meth_cg12580445 | 1.0000 | 0.9955 | 0.9991 | NA | 0.7636 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18276112 | 1.0000 | 0.9973 | 0.9991 | NA | 0.8205 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07565505 | 1.0000 | 0.9919 | 0.9985 | NA | 0.6867 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18456523 | 1.0000 | 0.9902 | 0.9982 | NA | 0.7065 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24797187 | 1.0000 | 0.9930 | 0.9957 | NA | 0.8706 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24865779 | 1.0000 | 0.9954 | 0.9982 | NA | 0.6904 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00680551 | 1.0000 | 0.9969 | 0.9984 | NA | 0.7186 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23557926 | 1.0000 | 0.9946 | 0.9987 | NA | 0.7505 | NA | not assessable | not assessable (no external cohort) |
| meth_cg13804575 | 1.0000 | 0.9969 | 0.9969 | NA | 0.7811 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22993195 | 1.0000 | 0.9951 | 0.9981 | NA | 0.6926 | NA | not assessable | not assessable (no external cohort) |
| meth_cg02983451 | 1.0000 | 0.9906 | 0.9955 | NA | 0.8184 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22081371 | 1.0000 | 0.9975 | 0.9985 | NA | 0.7242 | NA | not assessable | not assessable (no external cohort) |
| meth_cg13530938 | 1.0000 | 0.9956 | 0.9988 | NA | 0.7905 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00716257 | 1.0000 | 0.9938 | 0.9970 | NA | 0.7110 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22254463 | 1.0000 | 0.9898 | 0.9970 | NA | 0.6914 | NA | not assessable | not assessable (no external cohort) |
| meth_cg05130679 | 1.0000 | 0.9926 | 0.9966 | NA | 0.7221 | NA | not assessable | not assessable (no external cohort) |
| meth_cg15410675 | 1.0000 | 0.9937 | 0.9988 | NA | 0.7277 | NA | not assessable | not assessable (no external cohort) |
| meth_cg09439920 | 1.0000 | 0.9963 | 0.9988 | NA | 0.7608 | NA | not assessable | not assessable (no external cohort) |
| meth_cg17356733 | 1.0000 | 0.9949 | 0.9983 | NA | 0.7107 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04682802 | 1.0000 | 0.9957 | 0.9985 | NA | 0.7193 | NA | not assessable | not assessable (no external cohort) |
| meth_cg11564601 | 1.0000 | 0.9957 | 0.9989 | NA | 0.7428 | NA | not assessable | not assessable (no external cohort) |
| meth_cg10572943 | 1.0000 | 0.9956 | 0.9979 | NA | 0.6984 | NA | not assessable | not assessable (no external cohort) |
| meth_cg09866743 | 1.0000 | 0.9944 | 0.9988 | NA | 0.7498 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00857907 | 1.0000 | 0.9905 | 0.9976 | NA | 0.6777 | NA | not assessable | not assessable (no external cohort) |
| meth_cg27614178 | 1.0000 | 0.9932 | 0.9980 | NA | 0.7119 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18173184 | 1.0000 | 0.9944 | 0.9985 | NA | 0.7106 | NA | not assessable | not assessable (no external cohort) |
| meth_cg26353296 | 1.0000 | 0.9892 | 0.9983 | NA | 0.6974 | NA | not assessable | not assessable (no external cohort) |
| meth_cg20797905 | 1.0000 | 0.9955 | 0.9984 | NA | 0.7083 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07254421 | 1.0000 | 0.9907 | 0.9982 | NA | 0.6877 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07965640 | 1.0000 | 0.9930 | 0.9986 | NA | 0.7412 | NA | not assessable | not assessable (no external cohort) |
| meth_cg26354493 | 1.0000 | 0.9947 | 0.9982 | NA | 0.7282 | NA | not assessable | not assessable (no external cohort) |
| meth_cg26348348 | 1.0000 | 0.9974 | 0.9990 | NA | 0.7912 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18208842 | 1.0000 | 0.9946 | 0.9982 | NA | 0.7016 | NA | not assessable | not assessable (no external cohort) |
| meth_cg17415451 | 1.0000 | 0.9913 | 0.9985 | NA | 0.7172 | NA | not assessable | not assessable (no external cohort) |
| meth_cg08390172 | 1.0000 | 0.9901 | 0.9985 | NA | 0.7100 | NA | not assessable | not assessable (no external cohort) |
| meth_cg27628340 | 1.0000 | 0.9955 | 0.9989 | NA | 0.7732 | NA | not assessable | not assessable (no external cohort) |
| meth_cg20780180 | 1.0000 | 0.9963 | 0.9978 | NA | 0.7140 | NA | not assessable | not assessable (no external cohort) |
| meth_cg20391567 | 1.0000 | 0.9960 | 0.9986 | NA | 0.7625 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18450254 | 1.0000 | 0.9904 | 0.9975 | NA | 0.6831 | NA | not assessable | not assessable (no external cohort) |
| meth_cg10928544 | 1.0000 | 0.9920 | 0.9974 | NA | 0.7190 | NA | not assessable | not assessable (no external cohort) |
| meth_cg20312012 | 1.0000 | 0.9923 | 0.9960 | NA | 0.7495 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04942251 | 1.0000 | 0.9959 | 0.9983 | NA | 0.7279 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07937427 | 1.0000 | 0.9921 | 0.9982 | NA | 0.6870 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00314029 | 1.0000 | 0.9966 | 0.9986 | NA | 0.7399 | NA | not assessable | not assessable (no external cohort) |
| meth_cg02085507 | 1.0000 | 0.9954 | 0.9989 | NA | 0.7552 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18121066 | 1.0000 | 0.9973 | 0.9991 | NA | 0.8046 | NA | not assessable | not assessable (no external cohort) |
| meth_cg16493071 | 1.0000 | 0.9939 | 0.9960 | NA | 0.9086 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00800229 | 1.0000 | 0.9957 | 0.9985 | NA | 0.7306 | NA | not assessable | not assessable (no external cohort) |
| meth_cg12114888 | 1.0000 | 0.9969 | 0.9990 | NA | 0.8053 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23198707 | 1.0000 | 0.9956 | 0.9978 | NA | 0.7057 | NA | not assessable | not assessable (no external cohort) |
| meth_cg01623261 | 1.0000 | 0.9913 | 0.9974 | NA | 0.6875 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22728904 | 1.0000 | 0.9957 | 0.9985 | NA | 0.7240 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24162579 | 1.0000 | 0.9939 | 0.9982 | NA | 0.7421 | NA | not assessable | not assessable (no external cohort) |
| meth_cg08356637 | 1.0000 | 0.9949 | 0.9987 | NA | 0.7308 | NA | not assessable | not assessable (no external cohort) |
| meth_cg19516279 | 1.0000 | 0.9932 | 0.9987 | NA | 0.7350 | NA | not assessable | not assessable (no external cohort) |
| meth_cg01113811 | 1.0000 | 0.9949 | 0.9986 | NA | 0.7500 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00589541 | 1.0000 | 0.9958 | 0.9987 | NA | 0.7726 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07100000 | 1.0000 | 0.9962 | 0.9978 | NA | 0.7112 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24810735 | 1.0000 | 0.9877 | 0.9979 | NA | 0.7003 | NA | not assessable | not assessable (no external cohort) |
| meth_cg08481075 | 1.0000 | 0.9985 | 0.9968 | NA | 0.8104 | NA | not assessable | not assessable (no external cohort) |
| meth_cg25964728 | 1.0000 | 0.9967 | 0.9990 | NA | 0.7692 | NA | not assessable | not assessable (no external cohort) |
| meth_cg05890377 | 1.0000 | 0.9956 | 0.9988 | NA | 0.7639 | NA | not assessable | not assessable (no external cohort) |
| meth_cg25932869 | 1.0000 | 0.9940 | 0.9975 | NA | 0.7417 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24432768 | 1.0000 | 0.9923 | 0.9983 | NA | 0.7042 | NA | not assessable | not assessable (no external cohort) |
| meth_cg03930313 | 1.0000 | 0.9975 | 0.9966 | NA | 0.7498 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07570723 | 1.0000 | 0.9958 | 0.9989 | NA | 0.7377 | NA | not assessable | not assessable (no external cohort) |
| meth_cg02783887 | 1.0000 | 0.9932 | 0.9985 | NA | 0.7479 | NA | not assessable | not assessable (no external cohort) |
| meth_cg02863179 | 1.0000 | 0.9936 | 0.9986 | NA | 0.7428 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22753668 | 1.0000 | 0.9932 | 0.9983 | NA | 0.7396 | NA | not assessable | not assessable (no external cohort) |
| meth_cg06133110 | 1.0000 | 0.9905 | 0.9982 | NA | 0.6920 | NA | not assessable | not assessable (no external cohort) |
| meth_cg02013841 | 1.0000 | 0.9953 | 0.9984 | NA | 0.7484 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04850286 | 1.0000 | 0.9963 | 0.9983 | NA | 0.7298 | NA | not assessable | not assessable (no external cohort) |
| meth_cg02058870 | 1.0000 | 0.9934 | 0.9987 | NA | 0.7384 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23004031 | 1.0000 | 0.9953 | 0.9983 | NA | 0.7229 | NA | not assessable | not assessable (no external cohort) |
| meth_cg08613144 | 1.0000 | 0.9977 | 0.9991 | NA | 0.8087 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18261069 | 1.0000 | 0.9941 | 0.9985 | NA | 0.6949 | NA | not assessable | not assessable (no external cohort) |
| meth_cg10321156 | 1.0000 | 0.9961 | 0.9987 | NA | 0.7636 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00322927 | 1.0000 | 0.9935 | 0.9975 | NA | 0.6970 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24156393 | 1.0000 | 0.9943 | 0.9985 | NA | 0.7425 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07065111 | 1.0000 | 0.9955 | 0.9989 | NA | 0.7669 | NA | not assessable | not assessable (no external cohort) |
| meth_cg14377594 | 1.0000 | 0.9938 | 0.9988 | NA | 0.7365 | NA | not assessable | not assessable (no external cohort) |
| meth_cg25976440 | 1.0000 | 0.9975 | 0.9989 | NA | 0.7924 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00970396 | 1.0000 | 0.9970 | 0.9984 | NA | 0.7719 | NA | not assessable | not assessable (no external cohort) |
| meth_cg21143896 | 1.0000 | 0.9967 | 0.9989 | NA | 0.7588 | NA | not assessable | not assessable (no external cohort) |
| meth_cg03478739 | 1.0000 | 0.9951 | 0.9986 | NA | 0.7673 | NA | not assessable | not assessable (no external cohort) |
| meth_cg14429919 | 1.0000 | 0.9923 | 0.9982 | NA | 0.6978 | NA | not assessable | not assessable (no external cohort) |
| meth_cg16434331 | 1.0000 | 0.9910 | 0.9983 | NA | 0.7227 | NA | not assessable | not assessable (no external cohort) |
| meth_cg18420599 | 1.0000 | 0.9885 | 0.9983 | NA | 0.7146 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00081417 | 1.0000 | 0.9926 | 0.9984 | NA | 0.7422 | NA | not assessable | not assessable (no external cohort) |
| meth_cg11228682 | 1.0000 | 0.9905 | 0.9986 | NA | 0.7167 | NA | not assessable | not assessable (no external cohort) |
| meth_cg13045134 | 1.0000 | 0.9965 | 0.9979 | NA | 0.7000 | NA | not assessable | not assessable (no external cohort) |
| meth_cg09977969 | 1.0000 | 0.9932 | 0.9986 | NA | 0.7114 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00744866 | 1.0000 | 0.9931 | 0.9983 | NA | 0.6979 | NA | not assessable | not assessable (no external cohort) |
| meth_cg19021076 | 1.0000 | 0.9914 | 0.9986 | NA | 0.7128 | NA | not assessable | not assessable (no external cohort) |
| meth_cg15626285 | 1.0000 | 0.9950 | 0.9980 | NA | 0.7095 | NA | not assessable | not assessable (no external cohort) |
| meth_cg07863022 | 1.0000 | 0.9955 | 0.9958 | NA | 0.7830 | NA | not assessable | not assessable (no external cohort) |
| meth_cg12905085 | 1.0000 | 0.9921 | 0.9987 | NA | 0.7284 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00875805 | 1.0000 | 0.9937 | 0.9984 | NA | 0.7414 | NA | not assessable | not assessable (no external cohort) |
| meth_cg20378070 | 1.0000 | 0.9951 | 0.9990 | NA | 0.7585 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04085768 | 1.0000 | 0.9885 | 0.9977 | NA | 0.6794 | NA | not assessable | not assessable (no external cohort) |
| meth_cg11686819 | 1.0000 | 0.9934 | 0.9982 | NA | 0.7024 | NA | not assessable | not assessable (no external cohort) |
| meth_cg16562730 | 1.0000 | 0.9971 | 0.9988 | NA | 0.7596 | NA | not assessable | not assessable (no external cohort) |
| meth_cg15716680 | 1.0000 | 0.9963 | 0.9982 | NA | 0.7102 | NA | not assessable | not assessable (no external cohort) |
| meth_cg00768487 | 1.0000 | 0.9908 | 0.9979 | NA | 0.6830 | NA | not assessable | not assessable (no external cohort) |
| meth_cg21730858 | 1.0000 | 0.9958 | 0.9988 | NA | 0.7639 | NA | not assessable | not assessable (no external cohort) |
| meth_cg14826711 | 1.0000 | 0.9931 | 0.9984 | NA | 0.7236 | NA | not assessable | not assessable (no external cohort) |
| meth_cg10558233 | 1.0000 | 0.9969 | 0.9984 | NA | 0.7271 | NA | not assessable | not assessable (no external cohort) |
| meth_cg12444411 | 1.0000 | 0.9969 | 0.9991 | NA | 0.8019 | NA | not assessable | not assessable (no external cohort) |
| meth_cg09825309 | 1.0000 | 0.9936 | 0.9976 | NA | 0.6880 | NA | not assessable | not assessable (no external cohort) |
| meth_cg24453118 | 1.0000 | 0.9982 | 0.9980 | NA | 0.7281 | NA | not assessable | not assessable (no external cohort) |
| meth_cg03159396 | 1.0000 | 0.9945 | 0.9983 | NA | 0.7212 | NA | not assessable | not assessable (no external cohort) |
| meth_cg11574975 | 1.0000 | 0.9979 | 0.9988 | NA | 0.7700 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23285459 | 1.0000 | 0.9960 | 0.9991 | NA | 0.7861 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22193912 | 1.0000 | 0.9940 | 0.9987 | NA | 0.7505 | NA | not assessable | not assessable (no external cohort) |
| meth_cg17627629 | 1.0000 | 0.9977 | 0.9985 | NA | 0.7217 | NA | not assessable | not assessable (no external cohort) |
| meth_cg15925365 | 1.0000 | 0.9976 | 0.9970 | NA | 0.7477 | NA | not assessable | not assessable (no external cohort) |
| meth_cg15200418 | 1.0000 | 0.9971 | 0.9992 | NA | 0.7910 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04415176 | 1.0000 | 0.9928 | 0.9988 | NA | 0.7378 | NA | not assessable | not assessable (no external cohort) |
| meth_cg25296374 | 1.0000 | 0.9948 | 0.9987 | NA | 0.7452 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23076537 | 1.0000 | 0.9948 | 0.9983 | NA | 0.7281 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22272888 | 1.0000 | 0.9934 | 0.9971 | NA | 0.7061 | NA | not assessable | not assessable (no external cohort) |
| meth_cg11283983 | 1.0000 | 0.9941 | 0.9990 | NA | 0.7349 | NA | not assessable | not assessable (no external cohort) |
| meth_cg01308827 | 1.0000 | 0.9904 | 0.9982 | NA | 0.6826 | NA | not assessable | not assessable (no external cohort) |
| meth_cg25373630 | 1.0000 | 0.9958 | 0.9985 | NA | 0.7237 | NA | not assessable | not assessable (no external cohort) |
| meth_cg23907053 | 1.0000 | 0.9939 | 0.9987 | NA | 0.7375 | NA | not assessable | not assessable (no external cohort) |
| meth_cg04505435 | 1.0000 | 0.9984 | 0.9988 | NA | 0.7674 | NA | not assessable | not assessable (no external cohort) |
| meth_cg22212414 | 1.0000 | 0.9947 | 0.9983 | NA | 0.6969 | NA | not assessable | not assessable (no external cohort) |
| meth_cg10201671 | 1.0000 | 0.9966 | 0.9985 | NA | 0.7373 | NA | not assessable | not assessable (no external cohort) |
| meth_cg02976723 | 1.0000 | 0.9958 | 0.9970 | NA | 0.7362 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__gender.demographic_female | 1.0000 | 0.9948 | 0.9971 | NA | 0.8382 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__gender.demographic_male | 1.0000 | 0.9948 | 0.9971 | NA | 0.8382 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__race.demographic_american indian or alaska native | NA | NA | NA | NA | NA | NA | not assessable | not assessable (excluded from Cox model: separation) |
| clin_categorical__race.demographic_asian | 1.0000 | 0.9946 | 0.9980 | NA | 0.8337 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__race.demographic_black or african american | 1.0000 | 0.9876 | 0.9973 | NA | 0.6803 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__race.demographic_not reported | 1.0000 | 0.9979 | 0.9976 | NA | 0.8439 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__race.demographic_white | 1.0000 | 0.9905 | 0.9982 | NA | 0.6906 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__ethnicity.demographic_hispanic or latino | 1.0000 | 0.9945 | 0.9973 | NA | 0.7693 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__ethnicity.demographic_not hispanic or latino | 1.0000 | 0.9935 | 0.9929 | NA | 0.8346 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__ethnicity.demographic_not reported | 1.0000 | 0.9962 | 0.9957 | NA | 0.7265 | NA | not assessable | not assessable (no external cohort) |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage X | 1.0000 | 0.9962 | 0.9969 | NA | 0.8399 | NA | not assessable | not assessable (no external cohort) |

## Tekrar üretim

`python experiments/06_brca_real.py --replicates 50 --seed 20260917`

Bağımlılıklar: paket bağımlılıklarına ek lifelines, pyarrow, joblib, requests. Ham indirme betiği `experiments/06_download_metabric.py`. Argümanlar, sürümler, kaynak ve girdi SHA256: results/phase3_brca/metadata.json. Hasta/sample kimlikleri yalnız ham girdilerde kullanılır; çıktı tablolarında/loglarda yoktur. RSF/DeepSurv bu koşuya dahil değildir.
