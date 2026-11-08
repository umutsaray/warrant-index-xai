# Faz 21 ADIM 2 ön-kayıt — düzeltilmiş v3 boru hattı

**Statü:** sonuç görülmeden önce kaydedilmiş bağlayıcı yöntem planı. Bu dosyanın commit'i tamamlanmadan hiçbir v3 model fit'i, bootstrap, null-kalibrasyon veya sertifikasyon çalıştırılmayacaktır. Bu ön-kayıt v1'in dondurulmuş Tablo 3–9 sonuçlarını değiştirmez.

## 1. Sabit dış bölünme ve iç ayrım

- Dış test: Faz 7'de performans sonucu görülmeden seçilmiş hasta-gruplu fold 0; **168 satır**. Dış split seed'i **42** ve hasta anahtarı TCGA barkodunun ilk 12 karakteridir.
- Dış test, panel seçimi, eksik değer doldurma, ölçekleme, özellik sıralama, model uyarlama ve erken durdurma girdilerinin hiçbirine girmez.
- Dış eğitimin içinde `StratifiedGroupKFold(n_splits=7, shuffle=True, random_state=4201)` kullanılacaktır; stratifikasyon etiketi event, grup ilk-12-karakter hasta kimliğidir. İlk üretilen fold yaklaşık %14.3 iç-validation olarak ayrılır.
- İç-validation'ın hiçbir özellik değeri veya sonlanımı özellik filtreleme/sıralama, imputasyon, ölçekleme ya da herhangi bir fit için kullanılmaz. Tek istisna DeepSurv erken durdurma/checkpoint seçimidir.
- Öğrenilen bütün ön-işleme parametreleri yalnız iç-eğitimde öğrenilir ve iç-validation/dış teste değişmeden uygulanır.
- CoxPH ve RSF'nin nihai referans fit'i dış eğitimin tamamına yapılabilir; ancak panel ve bütün öğrenilmiş ön-işleme parametreleri yukarıdaki iç-eğitimden dondurulur. DeepSurv aynı dondurulmuş parametreleri kullanır ve iç-validation'ı yalnız erken durdurmada görür.

## 2. METABRIC platform beyaz listesi

Kaynak tam matris aynı çalışma/profilden alınmıştır:

- study: `brca_metabric`
- molecular profile: `brca_metabric_mrna`
- profil adı: Illumina HT-12 v3 microarray, continuous log2 intensity
- tam profil örnek sayısı: 1,980; önceden tanımlı OS-uygun ekspresyon kohortu: **1,979**
- tam matris gen satırı: 20,603
- ham yerel dosya: git-ignore'lı `data/external/metabric_phase21/data_mrna_illumina_microarray.txt`; commit edilmeyecektir
- ham dosya SHA-256: `4470069455a4ed38ffed5d12513e468bbf5b9a57d66fa1162e7bd9a64889ab7a`

**METABRIC beyaz listesi**, bu profilde dönen, 1,979 hastanın tamamında eksiksiz ve `SD(ddof=0)>0` olan gen satırlarıdır. METABRIC sonlanımı yalnız sabit 1,979 kohortu tanımlayan önceki uygunluk kontrolü dışında sıralama, seçim veya eşik belirlemede kullanılmaz; beyaz liste platform kapsamı ve ölçüm varlığına dayanır. Beyaz liste **20,592** satırdır. Kanonik temsil, `(Hugo_Symbol, Entrez_Gene_Id)` ile sıralı UTF-8 `SYMBOL<TAB>ENTREZ<LF>` satırlarıdır; SHA-256:

`3a7a633536a580065c1f93e766782b479ea81a113e3326f6969e34bdea9fe64a`

Yerel doğrulama kaydı `data/external/metabric_phase21/full_profile_whitelist.json` altında tutulur ve ham hasta-düzeyi verinin yanında git tarafından izlenmez.

## 3. v3 tarama ve model paneli

Her omik blokta (RNA, miRNA, methylation) yalnız iç-eğitim satırlarıyla aynı kural uygulanır:

1. en az iç-eğitim satırlarının %10'unda nonzero değer;
2. nonzero satırlarda en az 10 olay;
3. tek değişkenli Cox likelihood-ratio (LR) p-değeri; Breslow tie işlemi;
4. artan LR p-değeri, eşitlikte özellik kimliğinin artan leksikografik sırası;
5. blok başına ilk 200 özellik.

Klinik 16 özellik aynı tanım ve sırayla korunur. LR kuralı sonuç görülmeden önce seçilmiştir: hedefli doğrulamada gerçek ortak RNA sayısı LR için 130, skor testi için 100'dür ve LR ilk-200'ünün 113'ü tam filtrelenmiş ailede BH q<0.05'tir. Bu gerekçe skor/Wald alternatiflerine sonradan geçiş izni vermez.

## 4. Denetim ailesi

- RNA denetim ailesi: yalnız iç-eğitimde yukarıdaki filtre+LR kuralıyla sıralanan ve METABRIC beyaz listesinde bulunan genler arasından ilk **60**.
- Klinik denetim ailesi: iki kohortta ortak yaş özelliği ile AJCC evre gösterge özellikleri.
- Gerekçe: RSF/DeepSurv explainer maliyeti ve BH aile yükü; 178 özellik, 27-feature v1 aileye göre yaklaşık 6.6 kat daha ağırdır. Aile boyutu sonuca göre genişletilmeyecek/daraltılmayacaktır.
- CoxPH birincil olarak kendi tüm model girdilerinde değerlendirilir; ayrıca modeller arası eşleştirme için aynı 60+klinik ailede ve aynı 12+12 sorguda raporlanır.

## 5. Modeller ve sabit denetim ayarları

Model mimarisi/hiperparametreleri v1 ile aynıdır: CoxPH `penalizer=0.1` ve dondurulmuş retry; RSF 500 ağaç, `min_samples_split=10`, `min_samples_leaf=5`, `max_features="sqrt"`; DeepSurv aynı blok mimarisi, dropout 0.3, AdamW lr=0.001/wd=0.0001, en çok 300 epoch, patience 30.

- model seed: **42**
- patient-cluster bootstrap: **N=50**, seed **20260917**
- query seçimi: TCGA 12 + METABRIC 12, değiştirmeden ve yerine koymadan, seed **20264917**
- inferans bootstrap: **B=4,999**, seed **20260918**
- null-kalibrasyon: 2,000 dış tekrar × 999 iç bootstrap; `single_active` + `tie`; üç alpha; Holm ile 30 kontrol
- yöntem: `recursive_studentized`
- sabit eşikler: `tau=0.8`, `eta=0.02`, `a=b=0.5`, `u=3`
- WI: dondurulmuş min agregasyonu ve `G_claim=min(G,G_signed)`

RSF sampled inactive background k=16 ve 8/128 yol bütçesini; DeepSurv DeepExplainer k=16/128 background bütçesini kullanır. Null kapısı model başına geçmeden BH/BY sertifikasyonu çalıştırılmaz.

## 6. Önceden belirlenmiş çıktı şemaları

Tüm v3 çıktıları `phase21_reviewer1/v3/` altında yeni dosyalardır; v1 dosyalarının üzerine yazılmaz.

- `split_and_preprocessing.json`: dış/iç hasta hash'leri, satır/hasta/event sayıları, seed'ler, fit-scope ve girdi SHA-256'ları.
- `metabric_whitelist_summary.json`: yukarıdaki kaynak/hash/kapsam sayıları; hasta değerleri yok.
- `feature_screening.csv`: `block,feature,nonzero_n,nonzero_fraction,nonzero_events,lr_p,rank,selected_model_panel,metabric_whitelisted,selected_audit_family`.
- `model_metadata.json`: model ayarları, seed, checkpoint/fit uyarıları, C-index ve artefakt hash'leri.
- model başına `attributions_n50.npz` (yerel, git-ignore), `bootstrap_summary.json`, `feature_report.csv`, `null_gate.csv`, `feature_certification.csv`.
- `feature_report.csv`: `feature,block,F,S,U,G_magnitude,G_signed,G_claim,R,WI,bottleneck,point_adequate`.
- `feature_certification.csv`: bileşen p-değerleri, `p_joint,q_BH,q_BY,certified_BH,certified_BY`; ailedeki bütün özellikler yer alır.
- `matched_three_model.csv`: ortak aile/12+12 sorgu üzerinde üç modelin eşleştirilmiş sonuçları.

Hasta kimlikleri, ham METABRIC değerleri ve hasta-düzeyi tensörler izlenen çıktılara yazılmaz. Her çalıştırmadan önce/sonra `components.py` ve `aggregation.py` SHA-256'ları sırasıyla `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa` ve `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68` olmalıdır.

## 7. Kapılar

Bu plan **hangi sonucun görüleceğine bakmadan** uygulanacaktır. Bir v1 dondurulmuş dosyanın değişmesi, burada tanımlanmamış bir fit/remedy gereksinimi, hasta kesişimi, dış/validation sızıntısı, formül hash değişimi veya null-kapısı ihlali halinde sonraki aşama çalıştırılmaz ve durum açıkça raporlanır. Bu commit'ten sonra kullanıcı onayı olmadan hiçbir v3 model/denetim koşusu başlatılmayacaktır.
