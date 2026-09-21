# Faz 3 — TCGA-BRCA MANIFEST

Kaynak: `../Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/`

## Adım 0: dosyalar ve doğrulama

- `prepare_dataset.py`: mevcut, 16603 byte
- `impute_and_split.py`: mevcut, 15080 byte
- `select_features.py`: mevcut, 9984 byte
- `selected_features.json`: mevcut, 37117 byte
- `X_train_selected.parquet`: mevcut, 1606134 byte
- `X_test_selected.parquet`: mevcut, 526748 byte
- `X_train.parquet`: mevcut, 180258612 byte
- `X_test.parquet`: mevcut, 60278627 byte
- `y_train.parquet`: mevcut, 10032 byte
- `y_test.parquet`: mevcut, 3493 byte
- `coxph_baseline.joblib`: mevcut, 3392129 byte
- `train_baselines.py`: mevcut, 7290 byte
- train: X=(673, 616), y=(673, 2), event=106; beklenen 673/106: EŞLEŞİYOR
- train event/time mevcut; sıra uyumu=True; eksik X hücresi=0; time min=1.0.
- test: X=(169, 616), y=(169, 2), event=26; beklenen 169/26: EŞLEŞİYOR
- test event/time mevcut; sıra uyumu=True; eksik X hücresi=0; time min=1.0.

## SAPMA ve beklenmeyen yapı

- **İstenen 623 feature yerine 616 bulundu.** RNA=200, miRNA=200, methylation=200, clinical=16.
- Ham X_train ara matrisinde clinical=23; select_features.py AJCC alt evrelerini collapse ediyor. Seçilmiş panel ve selected_features.json 616 sütunla tutarlı; eksik 7 sütun uydurulmadı.
- İçte ayrıca başka analiz/kohort dizinleri var; istenen 673/169 ayrımıyla eşleşen kök dosyalar kullanıldı.
- Joblib sözlüğünde model/scaler/features var. Model tipi=CoxPHFitter, penalizer=0.1, feature=614; yükleme uyarı sayısı=0.
- Eğitim betiği: event veya censored grupta varyans<=1e-12 dışlama; StandardScaler train fit; penalizer=0.1, max_steps=50. Kaydedilmiş model kullanılır; bootstrap replikaları aynı ayarlarla yeniden eğitilir.
- Modelden dışlanan feature'lar: `clin_categorical__race.demographic_american indian or alaska native`, `mirna_hsa-mir-8073`.
- Klinik panel: yaş, cinsiyet, ırk, etnisite, Stage I/II/III/IV/X; grade/ER/PR/HER2 yok. Bunlar varmış gibi eklenmez.
- **Hasta çakışması:** train/testte 23 ortak hasta; train 628 benzersiz hasta. Audit'te ortak hastalar çıkarılır, kalan her hastada ilk örnek tutulur. Bootstrap hasta-kümesi bazındadır. Orijinal dosyalar değişmez.

## Tam seçilmiş feature listesi (hasta kimliği içermez)

- `rna_ENSG00000121351.8`
- `rna_ENSG00000229086.5`
- `rna_ENSG00000178522.14`
- `rna_ENSG00000252076.1`
- `rna_ENSG00000201919.1`
- `rna_ENSG00000234863.1`
- `rna_ENSG00000236485.1`
- `rna_ENSG00000213958.2`
- `rna_ENSG00000253940.1`
- `rna_ENSG00000207788.1`
- `rna_ENSG00000285569.2`
- `rna_ENSG00000230161.1`
- `rna_ENSG00000240320.1`
- `rna_ENSG00000286371.1`
- `rna_ENSG00000224250.1`
- `rna_ENSG00000201523.1`
- `rna_ENSG00000278687.1`
- `rna_ENSG00000199220.1`
- `rna_ENSG00000285684.1`
- `rna_ENSG00000249994.1`
- `rna_ENSG00000252978.1`
- `rna_ENSG00000227143.1`
- `rna_ENSG00000006116.4`
- `rna_ENSG00000278672.1`
- `rna_ENSG00000266448.1`
- `rna_ENSG00000224935.2`
- `rna_ENSG00000215368.5`
- `rna_ENSG00000259694.1`
- `rna_ENSG00000080910.13`
- `rna_ENSG00000286720.1`
- `rna_ENSG00000206891.1`
- `rna_ENSG00000212532.1`
- `rna_ENSG00000251264.1`
- `rna_ENSG00000237242.1`
- `rna_ENSG00000261340.1`
- `rna_ENSG00000133742.14`
- `rna_ENSG00000286103.1`
- `rna_ENSG00000113889.14`
- `rna_ENSG00000232395.1`
- `rna_ENSG00000260620.1`
- `rna_ENSG00000238077.1`
- `rna_ENSG00000260639.1`
- `rna_ENSG00000253776.1`
- `rna_ENSG00000270123.4`
- `rna_ENSG00000236656.1`
- `rna_ENSG00000229834.1`
- `rna_ENSG00000258803.1`
- `rna_ENSG00000200434.1`
- `rna_ENSG00000203876.10`
- `rna_ENSG00000253488.1`
- `rna_ENSG00000182591.6`
- `rna_ENSG00000260413.1`
- `rna_ENSG00000274244.1`
- `rna_ENSG00000243995.3`
- `rna_ENSG00000232728.1`
- `rna_ENSG00000201752.1`
- `rna_ENSG00000285990.1`
- `rna_ENSG00000241777.1`
- `rna_ENSG00000250923.1`
- `rna_ENSG00000213607.4`
- `rna_ENSG00000255618.2`
- `rna_ENSG00000254751.3`
- `rna_ENSG00000256484.1`
- `rna_ENSG00000249479.1`
- `rna_ENSG00000226695.1`
- `rna_ENSG00000224922.1`
- `rna_ENSG00000252891.1`
- `rna_ENSG00000169800.14`
- `rna_ENSG00000250391.1`
- `rna_ENSG00000258115.1`
- `rna_ENSG00000285831.1`
- `rna_ENSG00000198670.12`
- `rna_ENSG00000286252.1`
- `rna_ENSG00000251338.1`
- `rna_ENSG00000276118.1`
- `rna_ENSG00000207097.1`
- `rna_ENSG00000212383.1`
- `rna_ENSG00000279376.1`
- `rna_ENSG00000267125.2`
- `rna_ENSG00000226374.3`
- `rna_ENSG00000090402.8`
- `rna_ENSG00000167654.18`
- `rna_ENSG00000201412.1`
- `rna_ENSG00000285598.1`
- `rna_ENSG00000267603.1`
- `rna_ENSG00000233452.7`
- `rna_ENSG00000139874.6`
- `rna_ENSG00000248328.1`
- `rna_ENSG00000206871.1`
- `rna_ENSG00000261103.3`
- `rna_ENSG00000198099.9`
- `rna_ENSG00000169760.17`
- `rna_ENSG00000272855.1`
- `rna_ENSG00000225399.4`
- `rna_ENSG00000231426.7`
- `rna_ENSG00000231861.2`
- `rna_ENSG00000253346.1`
- `rna_ENSG00000239438.1`
- `rna_ENSG00000270069.1`
- `rna_ENSG00000248869.7`
- `rna_ENSG00000237560.7`
- `rna_ENSG00000174899.12`
- `rna_ENSG00000254403.1`
- `rna_ENSG00000176134.5`
- `rna_ENSG00000100911.16`
- `rna_ENSG00000113389.16`
- `rna_ENSG00000169432.18`
- `rna_ENSG00000157617.17`
- `rna_ENSG00000158122.12`
- `rna_ENSG00000232821.1`
- `rna_ENSG00000168917.9`
- `rna_ENSG00000270547.6`
- `rna_ENSG00000133800.9`
- `rna_ENSG00000077616.11`
- `rna_ENSG00000226674.11`
- `rna_ENSG00000233884.2`
- `rna_ENSG00000252969.1`
- `rna_ENSG00000215474.8`
- `rna_ENSG00000229237.2`
- `rna_ENSG00000201809.1`
- `rna_ENSG00000255162.1`
- `rna_ENSG00000286892.1`
- `rna_ENSG00000064042.18`
- `rna_ENSG00000138795.10`
- `rna_ENSG00000100908.14`
- `rna_ENSG00000147223.6`
- `rna_ENSG00000234723.1`
- `rna_ENSG00000266964.6`
- `rna_ENSG00000109158.11`
- `rna_ENSG00000060709.15`
- `rna_ENSG00000287616.1`
- `rna_ENSG00000279204.1`
- `rna_ENSG00000253177.2`
- `rna_ENSG00000004846.16`
- `rna_ENSG00000261472.1`
- `rna_ENSG00000023892.11`
- `rna_ENSG00000214575.9`
- `rna_ENSG00000231925.12`
- `rna_ENSG00000253105.6`
- `rna_ENSG00000183747.12`
- `rna_ENSG00000226622.6`
- `rna_ENSG00000172901.20`
- `rna_ENSG00000225131.2`
- `rna_ENSG00000118407.15`
- `rna_ENSG00000187531.14`
- `rna_ENSG00000108639.8`
- `rna_ENSG00000237870.6`
- `rna_ENSG00000230147.2`
- `rna_ENSG00000116663.11`
- `rna_ENSG00000124343.14`
- `rna_ENSG00000076043.11`
- `rna_ENSG00000134531.10`
- `rna_ENSG00000182173.13`
- `rna_ENSG00000267583.5`
- `rna_ENSG00000226592.1`
- `rna_ENSG00000172725.14`
- `rna_ENSG00000225056.1`
- `rna_ENSG00000240108.1`
- `rna_ENSG00000163632.14`
- `rna_ENSG00000233648.1`
- `rna_ENSG00000271776.1`
- `rna_ENSG00000187758.8`
- `rna_ENSG00000092010.15`
- `rna_ENSG00000104327.7`
- `rna_ENSG00000224365.1`
- `rna_ENSG00000204246.3`
- `rna_ENSG00000064961.19`
- `rna_ENSG00000229498.1`
- `rna_ENSG00000124440.16`
- `rna_ENSG00000205186.4`
- `rna_ENSG00000236466.2`
- `rna_ENSG00000286750.1`
- `rna_ENSG00000163145.13`
- `rna_ENSG00000181788.4`
- `rna_ENSG00000269124.1`
- `rna_ENSG00000246430.7`
- `rna_ENSG00000224609.8`
- `rna_ENSG00000287179.1`
- `rna_ENSG00000237202.1`
- `rna_ENSG00000227748.2`
- `rna_ENSG00000229950.1`
- `rna_ENSG00000249102.3`
- `rna_ENSG00000250993.1`
- `rna_ENSG00000254135.3`
- `rna_ENSG00000072657.9`
- `rna_ENSG00000121871.4`
- `rna_ENSG00000279499.1`
- `rna_ENSG00000280304.1`
- `rna_ENSG00000287785.1`
- `rna_ENSG00000169764.16`
- `rna_ENSG00000237417.1`
- `rna_ENSG00000140332.17`
- `rna_ENSG00000231672.8`
- `rna_ENSG00000163545.11`
- `rna_ENSG00000187242.6`
- `rna_ENSG00000250334.6`
- `rna_ENSG00000187486.5`
- `rna_ENSG00000095203.14`
- `rna_ENSG00000149557.14`
- `rna_ENSG00000122025.15`
- `mirna_hsa-mir-4533`
- `mirna_hsa-mir-548ac`
- `mirna_hsa-mir-4760`
- `mirna_hsa-mir-4493`
- `mirna_hsa-mir-640`
- `mirna_hsa-mir-4788`
- `mirna_hsa-mir-4776-1`
- `mirna_hsa-mir-22`
- `mirna_hsa-mir-6715a`
- `mirna_hsa-mir-452`
- `mirna_hsa-mir-891a`
- `mirna_hsa-mir-129-2`
- `mirna_hsa-mir-190b`
- `mirna_hsa-mir-874`
- `mirna_hsa-mir-224`
- `mirna_hsa-mir-219a-2`
- `mirna_hsa-mir-342`
- `mirna_hsa-mir-6883`
- `mirna_hsa-mir-133a-2`
- `mirna_hsa-mir-23c`
- `mirna_hsa-mir-618`
- `mirna_hsa-mir-1468`
- `mirna_hsa-mir-490`
- `mirna_hsa-mir-378a`
- `mirna_hsa-mir-486-1`
- `mirna_hsa-mir-141`
- `mirna_hsa-mir-486-2`
- `mirna_hsa-mir-6742`
- `mirna_hsa-mir-192`
- `mirna_hsa-mir-891b`
- `mirna_hsa-mir-221`
- `mirna_hsa-mir-193a`
- `mirna_hsa-mir-133a-1`
- `mirna_hsa-mir-200c`
- `mirna_hsa-mir-5094`
- `mirna_hsa-mir-378c`
- `mirna_hsa-mir-4435-2`
- `mirna_hsa-mir-4724`
- `mirna_hsa-mir-466`
- `mirna_hsa-mir-129-1`
- `mirna_hsa-mir-33b`
- `mirna_hsa-mir-3199-2`
- `mirna_hsa-mir-4458`
- `mirna_hsa-mir-4751`
- `mirna_hsa-mir-16-1`
- `mirna_hsa-mir-135a-2`
- `mirna_hsa-mir-5591`
- `mirna_hsa-mir-375`
- `mirna_hsa-mir-96`
- `mirna_hsa-mir-365b`
- `mirna_hsa-mir-664b`
- `mirna_hsa-mir-365a`
- `mirna_hsa-mir-4738`
- `mirna_hsa-mir-139`
- `mirna_hsa-mir-378d-2`
- `mirna_hsa-mir-16-2`
- `mirna_hsa-mir-2355`
- `mirna_hsa-mir-6798`
- `mirna_hsa-mir-4678`
- `mirna_hsa-mir-549a`
- `mirna_hsa-mir-222`
- `mirna_hsa-mir-4469`
- `mirna_hsa-mir-182`
- `mirna_hsa-mir-424`
- `mirna_hsa-mir-636`
- `mirna_hsa-mir-548at`
- `mirna_hsa-mir-3170`
- `mirna_hsa-mir-605`
- `mirna_hsa-mir-556`
- `mirna_hsa-mir-4732`
- `mirna_hsa-mir-892a`
- `mirna_hsa-mir-4528`
- `mirna_hsa-let-7b`
- `mirna_hsa-mir-4771-2`
- `mirna_hsa-mir-6874`
- `mirna_hsa-mir-431`
- `mirna_hsa-mir-191`
- `mirna_hsa-mir-1293`
- `mirna_hsa-mir-548o`
- `mirna_hsa-mir-3683`
- `mirna_hsa-mir-5683`
- `mirna_hsa-mir-4684`
- `mirna_hsa-mir-6780a`
- `mirna_hsa-mir-21`
- `mirna_hsa-mir-5682`
- `mirna_hsa-mir-585`
- `mirna_hsa-mir-5694`
- `mirna_hsa-mir-3688-1`
- `mirna_hsa-mir-4655`
- `mirna_hsa-mir-6796`
- `mirna_hsa-mir-28`
- `mirna_hsa-mir-4453`
- `mirna_hsa-mir-24-2`
- `mirna_hsa-mir-5191`
- `mirna_hsa-mir-147b`
- `mirna_hsa-mir-3907`
- `mirna_hsa-mir-144`
- `mirna_hsa-mir-24-1`
- `mirna_hsa-mir-4803`
- `mirna_hsa-mir-5691`
- `mirna_hsa-mir-548f-1`
- `mirna_hsa-mir-200a`
- `mirna_hsa-mir-1245b`
- `mirna_hsa-mir-3688-2`
- `mirna_hsa-mir-4532`
- `mirna_hsa-mir-5588`
- `mirna_hsa-mir-1225`
- `mirna_hsa-mir-335`
- `mirna_hsa-mir-4658`
- `mirna_hsa-mir-3661`
- `mirna_hsa-mir-3622a`
- `mirna_hsa-mir-146a`
- `mirna_hsa-mir-5100`
- `mirna_hsa-mir-662`
- `mirna_hsa-mir-4796`
- `mirna_hsa-mir-135a-1`
- `mirna_hsa-mir-574`
- `mirna_hsa-mir-3162`
- `mirna_hsa-mir-1258`
- `mirna_hsa-mir-145`
- `mirna_hsa-mir-6851`
- `mirna_hsa-mir-659`
- `mirna_hsa-mir-4670`
- `mirna_hsa-mir-5705`
- `mirna_hsa-mir-8073`
- `mirna_hsa-mir-6716`
- `mirna_hsa-mir-4714`
- `mirna_hsa-mir-5192`
- `mirna_hsa-mir-217`
- `mirna_hsa-mir-20b`
- `mirna_hsa-mir-6499`
- `mirna_hsa-mir-30a`
- `mirna_hsa-mir-142`
- `mirna_hsa-mir-429`
- `mirna_hsa-mir-6717`
- `mirna_hsa-mir-1257`
- `mirna_hsa-mir-6746`
- `mirna_hsa-mir-106a`
- `mirna_hsa-mir-6759`
- `mirna_hsa-mir-548f-4`
- `mirna_hsa-mir-6892`
- `mirna_hsa-mir-133b`
- `mirna_hsa-mir-4720`
- `mirna_hsa-mir-200b`
- `mirna_hsa-mir-3139`
- `mirna_hsa-mir-4520-2`
- `mirna_hsa-let-7i`
- `mirna_hsa-mir-148a`
- `mirna_hsa-mir-125b-1`
- `mirna_hsa-mir-541`
- `mirna_hsa-mir-6715b`
- `mirna_hsa-mir-320b-2`
- `mirna_hsa-mir-6830`
- `mirna_hsa-mir-922`
- `mirna_hsa-mir-600`
- `mirna_hsa-mir-125b-2`
- `mirna_hsa-mir-4668`
- `mirna_hsa-mir-1915`
- `mirna_hsa-mir-1-2`
- `mirna_hsa-mir-4450`
- `mirna_hsa-mir-5189`
- `mirna_hsa-mir-3619`
- `mirna_hsa-mir-1246`
- `mirna_hsa-mir-4716`
- `mirna_hsa-mir-143`
- `mirna_hsa-mir-4502`
- `mirna_hsa-mir-2114`
- `mirna_hsa-mir-4510`
- `mirna_hsa-mir-4641`
- `mirna_hsa-mir-6087`
- `mirna_hsa-mir-5693`
- `mirna_hsa-mir-548ai`
- `mirna_hsa-mir-194-2`
- `mirna_hsa-mir-326`
- `mirna_hsa-mir-6885`
- `mirna_hsa-mir-4435-1`
- `mirna_hsa-mir-6756`
- `mirna_hsa-mir-140`
- `mirna_hsa-mir-584`
- `mirna_hsa-mir-586`
- `mirna_hsa-mir-6511b-2`
- `mirna_hsa-mir-3115`
- `mirna_hsa-mir-6825`
- `mirna_hsa-mir-5186`
- `mirna_hsa-mir-183`
- `mirna_hsa-mir-5187`
- `mirna_hsa-mir-4452`
- `mirna_hsa-mir-451a`
- `mirna_hsa-mir-550a-2`
- `mirna_hsa-mir-6852`
- `mirna_hsa-mir-32`
- `mirna_hsa-mir-592`
- `mirna_hsa-mir-4444-2`
- `mirna_hsa-mir-6800`
- `mirna_hsa-mir-31`
- `mirna_hsa-mir-4681`
- `mirna_hsa-mir-181d`
- `mirna_hsa-mir-4770`
- `mirna_hsa-mir-3936`
- `mirna_hsa-mir-1294`
- `meth_cg23972735`
- `meth_cg09856092`
- `meth_cg09229620`
- `meth_cg03854238`
- `meth_cg01375382`
- `meth_cg06738031`
- `meth_cg13420364`
- `meth_cg02467858`
- `meth_cg12134633`
- `meth_cg18618429`
- `meth_cg21616785`
- `meth_cg07205203`
- `meth_cg22077197`
- `meth_cg09737197`
- `meth_cg03532263`
- `meth_cg15037823`
- `meth_cg14669379`
- `meth_cg01444716`
- `meth_cg10665891`
- `meth_cg00778807`
- `meth_cg04723534`
- `meth_cg26389380`
- `meth_cg16962463`
- `meth_cg16672562`
- `meth_cg09136878`
- `meth_cg01973456`
- `meth_cg15429127`
- `meth_cg26607031`
- `meth_cg21729992`
- `meth_cg07904452`
- `meth_cg19418318`
- `meth_cg07109453`
- `meth_cg17648080`
- `meth_cg27598107`
- `meth_cg18519308`
- `meth_cg19880942`
- `meth_cg14981189`
- `meth_cg23497383`
- `meth_cg00887153`
- `meth_cg00382999`
- `meth_cg10296205`
- `meth_cg27063372`
- `meth_cg27494647`
- `meth_cg23824801`
- `meth_cg04178787`
- `meth_cg14769786`
- `meth_cg22891070`
- `meth_cg26010218`
- `meth_cg06600725`
- `meth_cg08047907`
- `meth_cg08436419`
- `meth_cg25181170`
- `meth_cg18142262`
- `meth_cg24504442`
- `meth_cg00892742`
- `meth_cg27408285`
- `meth_cg11157127`
- `meth_cg22675660`
- `meth_cg19412295`
- `meth_cg24181266`
- `meth_cg08923594`
- `meth_cg04831594`
- `meth_cg24492202`
- `meth_cg17239057`
- `meth_cg16610086`
- `meth_cg08511651`
- `meth_cg10932018`
- `meth_cg03359666`
- `meth_cg17165836`
- `meth_cg04682791`
- `meth_cg19126300`
- `meth_cg12580445`
- `meth_cg18276112`
- `meth_cg07565505`
- `meth_cg18456523`
- `meth_cg24797187`
- `meth_cg24865779`
- `meth_cg00680551`
- `meth_cg23557926`
- `meth_cg13804575`
- `meth_cg22993195`
- `meth_cg02983451`
- `meth_cg22081371`
- `meth_cg13530938`
- `meth_cg00716257`
- `meth_cg22254463`
- `meth_cg05130679`
- `meth_cg15410675`
- `meth_cg09439920`
- `meth_cg17356733`
- `meth_cg04682802`
- `meth_cg11564601`
- `meth_cg10572943`
- `meth_cg09866743`
- `meth_cg00857907`
- `meth_cg27614178`
- `meth_cg18173184`
- `meth_cg26353296`
- `meth_cg20797905`
- `meth_cg07254421`
- `meth_cg07965640`
- `meth_cg26354493`
- `meth_cg26348348`
- `meth_cg18208842`
- `meth_cg17415451`
- `meth_cg08390172`
- `meth_cg27628340`
- `meth_cg20780180`
- `meth_cg20391567`
- `meth_cg18450254`
- `meth_cg10928544`
- `meth_cg20312012`
- `meth_cg04942251`
- `meth_cg07937427`
- `meth_cg00314029`
- `meth_cg02085507`
- `meth_cg18121066`
- `meth_cg16493071`
- `meth_cg00800229`
- `meth_cg12114888`
- `meth_cg23198707`
- `meth_cg01623261`
- `meth_cg22728904`
- `meth_cg24162579`
- `meth_cg08356637`
- `meth_cg19516279`
- `meth_cg01113811`
- `meth_cg00589541`
- `meth_cg07100000`
- `meth_cg24810735`
- `meth_cg08481075`
- `meth_cg25964728`
- `meth_cg05890377`
- `meth_cg25932869`
- `meth_cg24432768`
- `meth_cg03930313`
- `meth_cg07570723`
- `meth_cg02783887`
- `meth_cg02863179`
- `meth_cg22753668`
- `meth_cg06133110`
- `meth_cg02013841`
- `meth_cg04850286`
- `meth_cg02058870`
- `meth_cg23004031`
- `meth_cg08613144`
- `meth_cg18261069`
- `meth_cg10321156`
- `meth_cg00322927`
- `meth_cg24156393`
- `meth_cg07065111`
- `meth_cg14377594`
- `meth_cg25976440`
- `meth_cg00970396`
- `meth_cg21143896`
- `meth_cg03478739`
- `meth_cg14429919`
- `meth_cg16434331`
- `meth_cg18420599`
- `meth_cg00081417`
- `meth_cg11228682`
- `meth_cg13045134`
- `meth_cg09977969`
- `meth_cg00744866`
- `meth_cg19021076`
- `meth_cg15626285`
- `meth_cg07863022`
- `meth_cg12905085`
- `meth_cg00875805`
- `meth_cg20378070`
- `meth_cg04085768`
- `meth_cg11686819`
- `meth_cg16562730`
- `meth_cg15716680`
- `meth_cg00768487`
- `meth_cg21730858`
- `meth_cg14826711`
- `meth_cg10558233`
- `meth_cg12444411`
- `meth_cg09825309`
- `meth_cg24453118`
- `meth_cg03159396`
- `meth_cg11574975`
- `meth_cg23285459`
- `meth_cg22193912`
- `meth_cg17627629`
- `meth_cg15925365`
- `meth_cg15200418`
- `meth_cg04415176`
- `meth_cg25296374`
- `meth_cg23076537`
- `meth_cg22272888`
- `meth_cg11283983`
- `meth_cg01308827`
- `meth_cg25373630`
- `meth_cg23907053`
- `meth_cg04505435`
- `meth_cg22212414`
- `meth_cg10201671`
- `meth_cg02976723`
- `clin_numeric__age_at_index.demographic`
- `clin_categorical__gender.demographic_female`
- `clin_categorical__gender.demographic_male`
- `clin_categorical__race.demographic_american indian or alaska native`
- `clin_categorical__race.demographic_asian`
- `clin_categorical__race.demographic_black or african american`
- `clin_categorical__race.demographic_not reported`
- `clin_categorical__race.demographic_white`
- `clin_categorical__ethnicity.demographic_hispanic or latino`
- `clin_categorical__ethnicity.demographic_not hispanic or latino`
- `clin_categorical__ethnicity.demographic_not reported`
- `clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I`
- `clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II`
- `clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III`
- `clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV`
- `clin_categorical__ajcc_pathologic_stage.diagnoses_Stage X`

## Ek kohort ve yakınsama notları

- Kaynak primer-tümör-only değildir. Train örnek türü sayıları: 01=609, 11=60, 06=4; test: 01=153, 11=15, 06=1. Bunlar agregat sayılardır; hasta/sample ID yazılmaz.
- Bootstrap altıncı replikada meth_cg04831594 separation uyarısı verdi. Aynı örnekleme ve penalizer=0.1 ile uyarılı sütun dışlanır; kalan sayısal yakınsamama için solver üst sınırı 50'den 200'e, step_size 0.25'e geçirilir. Model hedefi/toleranslar değiştirilmez; fit denemeleri bootstrap_summary.json içinde tutulur.

## Faz 4 RSF pilot varlıkları

- Kaynak model: `../Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/rsf_baseline.joblib`.
- SHA256: `3ab8da5972812ec3411430419ec40862306316504b270a77a9b1a9ce96121574`.
- Model mevcut haliyle yeniden kullanıldı; ana RSF yeniden eğitilmedi.
- Panel/split uyumu: 616 model girdisi; train 673, test 169; hazırlanmış sütun sırası birebir eşleşiyor.
- Hyperparametreler: `n_estimators=500`, `min_samples_split=10`, `min_samples_leaf=5`, `max_features=sqrt`, `random_state=42`, kaynak eğitimde `n_jobs=-1`.
- Kaynak sonuç: test C-index 0.714795008912656; `baseline_results.json` ve `train_baselines.py` ile kayıtlı.
- Faz 4 N=2 pilot replikaları aynı model ayarları ve patient-cluster bootstrap ile ayrıca fit edilir; ana model artifact'ını değiştirmez.
- TreeSHAP uyumluluk testi: SHAP 0.52.0 / scikit-survival 0.28.0, `RandomSurvivalForest` için destek vermedi (`InvalidModelError`). Pilot model-agnostik permutation Shapley kullanır.
- Ortak dış destek 72 RNA + 5 klinik feature'dır. Pilot açıklama oyunu bu 77 feature'ı değiştirir; diğer 539 RSF girdisini TCGA train ortalamasında sabitler. Bu, tam 616-feature Shapley oyunu değildir.
- Faz 4.1 sağlamlık kontrolü sekiz feature, aynı N=2 ve 12+12 query ile yapıldı. On altı benzersiz TCGA train hastası üzerinden background-integrated sonuçta maksimum |ΔF|=0.003897, maksimum |ΔG_claim|=0.250446 bulundu; sabit background sağlam kabul edilmedi.
- Referans 128→256 kontrolünde maksimum |ΔF|=0.001367; candidate sabit olduğu için G_claim değişmedi. N=50 önerisi sampled background k=16 ve 8/128 bütçedir; tam koşu başlatılmadı.

## Faz 4.2 RSF N=20 sonucu

- Onaylanan tasarım: N=20 patient-cluster bootstrap, sampled background k=16, candidate/reference 8/128 ve aynı 12 TCGA + 12 METABRIC query.
- Ana `rsf_baseline.joblib` yeniden eğitilmedi; 20 bootstrap RSF modeli ayrıca fit edildi. Hiçbir fit uyarı üretmedi, replika/site atlanmadı.
- Bootstrap tahminlerinde `low_memory=True` yalnız çıktı saklama optimizasyonu olarak kullanıldı. Aynı seed/veride 500 ağacın yapısı birebir aynı; standart ve low-memory risk skorları maksimum `7.105427357601002e-15` farklıydı.
- Sonuç: 77/77 feature değerlendirildi, 23 feature WI>0.8 geçti. Bottleneck dağılımı R=72, G_claim=5.
- Aynı sampled-background tasarımında N=2 35/77, N=20 23/77 geçti; N=2 geçenlerin 22'si kaldı, 13'ü düştü, 1 yeni feature geçti.
- Ölçülen birikimli eğitim+açıklama süresi 2949.3 sn (49.2 dk). İlk yavaş site checkpoint'inden sonra low-memory completion invocation 2270.8 sn (37.8 dk) sürdü.
- Formül modüllerinin Faz 2.1 SHA256 değerleri değişmedi. N=50 ve tam 141+1466 query koşusu yapılmadı.

## Faz 5 sertifikasyon tasarımı ve pilotu

- `experiments/08_certification_pilot.py`, CoxPH N=50 ve RSF N=20 saklı replika artefaktlarından 5+5 feature için null-merkezli basic-bootstrap p-değeri pilotunu üretir; model fit veya explainer çağrısı yapmaz.
- `results/phase5_certification_pilot/inventory.json`, hangi replika-düzeyi girdilerin saklı veya deterministik yeniden kurulabilir olduğunu ve çıkarımın sabit query setlerine koşullu olduğunu kaydeder.
- `results/phase5_certification_pilot/pilot_pvalues.csv`, yalnız 10 feature'lık kalibrasyon kontrolüdür. BH/BY ve tam feature ailesi uygulaması kullanıcı onayına kadar başlatılmamıştır.
- `PHASE5_DESIGN_NOTE.md`, yöntem, sayısal pilot, sınırlamalar, formül SHA256 ve gizlilik/wall-clock kontrollerini raporlar.

## Faz 5.1 null kalibrasyon kapısı

- `experiments/09_null_calibration.py`, gerçek model/explainer çağırmadan analitik olarak theta=0.8 olan `single_active` ve `tie` null senaryolarını CoxPH N=50 ve RSF N=20 gürültü/bağımlılık yapılarıyla simüle eder.
- `results/phase51_null_calibration/type1_rates.csv`, 2 model × 2 senaryo × 5 bileşen × 3 alpha için 2000 dış simülasyon ve B=999 iç bootstrap ret oranlarını içerir. RSF satırları sabit 12+12 query ve saklı 8/128 attribution ortalamalarına koşullu olduğunu açıkça taşır.
- Kalibrasyon kapısı KALDI: RSF N=20 `single_active` senaryosunda F/U alpha=0.01 ve G_claim alpha=0.01/0.05 hücreleri Holm-düzeltilmiş binomial kontrolde anti-conservative bulundu. Tam BH/BY başlatılmadı.
- `PHASE5_1_NULL_CALIBRATION_REPORT.md`, kabul kararı, güven aralıkları, olası fakat uygulanmamış düzeltmeleri, SHA256, gizlilik ve wall-clock kontrollerini raporlar.

## Faz 5.2A CoxPH tam-aile sertifikasyonu

- `experiments/10_cox_certification.py`, Faz 3.2'de değerlendirilebilir bulunan 77 CoxPH feature'ın tamamında onaylanan B=4999 null-merkezli bootstrap bileşen p-değerlerini üretir; model veya explainer çağırmaz.
- `results/phase52_cox_certification/feature_certification.csv`, beş bileşen p-değeri, `p_joint=max_k p_jk`, BH/BY q-değerleri, α=0.05 kararları ve `WI>=0.8` nokta-tahmin karşılaştırmasını içerir.
- `PHASE5_2_COX_CERTIFICATION_REPORT.md`, m=77 tam ailede 22 nokta-yeterli feature'a karşı 8 BH ve 8 BY sertifikalı feature bulunduğunu raporlar. RSF sertifikasyonu bu koşuya dahil değildir.

## Faz 5.2B RSF remedy tasarımı ve kalibrasyonu

- `experiments/11_rsf_remedy_calibration.py`, recursive basic-IUT, F/U/G yapraklarında studentized recursive-IUT ve doğrudan m-out-of-n root testini aynı 2000×999 null bataryasında karşılaştırır; gerçek model veya explainer çağırmaz.
- `results/phase52_rsf_remedy/candidate_type1_rates.csv`, aday başına 30 RSF kontrolünü ve Holm-düzeltilmiş gate kararlarını içerir. Recursive basic ve m-out-of-n gate'i geçmedi; recursive studentized aday 30/30 kontrolde anti-conservative ihlal göstermedi.
- `results/phase52_rsf_remedy/actual_5feature_pilot.csv`, geçen yöntemin beş gerçek RSF feature üzerindeki pilotudur. Bütün satırlar sabit 12+12 query ve saklı 8/128 attribution ortalamalarına koşullu olduğunu taşır.
- `PHASE5_2_RSF_REMEDY_REPORT.md`, aday karşılaştırmasını ve öneriyi raporlar. RSF tam-aile ADIM 2 çalıştırılmadı; kullanıcı onayı beklenir.

## Faz 5.3 RSF sertifikasyonu ve CoxPH tutarlılığı

- `experiments/12_rsf_certification.py`, onaylanan recursive_studentized yöntemi RSF'nin 77 feature tam ailesine uygular. `results/phase53_rsf_certification/feature_certification.csv` içinde 23 nokta-yeterli feature'a karşı 15 BH ve 12 BY sertifikalı feature vardır.
- `PHASE5_3_RSF_CERTIFICATION_REPORT.md` ve bütün RSF tablo/metadata çıktıları sabit 12+12 query ile saklı 8/128 attribution ortalamalarına koşulluluk uyarısını taşır.
- `experiments/13_cox_recursive_certification.py`, CoxPH'nin aynı 77 feature ailesini recursive_studentized ile çalıştırır ve Faz 5.2A basic-root sonucuyla eşler. İki yöntem aynı 8 feature'ı BH ve BY ile sertifikalamıştır.
- `results/phase53_cox_recursive/component_p_feature_differences.csv`, yalnız üç G_claim p-değerinin daha muhafazakârlaştığını; sertifika kümesinin değişmediğini gösterir.
- `PHASE5_3_REPORT.md`, iki kolun ortak özetini ve recursive_studentized birincil / CoxPH basic-root robustness önerisini içerir.

## Faz 5.4 BH hayatta-kalma farkı tanısı

- `experiments/14_bh_survival_diagnostic.py`, yalnız mevcut Faz 5.3 CoxPH ve RSF sertifikasyon CSV'lerini okur; yeni bootstrap, model fit'i veya explainer çağrısı yapmaz.
- `results/phase54_bh_survival_diagnostic/point_pass_feature_diagnostics.csv`, nokta-geçen 22 CoxPH ve 23 RSF feature için WI/darboğaz marjinlerini, nokta ve inferans darboğazlarını, p/q değerlerini içerir.
- `margin_band_summary.csv` ve `p_vs_margin.png`, aynı marjin bantlarındaki p-değeri ilişkisini karşılaştırır. 0–0.03 bandında ortalama darboğaz p-değeri CoxPH 0.2644, RSF 0.5237'dir; iki modelde de bu banttan BH sertifikası çıkmamıştır.
- `PHASE5_4_BH_SURVIVAL_DIAGNOSTIC.md`, %36/%65 farkının ana nedenini marjin dağılımı olarak tanılar: marjini 0.03'ten büyük feature sayısı CoxPH'de 9, RSF'de 16'dır. Her iki modelin nokta-geçip BH'yi kaçıranlarında nokta darboğazı R'dir.
- Formül ve kaynak sertifikasyon dosyalarının SHA256 değerleri `metadata.json` içinde doğrulanır; kimlik-sızıntısı kontrolü geçer.

## Faz 5.5 R sıralama-evreni boyutu tanısı

- `experiments/15_rank_universe_diagnostic.py`, saklı CoxPH N=50 katsayıları ve sabit 141 satırlık TCGA audit kümesinden R'yi tanı amaçlı 77-feature evreninde yeniden hesaplar. Yeni bootstrap, model fit/tahmin, explainer veya BH/BY çalıştırmaz.
- `results/phase55_rank_universe_diagnostic/all_77_feature_comparison.csv`, özgün `R_614` ile tanısal `R_77` değerlerini ve gerçek RSF R/WI değerlerini aynı 77 feature üzerinde eşler. Resmî sertifikasyon tablolarına yazılmaz.
- `cox_original_22_comparison.csv`, özgün CoxPH nokta-geçen 22 feature için istenen ayrıntılı karşılaştırmadır. Ortalama R değişimi +0.00283'tür; `WI−0.8>0.03` sayısı 9 olarak kalır (RSF: 16).
- Hipotez desteklenmedi: aynı 22 feature'ın R marjinlerinde RSF'ye Wasserstein uzaklığı 0.05857'den 0.06076'ya çıktı; tüm 77 feature'da da 0.02174'ten 0.02381'e çıktı. Karşı-olgusal CoxPH nokta-geçen sayısı 22'den 23'e yükselse de marjin dağılımı RSF'ye yakınsamadı.
- `PHASE5_5_RANK_UNIVERSE_DIAGNOSTIC.md`, tablo, dağılım özeti, grafik, kapsam sınırları, SHA256 ve wall-clock doğrulamasını içerir.

## Faz 5.6 RSF N=50 ve FDR sertifikasyonu

- `experiments/16_rsf_n50.py`, Faz 4.2'nin sampled-background k=16, 8/128 bütçe ve sabit 12+12 query tasarımında RSF dizisini N=50'ye çıkarır. Doğrulanmış ilk 20 attribution korunur; aynı seed akışındaki 30 yeni bootstrap modeli fit edilir. Ana `rsf_baseline.joblib` değişmez.
- `results/phase56_rsf_n50/feature_report.csv`, 77/77 feature için N=50 birincil F/S/U/G_claim/R/WI sonuçlarını ve korunan N=20 kolonlarını içerir. 20 feature WI≥0.8 geçer; bottleneck dağılımı R=72, G_claim=5'tir.
- `bootstrap_summary.json`, 50 replikayı ve iki site kaydını eksiksiz içerir; fit uyarısı veya atlanan/değiştirilen replika yoktur. Yeni 30 replika 987.9 sn, N=50 toplam birikimli hesap 3937.2 sn sürmüştür.
- `experiments/17_rsf_n50_null_gate.py`, `recursive_studentized` yöntemini N=50 nuisance yapısıyla aynı 2000×999, iki senaryo, üç alpha ve 30 Holm kontrolünde yeniden kalibre eder. Gate 0/30 ihlalle geçmiştir.
- `experiments/18_rsf_n50_certification.py`, geçen yöntemle 77-feature tam ailede BH/BY uygular. N=50 sonucu: 20 nokta-geçen, 10 BH ve 10 BY sertifikalı feature.
- `experiments/19_phase56_comparison.py` ve `PHASE5_6_REPORT.md`, N=20→50 olgunlaşma hipotezini kısmen destekler: eski 23 kazananın 19'unda R düşer, fakat tüm 77 panelde ortalama değişim -0.0006 ve Wilcoxon p=0.743'tür. Yüksek-marjin payı %65.2'den %45.0'a iner, CoxPH'nin %9.1 düzeyine tam yakınsamaz.
- Bütün RSF inferans tabloları sabit 12+12 query ve saklı 8/128 attribution ortalamalarına koşulluluk uyarısını taşır; tekil permütasyon-yolu belirsizliği dahil değildir.

## Faz 6 DeepSurv keşif ve N=2 tasarım pilotu

- Kardeş projedeki hazır `explainable_dl_model.pt` doğrulandı: aynı 673/169 split, aynı 616-feature paneli, test C-index 0.745692 ve dört bloklu DeepSurv mimarisi kullanılır. Ana model yeniden eğitilmedi.
- `experiments/20_deepsurv_pilot.py`, gerçek checkpoint üzerinde DeepExplainer ve GradientExplainer desteğini test eder. İkisi de desteklenir; DeepExplainer additivity hatası 1.77e-07'dir. GradientExplainer, 8/128 `nsamples` bütçesini doğrudan temsil ettiği için approximate-F pilotunda kullanılır.
- `results/phase6_deepsurv_pilot/feature_report.csv`, N=2 ve sabit 12+12 query koşulunda 77/77 feature'ı raporlar: 15 feature WI≥0.8, bottleneck dağılımı R=75 ve G_claim=2'dir.
- Sekiz feature'lık sampled-background kontrolünde k=8→16 maksimum |ΔG_claim| 0.003799; reference 128→256 maksimum |ΔF| 0.000098 bulundu. Sampled k=16 ve reference=128 sağlık kontrolleri kararlıdır.
- Candidate=8 GradientExplainer için maksimum additivity hatası 0.439 risk-skoru birimine çıktığından tasarım kapısı `revision_required_before_N20` olarak kapatıldı. N=20/N=50, null kalibrasyon ve FDR başlatılmadı.
- `PHASE6_DESIGN_NOTE.md`, model keşfini, explainer testini, duyarlılık tablolarını, N=2 agregatlarını, süre projeksiyonlarını, SHA256 ve gizlilik kontrollerini içerir.

## Faz 6.1 DeepSurv explainer/bütçe revizyonu

- `experiments/21_deepsurv_explainer_revision.py`, aynı N=2 replikalar ve sekiz sağlık-kontrol feature'ında GradientExplainer 8/128→512/4096 bütçe eğrisini ve DeepExplainer nested TCGA background k16/k64/k128 tasarımını karşılaştırır.
- Additivity kapısı, maksimum reconstruction hatasının aynı çağrıdaki ampirik risk-skoru aralığına oranı ≤%5 olarak önceden tanımlandı. GradientExplainer ilk kez 256/2048'de geçti; DeepExplainer k16/k64 ve k16/k128 yaklaşık `3e-7` göreli hatayla geçti.
- Frozen `a=0.5` değiştirilmeden uygulanan 0/0.10/0.25/0.50 kontrollü candidate bozulma eğrisi iki seçenekte de monoton F düşüşü verdi. +0.25 bozulmada sekiz feature'ın en az yedisi F<0.8 oldu.
- Önerilen tasarım DeepExplainer candidate k=16 / reference k=128'dir. Yeni N=2 pilotta 77/77 feature değerlendirildi; 11 feature WI≥0.8 geçti ve bottleneck dağılımı R=75, G_claim=2 oldu.
- `PHASE6_1_EXPLAINER_REVISION.md` karşılaştırma, karar ölçütü, yeni pilot, wall-clock, SHA256 ve gizlilik sonuçlarını raporlar. N≥20, null kalibrasyon ve FDR başlatılmadı.

## Faz 6.2 DeepSurv N=50, null gate ve sertifikasyon

- `experiments/22_deepsurv_n50.py`, onaylanan DeepExplainer k16/k128 ve sabit 12+12 query tasarımını 50 patient-cluster bootstrap replikasına çıkarır. İlk N=2 attribution korunur; aynı seed akışında 48 yeni model fit edilir. Ana model değişmez.
- N=50 sonucu 77/77 feature, 8 nokta-yeterli feature ve R=77 bottleneck dağılımıdır. N=2 kazananlarının 1'i korundu, 10'u düştü ve 7 yeni feature geçti.
- SHAP'ın mutlak 0.01 iç additivity assertion'ı iki reference çağrısında tetiklendi. Replikalar atlanmadı; aynı attributionlar `check_additivity=False` ile yeniden hesaplanıp önceden tanımlı göreli ≤%5 kapısından geçirildi. En yüksek reference göreli hata %0.8075'tir.
- `experiments/23_deepsurv_n50_null_gate.py`, DeepSurv N=50 nuisance yapısında 2000×999 recursive_studentized bataryasını çalıştırır. Gate 0/30 Holm ihlaliyle geçti.
- `experiments/24_deepsurv_n50_certification.py`, geçen yöntemle tam 77-feature ailede B=4999 IUT + BH/BY uygular. Sekiz nokta-geçene karşı BH=0 ve BY=0 sertifika bulundu.
- `experiments/25_phase62_three_model_comparison.py`, CoxPH 22/8/8, RSF 20/10/10 ve DeepSurv 8/0/0 sonuçlarını karşılaştırır. DeepSurv nokta-geçen marjinleri düşük, 0.05–0.10 bandı payı %0 ve dağılımı CoxPH'ye daha yakındır.
- `PHASE6_2_REPORT.md`, winner's-curse sayısal desteğini ve N=2 üyeliğinin yüksek oynaklığını, darboğazları, adım bazlı wall-clock, SHA256 ve gizlilik kontrollerini raporlar.

## Makale ilk taslağı ve sayısal tutarlılık denetimi

- `RESULTS_LEDGER.md`, Methods–Availability taslağındaki bütün ampirik sayıların kaynak dosyasını ve dondurulmuş formül bloklarını kaydeder.
- `MANUSCRIPT_METHODS_TO_AVAILABILITY.md`, Abstract/Introduction/Related Work içermeyen ilk olgusal makale taslağıdır.
- `experiments/26_manuscript_consistency.py`, kaynak JSON/CSV dosyalarını yeniden okuyarak başlıca bütün sayıları, formül bloklarını, bölüm kapsamını ve kimlik-sızıntısı örüntülerini denetler; model veya explainer çağırmaz.
- `results/manuscript_consistency/consistency.json` ve `MANUSCRIPT_CONSISTENCY_REPORT.md`, 73/73 başarılı kontrolü, sıfır uyuşmazlığı ve ölçülen otomatik wall-clock sürelerini kaydeder.
- İstenen `section2-warrant-index-taslak.md` bulunamadığı için formül kaynağı olarak korunmuş Faz 2 spesifikasyonu kullanılmıştır; bu açık kaynak-provenans sınırı gönderim öncesinde kapatılmalıdır.

## BMC şekil ve tablo paketi

- `experiments/27_figures_tables.py`, yalnız dondurulmuş CSV/JSON sonuçlarını okuyarak dört yayın şekli ve iki `booktabs` LaTeX tablo snippet'i üretir; model, explainer veya bootstrap çağırmaz.
- `manuscript_assets/figure_1`–`figure_4`, her şeklin 600 dpi PNG, gömülü TrueType fontlu vektörel PDF ve Type 42 fontlu EPS sürümleridir. Figür 1/2/4 genişliği 170 mm, Figür 3 genişliği 85 mm'dir.
- `manuscript_assets/table_1.tex`, üç-model nokta/BH/BY karşılaştırmasıdır. `table_2.tex`, CoxPH 0/0/0/1/76, RSF 0/0/0/5/72 ve DeepSurv 0/0/0/0/77 F/S/U/G_claim/R darboğaz sayılarını verir.
- `FIGURES_TABLES_PROVENANCE.md`, her görsel ve tablonun kaynak artefaktını, dönüşümünü, çapraz kontrolünü ve biçim doğrulamasını kaydeder. Feature-düzeyi gerçek WI marjinleri mevcut olduğundan koşullu Figür 4 üretildi.

## Hasta-düzeyi split sızıntısı teşhisi

- `experiments/28_leakage_diagnostic.py`, üç saklı modeli yeniden eğitmeden orijinal 673/169 dosyadaki hasta kesişimini ve 146-satırlık disjoint test alt-kümesi C-index'lerini hesaplar.
- `LEAKAGE_DIAGNOSTIC.md`, sample-level split'in 23 hastayı iki tarafa taşıdığını; ortak hastaların train/test'te 25/23 satıra karşılık geldiğini ve model eğitiminden sonraki audit filtresinin yalnız XAI query havuzunu temizlediğini tam yol/satır kanıtlarıyla raporlar.
- `results/leakage_diagnostic/diagnostic.json`, orijinal C-index'lerin birebir yeniden üretimini ve disjoint test sonuçlarını saklar: CoxPH 0.700508, RSF 0.729272, DeepSurv 0.763959. Hiçbir gerçek hasta/örnek kimliği dışa aktarılmaz.

## Faz 7.0 hasta-gruplu split ve baseline yeniden eğitimi

- Kardeş projedeki `impute_and_split.py`, ilk 12 karakterlik hasta gruplarıyla `StratifiedGroupKFold(n_splits=5, random_state=42)` kullanır. Yeni split 674/168 satır, 616/153 benzersiz hasta ve sıfır hasta kesişimi üretir.
- `select_features.py` değişmeden yalnız yeni train istatistikleriyle 200 RNA + 200 miRNA + 200 methylation seçer; sabit klinik kategori şeması 16 klinik sütunu ve toplam 616 girdiyi korur.
- `results/phase7_retrain/` içinde ayrı CoxPH, RSF ve DeepSurv artefaktları vardır. Yeni C-index'ler sırasıyla 0.631128, 0.538521 ve 0.496498'dir; eski artefaktların 7/7 SHA256 kontrolü geçmiştir.
- `SPLIT_REDESIGN.md` split gerekçesini ve train-only feature-selection denetimini; `MODEL_RETRAIN_COMPARISON.md` eski/yeni skorları, yakınsama uyarısını, DeepSurv epoch değişimini ve artefakt yollarını raporlar.
- `experiments/29_phase7_retrain_validation.py` ile `results/phase7_retrain/validation.json`, split/panel/ayar/hash kontrollerini doğrular. Faz 7.1 audit pipeline'ı çalıştırılmamıştır.

## Faz 7.0b C-index çöküşü tanısı

- `experiments/30_phase70b_fold_sweep.py`, beş patient-grouped fold'un tamamında aynı preprocessing, train-only seçim ve CoxPH/RSF zincirini çalıştırır. RSF C-index aralığı 0.538521–0.789609'dur; Faz 7.0 fold 0 beş fold içinde en düşüktür.
- `experiments/31_phase70b_cox_diagnostic.py` ve `33_phase70b_cox_fold_retry.py`, uyarılı Cox fit'lerini aynı `penalizer=0.1` ile `max_steps=200, step_size=0.25` kullanarak tanısal olarak yeniden yakınsatır. Fold 0 skoru 0.631128→0.667704, fold 2 skoru 0.597101→0.681159 olur ve iki uyarı da kaybolur.
- Eski/yeni 616 panelde 333 ortak özellik vardır; RNA ortaklığı yalnız 43/200'dür. Beş grouped fold'un ikili toplam ortaklığı 277–328 arasında değişir.
- `experiments/32_phase70b_deepsurv_diagnostic.py`, epoch 38 validation checkpoint'inin test C-index'ini 0.496498 olarak yeniden üretir. Epoch 8'de gözlenen test maksimumu 0.557198'dir; hasta-gruplu iç validation da dış test uçurumunu gidermez.
- `PHASE7_COLLAPSE_DIAGNOSTIC.md` net sonucu fold/panel oynaklığı + Cox yakınsamama + DeepSurv overfitting/checkpoint uyumsuzluğu kombinasyonu olarak raporlar. Audit pipeline'ı çalıştırılmamış, dondurulmuş formüller değiştirilmemiştir.

## Faz 7.1 referans modeller ve beş-fold CV

- `experiments/35_phase71_reference_models.py`, fold-0 CoxPH'yi aynı `penalizer=0.1` ve doğrulanmış `max_steps=200, step_size=0.25` retry ile uyarısız yeniden fit eder. Resmi `coxph_baseline_v2.joblib` C-index'i 0.667704'tür; v1 korunur, RSF/DeepSurv hashleri değişmez.
- `PHASE7_REFERENCE_MODELS.md` Faz 7.2 yönlendirmesini sabitler: CoxPH-v2, mevcut fold-0 RSF ve mevcut fold-0 DeepSurv.
- `experiments/36_phase71_cv.py`, dondurulmuş DeepSurv'ü eksik grouped fold 1–4'te eğitir; fold 0 resmi checkpoint'i yeniden kullanır. CoxPH uyarılı fold 0/2 retry sonuçlarını ve değişmemiş RSF sonuçlarını birleştirir.
- Beş-fold ortalama [min, maks] C-index: CoxPH 0.672006 [0.650173, 0.705882], RSF 0.663116 [0.538521, 0.789609], DeepSurv 0.612080 [0.496498, 0.719193].
- `PHASE7_CV_SUMMARY.md` fold tablosunu, frozen ayarları ve süreleri; `RESULTS_LEDGER_PHASE7.md` kirli split, audit referansı ve CV genelleme tahminini ayrı kategorilerde kaydeder. Faz 7.2 audit pipeline'ı çalıştırılmamıştır.

## Faz 7.2 temiz referans modellerde tam XAI audit

- `experiments/37_phase72_panel.py` ve `38_phase72_download_metabric.py`, temiz fold-0 test havuzunu güncel train-only feature paneliyle yeniden eşler. Yeni ortak aile 27 feature, TCGA audit havuzu 153 benzersiz hasta ve tam METABRIC havuzu 1466 hastadır; RSF/DeepSurv query seti deterministik 12+12 satırdır.
- `experiments/39_phase72_cox.py`, `40_phase72_rsf.py` ve `41_phase72_deepsurv.py`, dondurulmuş model ve açıklayıcı tasarımlarıyla T=50 patient-cluster audit üretir. Nokta WI≥0.8 sayıları CoxPH/RSF/DeepSurv için 9/17/0'dır; DeepSurv göreli additivity hatası %5 sınırının çok altındadır.
- `experiments/42_phase72_calibration.py`, her modelin yeni N=50 nuisance yapısında 2000×999 `single_active`+`tie` bataryasını çalıştırır. Model başına 30 Holm kontrolünde ihlal yoktur.
- `experiments/43_phase72_certification.py`, gate sonrasında m=27 tam aileye B=4999 `recursive_studentized` BH/BY uygular. CoxPH 8/8, RSF 12/12, DeepSurv 0/0 feature sertifiye eder.
- `PHASE72_AUDIT_PANEL.md`, `PHASE72_COXPH.md`, `PHASE72_RSF.md`, `PHASE72_DEEPSURV.md` ve `PHASE72_CALIBRATION.md` adım raporlarıdır. `RESULTS_LEDGER_PHASE72.md`, eski kirli split ile temiz fold-0 C-index ve nokta/BH/BY sonuçlarını farklı aile boyutlarını açıkça göstererek yan yana tutar.
- `tests/test_phase72_artifacts.py`, panel/query şekillerini, model aile hizasını, formül SHA256 değerlerini, DeepSurv additivity kapısını ve sertifikasyonun ortak null gate'e bağlılığını denetler.

## Faz 7.3 özellik-düzeyi raporlama

- `experiments/45_phase73_reporting.py`, yalnız dondurulmuş Faz 7.2 certification CSV ve attribution NPZ artefaktlarını okur; model, explainer, bootstrap, kalibrasyon veya sertifikasyon çağrısı yapmaz.
- `PHASE73_COXPH_FULLTABLE.csv/.md`, `PHASE73_RSF_FULLTABLE.csv/.md` ve `PHASE73_DEEPSURV_FULLTABLE.csv/.md`, 27 özelliğin sembol/kanonik kimlik, tip, F/S/U/G_claim/R/WI, point kararı, p_joint, BH/BY q ve kararlarını taşır.
- `PHASE73_INTERSECTION.md`, CoxPH/RSF BH ve BY kümelerinde 2 CoxPH-only, 6 ortak, 6 RSF-only ve 14 birleşim sonucunu kaydeder; BH ve BY kümeleri aynıdır.
- `PHASE73_FEATURE_TYPE_DISTRIBUTION.md`, 22 RNA + 1 yaş + 4 AJCC-stage-indicator aile yapısını ve model alt kümelerindeki tip paylarını raporlar.
- `PHASE73_DIRECTIONALITY.md` ile eşlik eden CSV'ler, saklı signed attribution ortalamalarını verir. CoxPH–RSF ortak altı özellik aynı 12 TCGA query üzerinde karşılaştırılmış; RSF'de beş ortak RNA özelliğinin ortalaması tam sıfır ve kalan klinik özelliğin işareti ters olduğundan aynı işaret 0/6'dır.
- `PHASE73_SUMMARY.md` dört adımı ve ölçülen süreleri özetler. `tests/test_phase73_reporting.py`, 27 satırı, dondurulmuş 9/8/8–17/12/12–0/0/0 sayılarını, 2/6/6/14 kesişimini ve yön tablosunun bütünlüğünü doğrular.

## Faz 7.4 τ ve η duyarlılığı

- `experiments/46_phase74_sensitivity.py`, Faz 7.2'nin saklı T=50 attribution tensörlerinden aynı B=4999 inferential resample akışını yeniden üretir; yeni patient-cluster model fit'i veya explainer çağrısı yapmaz.
- `PHASE74_TAU_SENSITIVITY.csv/.md`, τ=0.70/0.75/0.80/0.85/0.90 için üç modelin point/BH/BY sayılarını verir. Referans τ=0.80 satırı Faz 7.2 feature-düzeyi component p-değerleri ve kararlarıyla tam eşleşir.
- DeepSurv point/BH/BY sayıları τ=0.70'te 21/6/6, τ=0.75'te 6/1/1, τ=0.80'de 0/0/0'dır; ana 0/0/0 sonucu aşağı yönlü τ değişimine duyarlıdır.
- `PHASE74_ETA_SENSITIVITY.csv/.md` ve `PHASE74_ETA_FEATURE_VALUES.csv`, η=0.01/0.02/0.05/0.10 R dağılımlarını ve darboğazları verir. η=0.10'da R-darboğaz CoxPH/RSF/DeepSurv için 4/25/25'tir; geniş toleransta dayanıklılık nonlinear modellerde sürer, CoxPH'de azalır.
- `PHASE74_SUMMARY.md` fizibilite envanterini ve Supplementary-hazır iki tabloyu içerir. `tests/test_phase74_sensitivity.py`, referans reprodüksiyonunu, τ monotonluğunu, η darboğazlarını ve hiçbir yeni model/explainer/bootstrap fit'i yapılmadığını doğrular.

## Faz 7.5 çoklu-fold audit sağlamlığı

- `PHASE75_FOLD_SELECTION.md`, RSF ve DeepSurv için beş-fold maksimumu olan fold 1'i ek audit olarak seçer; ortak medyan fold bulunmadığından isteğe bağlı ikinci audit çalıştırılmamıştır.
- `experiments/47_phase75_fold1_download_metabric.py`–`53_phase75_fold1_certification.py`, fold 1'e özgü train-only paneli ve 80-feature TCGA–METABRIC ailesini kurar; RSF ve DeepSurv N=50 audit, taze 2000×999 gate ve B=4999 sertifikasyonu tamamlar.
- Fold 1 RSF point/BH/BY sayıları 10/7/7, DeepSurv sayıları 7/0/0'dır. İki model de 0/30 Holm ihlaliyle kalibrasyon kapısını geçmiştir; R darboğazı sırasıyla 76/80 ve 80/80'dir.
- CoxPH deterministik bootstrap replika 15, `max_steps=50` ve dondurulmuş `max_steps=200, step_size=0.25` retry sonrasında yakınsamamıştır. Replika atlanmamış veya değiştirilmemiş; CoxPH gate ve sertifikasyonu NA bırakılmıştır. Tam tanı `results/phase75_fold1_cox/failed_replicate_15.json` içindedir.
- `PHASE75_MULTIFOLD_COMPARISON.md/.csv` fold 0 ve fold 1'i yan yana verir; `PHASE75_SUMMARY.md` süre, kapsam, koşulluluk ve bütünlük kontrollerini özetler. `tests/test_phase75_artifacts.py` panel, gate, sayılar, CoxPH stop kaydı ve formül hashlerini doğrular.

## Faz 7.6 informative ve degenerate warrant sınıflandırması

- `experiments/54_phase76_warrant_classification.py`, yalnız saklı attribution tensörlerini ve fitted fold-0 RSF ağaçlarını okur; model fit, explainer, bootstrap, kalibrasyon veya sertifikasyon çalıştırmaz.
- `PHASE76_WARRANT_CLASSIFICATION_FOLD0.csv/.md`, önceden tanımlı exact-zero kuralıyla CoxPH-v2 için 8 informative/0 degenerate, RSF için 4 informative/8 degenerate sonucu verir ve beklenen sayıları assertion ile doğrular.
- `PHASE76_RSF_SPLIT_USAGE.md`, `feature_importances_` çağrısının scikit-survival tarafından `NotImplementedError` verdiğini ve 500 ağacın `tree_.feature` dizilerinden doğrudan split kullanımını kaydeder. Altı degenerate feature kullanılmamış; OR6C1 10, AJCC Stage IV 25 ağaçta kullanılmıştır.
- `PHASE76_WARRANT_CLASSIFICATION_FOLD1_RSF.csv/.md`, fold-1 RSF'nin yedi BH-sertifikalı özelliğini 5 informative/2 degenerate olarak sınıflandırır. `PHASE76_SUMMARY.md` süreleri ve değişmeyen formül hashlerini verir.
- `tests/test_phase76_warrant_classification.py`, exact-zero sınıflarını, beklenen sayıları, split-kullanım tutarsızlığını ve formül hashlerini doğrular.

## Faz 7.7 cancellation kontrolü ve tam supplementary tablo

- `experiments/55_phase77_cancellation_full_table.py`, saklı attribution tensörlerinden hücre düzeyi exact-zero, M_j, P_j ve Q_j metriklerini üretir; yeni model fit, explainer, bootstrap, kalibrasyon veya sertifikasyon çalıştırmaz.
- `PHASE77_CANCELLATION_CHECK.csv/.md`, fold-0 sekiz ve fold-1 iki degenerate RSF fold-feature sonucu için 10/10 exact-zero assertion'ını doğrular: her birinde 600/600 hücre sıfır ve M_j=P_j=Q_j=0'dır. CoxPH-v2'nin sekiz ve fold-0 RSF'nin dört informative karşılaştırma satırı cellwise-zero değildir.
- `PHASE77_FULL_FEATURE_TABLE_FOLD0.csv/.md`, fold-0'ın 27 ortak özelliğini üç model için 81 satırda attribution özetleri, F/S/U/G/G_signed/R/WI, point ve BH/BY kararları, sertifikalı-satır etiketi ve RSF 500-ağaç split kullanım sayısıyla verir.
- `PHASE77_SUMMARY.md` adım sürelerini ve değişmeyen formül hashlerini kaydeder. `tests/test_phase77_cancellation_full_table.py`, 10/10 exact-zero sonucunu, informative karşılaştırmaları, 81-satırlık şemayı ve frozen kararları doğrular.

## Faz 8 ek sağlamlık ve public GitHub deposu

- `experiments/56_phase8_robustness.py`, saklı query matrisleri ve CoxPH N=50 attribution/log artefaktlarından sorgu-varyansı ve retry/direct diagnostiklerini üretir; yeni model fit, explainer, bootstrap, kalibrasyon veya sertifikasyon çalıştırmaz.
- `PHASE8_QUERY_VARIANCE_CHECK.csv/.md`, dört degenerate fold-feature çiftinin tamamında 24-sorgu varyansının sıfırdan büyük olduğunu gösterir; sıfır-varyans hipotezi 0/4 doğrulanmış, 4/4 reddedilmiştir.
- `PHASE8_COXPH_RETRY_BIAS_CHECK.csv/.md`, 28 retry ve 22 direct replikayı karşılaştırır. Beş feature nominal `|fark|>1.96×birleşik SE` eşiğini aşar; bu tanımlayıcı flag mevcut BH/BY kararlarını değiştirmez.
- Public MIT deposu `https://github.com/umutsaray/warrant-index-survival-xai` adresindedir. `data/`, `results/`, patient-level parquet, model/checkpoint ve attribution tensörleri izleme dışıdır. `PHASE8_GITHUB_REPO_SUMMARY.md` yayın envanterini kaydeder.
- `PHASE8_SUMMARY.md` üç adımı ve sürelerini; `tests/test_phase8_robustness.py` dört hipotez sonucunu, 28/22 partition'ı, beş nominal flag'i ve formül hashlerini doğrular.
