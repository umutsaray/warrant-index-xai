# Faz 10 — Sorgu yeniden-örnekleme duyarlılığı

Fold-0 dondurulmuş 27-feature panelinde RSF ve DeepSurv, orijinal 12+12 sorguya ek olarak 24+24 ve 48+48 boyutlarında üçer bağımsız deterministik seed ile yeniden değerlendirildi. Bu nokta-duyarlılık analizidir; yeni 4.999-draw sertifikasyon yapılmadı ve resmi kararlar değiştirilmedi.

## Nokta-adequacy kararlılığı

| model | query_per_cohort | min | max | mean |
| --- | --- | --- | --- | --- |
| DeepSurv | 12 | 0 | 0 | 0.0 |
| DeepSurv | 24 | 0 | 0 | 0.0 |
| DeepSurv | 48 | 0 | 0 | 0.0 |
| RSF | 12 | 17 | 17 | 17.0 |
| RSF | 24 | 11 | 15 | 12.666666666666666 |
| RSF | 48 | 12 | 18 | 14.333333333333334 |

## En yüksek WI aralığına sahip özellikler

| model | symbol | WI_min | WI_max | WI_range | WI_std | bottleneck_count | point_status_count |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RSF | AJCC Stage IV | 0.7099686028257446 | 0.9271742543171116 | 0.21720565149136706 | 0.07490281043541092 | 1 | 2 |
| RSF | ACCSL | 0.7096546310832016 | 0.9245525902668762 | 0.2148979591836746 | 0.07173686686940989 | 1 | 2 |
| RSF | MRGPRX1 | 0.657662160058314 | 0.8615123065270627 | 0.20385014646874866 | 0.07435009782220724 | 2 | 2 |
| RSF | Age at index | 0.6730840064257744 | 0.8454945054945058 | 0.1724104990687314 | 0.05767095753284441 | 2 | 2 |
| RSF | STXBP5-AS1 | 0.7320328011363697 | 0.8857896226615112 | 0.15375682152514147 | 0.05259633519965186 | 2 | 2 |
| RSF | OR13C3 | 0.8042857142857145 | 0.9278335949764521 | 0.12354788069073763 | 0.04649157071284024 | 1 | 1 |
| RSF | NPAP1 | 0.6971271585557295 | 0.8122605965463108 | 0.11513343799058129 | 0.03993806279972511 | 1 | 2 |
| RSF | OR6C1 | 0.7970643642072217 | 0.9037519623233912 | 0.10668759811616946 | 0.03748145412534604 | 1 | 2 |
| RSF | LHFPL3 | 0.699591836734694 | 0.8048822605965466 | 0.1052904238618525 | 0.03812603268834663 | 1 | 2 |
| DeepSurv | KCNH5 | 0.6603296703296707 | 0.7608634222919942 | 0.1005337519623235 | 0.038266695815118906 | 1 | 1 |
| RSF | KRTAP26-1 | 0.6789010989010978 | 0.7779434850863417 | 0.09904238618524386 | 0.033675813551833515 | 1 | 1 |
| DeepSurv | CACNG3 | 0.6882888540031395 | 0.7680219780219787 | 0.07973312401883914 | 0.028222961385590283 | 1 | 1 |
| DeepSurv | OR4C46 | 0.6874725274725282 | 0.7617896389324962 | 0.074317111459968 | 0.02539573353360771 | 1 | 1 |
| DeepSurv | RBMY1F | 0.6969073783359512 | 0.770156985871272 | 0.07324960753532084 | 0.024134030002208037 | 1 | 1 |
| DeepSurv | NPAP1 | 0.6273155416012564 | 0.6998430141287288 | 0.07252747252747249 | 0.02751906203130902 | 1 | 1 |
| DeepSurv | LHFPL3 | 0.6597174254317106 | 0.7300000000000005 | 0.0702825745682899 | 0.030086863687383284 | 1 | 1 |
| DeepSurv | OR6C1 | 0.674379905808477 | 0.7355416012558875 | 0.06116169544741057 | 0.021156612261268713 | 1 | 1 |
| DeepSurv | OR13C3 | 0.6810204081632647 | 0.737896389324961 | 0.05687598116169634 | 0.020781447923987187 | 1 | 1 |
| DeepSurv | FBXO40 | 0.6999215070643647 | 0.7524803767660924 | 0.05255886970172763 | 0.019517643702547195 | 1 | 1 |
| DeepSurv | KRTAP11-1 | 0.6517896389324962 | 0.7011930926216643 | 0.0494034536891681 | 0.018791282796814864 | 1 | 1 |

## RSF null-attribution warrant kontrolü

| feature | symbol | query_per_cohort | source_max_abs | external_max_abs | max_abs_attribution | mean_abs_attribution | WI_mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | AJCC Stage IV | 12 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9271742543171116 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | AJCC Stage IV | 24 | 0.0 | 0.0 | 0.0 | 0.0 | 0.8987702773417062 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | AJCC Stage IV | 48 | 6.8521503426580725 | 0.0 | 6.8521503426580725 | 0.005314033070748478 | 0.8391993720565147 |
| rna_ENSG00000169800.14 | RBMY1F | 12 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9271742543171116 |
| rna_ENSG00000169800.14 | RBMY1F | 24 | 0.0 | 0.0 | 0.0 | 0.0 | 0.8987702773417062 |
| rna_ENSG00000169800.14 | RBMY1F | 48 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9023600209314497 |
| rna_ENSG00000185926.1 | OR4C46 | 12 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9271742543171116 |
| rna_ENSG00000185926.1 | OR4C46 | 24 | 0.5131742536976911 | 0.404998980880231 | 0.5131742536976911 | 0.002532911508082596 | 0.8695081109366828 |
| rna_ENSG00000185926.1 | OR4C46 | 48 | 0.5131742536976911 | 0.404998980880231 | 0.5131742536976911 | 0.0017695415397411785 | 0.882637362637363 |
| rna_ENSG00000186965.5 | KRTAP19-2 | 12 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9271742543171116 |
| rna_ENSG00000186965.5 | KRTAP19-2 | 24 | 0.0 | 0.0 | 0.0 | 0.0 | 0.8987702773417062 |
| rna_ENSG00000186965.5 | KRTAP19-2 | 48 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9023600209314497 |
| rna_ENSG00000204246.3 | OR13C3 | 12 | 0.0 | 0.18151789148351538 | 0.18151789148351538 | 0.0012538549942534291 | 0.9278335949764521 |
| rna_ENSG00000204246.3 | OR13C3 | 24 | 0.3941418021214891 | 0.26268020833333394 | 0.3941418021214891 | 0.004148743963429512 | 0.841904761904762 |
| rna_ENSG00000204246.3 | OR13C3 | 48 | 0.46948344398656894 | 0.26268020833333483 | 0.46948344398656894 | 0.002845112291462879 | 0.8872318158032447 |
| rna_ENSG00000205330.4 | OR6C1 | 12 | 0.0 | 0.12454761904762002 | 0.12454761904762002 | 0.002082137028769846 | 0.9037519623233912 |
| rna_ENSG00000205330.4 | OR6C1 | 24 | 1.2261412136474634 | 0.4323392113095238 | 1.2261412136474634 | 0.009269979939722943 | 0.8096075353218214 |
| rna_ENSG00000205330.4 | OR6C1 | 48 | 1.2261412136474634 | 0.548322821969697 | 1.2261412136474634 | 0.00874625846724226 | 0.8166718995290428 |
| rna_ENSG00000244362.3 | KRTAP19-7 | 12 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9271742543171116 |
| rna_ENSG00000244362.3 | KRTAP19-7 | 24 | 0.0 | 0.0 | 0.0 | 0.0 | 0.8987702773417062 |
| rna_ENSG00000244362.3 | KRTAP19-7 | 48 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9023600209314497 |
| rna_ENSG00000254834.4 | OR5M10 | 12 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9271742543171116 |
| rna_ENSG00000254834.4 | OR5M10 | 24 | 0.0 | 0.0 | 0.0 | 0.0 | 0.8987702773417062 |
| rna_ENSG00000254834.4 | OR5M10 | 48 | 0.0 | 0.0 | 0.0 | 0.0 | 0.9023600209314497 |

Null-warrant özellik/boyut hücrelerinin tam-sıfır attribution sayısı **15/24**. Ayrıntılı CSV tüm model × özellik × boyut × seed bileşenlerini içerir.

Faz 7.6/7.7 sınıflaması TCGA/source tarafındaki 12 sorguya dayanıyordu; bu kapsamda sekiz özelliğin sekizi başlangıçta sıfırdır. Daha büyük TCGA sorgularında yalnız **4/8** özellik (`KRTAP19-2, KRTAP19-7, OR5M10, RBMY1F`) bütün boyut ve seed'lerde sıfır kalmıştır. Bu nedenle sekizli null-attribution örüntüsü sorgu örneklemesine duyarlıdır. RSF'de **10/27** özelliğin point-adequacy durumu ve **5/27** özelliğin darboğazı en az bir konfigürasyonda değişmiştir; DeepSurv'de iki sayı da **0/27**'dir. RSF yönü aynı ölçekte kalmakla birlikte 12+12'deki 17 point-adequate sonuç, 24+24'te 11–15 ve 48+48'te 12–18 aralığına açılmıştır; bu nedenle özellik düzeyinde güçlü kararlılık iddiası desteklenmez.

Çalışma süreleri: RSF **8136.0 sn**, DeepSurv **2276.8 sn**. Sonuçlar aynı audit havuzları, aynı leakage dışlamaları, aynı bootstrap modellerinin deterministik yeniden-fit karşılıkları ve dondurulmuş explainer bütçelerine koşulludur. Formül dosyaları değişmedi.
