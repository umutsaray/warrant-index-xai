from pathlib import Path
import hashlib,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase16'
def test_phase16_outputs_and_hashes():
 assert len(pd.read_csv(OUT/'PHASE16_ADIM1_FULL_POOL_RANKING.csv'))==60660
 assert len(pd.read_csv(OUT/'PHASE16_ADIM1_FOCUS_GENES.csv'))==8
 assert len(pd.read_csv(OUT/'PHASE16_ADIM2_PLATFORM_SHIFT.csv'))==54
 assert len(pd.read_csv(OUT/'PHASE16_ADIM2_G_FOUR_FORMS.csv'))==108
 assert len(pd.read_csv(OUT/'PHASE16_ADIM3_R_NULL_AND_NEAR_BOUNDARY.csv'))==72
 h={'components.py':'7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa','aggregation.py':'10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68','certification.py':'76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487'}
 for f,v in h.items():assert hashlib.file_digest(open(ROOT/'warrantxai'/f,'rb'),'sha256').hexdigest()==v
def test_phase16_directional_results():
 x=pd.read_csv(OUT/'PHASE16_ADIM1_FULL_POOL_RANKING.csv');assert x[x.wald_rank<=200].zero_rate.gt(.9).sum()==166;assert x[x.lrt_rank<=200].zero_rate.gt(.9).sum()==14
 assert x[x.fold0_active].lrt_rank.le(200).sum()==1 and x[x.fold0_active].winsor5_rank.le(200).sum()==3
 g=pd.read_csv(OUT/'PHASE16_ADIM2_G_FOUR_FORMS.csv');n=g.groupby('form').point_adequate.sum();assert n['official_location_and_scale']==9 and n['location_only']==7 and n['scale_only']==3
 r=pd.read_csv(OUT/'PHASE16_ADIM3_R_NULL_AND_NEAR_BOUNDARY.csv');assert not ((r.target_center==.8)&(r.rate>r.nominal_alpha)).any()
