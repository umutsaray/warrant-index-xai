from pathlib import Path
import hashlib
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase15'
def test_phase15_outputs_and_frozen_hashes():
 assert len(pd.read_csv(OUT/'PHASE15_ADIM1_DEEPSURV_VAL_NOOVERLAP.csv'))==15
 assert len(pd.read_csv(OUT/'PHASE15_ADIM2_SPARSE_LEVERAGE.csv'))==200
 assert len(pd.read_csv(OUT/'PHASE15_ADIM3_METABRIC_MAPPING.csv'))==27
 assert len(pd.read_csv(OUT/'PHASE15_ADIM3_G_STANDARDIZATION.csv'))==81
 assert len(pd.read_csv(OUT/'PHASE15_ADIM4_DUPLICATE_SIM.csv'))==36
 assert len(pd.read_csv(OUT/'PHASE15_ADIM5_KRTAP11_1_RESIDUAL.csv'))==24
 expected={'components.py':'7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa','aggregation.py':'10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68','certification.py':'76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487'}
 for f,h in expected.items():assert hashlib.file_digest(open(ROOT/'warrantxai'/f,'rb'),'sha256').hexdigest()==h
def test_phase15_directional_findings():
 val=pd.read_csv(OUT/'PHASE15_ADIM1_DEEPSURV_VAL_NOOVERLAP.csv');assert val[val.subset.eq('no_fit_patient_overlap')].c_index.mean()>val[val.subset.eq('all')].c_index.mean()
 g=pd.read_csv(OUT/'PHASE15_ADIM3_G_STANDARDIZATION.csv');counts=g.groupby('method').point_adequate.sum().to_dict();assert counts['official_cohort_standardized']==9 and counts['tcga_train_scaler_only_no_cohort_standardization']==3
 k=pd.read_csv(OUT/'PHASE15_ADIM5_KRTAP11_1_RESIDUAL.csv');assert (k.candidate8==0).all() and (k.reference128_recomputed==0).all() and k.reference_matches_stored_atol_1e_12.all()
