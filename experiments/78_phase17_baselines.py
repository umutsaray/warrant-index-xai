#!/usr/bin/env python3
"""Phase 17 stored-tensor baselines B1/B2/B3."""
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase17';OUT.mkdir(parents=True,exist_ok=True)
def main():
 v2=pd.read_csv(OUT/'PHASE17_ADIM3_V2_CERT_PRIMARY.csv').set_index(['model','feature']);rows=[]
 for model,src,off in [('CoxPH-v2','phase72_cox','phase72_coxph_certification'),('RSF','phase72_rsf','phase72_rsf_certification'),('DeepSurv','phase72_deepsurv','phase72_deepsurv_certification')]:
  if model=='CoxPH-v2':
   with np.load(ROOT/'results/phase72_cox/attributions_n50.npz') as z:
    fs=z['features'].astype(str).tolist();features=z['active_features'].astype(str).tolist();ix=np.array([fs.index(f) for f in features]);imp=np.stack([np.abs(z['source_phi'][:,:,ix]).mean(1),np.abs(z['external_phi']).mean(1)])
  else:
   with np.load(ROOT/'results'/src/'attributions_n50.npz') as z:features=z['features'].astype(str).tolist();imp=np.abs(z['candidate']).mean(2)
  conservative=imp.min(0);rank=np.empty_like(conservative)
  for t in range(50):rank[t]=pd.Series(-conservative[t]).rank(method='average').to_numpy()
  official=pd.read_csv(ROOT/'results'/off/'feature_certification.csv').set_index('feature');
  for j,f in enumerate(features):
   ci=np.quantile(conservative[:,j],[.025,.975]);ri=np.quantile(rank[:,j],[.025,.975]);q=v2.loc[(model,f)];rows.append(dict(model=model,feature=f,B1_importance_point=conservative[:,j].mean(),B1_ci95_lower=ci[0],B1_ci95_upper=ci[1],B1_significant=ci[0]>0,B2_rank_median=np.median(rank[:,j]),B2_rank_ci95_lower=ri[0],B2_rank_ci95_upper=ri[1],B2_rank_stable_top10=ri[1]<=10,B3_v1_point=bool(official.loc[f,'point_WI_ge_0_8']),B3_v2_point=bool(q.point_adequate6),v1_BH=bool(official.loc[f,'certified_BH_alpha_0_05']),v1_BY=bool(official.loc[f,'certified_BY_alpha_0_05']),v2_BH=bool(q.certified_BH),v2_BY=bool(q.certified_BY),v2_N=q.N))
 out=pd.DataFrame(rows);labels=pd.read_csv(ROOT/'PHASE76_WARRANT_CLASSIFICATION_FOLD0.csv');labelmap=dict(zip(labels['Feature ID'],labels.Label));out['phase76_warrant_label']=out.feature.map(labelmap).fillna('not classified');out.to_csv(OUT/'PHASE17_ADIM4_BASELINES.csv',index=False)
if __name__=='__main__':main()
