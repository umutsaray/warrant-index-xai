"""Phase 21 v3 preregistered secondary 2x2 cells and query-seed diagnostic."""
from __future__ import annotations
import importlib.util,json,sys,time
from pathlib import Path
import joblib,numpy as np,pandas as pd,torch
ROOT=Path(__file__).resolve().parents[1];LOCAL=ROOT/'results/phase21_v3';OUT=ROOT/'phase21_reviewer1/v3';DATA=ROOT/'results/phase7_retrain';SIB=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer';sys.path[:0]=[str(ROOT),str(ROOT/'experiments'),str(SIB)]
from warrantxai.certification import certify
from warrantxai.components import generalizability,generalizability_signed
from _shared import CALIBRATION
COMP5=['F','S','U','G_claim','R']
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;assert s.loader;s.loader.exec_module(m);return m
def main():
 shared=load('12_rsf_certification.py','p21sec_shared');v2=load('77_phase17_v2_certification.py','p21sec_v2');rng=np.random.default_rng(20260918+725000);idx=rng.integers(0,50,size=(4999,50));models={'CoxPH-v3':'cox','RSF-v3':'rsf','DeepSurv-v3':'deepsurv'};rows=[];original={k:CALIBRATION[k] for k in ['a','b','u','eta']}
 try:
  for model,stem in models.items():
   with np.load(LOCAL/f'{stem}_attributions_n50.npz') as z:phi=z['candidate'];ref=z['reference'];features=z['features'].astype(str).tolist()
   imp=np.abs(phi[0]).mean(axis=1).mean(axis=0);sm=float(np.median(imp[imp>0]));sites=[phi[0],phi[1]];npoint,npv=v2.n_stats(sites,.25*sm,idx,shared,load('11_rsf_remedy_calibration.py',f'p21sec_remedy_{stem}'))
   primary=pd.read_csv(OUT/f'{model.lower().replace("-","_")}_feature_certification.csv').set_index('feature').loc[features]
   # absolute WI5/WI6
   p5=primary[[f'p_{c}' for c in COMP5]].to_numpy();z5=primary[COMP5].to_numpy()
   for index in ['WI5','WI6']:
    z=z5 if index=='WI5' else np.column_stack([z5,npoint]);p=p5 if index=='WI5' else np.column_stack([p5,npv]);bh=certify(p,.05,'bh');by=certify(p,.05,'by');names=COMP5 if index=='WI5' else COMP5+['N']
    for j,f in enumerate(features):rows.append({'model':model,'feature':f,'tolerance':'absolute_v1','index':index,'scale_s_m':sm,'WI':z[j].min(),'bottleneck':names[int(z[j].argmin())],'point_adequate':z[j].min()>=.8,'p_joint':bh.joint_p[j],'q_BH':bh.adjusted_p[j],'q_BY':by.adjusted_p[j],'certified_BH':bh.rejected[j],'certified_BY':by.rejected[j]})
   CALIBRATION.update(a=sm,b=sm,u=6*sm,eta=sm);points,pvals=shared.component_statistics(phi,ref,idx)
   for index in ['WI5','WI6']:
    z=points if index=='WI5' else np.column_stack([points,npoint]);p=pvals if index=='WI5' else np.column_stack([pvals,npv]);bh=certify(p,.05,'bh');by=certify(p,.05,'by');names=COMP5 if index=='WI5' else COMP5+['N']
    for j,f in enumerate(features):rows.append({'model':model,'feature':f,'tolerance':'scaled_v2_k1_rho1','index':index,'scale_s_m':sm,'WI':z[j].min(),'bottleneck':names[int(z[j].argmin())],'point_adequate':z[j].min()>=.8,'p_joint':bh.joint_p[j],'q_BH':bh.adjusted_p[j],'q_BY':by.adjusted_p[j],'certified_BH':bh.rejected[j],'certified_BY':by.rejected[j]})
 finally:CALIBRATION.update(original)
 f=pd.DataFrame(rows);f.to_csv(OUT/'secondary_fourcell_features.csv',index=False);s=f.groupby(['model','tolerance','index']).agg(family=('feature','size'),point=('point_adequate','sum'),BH=('certified_BH','sum'),BY=('certified_BY','sum')).reset_index();s.to_csv(OUT/'secondary_fourcell.csv',index=False);(OUT/'SECONDARY_FOURCELL.md').write_text('# Faz 21 v3 ikincil 2x2\n\nBirincil sonuç absolute_v1 × WI5; diğer hücreler ikincildir.\n\n'+s.to_string(index=False)+'\n',encoding='utf-8');print(s.to_string(index=False))
if __name__=='__main__':main()
