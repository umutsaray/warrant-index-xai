"""Phase 21 v3 fresh null gates and primary certification."""
from __future__ import annotations
import importlib.util,json,sys,time
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import beta,binomtest,norm
ROOT=Path(__file__).resolve().parents[1];LOCAL=ROOT/'results/phase21_v3';OUT=ROOT/'phase21_reviewer1/v3';sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
from warrantxai.certification import certify
COMP=('F','S','U','G_claim','R');SEED=20260918
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;assert s.loader;s.loader.exec_module(m);return m
def nuisance(phi,ref,pilot,null,seed,draws=399):
 rng=np.random.default_rng(seed+9100);idx=rng.integers(0,50,size=(draws,50));_,sampled=pilot.score_bootstraps(phi=phi,reference=ref,target_indices=np.arange(phi.shape[-1]),bootstrap_indices=idx);flat=sampled.reshape(-1,5);var=flat.std(0,ddof=1)>1e-12;corr=np.eye(5)
 if var.sum()>1:corr[np.ix_(var,var)]=np.corrcoef(flat[:,var],rowvar=False)
 corr=null.nearest_correlation(corr);amp=np.clip(flat.std(0,ddof=1)*np.sqrt(150),.03,.18);amp[[0,2]]=np.minimum(amp[[0,2]],.09);return corr,amp
def holm(p):
 order=np.argsort(p);out=np.empty(len(p));run=0
 for rank,i in enumerate(order):run=max(run,(len(p)-rank)*p[i]);out[i]=min(1,run)
 return out
def main():
 started=time.perf_counter();remedy=load('11_rsf_remedy_calibration.py','p21v3_remedy');null=load('09_null_calibration.py','p21v3_null');pilot=load('08_certification_pilot.py','p21v3_pilot');shared=load('12_rsf_certification.py','p21v3_shared');models={'CoxPH-v3':'cox','RSF-v3':'rsf','DeepSurv-v3':'deepsurv'};tensors={};gate_rows=[]
 for mi,(model,stem) in enumerate(models.items()):
  with np.load(LOCAL/f'{stem}_attributions_n50.npz') as z:phi=z['candidate'];ref=z['reference'];features=z['features'].astype(str).tolist()
  tensors[model]=(phi,ref,features);corr,amp=nuisance(phi,ref,pilot,null,SEED+mi);rng=np.random.default_rng(SEED+720000+mi*100000);rows=[]
  for scenario in ('single_active','tie'):
   pv=np.empty((2000,5))
   for sim in range(2000):
    normal=np.stack([rng.multivariate_normal(np.zeros(5),corr,size=50) for _ in range(2)],axis=1);uniform=norm.cdf(normal);idx=rng.integers(0,50,size=(999,50));pv[sim]=remedy.candidate_pvalues(uniform,amp,idx,scenario,null)['recursive_studentized']
   for j,c in enumerate(COMP):
    for alpha in (.01,.05,.10):
     n=int((pv[:,j]<=alpha).sum());rows.append({'model':model,'scenario':scenario,'component':c,'alpha':alpha,'rejections':n,'empirical_type1':n/2000,'binomial_p':float(binomtest(n,2000,alpha,alternative='greater').pvalue)})
  f=pd.DataFrame(rows);f['holm_p']=holm(f.binomial_p.to_numpy());f['anti_conservative']=f.holm_p<.05;gate_rows.append(f);print('gate simulated',model,flush=True)
 gate=pd.concat(gate_rows,ignore_index=True);gate.to_csv(OUT/'null_gate.csv',index=False);summary=gate.groupby('model').anti_conservative.sum().rename('holm_flags').reset_index();summary['gate_passed']=summary.holm_flags.eq(0);summary.to_csv(OUT/'null_gate_summary.csv',index=False);(OUT/'NULL_GATE.md').write_text('# Faz 21 v3 null kapısı\n\n'+summary.to_string(index=False)+'\n',encoding='utf-8')
 if not summary.gate_passed.all():print(summary);raise SystemExit(2)
 rng=np.random.default_rng(SEED+530000);idx=rng.integers(0,50,size=(4999,50));allrows=[]
 for model,(phi,ref,features) in tensors.items():
  points,pvals=shared.component_statistics(phi,ref,idx);bh=certify(pvals,.05,'bh');by=certify(pvals,.05,'by');rows=[]
  for j,f in enumerate(features):rows.append({'model':model,'feature':f,**{c:points[j,k] for k,c in enumerate(COMP)},'WI':points[j].min(),'bottleneck':COMP[int(points[j].argmin())],'point_adequate':points[j].min()>=.8,**{f'p_{c}':pvals[j,k] for k,c in enumerate(COMP)},'p_joint':bh.joint_p[j],'q_BH':bh.adjusted_p[j],'q_BY':by.adjusted_p[j],'certified_BH':bh.rejected[j],'certified_BY':by.rejected[j]})
  frame=pd.DataFrame(rows);frame.to_csv(OUT/f'{model.lower().replace("-","_")}_feature_certification.csv',index=False);allrows.append(frame);print('certified',model,flush=True)
 combined=pd.concat(allrows,ignore_index=True);combined.to_csv(OUT/'matched_three_model.csv',index=False);agg=combined.groupby('model').agg(family=('feature','size'),point=('point_adequate','sum'),BH=('certified_BH','sum'),BY=('certified_BY','sum')).reset_index();agg.to_csv(OUT/'certification_summary.csv',index=False);(OUT/'CERTIFICATION.md').write_text('# Faz 21 v3 birincil sertifikasyon\n\nBirincil sonuç v1 mutlak toleranslı WI5, recursive_studentized B=4,999.\n\n'+agg.to_string(index=False)+'\n',encoding='utf-8');print(agg.to_string(index=False));print('seconds',time.perf_counter()-started)
if __name__=='__main__':main()
