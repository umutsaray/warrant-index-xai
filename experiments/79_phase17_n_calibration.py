#!/usr/bin/env python3
"""Phase 17 reduced recursive-studentized null calibration for the new N leaf."""
from pathlib import Path
import importlib.util,sys
import numpy as np,pandas as pd
from scipy.stats import beta,binomtest,norm
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase17';OUT.mkdir(parents=True,exist_ok=True);sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
SCALES={'CoxPH-v2':.023934215222165295,'RSF':.01707069437342852,'DeepSurv':.012307440008774696};ALPHAS=(.01,.05,.10);OUTER=1000;INNER=499;T=50
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def leaves(model):
 sm=SCALES[model];m0=.25*sm
 if model=='CoxPH-v2':
  with np.load(ROOT/'results/phase72_cox/attributions_n50.npz') as z:fs=z['features'].astype(str).tolist();act=z['active_features'].astype(str).tolist();ix=np.array([fs.index(f) for f in act]);sites=[z['source_phi'][:,:,ix],z['external_phi']]
 else:
  with np.load(ROOT/'results'/('phase72_rsf' if model=='RSF' else 'phase72_deepsurv')/'attributions_n50.npz') as z:sites=[z['candidate'][0],z['candidate'][1]]
 return np.stack([np.minimum(1,np.abs(x).mean(1)/m0) for x in sites])
def holm(p):
 o=np.argsort(p);z=np.maximum.accumulate((len(p)-np.arange(len(p)))*np.asarray(p)[o]);r=np.empty_like(z);r[o]=np.minimum(z,1);return r
def main():
 remedy=load('11_rsf_remedy_calibration.py','p17nremedy');rows=[]
 for mi,model in enumerate(SCALES):
  real=leaves(model);amp=min(.09,max(.03,float(real.reshape(-1,real.shape[-1]).std(0,ddof=1).mean()*np.sqrt(3*T))));rng=np.random.default_rng(20260918+760000+mi*100000)
  for scenario in ['single_active','tie']:
   pv=[]
   for _ in range(OUTER):
    centers=np.array([.8,.9 if scenario=='single_active' else .8]);u=rng.random((T,2));vals=centers+amp*(2*u-1);idx=rng.integers(0,T,size=(INNER,T));point=vals.mean(0);sample=vals[idx].mean(1);se=vals.std(0,ddof=1)/np.sqrt(T);sse=vals[idx].std(1,ddof=1)/np.sqrt(T);pv.append(float(remedy.studentized_p(point,sample,se,sse).max()))
   for a in ALPHAS:
    k=int((np.asarray(pv)<=a).sum());rows.append(dict(model=model,component='N',scenario=scenario,nominal_alpha=a,rejections=k,outer=OUTER,inner=INNER,empirical_type1=k/OUTER,exact_95ci_lower=0 if k==0 else beta.ppf(.025,k,OUTER-k+1),exact_95ci_upper=1 if k==OUTER else beta.ppf(.975,k+1,OUTER-k),one_sided_binomial_p=binomtest(k,OUTER,a,alternative='greater').pvalue,amplitude=amp,label='reduced 1000x499 v2 N calibration; not official gate'))
 out=pd.DataFrame(rows);out['holm_p_within_model_6']=np.nan
 for model in SCALES:
  m=out.model.eq(model);out.loc[m,'holm_p_within_model_6']=holm(out.loc[m,'one_sided_binomial_p']);out['anti_conservative_holm']=out.holm_p_within_model_6<.05
 out.to_csv(OUT/'PHASE17_ADIM6_N_CALIBRATION.csv',index=False)
if __name__=='__main__':main()
