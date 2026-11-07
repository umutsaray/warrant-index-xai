#!/usr/bin/env python3
"""Phase 15: diagnostic alternative tests and leverage measures for selected RNA."""
from __future__ import annotations
import importlib.util,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import chi2,norm
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase15';OUT.mkdir(parents=True,exist_ok=True)
SRC=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer'/'brca-multiomics-survival/src/feature_selection/select_features.py'
s=importlib.util.spec_from_file_location('p15select',SRC);sel=importlib.util.module_from_spec(s);sys.modules[s.name]=sel;s.loader.exec_module(sel)
def frame(p):return pd.read_parquet(p).set_index('sample_key')
def null_score(x,t,e):
 order=np.argsort(-t,kind='stable');ts=t[order];xs=x[order];es=e[order].astype(bool);ge=np.empty(len(t),int);a=0
 while a<len(t):
  b=a
  while b+1<len(t) and ts[b+1]==ts[a]:b+=1
  ge[a:b+1]=b;a=b+1
 er=np.flatnonzero(es);rr=ge[er];s0=np.arange(1,len(t)+1,dtype=float)[rr];s1=np.cumsum(xs)[rr];s2=np.cumsum(xs*xs)[rr]
 u=xs[er].sum()-(s1/s0).sum();info=(s2/s0-(s1/s0)**2).sum();ll0=-np.log(s0).sum()
 return u,info,ll0
def optimized_lrt(x,t,e,ll0):
 order=np.argsort(-t,kind='stable');ts=t[order];xs=x[order];es=e[order].astype(bool);ge=np.empty(len(t),int);a=0
 while a<len(t):
  b=a
  while b+1<len(t) and ts[b+1]==ts[a]:b+=1
  ge[a:b+1]=b;a=b+1
 er=np.flatnonzero(es);rr=ge[er]
 def nll(beta):
  eta=np.clip(xs*beta,-50,50);return -(eta[er].sum()-np.log(np.cumsum(np.exp(eta))[rr]).sum())
 fit=minimize_scalar(nll,bounds=(-10,10),method='bounded',options={'xatol':1e-10});llhat=-fit.fun
 return float(chi2.sf(max(0,2*(llhat-ll0)),1)),float(fit.x),float(llhat)
def main():
 x=frame(ROOT/'results/phase7_retrain/X_train.parquet');y=frame(ROOT/'results/phase7_retrain/y_train.parquet').loc[x.index]
 chosen=json.loads((ROOT/'results/phase7_retrain/selected_features.json').read_text(encoding='utf-8'))['selected_by_block']['rna']
 saved=pd.read_csv(ROOT/'results/phase7_retrain/univariate_cox_rna.tsv',sep='\t').set_index('feature').loc[chosen]
 vals=x[chosen].to_numpy(float);t=y.time.to_numpy(float);e=y.event.to_numpy(int)
 b,se,p,ll=sel.univariate_cox_chunk(vals,t,e,50,1e-7)
 np.testing.assert_allclose(p,saved.p_value.to_numpy(),rtol=1e-6,atol=0)
 mu=vals.mean(0);sd=vals.std(0);z=(vals-mu)/sd;zw=np.clip(z,-5,5)
 bw,sew,pw,llw=sel.univariate_cox_chunk(zw,t,e,50,1e-7)
 rows=[]
 for j,f in enumerate(chosen):
  u,info,ll0=null_score(z[:,j],t,e);ps=float(chi2.sf(u*u/info,1));pl,lrt_beta,lrt_ll=optimized_lrt(z[:,j],t,e,ll0);imax=int(np.argmax(np.abs(z[:,j])));nz=np.flatnonzero(vals[:,j]!=0)
  rows.append(dict(feature=f,nonzero_count=len(nz),zero_rate=float(np.mean(vals[:,j]==0)),max_abs_z=float(abs(z[imax,j])),max_z=float(z[imax,j]),max_z_event=int(e[imax]),max_z_time=float(t[imax]),max_z_squared_share=float(z[imax,j]**2/np.sum(z[:,j]**2)),nonzero_events='|'.join(map(str,e[nz].tolist())),nonzero_times='|'.join(map(lambda q:f'{q:.8g}',t[nz])),wald_beta=b[j],wald_p=p[j],lrt_beta=lrt_beta,lrt_log_likelihood=lrt_ll,lrt_p=pl,score_p=ps,winsor5_wald_beta=bw[j],winsor5_wald_p=pw[j],official_p_saved=saved.iloc[j].p_value,official_rank_saved=saved.iloc[j]['rank']))
 out=pd.DataFrame(rows)
 for col,pcol in [('wald_rank_200','wald_p'),('lrt_rank_200','lrt_p'),('score_rank_200','score_p'),('winsor5_rank_200','winsor5_wald_p')]:out[col]=out[pcol].rank(method='first').astype(int)
 active=set(np.load(ROOT/'results/phase72_panel/query_matrices.npz',allow_pickle=False)['active_features'].astype(str));out['fold0_active']=out.feature.isin(active);out['retained_top200_all_diagnostic_criteria']=True
 out.to_csv(OUT/'PHASE15_ADIM2_SPARSE_LEVERAGE.csv',index=False)
if __name__=='__main__':main()
