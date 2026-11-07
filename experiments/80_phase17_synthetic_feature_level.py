#!/usr/bin/env python3
"""Phase 17 real-valued feature-level synthetic validation (OLS ridge, two cohorts)."""
from __future__ import annotations
import argparse,importlib.util,sys,time
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase17';OUT.mkdir(parents=True,exist_ok=True);sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
from _shared import CALIBRATION
from warrantxai.certification import certify
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def ridge(x,y,lam=.1):return np.linalg.solve(x.T@x+lam*np.eye(x.shape[1]),x.T@y)
def n_stats(phi,m0,idx,shared,remedy):
 leaves=np.minimum(1,np.abs(phi).mean(2)/m0);B,T=idx.shape;counts=np.zeros((B,T),np.int16);np.add.at(counts,(np.arange(B)[:,None],idx),1);point=leaves.mean(1);sample,sse=shared.resampled_means(leaves,counts);se=leaves.std(1,ddof=1)/np.sqrt(T);p=remedy.studentized_p(point.T.reshape(-1),np.moveaxis(sample,1,2).reshape(B,-1),se.T.reshape(-1),np.moveaxis(sse,1,2).reshape(B,-1)).reshape(phi.shape[-1],2);return point.min(0),p.max(1)
def stats(phi,idx,shared,remedy,a,b,u,eta,add_n=False,sm=None):
 old={k:CALIBRATION[k] for k in ['a','b','u','eta']};CALIBRATION.update(a=a,b=b,u=u,eta=eta)
 try:pnt,pv=shared.component_statistics(phi,phi,idx)
 finally:CALIBRATION.update(old)
 if add_n:
  n,npv=n_stats(phi,.25*sm,idx,shared,remedy);pnt=np.column_stack([pnt,n]);pv=np.column_stack([pv,npv])
 bh=certify(pv,.05,'bh');by=certify(pv,.05,'by');return pnt,pv,bh.rejected,by.rejected
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repeats',type=int,default=20);a=ap.parse_args();start=time.perf_counter();shared=load('12_rsf_certification.py','p17synshared');remedy=load('11_rsf_remedy_calibration.py','p17synremedy');rng=np.random.default_rng(20261108);n,p,T,B=700,30,50,999;classes=np.array(['A']*8+['B']*6+['C']*6+['D']*5+['E']*5);beta0=np.r_[np.full(8,.05),np.zeros(6),np.full(6,.008),np.full(5,.05),np.full(5,.03)];beta1=beta0.copy();beta1[20:25]=-.05;cov=np.full((p,p),.3);np.fill_diagonal(cov,1);rows=[]
 for rep in range(a.repeats):
  xs=[rng.multivariate_normal(np.zeros(p),cov,n) for _ in range(2)];ys=[xs[0]@beta0+rng.normal(0,.1,n),xs[1]@beta1+rng.normal(0,.1,n)];phi=np.empty((2,T,n,p))
  for d in range(2):
   for t in range(T):
    ix=rng.integers(0,n,n);coef=ridge(xs[d][ix],ys[d][ix]);phi[d,t]=coef*(xs[d]-xs[d].mean(0))
  idx=rng.integers(0,T,size=(B,T));v1,_,v1bh,v1by=stats(phi,idx,shared,remedy,.5,.5,3,.02);m=np.abs(phi[0]).mean(1).mean(0);sm=float(np.median(m[m>0]));v2,_,v2bh,v2by=stats(phi,idx,shared,remedy,sm,sm,6*sm,sm,True,sm)
  conservative=np.abs(phi).mean(2).min(0);ci=np.quantile(conservative,[.025,.975],axis=0);ranks=np.stack([pd.Series(-conservative[t]).rank(method='average').to_numpy() for t in range(T)]);rhi=np.quantile(ranks,.975,axis=0)
  for j in range(p):rows.append(dict(repeat=rep,feature=j,truth_class=classes[j],truth_certifiable=classes[j]=='A',ambiguous_E=classes[j]=='E',scale_s_m=sm,v1_WI=v1[j].min(),v1_BH=v1bh[j],v1_BY=v1by[j],v2_WI6=v2[j].min(),v2_N=v2[j,-1],v2_BH=v2bh[j],v2_BY=v2by[j],B1_ci_lower=ci[0,j],B1_ci_upper=ci[1,j],B1_significant=ci[0,j]>0,B2_rank_upper=rhi[j],B2_top10=rhi[j]<=10))
  print('repeat',rep+1,'/',a.repeats,flush=True)
 out=pd.DataFrame(rows);out.to_csv(OUT/'PHASE17_ADIM5_SYNTH.csv',index=False);print('elapsed',time.perf_counter()-start)
if __name__=='__main__':main()
