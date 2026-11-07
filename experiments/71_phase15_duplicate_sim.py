#!/usr/bin/env python3
"""Phase 15 reduced duplicate-pair simulation (diagnostic, not a calibration gate)."""
from __future__ import annotations
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import beta
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase15';OUT.mkdir(parents=True,exist_ok=True)
ALPHAS=(.01,.05,.10);T=50;OUTER=1000;INNER=499
def scores(values,coef,counts):
 loss=np.abs(values[:,None]-values[None,:]);total=np.einsum('bi,ij,bj->b',counts,loss,counts,optimize=True);official=1-coef*total/(T*(T-1));den=T*T-np.sum(counts*counts,axis=1);corrected=1-coef*total/den;point=1-coef*loss.sum()/(T*(T-1));return point,official,corrected
def pval(point,sampled):return (1+np.sum(sampled-point>=point-.8))/(len(sampled)+1)
def ci(k,n):return (0 if k==0 else beta.ppf(.025,k,n-k+1),1 if k==n else beta.ppf(.975,k+1,n-k))
def main():
 rng=np.random.default_rng(20261108);records=[];specs=[]
 for scenario in ['single_active','tie']:
  specs.append(('S',scenario,.8))
  for center in [.85,.90]:specs.append(('R',scenario,center))
 for component,scenario,center in specs:
  po=np.empty(OUTER);pc=np.empty(OUTER)
  for o in range(OUTER):
   values=rng.random((T,2)) if component=='S' else (rng.random((T,2))>.5).astype(float)
   idx=rng.integers(0,T,size=(INNER,T));counts=np.zeros((INNER,T),np.int16);np.add.at(counts,(np.arange(INNER)[:,None],idx),1)
   if component=='S':coefs=[.6,.3 if scenario=='single_active' else .6]
   else:
    active_coef=2*(1-center);coefs=[active_coef,.1 if scenario=='single_active' else active_coef]
   pp=[];oo=[];cc=[]
   for site in range(2):
    point,official,corrected=scores(values[:,site],coefs[site],counts);pp.append(point);oo.append(official);cc.append(corrected)
   po[o]=max(pval(pp[s],oo[s]) for s in range(2));pc[o]=max(pval(pp[s],cc[s]) for s in range(2))
  for denominator,p in [('official_duplicates_in_denominator',po),('distinct_original_pairs_only',pc)]:
   for a in ALPHAS:
    k=int(np.sum(p<=a));lo,hi=ci(k,OUTER);records.append(dict(component=component,null_or_alternative_scenario=scenario,target_center=center,denominator=denominator,nominal_alpha=a,rejections=k,outer=OUTER,inner=INNER,rate=k/OUTER,exact_95ci_lower=lo,exact_95ci_upper=hi,label='reduced diagnostic; not the official 2000x999 gate'))
 pd.DataFrame(records).to_csv(OUT/'PHASE15_ADIM4_DUPLICATE_SIM.csv',index=False)
if __name__=='__main__':main()
