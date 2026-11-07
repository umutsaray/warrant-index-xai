#!/usr/bin/env python3
"""Phase 16 reduced R pair simulation at the null and near-boundary alternatives."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import beta
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase16';OUT.mkdir(parents=True,exist_ok=True)
T=50;OUTER=1000;INNER=499;CENTERS=(.80,.81,.82,.83,.85,.90);ALPHAS=(.01,.05,.10)
def site(values,coef,counts):
 loss=(values[:,None]!=values[None,:]).astype(float);total=np.einsum('bi,ij,bj->b',counts,loss,counts,optimize=True);point=1-coef*loss.sum()/(T*(T-1));official=1-coef*total/(T*(T-1));distinct=1-coef*total/(T*T-np.sum(counts*counts,axis=1));return point,official,distinct
def p(point,sample):return (1+np.sum(sample-point>=point-.8))/(len(sample)+1)
def ci(k):return (0 if k==0 else beta.ppf(.025,k,OUTER-k+1),1 if k==OUTER else beta.ppf(.975,k+1,OUTER-k))
def main():
 rng=np.random.default_rng(20261108);store={(s,c,d):[] for s in ['single_active','tie'] for c in CENTERS for d in ['official_duplicates_in_denominator','distinct_original_pairs_only']}
 for scenario in ['single_active','tie']:
  for _ in range(OUTER):
   v=(rng.random((T,2))>.5);idx=rng.integers(0,T,size=(INNER,T));counts=np.zeros((INNER,T),np.int16);np.add.at(counts,(np.arange(INNER)[:,None],idx),1)
   for center in CENTERS:
    co=[2*(1-center),.2 if scenario=='single_active' else 2*(1-center)];sites=[site(v[:,d],co[d],counts) for d in range(2)]
    store[(scenario,center,'official_duplicates_in_denominator')].append(max(p(x[0],x[1]) for x in sites));store[(scenario,center,'distinct_original_pairs_only')].append(max(p(x[0],x[2]) for x in sites))
 rows=[]
 for (scenario,center,den),vals in store.items():
  vals=np.asarray(vals)
  for a in ALPHAS:
   k=int(np.sum(vals<=a));lo,hi=ci(k);rows.append(dict(component='R',scenario=scenario,target_center=center,interpretation='null' if center==.8 else 'fixed alternative',denominator=den,nominal_alpha=a,rejections=k,outer=OUTER,inner=INNER,rate=k/OUTER,exact_95ci_lower=lo,exact_95ci_upper=hi,label='reduced diagnostic; not official 2000x999 gate'))
 pd.DataFrame(rows).to_csv(OUT/'PHASE16_ADIM3_R_NULL_AND_NEAR_BOUNDARY.csv',index=False)
if __name__=='__main__':main()
