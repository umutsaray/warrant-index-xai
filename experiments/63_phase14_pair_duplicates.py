#!/usr/bin/env python3
"""Phase 14 duplicate-replicate pair diagnostic for S and R."""
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
from _shared import CALIBRATION
def load(name,file):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file); m=importlib.util.module_from_spec(s); sys.modules[name]=m;s.loader.exec_module(m);return m
pilot=load('p14pilot','08_certification_pilot.py')
def losses():
 out={}
 with np.load(ROOT/'results/phase72_cox/attributions_n50.npz') as z:
  src=z['source_phi']; fs=z['features'].astype(str).tolist(); active=z['active_features'].astype(str).tolist(); idx=np.array([fs.index(x) for x in active])
 out['CoxPH-v2']=(active,pilot.pair_losses(src[None],idx,.5),pilot.order_pair_losses(np.abs(src).mean(1)[None],idx,.02))
 for model,folder in [('RSF','phase72_rsf'),('DeepSurv','phase72_deepsurv')]:
  with np.load(ROOT/f'results/{folder}/attributions_n50.npz') as z: phi=z['candidate']; f=z['features'].astype(str).tolist()
  idx=np.arange(len(f)); out[model]=(f,pilot.pair_losses(phi,idx,.5),pilot.order_pair_losses(np.abs(phi).mean(2),idx,.02))
 return out
def main():
 B,T=4999,50; rng=np.random.default_rng(20261107); indices=rng.integers(0,T,size=(B,T)); counts=np.zeros((B,T),np.int16); np.add.at(counts,(np.arange(B)[:,None],indices),1)
 duplicate=(np.sum(counts*(counts-1),axis=1)/(T*(T-1)))
 rows=[]
 for model,(features,sl,rl) in losses().items():
  for component,loss in [('S',sl),('R',rl)]:
   total=np.einsum('bt,sptu,bu->bsp',counts,loss,counts,optimize=True)
   official=(1-total/(T*(T-1))).min(1)
   denom=(T*T-np.sum(counts*counts,axis=1))[:,None,None]
   distinct=(1-total/denom).min(1)
   point=(1-loss.sum(axis=(-2,-1))/(T*(T-1))).min(0)
   for j,f in enumerate(features):
    rows.append(dict(model=model,feature=f,component=component,point=point[j],official_bootstrap_mean=official[:,j].mean(),official_bootstrap_se=official[:,j].std(ddof=1),official_mean_minus_point=official[:,j].mean()-point[j],distinct_bootstrap_mean=distinct[:,j].mean(),distinct_bootstrap_se=distinct[:,j].std(ddof=1),distinct_mean_minus_point=distinct[:,j].mean()-point[j],distinct_minus_official_mean=(distinct[:,j]-official[:,j]).mean(),duplicate_pair_fraction_mean=duplicate.mean(),duplicate_pair_fraction_sd=duplicate.std(ddof=1),expected_duplicate_fraction=1/T))
 pd.DataFrame(rows).to_csv(ROOT/'PHASE14_ADIM2_PAIR_DUPLICATES.csv',index=False)
 print({'duplicate_mean':duplicate.mean(),'expected':1/T,'rows':len(rows)})
if __name__=='__main__':main()
