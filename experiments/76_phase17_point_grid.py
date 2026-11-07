#!/usr/bin/env python3
"""Phase 17 v2 point-estimate grid from frozen attribution tensors."""
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase17';OUT.mkdir(parents=True,exist_ok=True);sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
SCALES={'CoxPH-v2':0.023934215222165295,'RSF':0.01707069437342852,'DeepSurv':0.012307440008774696};KS=(.25,.5,1.,2.,4.);RHOS=(.5,1.,2.);KAPPAS=(.1,.25,.5)
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def n_scores(item,m0):
 if 'external' in item:
  src=item['phi'][0][:,:,item['targets']];sites=[src,item['external']]
 else:sites=[item['phi'][0],item['phi'][1]]
 leaves=[np.minimum(1,np.abs(x).mean(axis=1)/m0).mean(axis=0) for x in sites];return np.minimum(leaves[0],leaves[1])
def main():
 p12=load('60_phase12_scale_tie_query_audit.py','p17p12');data=p12.load_inputs();rows=[]
 for model,item in data.items():
  sm=SCALES[model]
  for k in KS:
   for rho in RHOS:
    z=p12.scores(item,k*sm,k*sm,6*k*sm,rho*sm)
    for kap in KAPPAS:
     n=n_scores(item,kap*sm)
     for j,f in enumerate(item['features']):
      vals5=z.loc[f,['F','S','U','G_claim','R']];vals6=np.r_[vals5.to_numpy(float),n[j]];names=['F','S','U','G_claim','R','N'];rows.append(dict(model=model,feature=f,scale_s_m=sm,k=k,rho=rho,kappa=kap,a=k*sm,b=k*sm,u=6*k*sm,eta=rho*sm,m0=kap*sm,F=z.loc[f,'F'],S=z.loc[f,'S'],U=z.loc[f,'U'],G=z.loc[f,'G'],G_signed=z.loc[f,'G_signed'],G_claim=z.loc[f,'G_claim'],R=z.loc[f,'R'],N=n[j],WI5=vals5.min(),WI6=vals6.min(),bottleneck6=names[int(vals6.argmin())],point_adequate5=vals5.min()>.8,point_adequate6=vals6.min()>.8,WI6_equals_R=bool(vals6.min()==z.loc[f,'R'])))
 pd.DataFrame(rows).to_csv(OUT/'PHASE17_ADIM2_POINT_GRID.csv',index=False)
if __name__=='__main__':main()
