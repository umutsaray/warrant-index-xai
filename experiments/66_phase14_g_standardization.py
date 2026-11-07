#!/usr/bin/env python3
"""Phase 14 Cox G sensitivity to METABRIC standardization."""
from __future__ import annotations
import importlib.util,sys,joblib
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from warrantxai.components import generalizability,generalizability_signed
def loadphase3():
 s=importlib.util.spec_from_file_location('p14phase3',ROOT/'experiments/06_brca_real.py');m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def grid(x,n):
 x=np.sort(x); ix=np.minimum((np.arange(n)+.5)*len(x)//n,len(x)-1).astype(int);return x[ix]
def calc(src,ext):
 gm=generalizability(np.stack([np.abs(src).mean((0,1)),np.abs(ext).mean((0,1))]),.5)
 ms,me=src.mean(0),ext.mean(0);n=max(len(ms),len(me)); coupled=np.stack([np.column_stack([grid(ms[:,j],n) for j in range(src.shape[-1])]),np.column_stack([grid(me[:,j],n) for j in range(src.shape[-1])])])
 gs=generalizability_signed(coupled,.5);return gm,gs,np.minimum(gm,gs)
def main():
 with np.load(ROOT/'results/phase72_cox/attributions_n50.npz') as z: source=z['source_phi'];official_ext=z['external_phi'];coef=z['coefficients'];features=z['features'].astype(str).tolist();active=z['active_features'].astype(str).tolist()
 idx=np.array([features.index(f) for f in active]);src=source[:,:,idx]
 p3=loadphase3();p3.RAW=ROOT/'data/external/metabric'; ext,mapping,info=p3.external_columns(list(joblib.load(ROOT/'results/phase7_retrain/rsf_baseline.joblib')['features']))
 missing=[f for f in active if f not in ext]
 if missing:
  pd.DataFrame([{'status':'yapılamadı: yerel ham METABRIC cache aktif 27 özelliğin tamamını yeniden eşlemiyor','missing_feature_count':len(missing),'missing_features':'|'.join(missing)}]).to_csv(ROOT/'PHASE14_ADIM6_G_STANDARDIZATION.csv',index=False)
  return
 raw=np.column_stack([ext[f] for f in active]);raw=raw[np.isfinite(raw).all(1)]
 cox=joblib.load(ROOT/'results/phase7_retrain/coxph_baseline_v2.joblib');sc=cox['scaler']; raw_scaled=np.column_stack([(raw[:,j]-sc.mean_[idx[j]])/sc.scale_[idx[j]] for j in range(len(active))]);alt=coef[:,idx,None]*raw_scaled.T[None];alt=np.transpose(alt,(0,2,1))
 report=pd.read_csv(ROOT/'results/phase72_cox/feature_report.csv').set_index('feature').loc[active];rows=[]
 for method,ex in [('official_cohort_standardized',official_ext),('raw_METABRIC_TCGA_train_scaler',alt),('official_rank_grid_control',official_ext)]:
  gm,gs,gc=calc(src,ex)
  for j,f in enumerate(active):
   other=report.loc[f,['F','S','U','R']].to_numpy(float); vals=np.r_[other[:3],gc[j],other[3]]
   rows.append(dict(feature=f,method=method,G=gm[j],G_signed=gs[j],G_claim=gc[j],WI=vals.min(),bottleneck=['F','S','U','G_claim','R'][vals.argmin()],point_adequate=vals.min()>.8))
 pd.DataFrame(rows).to_csv(ROOT/'PHASE14_ADIM6_G_STANDARDIZATION.csv',index=False)
if __name__=='__main__':main()
