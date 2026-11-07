#!/usr/bin/env python3
"""Phase 15: reconstruct the Phase 7.2 METABRIC mapping and Cox G sensitivity."""
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import joblib,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase15';OUT.mkdir(parents=True,exist_ok=True);sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
from warrantxai.components import generalizability,generalizability_signed
def loadp3():
 s=importlib.util.spec_from_file_location('p15p3',ROOT/'experiments/06_brca_real.py');m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def calc(source,external,coef,anchors,external_scaled,phase3):
 gm=generalizability(np.stack([np.abs(source).mean((0,1)),np.abs(external).mean((0,1))]),.5)
 mc=coef.mean(0);grid=max(source.shape[1],external.shape[1]);ss=np.column_stack([phase3.empirical_rank_grid(mc[j]*anchors[:,j],grid) for j in range(len(mc))]);es=np.column_stack([phase3.empirical_rank_grid(mc[j]*external_scaled[:,j],grid) for j in range(len(mc))]);gs=generalizability_signed(np.stack([ss,es]),.5);return gm,gs,np.minimum(gm,gs)
def main():
 phase3=loadp3();phase3.RAW=ROOT/'data/external/metabric_phase72';train=pd.read_parquet(ROOT/'results/phase7_retrain/X_train_selected.parquet').set_index('sample_key');ext,mapping,info=phase3.external_columns(list(train.columns))
 with np.load(ROOT/'results/phase72_panel/query_matrices.npz') as z:active=z['active_features'].astype(str).tolist();stored=z['external_full'];source_pool=z['source_full']
 raw=np.column_stack([ext[f] for f in active]);complete=np.isfinite(raw).all(1);raw=raw[complete];idx=np.array([train.columns.get_loc(f) for f in active]);mu=train.mean().to_numpy()[idx];sd=train.std(ddof=0).to_numpy()[idx];aligned=mu+(raw-raw.mean(0))/raw.std(0)*sd
 if not np.allclose(aligned,stored,atol=1e-8,rtol=0):
  q=np.unravel_index(np.abs(aligned-stored).argmax(),aligned.shape);raise RuntimeError(f'external_full mismatch max={np.abs(aligned-stored).max()} feature={active[q[1]]} row={q[0]}')
 mp=mapping.set_index('feature').reindex(active).reset_index();mp['available_in_phase72_cache']=[f in ext for f in active];mp['complete_rows']=len(raw);mp['phase14_failure_cause']='Phase 14 diagnostic pointed helper at data/external/metabric instead of data/external/metabric_phase72';mp.to_csv(OUT/'PHASE15_ADIM3_METABRIC_MAPPING.csv',index=False)
 with np.load(ROOT/'results/phase72_cox/attributions_n50.npz') as z:source=z['source_phi'];stored_ext=z['external_phi'];coef_all=z['coefficients'];features=z['features'].astype(str).tolist()
 ci=np.array([features.index(f) for f in active]);coef=coef_all[:,ci];source=source[:,:,ci];cox=joblib.load(ROOT/'results/phase7_retrain/coxph_baseline_v2.joblib');sc=cox['scaler'];anchors=(source_pool-sc.mean_[ci])/sc.scale_[ci]
 variants={
  'official_cohort_standardized':(aligned-sc.mean_[ci])/sc.scale_[ci],
  'tcga_train_scaler_only_no_cohort_standardization':(raw-sc.mean_[ci])/sc.scale_[ci],
  'official_recalculation_control':(aligned-sc.mean_[ci])/sc.scale_[ci]}
 report=pd.read_csv(ROOT/'results/phase72_cox/feature_report.csv').set_index('feature').loc[active];rows=[]
 for method,es in variants.items():
  ep=np.transpose(coef[:,:,None]*es.T[None],(0,2,1));gm,gs,gc=calc(source,ep,coef,anchors,es,phase3)
  if method=='official_recalculation_control':
   np.testing.assert_allclose(ep,stored_ext,atol=1e-12,rtol=0);np.testing.assert_allclose(gm,report.G_magnitude,atol=1e-12,rtol=0);np.testing.assert_allclose(gs,report.G_signed,atol=1e-12,rtol=0)
  for j,f in enumerate(active):
   z=np.array([report.loc[f,'F'],report.loc[f,'S'],report.loc[f,'U'],gc[j],report.loc[f,'R']]);rows.append(dict(feature=f,method=method,G_magnitude=gm[j],G_signed=gs[j],G_claim=gc[j],WI=z.min(),point_adequate=z.min()>.8,bottleneck=['F','S','U','G_claim','R'][z.argmin()],external_raw_mean=raw[:,j].mean(),external_raw_sd=raw[:,j].std()))
 pd.DataFrame(rows).to_csv(OUT/'PHASE15_ADIM3_G_STANDARDIZATION.csv',index=False)
if __name__=='__main__':main()
