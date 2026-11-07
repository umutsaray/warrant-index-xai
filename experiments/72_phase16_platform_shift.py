#!/usr/bin/env python3
"""Phase 16: METABRIC platform shift and location/scale decomposition for Cox G."""
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import joblib,numpy as np,pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase16';OUT.mkdir(parents=True,exist_ok=True);sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
from warrantxai.components import generalizability,generalizability_signed
def loadp3():
 s=importlib.util.spec_from_file_location('p16p3',ROOT/'experiments/06_brca_real.py');m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def metrics(source,coef,anchors,ez,p3):
 ep=np.transpose(coef[:,:,None]*ez.T[None],(0,2,1));gm=generalizability(np.stack([np.abs(source).mean((0,1)),np.abs(ep).mean((0,1))]),.5);mc=coef.mean(0);grid=max(source.shape[1],len(ez));ss=np.column_stack([p3.empirical_rank_grid(mc[j]*anchors[:,j],grid) for j in range(len(mc))]);es=np.column_stack([p3.empirical_rank_grid(mc[j]*ez[:,j],grid) for j in range(len(mc))]);gs=generalizability_signed(np.stack([ss,es]),.5);return gm,gs,np.minimum(gm,gs)
def main():
 p3=loadp3();p3.RAW=ROOT/'data/external/metabric_phase72';train=pd.read_parquet(ROOT/'results/phase7_retrain/X_train_selected.parquet').set_index('sample_key');ext,_,_=p3.external_columns(list(train.columns))
 with np.load(ROOT/'results/phase72_panel/query_matrices.npz') as z:active=z['active_features'].astype(str).tolist();stored=z['external_full'];source_pool=z['source_full']
 raw_all=np.column_stack([ext[f] for f in active]);complete=np.isfinite(raw_all).all(1);raw=raw_all[complete];idx=np.array([train.columns.get_loc(f) for f in active]);tr=train[active].to_numpy(float);tmu=train[active].mean().to_numpy(float);tsd=train[active].std(ddof=0).to_numpy(float);tmin=tr.min(0);tmax=tr.max(0);mmu=raw.mean(0);msd=raw.std(0)
 shift=[]
 for scope,a in [('cox_complete_query_1466',raw),('all_metabric_available_per_feature',raw_all)]:
  for j,f in enumerate(active):
   v=a[:,j];v=v[np.isfinite(v)];shift.append(dict(feature=f,feature_type='RNA' if f.startswith('rna_') else ('age' if f.startswith('clin_numeric__age_') else 'AJCC'),scope=scope,n=len(v),tcga_train_mean=tmu[j],tcga_train_sd=tsd[j],tcga_train_min=tmin[j],tcga_train_max=tmax[j],metabric_raw_mean=v.mean(),metabric_raw_sd=v.std(),standardized_mean_shift=(v.mean()-tmu[j])/tsd[j],sd_ratio=v.std()/tsd[j],outside_tcga_train_range_fraction=np.mean((v<tmin[j])|(v>tmax[j]))))
 pd.DataFrame(shift).to_csv(OUT/'PHASE16_ADIM2_PLATFORM_SHIFT.csv',index=False)
 with np.load(ROOT/'results/phase72_cox/attributions_n50.npz') as z:source=z['source_phi'];coefall=z['coefficients'];features=z['features'].astype(str).tolist()
 ci=np.array([features.index(f) for f in active]);source=source[:,:,ci];coef=coefall[:,ci];cox=joblib.load(ROOT/'results/phase7_retrain/coxph_baseline_v2.joblib');sc=cox['scaler'];anchors=(source_pool-sc.mean_[ci])/sc.scale_[ci]
 aligned=tmu+(raw-mmu)/msd*tsd
 forms={'official_location_and_scale':(aligned-sc.mean_[ci])/sc.scale_[ci],'tcga_train_only':(raw-sc.mean_[ci])/sc.scale_[ci],'location_only':(raw-mmu)/sc.scale_[ci],'scale_only':(raw-sc.mean_[ci])/msd};rep=pd.read_csv(ROOT/'results/phase72_cox/feature_report.csv').set_index('feature').loc[active];rows=[]
 for form,ez in forms.items():
  gm,gs,gc=metrics(source,coef,anchors,ez,p3)
  for j,f in enumerate(active):
   z=np.array([rep.loc[f,'F'],rep.loc[f,'S'],rep.loc[f,'U'],gc[j],rep.loc[f,'R']]);rows.append(dict(feature=f,feature_type='RNA' if f.startswith('rna_') else ('age' if f.startswith('clin_numeric__age_') else 'AJCC'),form=form,G_magnitude=gm[j],G_signed=gs[j],G_claim=gc[j],WI=z.min(),point_adequate=z.min()>.8,bottleneck=['F','S','U','G_claim','R'][z.argmin()]))
 out=pd.DataFrame(rows);old=pd.read_csv(ROOT/'results/phase15/PHASE15_ADIM3_G_STANDARDIZATION.csv');omap={'official_cohort_standardized':'official_location_and_scale','tcga_train_scaler_only_no_cohort_standardization':'tcga_train_only'}
 for oldname,newname in omap.items():
  a=old[old.method.eq(oldname)].set_index('feature').loc[active];b=out[out.form.eq(newname)].set_index('feature').loc[active];np.testing.assert_allclose(b[['G_magnitude','G_signed','G_claim','WI']],a[['G_magnitude','G_signed','G_claim','WI']],atol=1e-12,rtol=0)
 out.to_csv(OUT/'PHASE16_ADIM2_G_FOUR_FORMS.csv',index=False)
if __name__=='__main__':main()
