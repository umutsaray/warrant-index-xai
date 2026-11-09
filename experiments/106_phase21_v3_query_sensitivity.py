"""Phase 21 v3 preregistered three-seed point-only query diagnostic."""
from __future__ import annotations
import importlib.util,json,sys
from pathlib import Path
import joblib,numpy as np,pandas as pd,torch
ROOT=Path(__file__).resolve().parents[1];LOCAL=ROOT/'results/phase21_v3';OUT=ROOT/'phase21_reviewer1/v3';DATA=ROOT/'results/phase7_retrain';SIB=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer';sys.path[:0]=[str(ROOT),str(ROOT/'experiments'),str(SIB)]
from warrantxai.components import generalizability,generalizability_signed
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;assert s.loader;s.loader.exec_module(m);return m
def g(src,ext):
 gm=generalizability(np.stack([np.abs(src).mean(0),np.abs(ext).mean(0)]),.5);gs=generalizability_signed(np.stack([np.sort(src,axis=0),np.sort(ext,axis=0)]),.5);return np.minimum(gm,gs)
def main():
 cache=SIB/'phase1_fold_cache';xr=pd.read_parquet(cache/'97cc92aa9eff44f1_test.parquet');yt=pd.read_parquet(DATA/'y_train.parquet').set_index('sample_key');ye=pd.read_parquet(DATA/'y_test.parquet').set_index('sample_key');xf=pd.read_parquet(cache/'97cc92aa9eff44f1_train.parquet');xtrain=pd.concat([xf,xr.loc[yt.index.intersection(xr.index)]]).loc[yt.index];xtest=xr.loc[ye.index]
 with np.load(LOCAL/'query_matrices.npz',allow_pickle=True) as z:active=z['active_features'].astype(str).tolist();extfull=z['external_full']
 up=xtest.loc[~xtest.index.astype(str).str[:12].duplicated()];pats=xtrain.index.astype(str).str[:12];uniq=xtrain.loc[~pats.duplicated()];br=np.random.default_rng(20260917+4200);first=np.sort(br.choice(len(uniq),16,False));rows=[]
 rh=load('07_rsf_pilot.py','p21q_rsf');rb=joblib.load(LOCAL/'rsf.joblib');rf=rb['features'];ri=np.array([rf.index(f) for f in active]);rbg=uniq[rf].iloc[first].to_numpy(np.float32);ram=xtrain[rf].mean().to_numpy()[ri]
 p6=load('20_deepsurv_pilot.py','p21q_p6');p61=load('21_deepsurv_explainer_revision.py','p21q_p61');p72=load('41_phase72_deepsurv.py','p21q_p72');ck=torch.load(LOCAL/'deepsurv.pt',map_location='cpu',weights_only=False);dm=p6.model_from_state(ck);di=np.array([ck['ordered_columns'].index(f) for f in active]);dbg=uniq[ck['ordered_columns']].iloc[first].to_numpy(np.float32)
 for model in ['RSF-v3','DeepSurv-v3']:
  base=pd.read_csv(OUT/f'{model.lower().replace("-","_")}_feature_certification.csv').set_index('feature').loc[active]
  for seed in [20264917,20264918,20264919]:
   rng=np.random.default_rng(seed);src=up[active].to_numpy()[np.sort(rng.choice(len(up),12,False))];ext=extfull[np.sort(rng.choice(len(extfull),12,False))]
   if model=='RSF-v3':
    a=rh.background_permutation_attributions(rb['model'],src,rbg,ram,ri,permutations=8,seed=seed+10000);b=rh.background_permutation_attributions(rb['model'],ext,rbg,ram,ri,permutations=8,seed=seed+10100)
   else:
    a=p72.deep_attributions_with_relative_gate(p61,p6,dm,ck,src,dbg,di,np.asarray(ck['means'],np.float32),np.asarray(ck['scales'],np.float32))[0];b=p72.deep_attributions_with_relative_gate(p61,p6,dm,ck,ext,dbg,di,np.asarray(ck['means'],np.float32),np.asarray(ck['scales'],np.float32))[0]
   gc=g(a,b)
   for j,f in enumerate(active):
    vals=np.array([base.loc[f,'F'],base.loc[f,'S'],base.loc[f,'U'],gc[j],base.loc[f,'R']]);rows.append({'model':model,'seed':seed,'feature':f,'G_claim_reference_model':gc[j],'diagnostic_WI':vals.min(),'point_adequate':vals.min()>=.8,'bottleneck':['F','S','U','G_claim','R'][int(vals.argmin())],'scope':'single reference model G; frozen N50 F/S/U/R; point-only'})
 result=pd.DataFrame(rows);result.to_csv(OUT/'query_seed_sensitivity_features.csv',index=False);summary=result.groupby(['model','seed']).point_adequate.sum().rename('point').reset_index();ranges=summary.groupby('model').point.agg(['min','max']).reset_index();summary=summary.merge(ranges,on='model');summary.to_csv(OUT/'query_seed_sensitivity.csv',index=False);(OUT/'QUERY_SEED_SENSITIVITY.md').write_text('# Faz 21 v3 sorgu-seed duyarlılığı\n\nBu ön-kayıtlı tanı yalnız nokta düzeyindedir. Bootstrap modelleri saklanmadığından üç seed tek referans model G hesabıyla, N=50 F/S/U/R sabit tutularak değerlendirilmiştir; resmi inferansı değiştirmez.\n\n'+summary.to_string(index=False)+'\n',encoding='utf-8');print(summary.to_string(index=False))
if __name__=='__main__':main()
