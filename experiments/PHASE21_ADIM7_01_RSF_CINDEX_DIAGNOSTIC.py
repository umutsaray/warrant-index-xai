"""Phase 21 ADIM 7.1: direction, bootstrap-CI, and seed diagnostics."""
from __future__ import annotations
import importlib.util,json,sys
from pathlib import Path
import joblib,numpy as np,pandas as pd,torch
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'results/phase7_retrain';LOCAL=ROOT/'results/phase21_v3';OUT=ROOT/'phase21_reviewer1';SIB=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer';sys.path[:0]=[str(SIB),str(ROOT)]
def deep_module():
 p=SIB/'train_explainable_dl.py';s=importlib.util.spec_from_file_location('p21a7deep',p);m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;assert s.loader;s.loader.exec_module(m);return m
def cidx(y,r):return float(concordance_index_censored(y.event.astype(bool).to_numpy(),y.time.to_numpy(),np.asarray(r))[0])
def cache_data():
 c=SIB/'phase1_fold_cache';xf=pd.read_parquet(c/'97cc92aa9eff44f1_train.parquet');xr=pd.read_parquet(c/'97cc92aa9eff44f1_test.parquet');yt=pd.read_parquet(DATA/'y_train.parquet').set_index('sample_key');ye=pd.read_parquet(DATA/'y_test.parquet').set_index('sample_key');xtr=pd.concat([xf,xr.loc[yt.index.intersection(xr.index)]]).loc[yt.index];return xf,xr.loc[yt.index.intersection(xr.index)],xtr,xr.loc[ye.index],yt,ye
def deep_risk(bundle_path,x):
 m=deep_module();ck=torch.load(bundle_path,map_location='cpu',weights_only=False);model=m.OmicsDeepSurv(ck['block_dims'],ck['dropout']);model.load_state_dict(ck['state_dict']);model.eval();a=(x[ck['ordered_columns']].to_numpy(np.float32)-ck['means'])/ck['scales'];return m.predict_risk(model,m.make_blocks(a,ck['slices'],torch.device('cpu')))
def patient_boot(y,risks,seed=20262107,b=1000):
 p=y.index.astype(str).str[:12];u=pd.unique(p);groups={q:np.flatnonzero(p==q) for q in u};rng=np.random.default_rng(seed);rows=[]
 for name,r in risks.items():
  vals=[]
  for _ in range(b):
   ix=np.concatenate([groups[u[j]] for j in rng.integers(0,len(u),len(u))])
   try:vals.append(float(concordance_index_censored(y.event.to_numpy()[ix].astype(bool),y.time.to_numpy()[ix],np.asarray(r)[ix])[0]))
   except ValueError:pass
  rows.append({'model':name,'point_cindex':cidx(y,r),'bootstrap_valid':len(vals),'ci95_lower':np.quantile(vals,.025),'ci95_upper':np.quantile(vals,.975),'bootstrap_seed':seed})
 return pd.DataFrame(rows)
def main():
 xf,xv,xtr,xte,yt,ye=cache_data();v3=joblib.load(LOCAL/'rsf.joblib');f=v3['features'];model=v3['model'];sets={'inner_train_577':xf,'inner_validation_97':xv,'outer_train_674':xtr,'outer_test_168':xte};ys={'inner_train_577':yt.loc[xf.index],'inner_validation_97':yt.loc[xv.index],'outer_train_674':yt,'outer_test_168':ye};rows=[]
 for scope,x in sets.items():
  r=model.predict(x[f].to_numpy(np.float32));rows.append({'model':'RSF-v3','scope':scope,'rows':len(x),'events':int(ys[scope].event.sum()),'cindex_risk':cidx(ys[scope],r),'cindex_negated_risk':cidx(ys[scope],-r),'risk_definition':'model.predict; higher=worse'})
 direction=pd.DataFrame(rows);direction.to_csv(OUT/'PHASE21_ADIM7_01_RSF_DIRECTION.csv',index=False)
 # v3 test risks
 cb=joblib.load(LOCAL/'coxph.joblib');cz=cb['scaler'].transform(xte[cb['features']]);cr=cb['model'].predict_partial_hazard(pd.DataFrame(cz,columns=cb['features'])).to_numpy();rr=model.predict(xte[f].to_numpy(np.float32));dr=deep_risk(LOCAL/'deepsurv.pt',xte)
 # v1 clean fold-0 risks
 x1tr=pd.read_parquet(DATA/'X_train_selected.parquet').set_index('sample_key');x1te=pd.read_parquet(DATA/'X_test_selected.parquet').set_index('sample_key');c1=joblib.load(DATA/'coxph_baseline_v2.joblib');z1=c1['scaler'].transform(x1te[c1['features']]);c1r=c1['model'].predict_partial_hazard(pd.DataFrame(z1,columns=c1['features'])).to_numpy();r1b=joblib.load(DATA/'rsf_baseline.joblib');r1=r1b['model'].predict(x1te[r1b['features']].to_numpy(np.float32));d1=deep_risk(DATA/'explainable_dl_model.pt',x1te)
 ci=patient_boot(ye,{'CoxPH-v3':cr,'RSF-v3':rr,'DeepSurv-v3':dr,'CoxPH-v1':c1r,'RSF-v1':r1,'DeepSurv-v1':d1});ci.to_csv(OUT/'PHASE21_ADIM7_02_CINDEX_BOOTSTRAP_CI.csv',index=False)
 surv=Surv.from_arrays(event=yt.event.astype(bool),time=yt.time);seedrows=[]
 for seed in (43,44,45):
  m=RandomSurvivalForest(n_estimators=500,min_samples_split=10,min_samples_leaf=5,max_features='sqrt',random_state=seed,n_jobs=-1,low_memory=True);m.fit(xtr[f].to_numpy(np.float32),surv);risk=m.predict(xte[f].to_numpy(np.float32));seedrows.append({'random_state':seed,'test_cindex_risk':cidx(ye,risk),'test_cindex_negated_risk':cidx(ye,-risk)})
 seeds=pd.DataFrame(seedrows);seeds.to_csv(OUT/'PHASE21_ADIM7_03_RSF_SEED_DIAGNOSTIC.csv',index=False)
 bug=not np.isclose(direction.loc[direction.scope=='outer_test_168','cindex_risk'].iloc[0],.38443579766536967)
 report=f'''# Faz 21 ADIM 7.1 — RSF-v3 C-index doğrulaması

## Yön
```
{direction.to_string(index=False)}
```

## Hasta-bootstrap %95 GA
```
{ci.to_string(index=False)}
```

## Alternatif RSF seed'leri
```
{seeds.to_string(index=False)}
```

Kod yolu: ortak skor fonksiyonu bu dosyada satır 12'de `concordance_index_censored(event,time,risk)` çağrısıdır. v3 ve v1 riskleri sırasıyla satır 33 ve 35'te üretilip satır 36'da aynı fonksiyona verilir. İlk v3 Cox/RSF yolu `experiments/101_phase21_v3_stage1.py:194-202`; DeepSurv yolu aynı dosyada `train_deep` içindeki satır 139'dur. RSF için her iki sürümde de `model.predict` riski kullanılır.

Yön/kod uyuşmazlığı: {'EVET' if bug else 'HAYIR'}.
''';(OUT/'PHASE21_ADIM7_01_RSF_CINDEX_DIAGNOSTIC.md').write_text(report,encoding='utf-8');print(report)
 if bug:raise SystemExit(2)
if __name__=='__main__':main()
