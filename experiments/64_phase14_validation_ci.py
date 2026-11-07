#!/usr/bin/env python3
"""Phase 14 DeepSurv validation audit and stored-model test C-index CIs."""
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np,pandas as pd,torch
from sklearn.model_selection import train_test_split
from sksurv.metrics import concordance_index_censored
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer'/'train_explainable_dl.py'
s=importlib.util.spec_from_file_location('p14deep',SRC); deep=importlib.util.module_from_spec(s);sys.modules[s.name]=deep;s.loader.exec_module(deep)
def folder(f): return ROOT/'results/phase7_retrain' if f==0 else ROOT/f'results/phase70b_collapse_diagnostic/fold_{f}'
def frame(p): return pd.read_parquet(p).set_index('sample_key')
def boot_ci(e,t,r,groups,B=2000):
 rng=np.random.default_rng(20261107); vals=[]; ug=np.unique(groups)
 loc={g:np.flatnonzero(groups==g) for g in ug}
 for _ in range(B):
  ix=np.concatenate([loc[g] for g in rng.choice(ug,len(ug),replace=True)])
  try: vals.append(concordance_index_censored(e[ix].astype(bool),t[ix],r[ix])[0])
  except ValueError: pass
 return np.quantile(vals,[.025,.975]),len(vals)
def main():
 cv=pd.read_csv(ROOT/'results/phase71_cv/cv_fold_scores.csv'); rows=[]; cis=[]
 for f in range(5):
  d=folder(f); ck=torch.load(d/'explainable_dl_model.pt',map_location='cpu',weights_only=False)
  xtr=frame(d/'X_train_selected.parquet')[ck['ordered_columns']]; ytr=frame(d/'y_train.parquet').loc[xtr.index]
  xte=frame(d/'X_test_selected.parquet')[ck['ordered_columns']]; yte=frame(d/'y_test.parquet').loc[xte.index]
  ix=np.arange(len(xtr)); fit,val=train_test_split(ix,test_size=.15,random_state=42,stratify=ytr.event.to_numpy())
  model=deep.OmicsDeepSurv(ck['block_dims'],dropout=ck['dropout']); model.load_state_dict(ck['state_dict']); model.eval()
  z=(xte.to_numpy(np.float32)-ck['means'])/ck['scales']; risk=deep.predict_risk(model,deep.make_blocks(z,ck['slices'],torch.device('cpu')))
  ci,n=boot_ci(yte.event.to_numpy(),yte.time.to_numpy(),risk,xte.index.astype(str).str[:12].to_numpy())
  fp=set(xtr.iloc[fit].index.astype(str).str[:12]); vp=set(xtr.iloc[val].index.astype(str).str[:12]); tp=set(xte.index.astype(str).str[:12])
  common=dict(fold=f,train_rows=len(xtr),train_patients=len(fp|vp),train_events=int(ytr.event.sum()),test_rows=len(xte),test_patients=len(tp),test_events=int(yte.event.sum()),train_test_patient_overlap=len((fp|vp)&tp))
  for model in ['CoxPH','RSF']:
   score=float(cv.loc[cv.fold.eq(f),{'CoxPH':'coxph_converged_c_index','RSF':'rsf_test_c_index'}[model]].iloc[0])
   rows.append(dict(**common,model=model,fit_rows=np.nan,fit_patients=np.nan,fit_events=np.nan,val_rows=np.nan,val_patients=np.nan,val_events=np.nan,fit_val_patient_overlap=np.nan,best_validation_c_index=np.nan,test_c_index=score,best_epoch=np.nan,stopped_epoch=np.nan))
   cis.append(dict(fold=f,model=model,test_c_index=score,bootstrap_valid_draws=0,ci95_lower=np.nan,ci95_upper=np.nan,status='CI unavailable: fold-specific fitted predictions were not retained'))
  rows.append(dict(**common,model='DeepSurv',fit_rows=len(fit),fit_patients=len(fp),fit_events=int(ytr.iloc[fit].event.sum()),val_rows=len(val),val_patients=len(vp),val_events=int(ytr.iloc[val].event.sum()),fit_val_patient_overlap=len(fp&vp),best_validation_c_index=ck['best_validation_c_index'],test_c_index=ck['test_c_index'],best_epoch=ck['best_epoch'],stopped_epoch=ck['stopped_epoch']))
  cis.append(dict(fold=f,model='DeepSurv',test_c_index=ck['test_c_index'],bootstrap_valid_draws=n,ci95_lower=ci[0],ci95_upper=ci[1],status='patient-cluster bootstrap CI'))
 pd.DataFrame(rows).to_csv(ROOT/'PHASE14_ADIM4_DEEPSURV_VALIDATION.csv',index=False); pd.DataFrame(cis).to_csv(ROOT/'PHASE14_ADIM4_CINDEX_CI.csv',index=False)
if __name__=='__main__':main()
