#!/usr/bin/env python3
"""Phase 15: score frozen DeepSurv checkpoints on overlap-stratified validation rows."""
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np,pandas as pd,torch
from scipy.stats import norm
from sklearn.model_selection import train_test_split
from sksurv.metrics import concordance_index_censored
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/phase15';OUT.mkdir(parents=True,exist_ok=True)
SRC=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer'/'train_explainable_dl.py'
s=importlib.util.spec_from_file_location('p15deep',SRC);deep=importlib.util.module_from_spec(s);sys.modules[s.name]=deep;s.loader.exec_module(deep)
def folder(f):return ROOT/'results/phase7_retrain' if f==0 else ROOT/f'results/phase70b_collapse_diagnostic/fold_{f}'
def frame(p):return pd.read_parquet(p).set_index('sample_key')
def cindex(y,r):return float(concordance_index_censored(y.event.to_numpy(bool),y.time.to_numpy(float),r)[0])
def patient_boot(y,r,B=2000):
 groups=y.index.astype(str).str[:12].to_numpy();ug=np.unique(groups);loc={g:np.flatnonzero(groups==g) for g in ug};rng=np.random.default_rng(20261108);v=[]
 for _ in range(B):
  ix=np.concatenate([loc[g] for g in rng.choice(ug,len(ug),replace=True)])
  try:v.append(cindex(y.iloc[ix],r[ix]))
  except ValueError:pass
 return np.quantile(v,[.025,.975]),len(v)
def main():
 rows=[]
 for f in range(5):
  d=folder(f);ck=torch.load(d/'explainable_dl_model.pt',map_location='cpu',weights_only=False)
  x=frame(d/'X_train_selected.parquet')[ck['ordered_columns']];y=frame(d/'y_train.parquet').loc[x.index]
  xt=frame(d/'X_test_selected.parquet')[ck['ordered_columns']];yt=frame(d/'y_test.parquet').loc[xt.index]
  fi,vi=train_test_split(np.arange(len(x)),test_size=.15,random_state=42,stratify=y.event.to_numpy())
  m=deep.OmicsDeepSurv(ck['block_dims'],dropout=ck['dropout']);m.load_state_dict(ck['state_dict']);m.eval()
  def risk(a):
   z=(a.to_numpy(np.float32)-ck['means'])/ck['scales'];return deep.predict_risk(m,deep.make_blocks(z,ck['slices'],torch.device('cpu')))
  rv=risk(x.iloc[vi]);rt=risk(xt);yv=y.iloc[vi];fit_pat=set(x.iloc[fi].index.astype(str).str[:12]);vp=yv.index.astype(str).str[:12]
  non=~vp.isin(fit_pat);over=~non
  full=cindex(yv,rv);test=cindex(yt,rt)
  if abs(full-ck['best_validation_c_index'])>1e-6 or abs(test-ck['test_c_index'])>1e-6:raise RuntimeError(f'fold {f} checkpoint score mismatch: {full}, {test}')
  for subset,mask in [('all',np.ones(len(yv),bool)),('no_fit_patient_overlap',np.asarray(non)),('overlap_only',np.asarray(over))]:
   yy=yv.iloc[np.flatnonzero(mask)];rr=rv[mask];score=cindex(yy,rr) if len(yy)>1 and yy.event.sum()>0 else np.nan
   lo=hi=np.nan;draws=0
   if subset=='no_fit_patient_overlap' and np.isfinite(score):(lo,hi),draws=patient_boot(yy,rr)
   rows.append(dict(fold=f,subset=subset,rows=len(yy),patients=yy.index.astype(str).str[:12].nunique(),events=int(yy.event.sum()),c_index=score,ci95_lower=lo,ci95_upper=hi,bootstrap_draws=draws,checkpoint_val_c=ck['best_validation_c_index'],test_c_index=test,fit_val_overlap_patients=len(fit_pat&set(vp)),epoch_log_status='yapılamadı: epoch-bazlı eğitim günlüğü checkpoint içinde saklı değil'))
 pd.DataFrame(rows).to_csv(OUT/'PHASE15_ADIM1_DEEPSURV_VAL_NOOVERLAP.csv',index=False)
if __name__=='__main__':main()
