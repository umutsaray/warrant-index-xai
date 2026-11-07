#!/usr/bin/env python3
"""Phase 15: reconstruct the two KRTAP11-1 residual RSF replicates."""
from __future__ import annotations
import importlib.util,sys,time
from pathlib import Path
import joblib,numpy as np,pandas as pd
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase15';OUT.mkdir(parents=True,exist_ok=True);DATA=ROOT/'results/phase7_retrain';SEED=20260917
def loadpilot():
 s=importlib.util.spec_from_file_location('p15pilot',ROOT/'experiments/07_rsf_pilot.py');m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def frame(n):return pd.read_parquet(DATA/n).set_index('sample_key')
def main():
 focus='rna_ENSG00000182591.6';targets={3,5};pilot=loadpilot();train=frame('X_train_selected.parquet');y=frame('y_train.parquet').loc[train.index];features=joblib.load(DATA/'rsf_baseline.joblib')['features']
 with np.load(ROOT/'results/phase72_panel/query_matrices.npz') as z:active=z['active_features'].astype(str).tolist();query=z['external_query12']
 with np.load(ROOT/'results/phase72_rsf/attributions_n50.npz') as z:stored_c=z['candidate'];stored_r=z['reference']
 ai=np.array([features.index(f) for f in active]);j=active.index(focus);am=train.mean().to_numpy()[ai];patients=train.index.astype(str).str[:12];unique=train.loc[~patients.duplicated()]
 rng=np.random.default_rng(SEED+4200);bg=unique.iloc[np.sort(rng.choice(len(unique),16,replace=False))].to_numpy(np.float32);bg[:,ai]=am.astype(np.float32)
 groups=[np.flatnonzero(patients==p) for p in patients.unique()];brng=np.random.default_rng(SEED);indices={}
 for rep in range(6):
  ix=np.concatenate([groups[k] for k in brng.integers(0,len(groups),len(groups))])
  if rep in targets:indices[rep]=ix
 surv=Surv.from_arrays(y.event.astype(bool),y.time.to_numpy());values=train.to_numpy(np.float32);rows=[]
 for rep in sorted(targets):
  m=RandomSurvivalForest(n_estimators=500,min_samples_split=10,min_samples_leaf=5,max_features='sqrt',random_state=SEED+rep,n_jobs=-1,low_memory=True).fit(values[indices[rep]],surv[indices[rep]])
  ref=pilot.background_permutation_attributions(m,query,bg,am,ai,permutations=128,seed=SEED+20000+100+rep)
  gj=ai[j];tree_thresholds=[tree.tree_.threshold[tree.tree_.feature==gj] for tree in m.estimators_];base=am[j]
  for qi,qv in enumerate(query[:,j]):
   per_tree=[int(np.any((qv<=th)!=(base<=th))) for th in tree_thresholds]
   rows.append(dict(replicate=rep,query=qi,feature=focus,base_value=base,query_value=qv,trees_with_separating_split=sum(per_tree),total_feature_split_nodes=sum(map(len,tree_thresholds)),candidate8=stored_c[1,rep,qi,j],reference128_recomputed=ref[qi,j],reference128_stored=stored_r[1,rep,qi,j],reference_matches_stored_atol_1e_12=abs(ref[qi,j]-stored_r[1,rep,qi,j])<=1e-12))
 pd.DataFrame(rows).to_csv(OUT/'PHASE15_ADIM5_KRTAP11_1_RESIDUAL.csv',index=False)
if __name__=='__main__':main()
