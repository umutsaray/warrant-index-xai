#!/usr/bin/env python3
"""Phase 14 deterministic reconstruction of frozen RSF bootstrap forests."""
from __future__ import annotations
import argparse,importlib.util,sys,time
from pathlib import Path
import joblib,numpy as np,pandas as pd
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv
ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'results/phase7_retrain'; SEED=20260917
def loadpilot():
 s=importlib.util.spec_from_file_location('p14rsfpilot',ROOT/'experiments/07_rsf_pilot.py');m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def frame(n):return pd.read_parquet(DATA/n).set_index('sample_key')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--replicates',type=int,default=3);a=ap.parse_args(); pilot=loadpilot(); start=time.perf_counter()
 train=frame('X_train_selected.parquet'); y=frame('y_train.parquet').loc[train.index]; saved=joblib.load(DATA/'rsf_baseline.joblib');features=saved['features']
 with np.load(ROOT/'results/phase72_panel/query_matrices.npz') as z: active=z['active_features'].astype(str).tolist();queries=[z['source_query12'],z['external_query12']]
 with np.load(ROOT/'results/phase72_rsf/attributions_n50.npz') as z: stored=z['candidate']
 ai=np.array([features.index(f) for f in active]); am=train.mean().to_numpy()[ai]; patients=train.index.astype(str).str[:12]; unique=train.loc[~patients.duplicated()]
 rng=np.random.default_rng(SEED+4200); sel=np.sort(rng.choice(len(unique),16,replace=False)); bg=unique.iloc[sel].to_numpy(np.float32);bg[:,ai]=am.astype(np.float32)
 groups=[np.flatnonzero(patients==p) for p in patients.unique()]; brng=np.random.default_rng(SEED); values=train.to_numpy(np.float32);surv=Surv.from_arrays(y.event.astype(bool),y.time.to_numpy())
 repro=[];trav=[]
 for rep in range(a.replicates):
  ix=np.concatenate([groups[i] for i in brng.integers(0,len(groups),len(groups))]); m=RandomSurvivalForest(n_estimators=500,min_samples_split=10,min_samples_leaf=5,max_features='sqrt',random_state=SEED+rep,n_jobs=-1,low_memory=True);m.fit(values[ix],surv[ix])
  for site in range(2):
   got=pilot.background_permutation_attributions(m,queries[site],bg,am,ai,permutations=8,seed=SEED+10000+100*site+rep); diff=np.abs(got-stored[site,rep])
   repro.append(dict(replicate=rep,site=site,max_abs_difference=diff.max(),mean_abs_difference=diff.mean(),matching_atol_1e_12=bool(np.all(diff<=1e-12)),worst_feature=active[np.unravel_index(diff.argmax(),diff.shape)[1]],worst_query=np.unravel_index(diff.argmax(),diff.shape)[0]))
  for j,(f,gj) in enumerate(zip(active,ai)):
   th=[]
   for tree in m.estimators_: th.extend(tree.tree_.threshold[tree.tree_.feature==gj])
   th=np.asarray(th); base=am[j]
   for site,name in enumerate(['TCGA','METABRIC']):
    q=queries[site][:,j]; separated=np.array([np.sum((v<=th)!=(base<=th)) for v in q]) if len(th) else np.zeros(12,int)
    trav.append(dict(replicate=rep,site=name,feature=f,split_nodes=len(th),queries_with_separating_split=int((separated>0).sum()),separating_nodes_total=int(separated.sum()),separating_ratio_mean=float(np.mean(separated/len(th))) if len(th) else 0,max_abs_stored_candidate=float(np.abs(stored[site,rep,:,j]).max())))
  print('rep',rep,'done',flush=True)
 pd.DataFrame(repro).to_csv(ROOT/'PHASE14_ADIM7_BOOTSTRAP_REPRO.csv',index=False);pd.DataFrame(trav).to_csv(ROOT/'PHASE14_ADIM7_BOOTSTRAP_TRAVERSAL.csv',index=False)
 print('elapsed',time.perf_counter()-start,'all_match',all(x['matching_atol_1e_12'] for x in repro))
if __name__=='__main__':main()
