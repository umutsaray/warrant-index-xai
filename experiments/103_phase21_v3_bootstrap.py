"""Phase 21 v3 checkpointed N=50 attribution bootstrap."""
from __future__ import annotations
import argparse, importlib.util, json, sys, time, warnings
from pathlib import Path
import joblib, numpy as np, pandas as pd, torch
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceWarning
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'results/phase7_retrain';LOCAL=ROOT/'results/phase21_v3';OUT=ROOT/'phase21_reviewer1/v3';SIB=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer';SEED=20260917
sys.path[:0]=[str(ROOT),str(ROOT/'experiments'),str(SIB)]
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;assert s.loader;s.loader.exec_module(m);return m
def data():
 cache=SIB/'phase1_fold_cache';xf=pd.read_parquet(cache/'97cc92aa9eff44f1_train.parquet');xr=pd.read_parquet(cache/'97cc92aa9eff44f1_test.parquet');y=pd.read_parquet(DATA/'y_train.parquet').set_index('sample_key');x=pd.concat([xf,xr.loc[y.index.intersection(xr.index)]]).loc[y.index];return xf,x,y
def setup():
 xf,x,y=data();q=np.load(LOCAL/'query_matrices.npz',allow_pickle=True);active=q['active_features'].astype(str).tolist();queries=[q['source_query12'],q['external_query12']];features=joblib.load(LOCAL/'rsf.joblib')['features'];x=x[features]
 pats=x.index.astype(str).str[:12];groups=[np.flatnonzero(pats==p) for p in pats.unique()];rng=np.random.default_rng(SEED);indices=[np.concatenate([groups[i] for i in rng.integers(0,len(groups),len(groups))]) for _ in range(50)]
 uniq=x.loc[~pats.duplicated()];br=np.random.default_rng(SEED+4200);first=np.sort(br.choice(len(uniq),16,False));rem=np.setdiff1d(np.arange(len(uniq)),first);extra=np.sort(br.choice(rem,112,False));nested=np.r_[first,extra];ai=np.array([features.index(f) for f in active]);return x,y,active,queries,indices,uniq,nested,ai
def run_rsf(start_rep=0, stop_rep=50):
 x,y,active,queries,indices,uniq,nested,ai=setup();rh=load('07_rsf_pilot.py','p21v3_rsf');cp=LOCAL/'checkpoints_rsf';cp.mkdir(exist_ok=True);am=x.mean().to_numpy()[ai];bg=uniq.iloc[nested[:16]].to_numpy(np.float32);surv=Surv.from_arrays(event=y.event.astype(bool),time=y.time);records=[]
 for rep in range(start_rep,stop_rep):
  ix=indices[rep]
  p=cp/f'rep{rep:02d}.npz';m=p.with_suffix('.json')
  if p.exists() and m.exists():records.append(json.loads(m.read_text()));print('RSF checkpoint',rep+1,flush=True);continue
  start=time.perf_counter();model=RandomSurvivalForest(n_estimators=500,min_samples_split=10,min_samples_leaf=5,max_features='sqrt',random_state=SEED+rep,n_jobs=6,low_memory=True);model.fit(x.to_numpy(np.float32)[ix],surv[ix]);fit=time.perf_counter()-start;c=[];r=[];t=time.perf_counter()
  for site,q in enumerate(queries):
   c.append(rh.background_permutation_attributions(model,q,bg,am,ai,permutations=8,seed=SEED+10000+100*site+rep));r.append(rh.background_permutation_attributions(model,q,bg,am,ai,permutations=128,seed=SEED+20000+100*site+rep))
  rec={'replicate':rep,'fit_seconds':fit,'attribution_seconds':time.perf_counter()-t,'rows':len(ix),'events':int(y.event.to_numpy()[ix].sum()),'warnings':[]};np.savez_compressed(p,candidate=np.stack(c),reference=np.stack(r));m.write_text(json.dumps(rec,indent=2));records.append(rec);print('RSF',rep+1,'/50',rec,flush=True)
 if all((cp/f'rep{rep:02d}.npz').exists() and (cp/f'rep{rep:02d}.json').exists() for rep in range(50)):
  cand=[];ref=[];allrecords=[]
  for rep in range(50):
   with np.load(cp/f'rep{rep:02d}.npz') as z:cand.append(z['candidate']);ref.append(z['reference'])
   allrecords.append(json.loads((cp/f'rep{rep:02d}.json').read_text()))
  np.savez_compressed(LOCAL/'rsf_attributions_n50.npz',candidate=np.stack(cand,axis=1),reference=np.stack(ref,axis=1),features=np.asarray(active));(OUT/'RSF_BOOTSTRAP.md').write_text('# RSF-v3 N=50\n\n50/50 checkpoint tamamlandı.\n');(OUT/'rsf_bootstrap_summary.json').write_text(json.dumps(allrecords,indent=2))
def run_deep():
 x,y,active,queries,indices,uniq,nested,ai=setup();p6=load('20_deepsurv_pilot.py','p21v3_p6');p61=load('21_deepsurv_explainer_revision.py','p21v3_p61');p72=load('41_phase72_deepsurv.py','p21v3_p72');ck=torch.load(LOCAL/'deepsurv.pt',map_location='cpu',weights_only=False);x=x[ck['ordered_columns']];raw=x.to_numpy(np.float32);cp=LOCAL/'checkpoints_deep';cp.mkdir(exist_ok=True);bgs={k:uniq[ck['ordered_columns']].iloc[nested[:k]].to_numpy(np.float32) for k in (16,128)};records=[]
 for rep,ix in enumerate(indices):
  p=cp/f'rep{rep:02d}.npz';m=p.with_suffix('.json')
  if p.exists() and m.exists():records.append(json.loads(m.read_text()));print('Deep checkpoint',rep+1,flush=True);continue
  model,means,scales,trainrec=p6.train_bootstrap_model(ck,raw,y.event.to_numpy(bool),y.time.to_numpy(np.float32),ix,rep);c=[];r=[];sites=[]
  for site,q in enumerate(queries):
   cv=p72.deep_attributions_with_relative_gate(p61,p6,model,ck,q,bgs[16],ai,means,scales);rv=p72.deep_attributions_with_relative_gate(p61,p6,model,ck,q,bgs[128],ai,means,scales);c.append(cv[0]);r.append(rv[0]);sites.append({'site':site,'candidate_error':cv[1],'reference_error':rv[1],'candidate_seconds':cv[3],'reference_seconds':rv[3]})
  rec={'replicate':rep,'training':trainrec,'sites':sites};np.savez_compressed(p,candidate=np.stack(c),reference=np.stack(r));m.write_text(json.dumps(rec,indent=2));records.append(rec);print('Deep',rep+1,'/50',flush=True)
 cand=[];ref=[]
 for rep in range(50):
  with np.load(cp/f'rep{rep:02d}.npz') as z:cand.append(z['candidate']);ref.append(z['reference'])
 np.savez_compressed(LOCAL/'deepsurv_attributions_n50.npz',candidate=np.stack(cand,axis=1),reference=np.stack(ref,axis=1),features=np.asarray(active));(OUT/'DEEPSURV_BOOTSTRAP.md').write_text('# DeepSurv-v3 N=50\n\n50/50 checkpoint tamamlandı.\n');(OUT/'deepsurv_bootstrap_summary.json').write_text(json.dumps(records,indent=2))
def run_cox():
 x,y,active,queries,indices,uniq,nested,ai=setup();bundle=joblib.load(LOCAL/'coxph.joblib');features=bundle['features'];sc=bundle['scaler'];a=sc.transform(x[features]);cp=LOCAL/'checkpoints_cox';cp.mkdir(exist_ok=True);coefs=[];records=[]
 for rep,ix in enumerate(indices):
  p=cp/f'rep{rep:02d}.npz';m=p.with_suffix('.json')
  if p.exists() and m.exists():
   with np.load(p) as z:coefs.append(z['coefficients']);records.append(json.loads(m.read_text()));continue
  vals=a[ix];ev=y.event.to_numpy()[ix].astype(bool);valid=(vals[ev].var(0)>1e-12)&(vals[~ev].var(0)>1e-12);frame=pd.DataFrame(vals[:,valid],columns=np.asarray(features)[valid]);frame['time']=y.time.to_numpy()[ix];frame['event']=ev;model=CoxPHFitter(penalizer=.1)
  with warnings.catch_warnings(record=True) as caught:warnings.simplefilter('always');model.fit(frame,'time','event',fit_options={'max_steps':200,'step_size':.25})
  conv=[w for w in caught if issubclass(w.category,ConvergenceWarning)]
  if conv:raise RuntimeError(f'Cox rep {rep} convergence warning')
  coef=np.zeros(len(features));coef[valid]=model.params_.to_numpy();rec={'replicate':rep,'excluded':np.asarray(features)[~valid].tolist(),'warnings':[]};np.savez_compressed(p,coefficients=coef);m.write_text(json.dumps(rec,indent=2));coefs.append(coef);records.append(rec);print('Cox',rep+1,'/50',flush=True)
 coefs=np.asarray(coefs);qs=[(q-sc.mean_[ai])/sc.scale_[ai] for q in queries];phi=np.stack([coefs[:,None,ai]*q[None] for q in qs]);np.savez_compressed(LOCAL/'cox_attributions_n50.npz',candidate=phi,reference=phi,features=np.asarray(active),coefficients=coefs);(OUT/'COX_BOOTSTRAP.md').write_text('# CoxPH-v3 N=50\n\n50/50 checkpoint tamamlandı.\n');(OUT/'cox_bootstrap_summary.json').write_text(json.dumps(records,indent=2))
def main():
 p=argparse.ArgumentParser();p.add_argument('--model',choices=['cox','rsf','deep'],required=True);p.add_argument('--start',type=int,default=0);p.add_argument('--stop',type=int,default=50);a=p.parse_args()
 if a.model=='rsf':run_rsf(a.start,a.stop)
 else:{'cox':run_cox,'deep':run_deep}[a.model]()
if __name__=='__main__':main()
