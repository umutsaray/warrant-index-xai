"""Phase 21 v3 resource gate: build 65-feature queries and time one replicate/model."""
from __future__ import annotations
import importlib.util, json, sys, time, warnings
from pathlib import Path
import joblib, numpy as np, pandas as pd, torch
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv

ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'results/phase7_retrain'; LOCAL=ROOT/'results/phase21_v3'; OUT=ROOT/'phase21_reviewer1/v3'
SIB=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer'; SEED=20260917
sys.path[:0]=[str(ROOT),str(ROOT/'experiments'),str(SIB)]
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;assert s.loader;s.loader.exec_module(m);return m
def matrices():
 sig='97cc92aa9eff44f1'; cache=SIB/'phase1_fold_cache'; xf=pd.read_parquet(cache/f'{sig}_train.parquet'); xr=pd.read_parquet(cache/f'{sig}_test.parquet')
 yt=pd.read_parquet(DATA/'y_train.parquet').set_index('sample_key'); ye=pd.read_parquet(DATA/'y_test.parquet').set_index('sample_key')
 val=set(json.loads((OUT/'stage1_metadata.json').read_text())['leakage'].values()) # force metadata presence
 # cache test is remainder in all-y order; select by keys.
 xtrain=pd.concat([xf,xr.loc[yt.index.intersection(xr.index)]]).loc[yt.index]; xtest=xr.loc[ye.index]
 return xf,xtrain,xtest,yt,ye
def external(active,xf):
 audit=json.loads((LOCAL/'audit_features.json').read_text()); fmap=audit['rna_symbols']; symbols={f:fmap[f] for f in active if f.startswith('rna_')}
 raw=ROOT/'data/external/metabric_phase21/data_mrna_illumina_microarray.txt'; wanted=set(symbols.values()); chunks=[]
 for c in pd.read_csv(raw,sep='\t',chunksize=256,low_memory=False):
  z=c[c.Hugo_Symbol.astype(str).isin(wanted)]
  if len(z):chunks.append(z)
 expr=pd.concat(chunks).drop_duplicates('Hugo_Symbol').set_index('Hugo_Symbol'); samples=expr.columns[2:].tolist() if 'Entrez_Gene_Id' in expr.columns else expr.columns[1:].tolist()
 # after set_index, columns are Entrez_Gene_Id + samples
 samples=[c for c in expr.columns if c!='Entrez_Gene_Id']
 old=ROOT/'data/external/metabric_phase72'; si=pd.DataFrame(json.loads((old/'samples.json').read_text())).set_index('sampleId')
 pr=pd.DataFrame(json.loads((old/'patient_clinical.json').read_text())).pivot(index='patientId',columns='clinicalAttributeId',values='value')
 sr=pd.DataFrame(json.loads((old/'sample_clinical.json').read_text())).pivot(index='sampleId',columns='clinicalAttributeId',values='value')
 linked=si.reindex(samples).patientId; tm=pd.to_numeric(pr.reindex(linked).OS_MONTHS,errors='coerce').set_axis(samples); st=pr.reindex(linked).OS_STATUS.set_axis(samples)
 keep=tm.notna()&tm.gt(0)&st.isin(['0:LIVING','1:DECEASED']); samples=list(np.asarray(samples)[keep.to_numpy()]); linked=si.reindex(samples).patientId
 cols=[]
 for f in active:
  if f.startswith('rna_'): cols.append(pd.to_numeric(expr.loc[symbols[f],samples],errors='coerce').to_numpy(float))
  elif 'age_at_index' in f: cols.append(pd.to_numeric(pr.reindex(linked).AGE_AT_DIAGNOSIS,errors='coerce').to_numpy(float))
  else:
   code={'Stage I':1,'Stage II':2,'Stage III':3,'Stage IV':4}[next(x for x in ['Stage I','Stage II','Stage III','Stage IV'] if f.endswith('_'+x))]
   v=pd.to_numeric(sr.reindex(samples).TUMOR_STAGE,errors='coerce').to_numpy(float);cols.append(np.where(np.isfinite(v),(v==code).astype(float),np.nan))
 arr=np.column_stack(cols); complete=np.isfinite(arr).all(1);arr=arr[complete]
 mu=xf[active].mean().to_numpy();sd=xf[active].std(ddof=0).to_numpy(); es=arr.std(0); 
 if np.any(es==0):raise RuntimeError('constant external audit feature')
 return mu+(arr-arr.mean(0))/es*sd
def main():
 started=time.perf_counter(); xf,xtrain,xtest,y,_=matrices(); active=json.loads((LOCAL/'audit_features.json').read_text())['features']; ext=external(active,xf)
 rng=np.random.default_rng(20264917); patients=xtest.index.astype(str).str[:12]; unique=xtest.loc[~patients.duplicated()]; source=unique[active].to_numpy()[np.sort(rng.choice(len(unique),12,False))]; external_q=ext[np.sort(rng.choice(len(ext),12,False))]
 np.savez_compressed(LOCAL/'query_matrices.npz',source_query12=source,external_query12=external_q,external_full=ext,active_features=np.asarray(active),model_features=np.asarray(xtrain.columns,dtype=object))
 train_pat=xtrain.index.astype(str).str[:12];groups=[np.flatnonzero(train_pat==p) for p in train_pat.unique()];br=np.random.default_rng(SEED);ix=np.concatenate([groups[i] for i in br.integers(0,len(groups),len(groups))]); surv=Surv.from_arrays(event=y.event.astype(bool),time=y.time)
 # deterministic sampled inactive background
 uniq=xtrain.loc[~train_pat.duplicated()]; bg_rng=np.random.default_rng(SEED+4200); sel=np.sort(bg_rng.choice(len(uniq),16,False)); bg=uniq.iloc[sel].to_numpy(np.float32); ai=np.array([xtrain.columns.get_loc(f) for f in active]); am=xtrain.mean().to_numpy()[ai]
 rh=load('07_rsf_pilot.py','p21gate_rsf'); t=time.perf_counter(); model=RandomSurvivalForest(n_estimators=500,min_samples_split=10,min_samples_leaf=5,max_features='sqrt',random_state=SEED,n_jobs=-1,low_memory=True);model.fit(xtrain.to_numpy(np.float32)[ix],surv[ix]);fit_r=time.perf_counter()-t
 t=time.perf_counter()
 for site,q in enumerate([source,external_q]):
  rh.background_permutation_attributions(model,q,bg,am,ai,permutations=8,seed=SEED+10000+100*site)
  rh.background_permutation_attributions(model,q,bg,am,ai,permutations=128,seed=SEED+20000+100*site)
 attr_r=time.perf_counter()-t
 # DeepSurv timing uses the frozen v3 architecture and the validated Phase-7.2 bootstrap trainer.
 p6=load('20_deepsurv_pilot.py','p21gate_p6');p61=load('21_deepsurv_explainer_revision.py','p21gate_p61');p72=load('41_phase72_deepsurv.py','p21gate_p72');ck=torch.load(LOCAL/'deepsurv.pt',map_location='cpu',weights_only=False); raw=xtrain[ck['ordered_columns']].to_numpy(np.float32)
 t=time.perf_counter();dm,means,scales,record=p6.train_bootstrap_model(ck,raw,y.event.to_numpy(bool),y.time.to_numpy(np.float32),ix,0);fit_d=time.perf_counter()-t
 db=xtrain[ck['ordered_columns']].loc[~train_pat.duplicated()].iloc[sel].to_numpy(np.float32);dai=np.array([ck['ordered_columns'].index(f) for f in active]);t=time.perf_counter()
 rem=np.setdiff1d(np.arange(len(uniq)),sel); ex=np.sort(bg_rng.choice(rem,112,False)); db128=xtrain[ck['ordered_columns']].loc[~train_pat.duplicated()].iloc[np.r_[sel,ex]].to_numpy(np.float32)
 for q in [source,external_q]:
  p72.deep_attributions_with_relative_gate(p61,p6,dm,ck,q,db,dai,means,scales)
  p72.deep_attributions_with_relative_gate(p61,p6,dm,ck,q,db128,dai,means,scales)
 attr_d=time.perf_counter()-t
 # Null gate is tensor-only; prior 65-feature calibration cost is conservatively bounded at 30 minutes/model.
 rows=[]
 for name,fit,attr in [('RSF-v3',fit_r,attr_r),('DeepSurv-v3',fit_d,attr_d)]:
  n50=50*(fit+attr); null=1800.; total=n50+null;rows.append({'model':name,'single_fit_seconds':fit,'single_attribution_seconds':attr,'projected_n50_seconds':n50,'projected_null_gate_seconds':null,'projected_total_seconds':total,'projected_total_hours':total/3600,'pass_24h':total<=86400})
 frame=pd.DataFrame(rows);frame.to_csv(OUT/'resource_gate.csv',index=False); payload={'family_features':len(active),'rows':rows,'gate_pass':bool(frame.pass_24h.all()),'wall_clock_seconds':time.perf_counter()-started};(OUT/'resource_gate.json').write_text(json.dumps(payload,indent=2));(OUT/'RESOURCE_GATE.md').write_text('# Faz 21 v3 kaynak kapısı\n\n'+frame.to_string(index=False)+'\n\n24 saat kapısı: '+('GEÇTİ' if payload['gate_pass'] else 'KALDI')+'\n',encoding='utf-8');print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
