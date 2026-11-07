#!/usr/bin/env python3
"""Phase 14 fast checks: R universe, calibration tables, and gene selection."""
from __future__ import annotations
import hashlib,json,re,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from warrantxai.components import rank_reproducibility
OUT=ROOT/'results/phase14'; OUT.mkdir(parents=True,exist_ok=True)
ETAS=(.01,.02,.05,.10); SIZES=(27,100,300,610)
def sha(p): return hashlib.file_digest(open(p,'rb'),'sha256').hexdigest()
def r_universe():
 with np.load(ROOT/'results/phase72_cox/attributions_n50.npz',allow_pickle=False) as z:
  src=z['source_phi']; feats=z['features'].astype(str); active=z['active_features'].astype(str)
 idx=np.array([list(feats).index(x) for x in active]); imp=np.abs(src).mean(1)[None]
 rep=pd.read_csv(ROOT/'results/phase72_cox/feature_report.csv').set_index('feature').loc[active]
 r0=rank_reproducibility(imp,.02)[idx]; np.testing.assert_allclose(r0,rep.R,atol=1e-12,rtol=0)
 rows=[]
 for eta in ETAS:
  for size in SIZES:
   seeds=range(5) if size in (100,300) else [0]
   for seed in seeds:
    if size==610: subset=np.arange(610)
    elif size==27: subset=idx
    else:
     rest=np.setdiff1d(np.arange(610),idx); rng=np.random.default_rng(20261107+size*10+seed)
     subset=np.r_[idx,np.sort(rng.choice(rest,size-27,replace=False))]
    all_r=rank_reproducibility(imp[:,:,subset],eta)
    r=all_r[idx] if size==610 else all_r[:27]
    for j,f in enumerate(active):
     other=rep.loc[f,['F','S','U','G_claim']].to_numpy(float); z=np.r_[other,r[j]]
     rows.append(dict(feature=f,eta=eta,universe=size,seed=seed,R=r[j],WI=z.min(),bottleneck=['F','S','U','G_claim','R'][z.argmin()],point_adequate=z.min()>.8))
 out=pd.DataFrame(rows); out.to_csv(ROOT/'PHASE14_ADIM1_R_UNIVERSE.csv',index=False); return out
def calibration():
 main=pd.read_csv(ROOT/'results/phase72_calibration/type1_rates.csv'); main.insert(0,'source','phase72')
 remedy=pd.read_csv(ROOT/'results/phase52_rsf_remedy/candidate_type1_rates.csv').rename(columns={'candidate':'method'}); remedy.insert(0,'source','phase52_rsf_n20')
 cols=['source','model','method','null_scenario','component','nominal_alpha','rejections','outer_simulations','empirical_type1','exact_95ci_lower','exact_95ci_upper','one_sided_binomial_p','anti_conservative_gate']
 out=pd.concat([main[cols],remedy[cols]],ignore_index=True); out.to_csv(ROOT/'PHASE14_ADIM5_CALIBRATION_TYPE1.csv',index=False); return out
def gene_selection():
 base=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer'
 panel0=np.load(ROOT/'results/phase72_panel/query_matrices.npz',allow_pickle=True)['active_features'].astype(str)
 panel1=np.load(ROOT/'results/phase75_fold1_panel/query_matrices.npz',allow_pickle=True)['active_features'].astype(str)
 selected=json.loads((ROOT/'results/phase7_retrain/selected_features.json').read_text(encoding='utf-8'))['selected_by_block']['rna']
 wanted=set(selected)|set([x for x in panel1 if x.startswith('rna_')])
 train=pd.read_parquet(ROOT/'results/phase7_retrain/X_train.parquet').set_index('sample_key'); test=pd.read_parquet(ROOT/'results/phase7_retrain/X_test.parquet').set_index('sample_key')
 male_col='clin_categorical__gender.demographic_male'; female_col='clin_categorical__gender.demographic_female'
 rank=pd.read_csv(ROOT/'results/phase7_retrain/univariate_cox_rna.tsv',sep='\t') if (ROOT/'results/phase7_retrain/univariate_cox_rna.tsv').exists() else pd.read_csv(base/'univariate_cox_rna.tsv',sep='\t')
 rank=rank.set_index('feature'); fmap=pd.read_csv(ROOT/'results/phase72_panel/feature_mapping.csv').set_index('feature')
 with np.load(ROOT/'results/phase72_panel/query_matrices.npz',allow_pickle=True) as z: ext=pd.DataFrame(z['external_full'],columns=z['active_features'].astype(str))
 rows=[]
 for f in sorted(wanted):
  tr=train[f] if f in train else pd.Series(dtype=float); te=test[f] if f in test else pd.Series(dtype=float)
  ensg=f.removeprefix('rna_').split('.')[0]; symbol=fmap.loc[f,'symbol'] if f in fmap.index and 'symbol' in fmap.columns else None
  chromosome='Y' if symbol=='RBMY1F' else 'not available in local annotation'
  e=ext[f] if f in ext else pd.Series(dtype=float); rr=rank.loc[f] if f in rank.index else {}
  tr_nz=tr.ne(0) if len(tr) else pd.Series(False,index=train.index); te_nz=te.ne(0) if len(te) else pd.Series(False,index=test.index)
  rows.append(dict(feature=f,ENSG=ensg,symbol=symbol,chromosome=chromosome,fold0_active=f in panel0,fold1_active=f in panel1,
   train_zero_rate=float((tr==0).mean()) if len(tr) else np.nan,train_min=tr.min(),train_median=tr.median(),train_max=tr.max(),train_nonzero=int((tr!=0).sum()) if len(tr) else 0,
   test_zero_rate=float((te==0).mean()) if len(te) else np.nan,test_min=te.min(),test_median=te.median(),test_max=te.max(),test_nonzero=int((te!=0).sum()) if len(te) else 0,
   metabric_available=f in ext,metabric_zero_rate=float((e==0).mean()) if len(e) else np.nan,metabric_min=e.min(),metabric_median=e.median(),metabric_max=e.max(),
   cox_coefficient=rr.get('coefficient',np.nan),cox_p=rr.get('p_value',np.nan),rna_rank=rr.get('rank',np.nan),
   train_nonzero_male=int(((train[male_col]>0)&tr_nz).sum()) if male_col in train else np.nan,
   train_nonzero_female=int(((train[female_col]>0)&tr_nz).sum()) if female_col in train else np.nan,
   test_nonzero_male=int(((test[male_col]>0)&te_nz).sum()) if male_col in test else np.nan,
   test_nonzero_female=int(((test[female_col]>0)&te_nz).sum()) if female_col in test else np.nan))
 out=pd.DataFrame(rows); out.to_csv(ROOT/'PHASE14_ADIM3_GENE_SELECTION.csv',index=False); return out,rank
def main():
 before={f:sha(ROOT/'warrantxai'/f) for f in ['components.py','aggregation.py','certification.py']}
 r=r_universe(); c=calibration(); g,rank=gene_selection(); after={f:sha(ROOT/'warrantxai'/f) for f in before}; assert before==after
 print(json.dumps({'r_rows':len(r),'calibration_rows':len(c),'gene_rows':len(g),'hashes':before},indent=2))
if __name__=='__main__': main()
