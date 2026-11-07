#!/usr/bin/env python3
"""Phase 16: diagnostic full RNA-pool Wald/LRT/score/winsorized rankings."""
from __future__ import annotations
import importlib.util,json,sys,time
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase16';OUT.mkdir(parents=True,exist_ok=True)
SRC=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer'/'brca-multiomics-survival/src/feature_selection/select_features.py'
s=importlib.util.spec_from_file_location('p16select',SRC);sel=importlib.util.module_from_spec(s);sys.modules[s.name]=sel;s.loader.exec_module(sel)
def frame(p):return pd.read_parquet(p).set_index('sample_key')
def risk_setup(t,e):
 order=np.argsort(-t,kind='stable');ts=t[order];es=e[order].astype(bool);ge=np.empty(len(t),int);a=0
 while a<len(t):
  b=a
  while b+1<len(t) and ts[b+1]==ts[a]:b+=1
  ge[a:b+1]=b;a=b+1
 er=np.flatnonzero(es);return order,er,ge[er]
def lrt_score(z,order,er,rr):
 x=z[order];event_sum=x[er].sum(0);s0=np.arange(1,len(x)+1,dtype=float)[rr];s1=np.cumsum(x,axis=0)[rr];s2=np.cumsum(x*x,axis=0)[rr];u0=event_sum-(s1/s0[:,None]).sum(0);i0=(s2/s0[:,None]-(s1/s0[:,None])**2).sum(0);score=chi2.sf(np.divide(u0*u0,i0,out=np.zeros_like(u0),where=i0>0),1);ll0=-np.log(s0).sum()
 lo=np.full(x.shape[1],-10.);hi=np.full(x.shape[1],10.)
 for _ in range(45):
  b=(lo+hi)/2;eta=np.clip(x*b,-50,50);w=np.exp(eta);a0=np.cumsum(w,axis=0)[rr];a1=np.cumsum(w*x,axis=0)[rr];u=event_sum-(a1/a0).sum(0);lo=np.where(u>0,b,lo);hi=np.where(u>0,hi,b)
 beta=(lo+hi)/2;eta=np.clip(x*beta,-50,50);ll=eta[er].sum(0)-np.log(np.cumsum(np.exp(eta),axis=0)[rr]).sum(0);lrt=chi2.sf(np.maximum(0,2*(ll-ll0)),1);return lrt,score,beta,ll
def ranks(df,pcol,tie=None):
 cols=[pcol]+([tie] if tie else []);asc=[True]+([False] if tie else []);order=df.sort_values(cols,ascending=asc,na_position='last',kind='stable').index;out=pd.Series(np.arange(1,len(df)+1),index=order);return out.reindex(df.index).to_numpy()
def bh(p):
 p=np.asarray(p);o=np.argsort(p);q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];r=np.empty_like(q);r[o]=np.minimum(q,1);return r
def main():
 start=time.perf_counter();x=frame(ROOT/'results/phase7_retrain/X_train.parquet');y=frame(ROOT/'results/phase7_retrain/y_train.parquet').loc[x.index];cols=[c for c in x if c.startswith('rna_')];saved=pd.read_csv(ROOT/'results/phase7_retrain/univariate_cox_rna.tsv',sep='\t').set_index('feature').loc[cols];t=y.time.to_numpy(float);e=y.event.to_numpy(int);order,er,rr=risk_setup(t,e);parts=[]
 for a in range(0,len(cols),512):
  names=cols[a:a+512];v=x[names].to_numpy(float);b,se,p,ll=sel.univariate_cox_chunk(v,t,e,50,1e-7);mu=v.mean(0);sd=v.std(0);valid=sd>0;z=np.zeros_like(v);z[:,valid]=(v[:,valid]-mu[valid])/sd[valid];lrt,score,lrtb,lrtll=lrt_score(z[:,valid],order,er,rr);lp=np.full(len(names),np.nan);sp=lp.copy();lb=lp.copy();lll=lp.copy();lp[valid]=lrt;sp[valid]=score;lb[valid]=lrtb;lll[valid]=lrtll;bw,sew,pw,llw=sel.univariate_cox_chunk(np.clip(z,-5,5),t,e,50,1e-7)
  parts.append(pd.DataFrame({'feature':names,'wald_beta':b,'wald_p':p,'wald_ll':ll,'lrt_beta':lb,'lrt_ll':lll,'lrt_p':lp,'score_p':sp,'winsor5_beta':bw,'winsor5_p':pw,'zero_rate':(v==0).mean(0),'nonzero_count':(v!=0).sum(0)}));print(min(a+512,len(cols)),len(cols),flush=True)
 out=pd.concat(parts,ignore_index=True);np.testing.assert_allclose(out.wald_p,saved.p_value,rtol=1e-6,atol=0,equal_nan=True);out['wald_rank']=ranks(out,'wald_p','wald_ll');out['lrt_rank']=ranks(out,'lrt_p','lrt_ll');out['score_rank']=ranks(out,'score_p');out['winsor5_rank']=ranks(out,'winsor5_p','wald_ll')
 chosen=json.loads((ROOT/'results/phase7_retrain/selected_features.json').read_text(encoding='utf-8'))['selected_by_block']['rna'];got=out.sort_values(['wald_p','wald_ll'],ascending=[True,False]).head(200).feature.tolist();assert got==chosen,'official first-200 set/order mismatch'
 out['official_selected_200']=out.feature.isin(chosen);active=set(np.load(ROOT/'results/phase72_panel/query_matrices.npz',allow_pickle=False)['active_features'].astype(str));out['fold0_active']=out.feature.isin(active);out.to_csv(OUT/'PHASE16_ADIM1_FULL_POOL_RANKING.csv',index=False)
 fmap=pd.read_csv(ROOT/'results/phase72_panel/feature_mapping.csv')[['feature','symbol']];focus={'RBMY1F','OR6C1','OR4C46','OR13C3','OR5M10','KRTAP19-2','KRTAP19-7','KRTAP11-1'};q=out.merge(fmap,on='feature',how='left');q=q[q.symbol.isin(focus)].copy();q['lrt_bh_q_full_pool']=bh(out.lrt_p.fillna(1))[q.index];q.to_csv(OUT/'PHASE16_ADIM1_FOCUS_GENES.csv',index=False);print('elapsed',time.perf_counter()-start)
if __name__=='__main__':main()
