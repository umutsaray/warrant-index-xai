#!/usr/bin/env python3
"""Phase 17 primary v2 six-component recursive-studentized certification."""
from __future__ import annotations
import importlib.util,sys,time
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase17';OUT.mkdir(parents=True,exist_ok=True);sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
from _shared import CALIBRATION
from warrantxai.certification import certify
SCALES={'CoxPH-v2':0.023934215222165295,'RSF':0.01707069437342852,'DeepSurv':0.012307440008774696};COMP=('F','S','U','G_claim','R','N')
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def n_stats(sites,m0,indices,shared,remedy):
 leaves=np.stack([np.minimum(1,np.abs(x).mean(axis=1)/m0) for x in sites]);t=leaves.shape[1];draws=len(indices);counts=np.zeros((draws,t),np.int16);np.add.at(counts,(np.arange(draws)[:,None],indices),1);point=leaves.mean(1);sample,seboot=shared.resampled_means(leaves,counts);se=leaves.std(1,ddof=1)/np.sqrt(t);p=remedy.studentized_p(point.T.reshape(-1),np.moveaxis(sample,1,2).reshape(draws,-1),se.T.reshape(-1),np.moveaxis(seboot,1,2).reshape(draws,-1)).reshape(leaves.shape[2],leaves.shape[0]);return point.min(0),p.max(1)
def main():
 started=time.perf_counter();cert=load('43_phase72_certification.py','p17v2cert');shared=load('12_rsf_certification.py','p17v2shared');remedy=load('11_rsf_remedy_calibration.py','p17v2remedy');rng=np.random.default_rng(cert.SEED+725_000);indices=rng.integers(0,50,size=(cert.DRAWS,50));original={k:CALIBRATION[k] for k in ['a','b','u','eta']};rows=[];compare=[]
 try:
  for model,src,off in [('CoxPH-v2','phase72_cox','phase72_coxph_certification'),('RSF','phase72_rsf','phase72_rsf_certification'),('DeepSurv','phase72_deepsurv','phase72_deepsurv_certification')]:
   sm=SCALES[model];CALIBRATION.update(a=sm,b=sm,u=6*sm,eta=sm)
   if model=='CoxPH-v2':
    points,pvals=cert.cox_statistics(indices)
    with np.load(ROOT/'results/phase72_cox/attributions_n50.npz') as z:
     fs=z['features'].astype(str).tolist();act=z['active_features'].astype(str).tolist();ix=np.array([fs.index(f) for f in act]);sites=[z['source_phi'][:,:,ix],z['external_phi']];features=act
   else:
    with np.load(ROOT/'results'/src/'attributions_n50.npz') as z:phi=z['candidate'];ref=z['reference'];features=z['features'].astype(str).tolist()
    points,pvals=shared.component_statistics(phi,ref,indices);sites=[phi[0],phi[1]]
   npnt,npv=n_stats(sites,.25*sm,indices,shared,remedy);p6=np.column_stack([pvals,npv]);z6=np.column_stack([points,npnt]);bh=certify(p6,.05,'bh');by=certify(p6,.05,'by');official=pd.read_csv(ROOT/'results'/off/'feature_certification.csv').set_index('feature')
   for j,f in enumerate(features):
    vals=z6[j];rows.append(dict(model=model,feature=f,scale_s_m=sm,k=1,rho=1,kappa=.25,a=sm,b=sm,u=6*sm,eta=sm,m0=.25*sm,**{c:vals[q] for q,c in enumerate(COMP)},WI6=vals.min(),bottleneck6=COMP[int(vals.argmin())],point_adequate6=vals.min()>.8,**{f'p_{c}':p6[j,q] for q,c in enumerate(COMP)},p_joint=bh.joint_p[j],q_BH=bh.adjusted_p[j],q_BY=by.adjusted_p[j],certified_BH=bh.rejected[j],certified_BY=by.rejected[j]))
    o=official.loc[f];compare.append(dict(model=model,feature=f,v1_point=bool(o.point_WI_ge_0_8),v1_BH=bool(o.certified_BH_alpha_0_05),v1_BY=bool(o.certified_BY_alpha_0_05),v2_point=bool(vals.min()>.8),v2_BH=bool(bh.rejected[j]),v2_BY=bool(by.rejected[j]),v1_WI=o.WI,v2_WI6=vals.min(),v2_N=npnt[j],change_BH='retained' if o.certified_BH_alpha_0_05 and bh.rejected[j] else ('dropped' if o.certified_BH_alpha_0_05 else ('new' if bh.rejected[j] else 'neither'))))
 finally:CALIBRATION.update(original)
 pd.DataFrame(rows).to_csv(OUT/'PHASE17_ADIM3_V2_CERT_PRIMARY.csv',index=False);pd.DataFrame(compare).to_csv(OUT/'PHASE17_ADIM3_V1_VS_V2_81.csv',index=False);print('elapsed',time.perf_counter()-started)
if __name__=='__main__':main()
