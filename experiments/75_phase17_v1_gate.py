#!/usr/bin/env python3
"""Phase 17 v1 reproduction gate; writes only new phase17 diagnostics."""
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/phase17';OUT.mkdir(parents=True,exist_ok=True);sys.path[:0]=[str(ROOT),str(ROOT/'experiments')]
from warrantxai.certification import certify
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def main():
 cert=load('43_phase72_certification.py','p17cert');shared=load('12_rsf_certification.py','p17shared');rng=np.random.default_rng(cert.SEED+725_000);indices=rng.integers(0,50,size=(cert.DRAWS,50));rows=[]
 cp,cv=cert.cox_statistics(indices);models=[('CoxPH-v2','phase72_cox','phase72_coxph_certification',cp,cv)]
 for model,src,off in [('RSF','phase72_rsf','phase72_rsf_certification'),('DeepSurv','phase72_deepsurv','phase72_deepsurv_certification')]:
  with np.load(ROOT/'results'/src/'attributions_n50.npz') as z:p,v=shared.component_statistics(z['candidate'],z['reference'],indices)
  models.append((model,src,off,p,v))
 for model,src,off,points,pvals in models:
  fr=pd.read_csv(ROOT/'results'/src/'feature_report.csv');official=pd.read_csv(ROOT/'results'/off/'feature_certification.csv');np.testing.assert_allclose(points,fr[['F','S','U','G_claim','R']],atol=1e-12,rtol=0);bh=certify(pvals,.05,'bh');by=certify(pvals,.05,'by');np.testing.assert_allclose(bh.joint_p,official.p_joint,atol=1e-12,rtol=0);np.testing.assert_allclose(bh.adjusted_p,official.q_BH,atol=1e-12,rtol=0);np.testing.assert_allclose(by.adjusted_p,official.q_BY,atol=1e-12,rtol=0);assert np.array_equal(bh.rejected,official.certified_BH_alpha_0_05) and np.array_equal(by.rejected,official.certified_BY_alpha_0_05)
  for j,f in enumerate(fr.feature):rows.append(dict(model=model,feature=f,F=points[j,0],S=points[j,1],U=points[j,2],G_claim=points[j,3],R=points[j,4],p_joint=bh.joint_p[j],q_BH=bh.adjusted_p[j],q_BY=by.adjusted_p[j],certified_BH=bh.rejected[j],certified_BY=by.rejected[j],official_exact_match=True))
 pd.DataFrame(rows).to_csv(OUT/'PHASE17_ADIM0_V1_REPRODUCTION.csv',index=False)
if __name__=='__main__':main()
