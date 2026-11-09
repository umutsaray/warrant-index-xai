"""Phase 21 ADIM 7.2 complementary tables and requested diagnostics."""
from __future__ import annotations
import importlib.util,json,sys
from pathlib import Path
import joblib,numpy as np,pandas as pd,torch
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'results/phase7_retrain';LOCAL=ROOT/'results/phase21_v3';V3=ROOT/'phase21_reviewer1/v3';OUT=ROOT/'phase21_reviewer1';SIB=ROOT.parent/'Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer';sys.path[:0]=[str(ROOT),str(ROOT/'experiments'),str(SIB)]
from warrantxai.components import generalizability,generalizability_signed
def load(file,name):
 s=importlib.util.spec_from_file_location(name,ROOT/'experiments'/file);m=importlib.util.module_from_spec(s);sys.modules[name]=m;assert s.loader;s.loader.exec_module(m);return m
def block(f):return 'RNA' if f.startswith('rna_') else 'clinical'
def g(a,b):
 gm=generalizability(np.stack([np.abs(a).mean(0),np.abs(b).mean(0)]),.5);gs=generalizability_signed(np.stack([np.sort(a,axis=0),np.sort(b,axis=0)]),.5);return gm,gs,np.minimum(gm,gs)
def external_raw(active,xf):
 audit=json.loads((LOCAL/'audit_features.json').read_text());symbols={f:audit['rna_symbols'][f] for f in active if f.startswith('rna_')};wanted=set(symbols.values());chunks=[];path=ROOT/'data/external/metabric_phase21/data_mrna_illumina_microarray.txt'
 for c in pd.read_csv(path,sep='\t',chunksize=256,low_memory=False):
  z=c[c.Hugo_Symbol.astype(str).isin(wanted)]
  if len(z):chunks.append(z)
 expr=pd.concat(chunks).drop_duplicates('Hugo_Symbol').set_index('Hugo_Symbol');samples=[c for c in expr.columns if c!='Entrez_Gene_Id'];old=ROOT/'data/external/metabric_phase72';si=pd.DataFrame(json.loads((old/'samples.json').read_text())).set_index('sampleId');pr=pd.DataFrame(json.loads((old/'patient_clinical.json').read_text())).pivot(index='patientId',columns='clinicalAttributeId',values='value');sr=pd.DataFrame(json.loads((old/'sample_clinical.json').read_text())).pivot(index='sampleId',columns='clinicalAttributeId',values='value');linked=si.reindex(samples).patientId;tm=pd.to_numeric(pr.reindex(linked).OS_MONTHS,errors='coerce').set_axis(samples);st=pr.reindex(linked).OS_STATUS.set_axis(samples);keep=tm.notna()&tm.gt(0)&st.isin(['0:LIVING','1:DECEASED']);samples=list(np.asarray(samples)[keep.to_numpy()]);linked=si.reindex(samples).patientId;cols=[]
 for f in active:
  if f.startswith('rna_'):cols.append(pd.to_numeric(expr.loc[symbols[f],samples],errors='coerce').to_numpy(float))
  elif 'age_at_index' in f:cols.append(pd.to_numeric(pr.reindex(linked).AGE_AT_DIAGNOSIS,errors='coerce').to_numpy(float))
  else:
   stage=next(x for x in ['Stage I','Stage II','Stage III','Stage IV'] if f.endswith('_'+x));code={'Stage I':1,'Stage II':2,'Stage III':3,'Stage IV':4}[stage];v=pd.to_numeric(sr.reindex(samples).TUMOR_STAGE,errors='coerce').to_numpy(float);cols.append(np.where(np.isfinite(v),(v==code).astype(float),np.nan))
 arr=np.column_stack(cols);arr=arr[np.isfinite(arr).all(1)];mu=xf[active].mean().to_numpy();sd=xf[active].std(ddof=0).to_numpy();aligned=mu+(arr-arr.mean(0))/arr.std(0)*sd;return arr,aligned
def epoch_curve(xf,xv,xt,yf,yv,yt,ck):
 m=load('20_deepsurv_pilot.py','p21a7deeptrain');from train_explainable_dl import OmicsDeepSurv,make_blocks,predict_risk,cox_negative_log_likelihood,cindex,set_seed
 set_seed(42);cols=ck['ordered_columns'];a=lambda x:(x[cols].to_numpy(np.float32)-ck['means'])/ck['scales'];af,av,at=a(xf),a(xv),a(xt);fb=make_blocks(af,ck['slices'],torch.device('cpu'));vb=make_blocks(av,ck['slices'],torch.device('cpu'));tb=make_blocks(at,ck['slices'],torch.device('cpu'));model=OmicsDeepSurv(ck['block_dims'],.3);opt=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001);tt=torch.tensor(yf.time.to_numpy(),dtype=torch.float32);ee=torch.tensor(yf.event.to_numpy(),dtype=torch.float32);rows=[];best=-np.inf;wait=0
 for epoch in range(1,301):
  model.train();opt.zero_grad();loss=cox_negative_log_likelihood(model(fb),tt,ee);loss.backward();opt.step();vr=predict_risk(model,vb);tr=predict_risk(model,tb);vs=cindex(yv.event.to_numpy(),yv.time.to_numpy(),vr);ts=cindex(yt.event.to_numpy(),yt.time.to_numpy(),tr);rows.append({'epoch':epoch,'train_loss':float(loss.detach()),'validation_cindex':vs,'test_cindex':ts})
  if vs>best+1e-6:best=vs;wait=0
  else:wait+=1
  if wait>=30:break
 return pd.DataFrame(rows)
def main():
 models={'CoxPH-v3':'cox','RSF-v3':'rsf','DeepSurv-v3':'deepsurv'};cert=[];heat=[]
 rb=joblib.load(LOCAL/'rsf.joblib');rf=rb['features']
 for model,stem in models.items():
  tab=pd.read_csv(V3/f'{model.lower().replace("-","_")}_feature_certification.csv');
  with np.load(LOCAL/f'{stem}_attributions_n50.npz') as z:phi=z['candidate'];features=z['features'].astype(str).tolist()
  for j,f in enumerate(features):
   vals=phi[0,:,:,j];M=float(np.abs(vals).mean());P=float(np.count_nonzero(vals)/vals.size);Q=float(np.abs(vals).max());row=tab.loc[tab.feature==f].iloc[0];heat.append({'model':model,'feature':f,'block':block(f),**{c:row[c] for c in ['F','S','U','G_claim','R','WI']}})
   if bool(row.certified_BH):
    trees=np.nan
    if model=='RSF-v3':ix=rf.index(f);trees=sum(np.any(t.tree_.feature==ix) for t in rb['model'].estimators_)
    cert.append({'model':model,'feature':f,'block':block(f),'WI':row.WI,'bottleneck':row.bottleneck,'M_j':M,'P_j':P,'Q_j':Q,'warrant_type':'informative' if Q>0 else 'unexercised','RSF_trees_with_split':trees})
 certified=pd.DataFrame(cert);certified.to_csv(OUT/'PHASE21_ADIM7_04_CERTIFIED_FEATURES.csv',index=False);heat=pd.DataFrame(heat);heat.to_csv(OUT/'PHASE21_ADIM7_05_HEATMAP_DATA.csv',index=False);bn=heat.assign(bottleneck=heat[['F','S','U','G_claim','R']].idxmin(axis=1)).groupby(['model','bottleneck']).size().rename('n').reset_index();bn.to_csv(OUT/'PHASE21_ADIM7_05_BOTTLENECK_DISTRIBUTION.csv',index=False)
 # Design statement and full 2x2 copy under ADIM7 names.
 primary=pd.read_csv(V3/'certification_summary.csv');design=pd.DataFrame([{'result':'primary_Cox_7_4_4','TCGA_queries':12,'METABRIC_queries':12,'R_universe':65,'family':65,'point':7,'BH':4,'BY':4},{'result':'matched_three_model_Cox','TCGA_queries':12,'METABRIC_queries':12,'R_universe':65,'family':65,'point':int(primary.loc[primary.model=='CoxPH-v3','point'].iloc[0]),'BH':int(primary.loc[primary.model=='CoxPH-v3','BH'].iloc[0]),'BY':int(primary.loc[primary.model=='CoxPH-v3','BY'].iloc[0])}]);design.to_csv(OUT/'PHASE21_ADIM7_06_COX_DESIGN.csv',index=False);pd.read_csv(V3/'secondary_fourcell.csv').to_csv(OUT/'PHASE21_ADIM7_07_FOURCELL.csv',index=False)
 # Standardization single-reference diagnostic.
 cache=SIB/'phase1_fold_cache';xf=pd.read_parquet(cache/'97cc92aa9eff44f1_train.parquet');xr=pd.read_parquet(cache/'97cc92aa9eff44f1_test.parquet');yt=pd.read_parquet(DATA/'y_train.parquet').set_index('sample_key');ye=pd.read_parquet(DATA/'y_test.parquet').set_index('sample_key');xtrain=pd.concat([xf,xr.loc[yt.index.intersection(xr.index)]]).loc[yt.index];xtest=xr.loc[ye.index]
 with np.load(LOCAL/'query_matrices.npz',allow_pickle=True) as z:active=z['active_features'].astype(str).tolist();source=z['source_query12'];official=z['external_query12']
 raw,aligned=external_raw(active,xf);rng=np.random.default_rng(20264917);rawq=raw[np.sort(rng.choice(len(raw),12,False))]
 pats=xtrain.index.astype(str).str[:12];uniq=xtrain.loc[~pats.duplicated()];br=np.random.default_rng(20260917+4200);first=np.sort(br.choice(len(uniq),16,False));rows=[];rh=load('07_rsf_pilot.py','p21a7rsf');ri=np.array([rf.index(f) for f in active]);rbg=uniq[rf].iloc[first].to_numpy(np.float32);ram=xtrain[rf].mean().to_numpy()[ri]
 p6=load('20_deepsurv_pilot.py','p21a7p6');p61=load('21_deepsurv_explainer_revision.py','p21a7p61');p72=load('41_phase72_deepsurv.py','p21a7p72');ck=torch.load(LOCAL/'deepsurv.pt',map_location='cpu',weights_only=False);dm=p6.model_from_state(ck);di=np.array([ck['ordered_columns'].index(f) for f in active]);dbg=uniq[ck['ordered_columns']].iloc[first].to_numpy(np.float32);cb=joblib.load(LOCAL/'coxph.joblib');ci=np.array([cb['features'].index(f) for f in active]);coef=cb['model'].params_.to_numpy()[ci]
 for model in models:
  base=pd.read_csv(V3/f'{model.lower().replace("-","_")}_feature_certification.csv').set_index('feature').loc[active]
  if model=='CoxPH-v3':src=coef* ((source-cb['scaler'].mean_[ci])/cb['scaler'].scale_[ci]);offs=[coef*((official-cb['scaler'].mean_[ci])/cb['scaler'].scale_[ci]),coef*((rawq-cb['scaler'].mean_[ci])/cb['scaler'].scale_[ci])]
  elif model=='RSF-v3':src=rh.background_permutation_attributions(rb['model'],source,rbg,ram,ri,permutations=8,seed=20274917);offs=[rh.background_permutation_attributions(rb['model'],q,rbg,ram,ri,permutations=8,seed=20275017) for q in [official,rawq]]
  else:src=p72.deep_attributions_with_relative_gate(p61,p6,dm,ck,source,dbg,di,ck['means'],ck['scales'])[0];offs=[p72.deep_attributions_with_relative_gate(p61,p6,dm,ck,q,dbg,di,ck['means'],ck['scales'])[0] for q in [official,rawq]]
  for method,ext in zip(['official_cohort_standardized','tcga_train_scaler_only_no_cohort_standardization'],offs):
   gm,gs,gc=g(src,ext);point=0
   for j,f in enumerate(active):
    vals=np.array([base.loc[f,'F'],base.loc[f,'S'],base.loc[f,'U'],gc[j],base.loc[f,'R']]);point+=vals.min()>=.8;rows.append({'model':model,'method':method,'feature':f,'G_claim':gc[j],'diagnostic_WI':vals.min(),'point_adequate':vals.min()>=.8})
 std=pd.DataFrame(rows);std.to_csv(OUT/'PHASE21_ADIM7_08_STANDARDIZATION_FEATURES.csv',index=False);stds=std.groupby(['model','method']).agg(point=('point_adequate','sum'),median_G_claim=('G_claim','median')).reset_index();stds.to_csv(OUT/'PHASE21_ADIM7_08_STANDARDIZATION_SUMMARY.csv',index=False)
 # DeepSurv epoch curve.
 xv=xr.loc[yt.index.intersection(xr.index)];curve=epoch_curve(xf,xv,xr.loc[ye.index],yt.loc[xf.index],yt.loc[xv.index],ye,ck);curve.to_csv(OUT/'PHASE21_ADIM7_09_DEEPSURV_EPOCH_CURVE.csv',index=False);rho=float(spearmanr(curve.validation_cindex,curve.test_cindex).statistic);pd.DataFrame([{'best_validation_epoch':int(curve.loc[curve.validation_cindex.idxmax(),'epoch']),'best_validation_cindex':curve.validation_cindex.max(),'test_cindex_at_best_validation':curve.loc[curve.validation_cindex.idxmax(),'test_cindex'],'validation_minus_test':curve.validation_cindex.max()-curve.loc[curve.validation_cindex.idxmax(),'test_cindex'],'spearman_validation_test':rho,'v1_reported_mean_gap':.264}]).to_csv(OUT/'PHASE21_ADIM7_09_DEEPSURV_EPOCH_SUMMARY.csv',index=False)
 # v1-v3 overlap.
 v1active=json.loads((ROOT/'results/phase72_panel/metadata.json').read_text())['active_features'];v3set=set(active);v1pairs=[]
 for model,path in [('CoxPH-v3',ROOT/'results/phase72_coxph_certification/feature_certification.csv'),('RSF-v3',ROOT/'results/phase72_rsf_certification/feature_certification.csv'),('DeepSurv-v3',ROOT/'results/phase72_deepsurv_certification/feature_certification.csv')]:
  old=pd.read_csv(path);col='certified_BH_alpha_0_05';new=pd.read_csv(V3/f'{model.lower().replace("-","_")}_feature_certification.csv').set_index('feature')
  for f in old.loc[old[col].astype(bool),'feature']:v1pairs.append({'model':model,'feature':f,'feature_in_v3_family':f in v3set,'certified_in_v3':bool(f in new.index and new.loc[f,'certified_BH'])})
 pairs=pd.DataFrame(v1pairs);pairs.to_csv(OUT/'PHASE21_ADIM7_10_V1_V3_CERTIFIED_PAIR_OVERLAP.csv',index=False);over=pd.DataFrame([{'v1_features':27,'v1_features_in_v3_family':len(set(v1active)&v3set),'v1_certified_model_feature_pairs':len(pairs),'pairs_feature_in_v3_family':int(pairs.feature_in_v3_family.sum()),'pairs_certified_in_v3':int(pairs.certified_in_v3.sum())}]);over.to_csv(OUT/'PHASE21_ADIM7_10_V1_V3_OVERLAP_SUMMARY.csv',index=False)
 unex=int(((certified.model=='RSF-v3')&(certified.warrant_type=='unexercised')).sum());report=f'''# Faz 21 ADIM 7.2 — tamamlayıcı sayılar

RSF-v3 BH sertifikalı 14 satırda unexercised (Q_j=0): **{unex}**.

Bağlayıcı dağılım:\n```\n{bn.to_string(index=False)}\n```

Cox 7/4/4 tasarımı: 12 TCGA + 12 METABRIC sorgu, 65-feature R evreni ve 65 aile; `matched_three_model.csv` ile aynıdır, tam 153/1,466 veya 610-girdi sonucu değildir.

Standardizasyon:\n```\n{stds.to_string(index=False)}\n```

DeepSurv validation-test epoch Spearman: {rho}.\n\nV1-v3 örtüşme:\n```\n{over.to_string(index=False)}\n```
''';(OUT/'PHASE21_ADIM7_02_COMPLEMENTARY.md').write_text(report,encoding='utf-8');print(report)
if __name__=='__main__':main()
