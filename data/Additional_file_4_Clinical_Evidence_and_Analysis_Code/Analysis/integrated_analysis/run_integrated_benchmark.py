#!/usr/bin/env python3
"""Diagnostic reconstruction of the uncorrected legacy virtual-KO rank pipeline.
The fusion output of this script is NOT the primary revised analysis. Numerical-zero
runs are flagged and must be excluded before interpreting perturbation enrichment.
No patient-risk labels, candidate-specific tuning, or outcome-directed exclusions.
"""
from pathlib import Path
import argparse, hashlib, json, platform, warnings
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, beta, norm

parser=argparse.ArgumentParser()
parser.add_argument('--input-root',type=Path,required=True)
parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'uncorrected_vko_diagnostic')
args=parser.parse_args(); ROOT=args.input_root; OUT=args.output_dir; OUT.mkdir(parents=True,exist_ok=True)
SD=ROOT/'Source_Data'; RNG=np.random.default_rng(20260919); NNULL=4000; NFUSION=20000; NWEIGHT=5000
manifest=[]
def read(rel,**kw):
 p=SD/rel; manifest.append({'path':'Source_Data/'+rel,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
 return pd.read_csv(p,**kw)
def save(x,name):x.to_csv(OUT/name,index=False)
def bh(p):
 p=np.asarray(p,float);z=np.full_like(p,np.nan);m=np.isfinite(p);v=p[m];o=np.argsort(v);q=np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];r=np.empty_like(q);r[o]=np.minimum(q,1);z[m]=r;return z
def corr(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float);m=np.isfinite(a)&np.isfinite(b)
 return float(spearmanr(a[m],b[m]).statistic) if m.sum()>2 and np.std(a[m])>0 and np.std(b[m])>0 else np.nan
def score_rank(x):
 x=np.asarray(x,float);o=np.full(x.shape,np.nan);m=np.isfinite(x);o[m]=(rankdata(x[m])-0.5)/m.sum();return o
def rank_high(x):return rankdata(-np.asarray(x),method='average')
def top_set(x,k=5):return set(np.flatnonzero(rank_high(x)<=k))
def jac(a,b):return len(a&b)/len(a|b) if len(a|b) else np.nan

old=read('Integration/candidate_consensus_ranking.csv')
tar=old.loc[old.candidate_class.eq('exposure_proximal_target')].sort_values('gene_symbol').reset_index(drop=True)
GENES=tar.gene_symbol.tolist(); N=len(GENES); assert N==17
prot=read('Human_DILI_Proteomics/candidate_evidence_table.csv'); ANCHORS=prot.protein.tolist()
vko=read('Virtual_Knockout/virtual_knockout_gene_results.csv.gz'); runmax=vko.groupby(['ko_label','compartment','seed']).distance.max().to_dict(); sc=read('Liver_Single_Cell/target_expression_by_compartment.csv')
# Fixed primary: magnitude of continuous etiology contrast, without FDR prefiltering.
weight_rows=[]
for cohort in ['discovery','confirmatory']:
 for contrast in ['DO_vs_NDO','DO_vs_HV']:
  raw=prot[f'{cohort}_{contrast}_cliffs_delta_group1_minus_group2'].abs().to_numpy()
  for g,w in zip(ANCHORS,raw):weight_rows.append({'cohort':cohort,'contrast':contrast,'protein':g,'absolute_cliffs_delta':w,'normalized_weight':w/np.nansum(raw)})
weights=pd.DataFrame(weight_rows);save(weights,'clinical_projection_weights.csv')
weight_maps={(c,t):dict(zip(z.protein,z.absolute_cliffs_delta)) for (c,t),z in weights.groupby(['cohort','contrast'])}
weight_maps[('unweighted','all')]=dict.fromkeys(ANCHORS,1.)
COMP=['Hepatocyte','Endothelial','Macrophage']; seedids=sorted(vko.seed.unique())
# Each candidate is evaluated in every available liver compartment; no max selection.
# Expression-matched null pools conditioned on run gene universe and KO removal.
runrows=[]; edge_rows=[]; nullrows=[]; leaveprotein_rows=[]
for comp in COMP:
 mat=read(f'Liver_Single_Cell/{comp.lower()}_network_matrix.csv.gz',index_col=0)
 det=(mat>0).mean(axis=1);avg=mat.mean(axis=1)
 feat=pd.DataFrame({'mean':np.log1p(avg),'detection':np.log(det.clip(1e-4,1-1e-4)/(1-det.clip(1e-4,1-1e-4)))})
 feat=(feat-feat.mean())/feat.std(ddof=0)
 sub=vko.loc[(vko.compartment.eq(comp)) & (vko.ko_class.eq('drug_target')) & vko.ko_label.isin(GENES)]
 for (seed,target),run in sub.groupby(['seed','ko_label'],sort=True):
  run=run.loc[~run.is_knocked_gene.astype(bool)].drop_duplicates('gene').set_index('gene')
  ag=[g for g in ANCHORS if g in run.index and g in feat.index]
  bg=sorted(set(run.index)&set(feat.index)-set(ANCHORS)-{target})
  assert len(ag)>=4 and len(bg)>=50
  vals=run.loc[ag,'distance_percentile'].to_numpy(float)
  targetdet=float(sc.loc[sc.gene_symbol.eq(target)&sc.compartment.eq(comp),'mean_detection_fraction'].iloc[0])
  # Nearest 30 within the run; all anchors excluded from background.
  pool=[]
  for g in ag:
   ds=((feat.loc[bg]-feat.loc[g])**2).sum(axis=1).sort_values(kind='stable')
   pool.append(ds.index[:min(30,len(ds))].tolist())
  # The same matched gene draws are reused for each clinical-weight specification.
  draws=np.column_stack([RNG.choice(p,size=NNULL,replace=True) for p in pool])
  # Distinct draws sampled independently per anchor; duplicate-background draws allowed.
  nv=np.column_stack([run.loc[draws[:,j],'distance_percentile'].to_numpy() for j in range(len(ag))])
  for key,wmap in weight_maps.items():
   w=np.array([wmap[g] for g in ag]); w=w/w.sum(); obs=float(vals@w); nul=nv@w
   pp=(1+np.sum(nul>=obs))/(NNULL+1); z=(obs-nul.mean())/nul.std(ddof=1)
   runrows.append({'target':target,'compartment':comp,'seed':seed,'cohort':key[0],'contrast':key[1], 'n_proteins_available':len(ag),'clinical_weight_coverage':sum(wmap[g] for g in ag)/sum(wmap.values()),'target_detection_fraction':targetdet,'weighted_displacement':obs,'matched_null_mean':float(nul.mean()),'matched_null_sd':float(nul.std(ddof=1)),'matched_z':float(z),'matched_empirical_p':pp,'matched_upper_percentile':1-pp,'n_null_draws':NNULL,'maximum_displacement':float(runmax[(target,comp,seed)]),'numerically_nonzero':bool(runmax[(target,comp,seed)]>1e-12)})
   if key==('confirmatory','DO_vs_NDO'):
    for g,d,ww in zip(ag,vals,w):edge_rows.append({'target':target,'compartment':comp,'seed':seed,'protein':g,'unsigned_displacement_percentile':d,'conditional_clinical_weight':ww,'weighted_contribution':d*ww,'target_detection_fraction':targetdet})
    for omit in ANCHORS:
     wo=np.array([0 if g==omit else wmap[g] for g in ag]);wo=wo/wo.sum()
     leaveprotein_rows.append({'target':target,'compartment':comp,'seed':seed,'omitted_protein':omit,'weighted_displacement':float(vals@wo),'target_detection_fraction':targetdet})
runs=pd.DataFrame(runrows);runs['matched_q_within_specification']=runs.groupby(['cohort','contrast']).matched_empirical_p.transform(bh)
save(runs,'target_protein_projection_by_seed.csv');save(pd.DataFrame(edge_rows),'target_protein_bridge_edges.csv')
save(pd.DataFrame(leaveprotein_rows),'leave_one_protein_projection_by_seed.csv')

# Pool within target over measured compartments, using donor-equal target detection weights.
def aggregate_run(z):
 bycomp=z.groupby(['target','cohort','contrast','compartment']).agg(weighted_displacement=('weighted_displacement','median'),matched_z=('matched_z','median'),detection=('target_detection_fraction','first'),n_seeds=('seed','nunique')).reset_index()
 rows=[]
 for key,part in bycomp.groupby(['target','cohort','contrast']):
  w=part.detection.to_numpy(); w=w/w.sum() if w.sum()>0 else np.repeat(1/len(w),len(w))
  rows.append(dict(zip(['target','cohort','contrast'],key))|{'bridge_score':float(part.weighted_displacement@w),'bridge_matched_z_descriptive':float(part.matched_z@w),'n_available_compartments':len(part),'minimum_seeds':int(part.n_seeds.min()),'sum_target_detection_across_available_compartments':float(part.detection.sum())})
 return pd.DataFrame(rows),bycomp
agg,bycomp=aggregate_run(runs);save(agg,'target_projection_summary.csv');save(bycomp,'target_projection_by_compartment.csv')
primary=agg.loc[agg.cohort.eq('confirmatory')&agg.contrast.eq('DO_vs_NDO')].set_index('target').reindex(GENES)
assert primary.bridge_score.notna().all()
# Primary revised feature matrix; three inherited domains remain numerically untouched.
Xraw=pd.DataFrame({'target':GENES,'pharmacology':tar.pharmacology,'clinical_network_bridge':primary.bridge_score.to_numpy(),'HepG2_perturbation':tar.HepG2_perturbation,'STRING_context':tar.STRING_context})
FAMS=Xraw.columns[1:].tolist();X=np.column_stack([score_rank(Xraw[c]) for c in FAMS]);M=np.isfinite(X)
MISSING_FLOOR=0.5/N
Xcon=np.where(M,X,MISSING_FLOOR)
save(Xraw,'integrated_source_features.csv');scoretbl=pd.DataFrame(X,columns=FAMS);scoretbl.insert(0,'target',GENES);save(scoretbl,'integrated_family_percentile_scores.csv')
legacyraw=tar[['pharmacology','liver_context','HepG2_perturbation','STRING_context']].to_numpy()
legacypercent=np.column_stack([score_rank(legacyraw[:,j]) for j in range(4)])

def fuse(x,method,w=None,missing='floor'):
 x=np.asarray(x,float); mask=np.isfinite(x)
 if w is None:w=np.ones(x.shape[-1])
 w=np.asarray(w,float)
 if missing=='floor':xx=np.where(mask,x,MISSING_FLOOR);ww=np.broadcast_to(w,x.shape)
 elif missing=='neutral':xx=np.where(mask,x,0.5);ww=np.broadcast_to(w,x.shape)
 else:xx=np.where(mask,x,0.5);ww=np.broadcast_to(w,x.shape)*mask
 ww=ww/ww.sum(axis=-1,keepdims=True)
 if method=='mean':return np.sum(xx*ww,axis=-1)
 if method=='geometric':return np.exp(np.sum(np.log(np.clip(xx,1e-8,1))*ww,axis=-1))
 if method=='maximin':return np.min(np.where(ww>0,xx,np.inf),axis=-1)
 if method=='median':return np.median(xx,axis=-1)+1e-7*np.sum(xx*ww,axis=-1)
 if method=='RRA':
  # Upper-percentile scores transformed to lower-is-better rank fractions.
  p=np.sort(1-xx,axis=-1); n=p.shape[-1];ks=np.arange(1,n+1)
  raw=np.min(beta.cdf(p,ks,n-ks+1),axis=-1)
  return -np.log10(np.maximum(raw,1e-300))
 raise ValueError(method)
METHODS=['median','mean','geometric','maximin','RRA']
full={m:fuse(X,m) for m in METHODS}; records=[]
for m,s in full.items():
 for g,v,rr in zip(GENES,s,rank_high(s)):records.append({'target':g,'algorithm':m,'score':v,'rank':rr})
for g,r in zip(GENES,tar.consensus_rank):records.append({'target':g,'algorithm':'legacy_consensus','score':-float(r),'rank':r})
rankings=pd.DataFrame(records); save(rankings,'algorithm_candidate_rankings.csv')
# Candidate-identification null: independent source rank permutations preserve each source margin.
# This is a competitive diagnostic conditional on these 17 candidates, not a causal null.
nullX=np.stack([np.column_stack([RNG.permutation(Xcon[:,j]) for j in range(4)]) for _ in range(NFUSION)])
nullr=[]
for method in METHODS:
 ns=fuse(nullX,method)
 obs=full[method];p=(1+(ns>=obs[None,:]).sum(axis=0))/(NFUSION+1)
 for g,o,v,q in zip(GENES,obs,p,bh(p)):nullr.append({'target':g,'algorithm':method,'observed_score':o,'permutation_p':v,'bh_q':q,'permutations':NFUSION,'null_interpretation':'competitive rank coincidence; family-exchangeability sensitivity'})
save(pd.DataFrame(nullr),'algorithm_permutation_calibration.csv')
# Held-source concordance: re-fit fusion with source excluded; score against that source only.
held=[];losor=[]
for j,source in enumerate(FAMS):
 train=np.delete(X,j,axis=1);test=X[:,j];mask=np.isfinite(test)
 for method in METHODS:
  ss=fuse(train,method);rh=corr(ss,test); rr=rank_high(ss)
  for g,r in zip(GENES,rr):losor.append({'target':g,'algorithm':method,'omitted_source':source,'rank':r})
  perm=RNG.permutation(np.tile(np.arange(N),1)) # no fitted labels; RNG stream documented.
  pn=np.array([corr(ss[mask],RNG.permutation(test[mask])) for _ in range(1999)])
  held.append({'algorithm':method,'held_source':source,'n_candidates':int(mask.sum()),'spearman_rho':rh,'two_sided_permutation_p':(1+(abs(pn)>=abs(rh)).sum())/(len(pn)+1),'top5_jaccard_vs_full':jac(top_set(ss),top_set(full[method])),'rank_spearman_vs_full':corr(ss,full[method]),'evaluation':'source agreement, not clinical discrimination'})
# The genuine legacy algorithm undergoes same source withholding on its original features.
for j,source in enumerate(FAMS):
 train=np.delete(legacypercent,j,axis=1);test=legacypercent[:,j];ss=fuse(train,'median');mask=np.isfinite(test)
 held.append({'algorithm':'legacy_consensus','held_source':source,'n_candidates':int(mask.sum()),'spearman_rho':corr(ss,test),'two_sided_permutation_p':np.nan,'top5_jaccard_vs_full':jac(top_set(ss),top_set(-tar.consensus_rank.to_numpy())),'rank_spearman_vs_full':corr(ss,-tar.consensus_rank.to_numpy()),'evaluation':'legacy feature matrix; its liver source differs from revised matrix'})
held=pd.DataFrame(held);save(held,'algorithm_held_source_evaluation.csv');save(pd.DataFrame(losor),'algorithm_leave_one_family_out_ranks.csv')
summ=held.groupby('algorithm').agg(mean_held_source_rho=('spearman_rho','mean'),min_held_source_rho=('spearman_rho','min'),max_held_source_rho=('spearman_rho','max'),mean_loso_rank_rho=('rank_spearman_vs_full','mean'),minimum_loso_rank_rho=('rank_spearman_vs_full','min'),mean_loso_top5_jaccard=('top5_jaccard_vs_full','mean')).reset_index()
save(summ,'algorithm_benchmark_summary.csv')
# Weight perturbation: 5000 symmetric Dirichlet weights, no optimization or best draw selection.
weightdraw=RNG.dirichlet(np.ones(4)*4,size=NWEIGHT)
ws=[]
for method in ['mean','geometric']:
 scores=np.stack([fuse(X,method,w=w) for w in weightdraw]);rr=np.apply_along_axis(rank_high,1,scores)
 for i,g in enumerate(GENES):ws.append({'target':g,'algorithm':method,'n_draws':NWEIGHT,'rank_median':np.median(rr[:,i]),'rank_q025':np.quantile(rr[:,i],.025),'rank_q975':np.quantile(rr[:,i],.975),'top5_frequency':np.mean(rr[:,i]<=5),'first_rank_frequency':np.mean(rr[:,i]==1)})
save(pd.DataFrame(ws),'algorithm_weight_sensitivity.csv')
# Masking mechanisms: observed-only, neutral, pessimistic floor; all genes retained.
miss=[]
for mode in ['floor','neutral','observed_only']:
 for method in METHODS:
  s=fuse(X,method,missing=mode)
  for g,rr in zip(GENES,rank_high(s)):miss.append({'target':g,'algorithm':method,'missing_rule':mode,'rank':rr})
save(pd.DataFrame(miss),'algorithm_missingness_sensitivity.csv')
# Replace clinical-network feature only for coherent sensitivity specifications.
specs=[]
for (cohort,contrast),z in agg.groupby(['cohort','contrast']):
 s=z.set_index('target').reindex(GENES).bridge_score.to_numpy();xx=X.copy();xx[:,1]=score_rank(s)
 for method in METHODS:
  for g,rr in zip(GENES,rank_high(fuse(xx,method))):specs.append({'target':g,'algorithm':method,'specification':cohort+'_'+contrast,'rank':rr})
for seed in seedids:
 a,_=aggregate_run(runs.loc[runs.seed.eq(seed)]);z=a.loc[a.cohort.eq('confirmatory')&a.contrast.eq('DO_vs_NDO')].set_index('target').reindex(GENES)
 xx=X.copy();xx[:,1]=score_rank(z.bridge_score.to_numpy())
 for method in METHODS:
  for g,rr in zip(GENES,rank_high(fuse(xx,method))):specs.append({'target':g,'algorithm':method,'specification':'network_seed_'+str(seed),'rank':rr})
# Leave-one-protein-out, then same compartment/seed aggregation, retains all omitted proteins.
lpo=pd.DataFrame(leaveprotein_rows)
for omit,z in lpo.groupby('omitted_protein'):
 by=z.groupby(['target','compartment']).agg(v=('weighted_displacement','median'),d=('target_detection_fraction','first')).reset_index()
 ss={g:float(p.v@p.d/p.d.sum()) for g,p in by.groupby('target')}
 xx=X.copy();xx[:,1]=score_rank([ss[g] for g in GENES])
 for method in METHODS:
  for g,rr in zip(GENES,rank_high(fuse(xx,method))):specs.append({'target':g,'algorithm':method,'specification':'omit_protein_'+omit,'rank':rr})
specs=pd.DataFrame(specs);save(specs,'algorithm_clinical_seed_protein_sensitivity.csv')
# Stable target set is a transparently descriptive budget-limited shortlist, not significance.
# Primary geometric shortlist top5; intersection with baseline and mean is shown unfiltered.
targetsummary=tar[['gene_symbol','consensus_rank','max_pchembl','n_documents','is_landmark','is_bing']].rename(columns={'gene_symbol':'target','consensus_rank':'legacy_consensus_rank'})
targetsummary=targetsummary.merge(rankings.pivot(index='target',columns='algorithm',values='rank').reset_index().drop(columns='legacy_consensus'),on='target')
targetsummary=targetsummary.merge(primary[['bridge_score','bridge_matched_z_descriptive','n_available_compartments']].reset_index(),on='target')
wss=pd.DataFrame(ws).loc[lambda d:d.algorithm.eq('geometric')].drop(columns='algorithm')
targetsummary=targetsummary.merge(wss,on='target')
for label,filt in [('clinical',specs.specification.str.startswith(('discovery','confirmatory','unweighted'))),('seed',specs.specification.str.startswith('network_seed')),('protein',specs.specification.str.startswith('omit_protein'))]:
 z=specs.loc[filt&specs.algorithm.eq('geometric')].groupby('target')['rank'].agg(['min','median','max']).add_prefix(label+'_rank_').reset_index();targetsummary=targetsummary.merge(z,on='target')
targetsummary['primary_top5']=targetsummary.geometric<=5
# Bridge scores are not tested by aggregate normal approximation.
save(targetsummary.sort_values('geometric'),'candidate_integrated_summary.csv')
# Cross-cohort agreement in projected target scores does not equate to independent replication.
compare=[]
for contrast in ['DO_vs_NDO','DO_vs_HV']:
 wide=agg.loc[agg.contrast.eq(contrast)].pivot(index='target',columns='cohort',values='bridge_score')
 compare.append({'comparison':'discovery_vs_confirmatory_projection_'+contrast,'n':len(wide),'spearman_rho':corr(wide.discovery,wide.confirmatory)})
for method in METHODS:
 compare.append({'comparison':'legacy_vs_'+method,'n':N,'spearman_rho':corr(-tar.consensus_rank,full[method])})
save(pd.DataFrame(compare),'cross_algorithm_projection_concordance.csv')
# Evidence changes selected by workflow definition, never by agreement with favored genes.
save(pd.DataFrame(manifest).drop_duplicates('path'),'input_manifest.csv')
config={'candidate_universe':GENES,'n_candidates':N,'n_human_proteins':len(ANCHORS),'vko_gene_rows':len(vko),'primary_clinical_weight':'absolute confirmatory DO-versus-NDO Cliffs delta, no significance prefilter','primary_compartment_pooling':'target detection-weighted average across available compartments; seed median within compartment','matching_features':['log1p mean normalized expression','logit detection fraction'],'matched_neighbors':30,'matched_null_draws':NNULL,'rank_null_draws':NFUSION,'weight_draws':NWEIGHT,'weight_dirichlet_alpha':[4,4,4,4],'seed':20260919,'missing_floor':MISSING_FLOOR,'primary_fusion':'geometric','primary_fusion_rationale':'balanced multiplicative support across fixed source families; not chosen by held-source performance','algorithms':METHODS,'python':platform.python_version(),'limitations':['All vKO effects are unsigned topology displacement','Clinical protein source is multi-drug DILI, not nintedanib-specific','RRA null is an exchangeability diagnostic; correlated source errors may invalidate nominal calibration','No candidate-specific clinical truth labels, hence no target-prediction AUROC or clinical predictive performance','Rank robustness measures computational reproducibility, not causal target confirmation']}
(OUT/'analysis_configuration.json').write_text(json.dumps(config,indent=2))
print(summ.to_string(index=False));print(targetsummary.sort_values('geometric')[['target','legacy_consensus_rank','geometric','bridge_score','bridge_matched_z_descriptive','top5_frequency','seed_rank_min','seed_rank_max']].to_string(index=False))
