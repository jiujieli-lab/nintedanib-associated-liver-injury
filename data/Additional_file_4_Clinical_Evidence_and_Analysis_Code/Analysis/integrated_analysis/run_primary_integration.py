#!/usr/bin/env python3
"""Source-aware continuous integration of donor-balanced diffusion with public data.
Primary algorithm fixed before outcome inspection; all 17 pharmacology candidates retained.
"""
from pathlib import Path
import argparse,hashlib,json,platform
import numpy as np,pandas as pd
from scipy.stats import rankdata,beta,spearmanr
p=argparse.ArgumentParser();p.add_argument('--input-root',type=Path,required=True);p.add_argument('--network-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'primary');a=p.parse_args()
O=a.output_dir;O.mkdir(parents=True,exist_ok=True);RNG=np.random.default_rng(20260919);NNULL=20000;NW=5000
manifest=[]
def read(path):
 manifest.append({'filename':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size});return pd.read_csv(path)
def save(x,f):x.to_csv(O/f,index=False)
def bh(p):
 p=np.asarray(p);o=np.argsort(p);q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];z=np.empty_like(q);z[o]=np.minimum(1,q);return z
def rankscore(x):
 x=np.asarray(x,float);m=np.isfinite(x);z=np.full_like(x,np.nan);z[m]=(rankdata(x[m])-.5)/m.sum();return z
def rankhigh(x):return rankdata(-np.asarray(x))
def corr(x,y):
 x=np.asarray(x);y=np.asarray(y);m=np.isfinite(x)&np.isfinite(y)
 return float(spearmanr(x[m],y[m]).statistic) if m.sum()>2 else np.nan
def top(x):return set(np.flatnonzero(rankhigh(x)<=5))
def jac(x,y):return len(x&y)/len(x|y) if len(x|y) else np.nan
old=read(a.input_root/'Source_Data/Integration/candidate_consensus_ranking.csv');old=old.loc[old.candidate_class.eq('exposure_proximal_target')].sort_values('gene_symbol').reset_index(drop=True)
G=old.gene_symbol.tolist();N=len(G);assert N==17
net=read(a.network_dir/'candidate_network_summary.csv');spec=read(a.network_dir/'network_specification_scores.csv');lodo=read(a.network_dir/'network_donor_lodo_scores.csv')
donor_detection=read(a.input_root/'Source_Data/Liver_Single_Cell/donor_level_target_expression.csv')
# Aggregate matched z only descriptively. Matched z is never assigned an aggregate z-test.
detection=net[['gene_symbol','compartment','donor_mean_detection']].copy()
def pool(z):
 if 'donor_mean_detection' not in z:z=z.merge(detection,on=['gene_symbol','compartment'],validate='many_to_one')
 rows=[]
 for g in G:
  q=z.loc[z.gene_symbol.eq(g)].dropna(subset=['matched_null_z','donor_mean_detection']);d=q.donor_mean_detection.to_numpy();ww=d/d.sum() if d.sum()>0 else np.ones(len(d))/max(1,len(d))
  rows.append({'target':g,'detection_weighted_network_z':float(q.matched_null_z@ww) if len(q) else np.nan,'n_available_compartments':len(q),'maximum_detection_fraction':float(d.max()) if len(d) else np.nan,'sum_compartment_detection':float(d.sum())})
 return pd.DataFrame(rows)
network=pool(net);save(network,'network_target_aggregate.csv')
F=['pharmacology','clinical_network','HepG2_perturbation','STRING_context']
raw=pd.DataFrame({'target':G,'pharmacology':old.pharmacology,'clinical_network':network.detection_weighted_network_z,'HepG2_perturbation':old.HepG2_perturbation,'STRING_context':old.STRING_context})
X=np.column_stack([rankscore(raw[c]) for c in F]);floor=.5/N;METHOD=['median','mean','geometric','maximin','RRA']
def fuse(x,m,w=None,missing='floor'):
 x=np.asarray(x,float);ms=np.isfinite(x);ww=np.ones(x.shape[-1]) if w is None else np.asarray(w)
 xx=np.where(ms,x,floor if missing=='floor' else .5);ww=np.broadcast_to(ww,x.shape)*(ms if missing=='observed_only' else 1);ww=ww/ww.sum(axis=-1,keepdims=True)
 if m=='mean':return np.sum(xx*ww,axis=-1)
 if m=='geometric':return np.exp(np.sum(np.log(xx)*ww,axis=-1))
 if m=='median':return (np.nanmedian(np.where(ms,xx,np.nan),axis=-1) if missing=='observed_only' else np.median(xx,axis=-1))+1e-7*np.sum(xx*ww,axis=-1)
 if m=='maximin':return np.min(np.where(ww>0,xx,np.inf),axis=-1)
 if m=='RRA':
  if missing=='observed_only':
   result=np.empty(x.shape[:-1])
   for idx in np.ndindex(x.shape[:-1]):
    vv=np.sort(1-x[idx][ms[idx]]);n=len(vv);result[idx]=-np.log10(np.clip(np.min(beta.cdf(vv,np.arange(1,n+1),np.arange(n,0,-1))),1e-300,1))
   return result
  xx=np.sort(1-xx,axis=-1);n=xx.shape[-1];return -np.log10(np.clip(np.min(beta.cdf(xx,np.arange(1,n+1),np.arange(n,0,-1)),axis=-1),1e-300,1))
full={m:fuse(X,m) for m in METHOD}
rows=[]
for m in METHOD:
 for g,s,r in zip(G,full[m],rankhigh(full[m])):rows.append({'target':g,'algorithm':m,'score':s,'rank':r})
for g,r in zip(G,old.consensus_rank):rows.append({'target':g,'algorithm':'legacy_consensus','score':-r,'rank':r})
rankings=pd.DataFrame(rows);save(rankings,'algorithm_candidate_rankings.csv');save(raw,'integrated_source_features.csv')
scores=pd.DataFrame(X,columns=F);scores.insert(0,'target',G);save(scores,'integrated_family_percentile_scores.csv')
# Family-aware competitive permutation: clinical-network and STRING shuffled jointly,
# preserving their empirical dependence; P and HepG2 are separate blocks. Missingness fixed.
B=np.empty((NNULL,N,4));mask=np.isfinite(X)
for b in range(NNULL):
 xx=X.copy()
 for j in [0,2]:
  idx=np.flatnonzero(mask[:,j]);xx[idx,j]=X[RNG.permutation(idx),j]
 for pattern in set(map(tuple,mask[:,[1,3]])):
  idx=np.flatnonzero((mask[:,[1,3]]==np.array(pattern)).all(axis=1));xx[np.ix_(idx,[1,3])]=X[np.ix_(RNG.permutation(idx),[1,3])]
 B[b]=xx
null=[]
for m in METHOD:
 ss=fuse(B,m);p=(1+(ss>=full[m]).sum(axis=0))/(NNULL+1)
 for g,pp,qq in zip(G,p,bh(p)):null.append({'target':g,'algorithm':m,'competitive_permutation_p':pp,'bh_q':qq,'n_permutations':NNULL,'calibration':'network pair block permutation; candidate missingness preserved'})
save(pd.DataFrame(null),'algorithm_block_permutation_calibration.csv')
# Feature-family omission is not an independent patient validation. The entire
# clinical-network composite is omitted together so no clinical weights remain in that fold.
held=[];loso=[]
def heldp(x,y):
 mm=np.isfinite(x)&np.isfinite(y);xx=rankdata(x[mm]);yy=rankdata(y[mm]);xx-=xx.mean();yy-=yy.mean();den=np.sqrt((xx@xx)*(yy@yy));obs=xx@yy/den
 perms=np.vstack([RNG.permutation(yy) for _ in range(4999)]);v=perms@xx/den
 return float(obs),float((1+(abs(v)>=abs(obs)-1e-12).sum())/5000)
for j,f in enumerate(F):
 xx=np.delete(X,j,axis=1)
 for m in METHOD:
  ss=fuse(xx,m);rh,pv=heldp(ss,X[:,j]);rr=rankhigh(ss)
  held.append({'algorithm':m,'omitted_family':f,'held_feature_spearman':rh,'permutation_p':pv,'rank_spearman_vs_full':corr(ss,full[m]),'top5_jaccard_vs_full':jac(top(ss),top(full[m])),'n_held_values':int(np.isfinite(X[:,j]).sum())})
  for g,r in zip(G,rr):loso.append({'target':g,'algorithm':m,'omitted_family':f,'rank':r})
# Legacy reconstructed as median of descending ranks with average-rank tiebreaker,
# including bottom-rank missingness exactly as its source code.
LEG=old[['pharmacology','liver_context','HepG2_perturbation','STRING_context']]
def legacy(z):
 r=z.rank(ascending=False,method='average',na_option='bottom');tu=list(zip(r.median(axis=1),r.mean(axis=1)));return -pd.Series(tu).rank(method='dense').to_numpy()
assert np.array_equal(-legacy(LEG),old.consensus_rank.to_numpy())
for j,f in enumerate(F):
 ss=legacy(LEG.drop(columns=LEG.columns[j]));test=LEG.iloc[:,j].to_numpy();rh,pv=heldp(ss,test)
 held.append({'algorithm':'legacy_consensus','omitted_family':f,'held_feature_spearman':rh,'permutation_p':pv,'rank_spearman_vs_full':corr(ss,-old.consensus_rank),'top5_jaccard_vs_full':jac(top(ss),top(-old.consensus_rank)),'n_held_values':int(np.isfinite(test).sum())})
held=pd.DataFrame(held);save(held,'algorithm_family_omission_evaluation.csv');save(pd.DataFrame(loso),'algorithm_family_omission_ranks.csv')
summary=held.groupby('algorithm').agg(mean_held_feature_rho=('held_feature_spearman','mean'),minimum_held_feature_rho=('held_feature_spearman','min'),maximum_held_feature_rho=('held_feature_spearman','max'),mean_omission_rank_rho=('rank_spearman_vs_full','mean'),minimum_omission_rank_rho=('rank_spearman_vs_full','min'),mean_omission_top5_jaccard=('top5_jaccard_vs_full','mean')).reset_index()
save(summary,'algorithm_benchmark_summary.csv')
# Hyperparameters, clinical seeds, donors: recompute pooled network before fusion.
sensitivity=[];poolrows=[]
for (sp,weight),z in spec.groupby(['specification','weight_profile']):
 agg=pool(z);xx=X.copy();xx[:,1]=rankscore(agg.detection_weighted_network_z)
 for rr in agg.to_dict('records'):poolrows.append(rr|{'specification':sp,'weight_profile':weight})
 for m in METHOD:
  for g,r in zip(G,rankhigh(fuse(xx,m))):sensitivity.append({'target':g,'algorithm':m,'sensitivity_type':'network_specification','specification':sp+'__'+weight,'rank':r})
for (donor,weight),z in lodo.groupby(['omitted_donor','weight_profile']):
 train_detection=donor_detection.loc[~donor_detection.donor.eq(donor)].groupby(['gene_symbol','compartment'],as_index=False).detection_fraction.mean().rename(columns={'detection_fraction':'donor_mean_detection'})
 z=z.merge(train_detection,on=['gene_symbol','compartment'],validate='many_to_one')
 agg=pool(z);xx=X.copy();xx[:,1]=rankscore(agg.detection_weighted_network_z)
 for m in METHOD:
  for g,r in zip(G,rankhigh(fuse(xx,m))):sensitivity.append({'target':g,'algorithm':m,'sensitivity_type':'donor_omission','specification':str(donor)+'__'+weight,'rank':r})
# Missingness sensitivity applies to all 17 genes and does not reinterpret absence as zero biology.
for rule in ['floor','neutral','observed_only']:
 for m in METHOD:
  for g,r in zip(G,rankhigh(fuse(X,m,missing=rule))):sensitivity.append({'target':g,'algorithm':m,'sensitivity_type':'missingness','specification':rule,'rank':r})
# Measured-only perturbation sensitivity; inferred transcript values masked together.
xmeasure=X.copy();xmeasure[~old.is_landmark.eq(1).to_numpy(),2]=np.nan
for m in METHOD:
 for g,r in zip(G,rankhigh(fuse(xmeasure,m))):sensitivity.append({'target':g,'algorithm':m,'sensitivity_type':'measurement','specification':'HepG2_landmarks_only','rank':r})
sens=pd.DataFrame(sensitivity);save(sens,'algorithm_specification_donor_missingness_sensitivity.csv');save(pd.DataFrame(poolrows),'network_specification_target_aggregate.csv')
ws=[];draws=RNG.dirichlet([4,4,4,4],size=NW)
for m in ['mean','geometric']:
 rr=np.vstack([rankhigh(fuse(X,m,w)) for w in draws])
 for i,g in enumerate(G):ws.append({'target':g,'algorithm':m,'n_weight_draws':NW,'rank_median':np.median(rr[:,i]),'rank_q025':np.quantile(rr[:,i],.025),'rank_q975':np.quantile(rr[:,i],.975),'top5_frequency':np.mean(rr[:,i]<=5),'first_rank_frequency':np.mean(rr[:,i]==1)})
ws=pd.DataFrame(ws);save(ws,'algorithm_weight_sensitivity.csv')
result=old[['gene_symbol','consensus_rank','max_pchembl','n_documents','is_landmark','is_bing']].rename(columns={'gene_symbol':'target','consensus_rank':'legacy_consensus_rank'})
result=result.merge(rankings.pivot(index='target',columns='algorithm',values='rank').reset_index().drop(columns='legacy_consensus'),on='target').merge(network,on='target').merge(ws.loc[ws.algorithm.eq('geometric')].drop(columns='algorithm'),on='target')
for ty in sens.sensitivity_type.unique():
 z=sens.loc[sens.algorithm.eq('geometric')&sens.sensitivity_type.eq(ty)].groupby('target')['rank'].agg(['min','median','max']).add_prefix(ty+'_rank_').reset_index();result=result.merge(z,on='target')
result['primary_experimental_priority_top5']=result.geometric<=5
result['causal_confirmation']='not established';result['primary_score_inference']='relative experimental priority; not DILI mediation probability'
save(result.sort_values('geometric'),'candidate_integrated_summary.csv')
comparison=[]
for m in METHOD:comparison.append({'comparison':'legacy_vs_'+m,'n':N,'spearman_rank_rho':corr(-old.consensus_rank,full[m]),'top5_jaccard':jac(top(-old.consensus_rank),top(full[m]))})
save(pd.DataFrame(comparison),'algorithm_pairwise_comparison.csv')
# Sensitivity rank correlations are computed against fixed primary, never used to select its settings.
ss=[]
for (m,ty,sp),z in sens.groupby(['algorithm','sensitivity_type','specification']):
 rr=z.set_index('target').reindex(G)['rank'].to_numpy();fr=rankhigh(full[m]);ss.append({'algorithm':m,'sensitivity_type':ty,'specification':sp,'rank_rho_vs_primary':corr(rr,fr),'top5_jaccard_vs_primary':jac(set(np.flatnonzero(rr<=5)),set(np.flatnonzero(fr<=5)))})
save(pd.DataFrame(ss),'algorithm_sensitivity_metrics.csv')
save(pd.DataFrame(manifest),'input_manifest.csv')
config={'primary_algorithm':'geometric percentile fusion','primary_network':'donor-balanced LedoitWolf positive k20 restart0.5, discovery DOHV absolute Cliff weights','target_universe':G,'n_targets':N,'n_human_proteins':13,'no_candidate_specific_tuning':True,'pooling':'target-detection-weighted matched z across available compartments, descriptive score only','family_columns':F,'competitive_permutations':NNULL,'permutation_blocks':[['pharmacology'],['HepG2_perturbation'],['clinical_network','STRING_context']],'permutation_missingness':'each candidate measurement mask remains fixed','weight_draws':NW,'Dirichlet_alpha':[4,4,4,4],'random_seed':20260919,'missing_floor':floor,'primary_selection':'geometric rule chosen for balanced support, not benchmark performance','donor_omission_pooling':'target detection weights recalculated from retained donors only','evaluation':'component omission, donor robustness and ranking stability; no patient-level predictive validation','python':platform.python_version()}
(O/'analysis_configuration.json').write_text(json.dumps(config,indent=2))
print(summary.to_string(index=False));print(result.sort_values('geometric')[['target','legacy_consensus_rank','geometric','detection_weighted_network_z','top5_frequency','donor_omission_rank_min','donor_omission_rank_max']].to_string(index=False))
