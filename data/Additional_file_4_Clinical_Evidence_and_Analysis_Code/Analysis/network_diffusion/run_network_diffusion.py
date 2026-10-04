#!/usr/bin/env python3
"""Donor-blocked clinical-protein-seeded graph diffusion. No outcome-guided tuning.
Usage: python run_network_diffusion.py --package /path/to/recovered/package --output /path/to/output
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse, hashlib, json, time, platform, itertools
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.linalg import lu_factor, lu_solve
from scipy.stats import false_discovery_control, spearmanr, rankdata
from sklearn.covariance import LedoitWolf
from threadpoolctl import threadpool_limits

CANDIDATES=['FLT4','PDGFRA','FGFR1','KDR','JAK1','LYN','FLT1','FGFR2','FGFR3','ABL1','CDK4','FLT3','FGFR4','TYK2','AXL','LRRK2','MET']
COMPARTMENTS=['Hepatocyte','Endothelial','Macrophage']
ANCHORS=['ACO1','ASS1','FAH','CPS1','ALDOB','HPD','OTC','DMGDH','GSTA1','FBP1','PCK2','CES1','LECT2']
N_NULL=1999
SEED=20260919
GRID=list(itertools.product(['positive','absolute'],[10,20,40],[0.3,0.5,0.7]))
PRIMARY=('positive',20,0.5)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bh(x):
 a=np.asarray(x,float); ok=np.isfinite(a); out=np.full(a.shape,np.nan)
 if ok.any():out[ok]=false_discovery_control(a[ok])
 return out

def weighted_nulls(X,donor,genes,clinical_weights,degree=None,pool_size=50):
 # Equal-donor means and detection fractions, fixed across all network specifications.
 ds=[]; ms=[]
 for d in np.unique(donor):
  v=X[donor==d];ds.append((v>0).mean(0));ms.append(v.mean(0))
 det=np.mean(ds,0);mean=np.mean(ms,0)
 feat=np.c_[np.log1p(mean),np.log((det+0.005)/(1-det+0.005))]
 if degree is not None:feat=np.c_[feat,np.log1p(degree)]
 feat=(feat-feat.mean(0))/np.maximum(feat.std(0),1e-8)
 bg=np.array([i for i,g in enumerate(genes) if g not in ANCHORS and g not in CANDIDATES])
 pools={}
 for g in ANCHORS:
  if g in genes:
   idx=genes.index(g);dist=((feat[bg]-feat[idx])**2).sum(1)
   pools[g]=bg[np.argsort(dist,kind='stable')[:pool_size]]
 rng=np.random.default_rng(SEED)
 null_indices=[]
 present=[g for g in ANCHORS if g in genes]
 for b in range(N_NULL):
  chosen=[]
  for g in present:
   opts=np.setdiff1d(pools[g],chosen,assume_unique=False)
   chosen.append(int(rng.choice(opts)))
  null_indices.append(chosen)
 null_indices=np.asarray(null_indices)
 matrices={}; weights_table=[]
 for label,w in clinical_weights.items():
  wv=np.array([w[g] for g in present]);wv=wv/wv.sum()
  v=np.zeros(len(genes));v[[genes.index(g) for g in present]]=wv
  null=np.zeros((len(genes),N_NULL))
  for j in range(len(present)):null[null_indices[:,j],np.arange(N_NULL)]+=wv[j]
  matrices[label]=(v,null)
  for g,ww in zip(present,wv):weights_table.append(dict(weight_profile=label,gene_symbol=g,weight=ww))
 pool_rows=[]
 for g,p in pools.items():
  for i in p:pool_rows.append(dict(anchor_gene=g,matched_gene=genes[i],anchor_mean_expression=mean[genes.index(g)],matched_mean_expression=mean[i],anchor_detection_fraction=det[genes.index(g)],matched_detection_fraction=det[i]))
 return matrices,mean,det,pool_rows,weights_table

def donor_covariances(X,donors):
 covs={};rawcovs={};qc=[]
 for d in np.unique(donors):
  x=X[donors==d]; sd=x.std(0)
  z=(x-x.mean(0))/np.where(sd>0,sd,1)
  lw=LedoitWolf(assume_centered=True).fit(z)
  c=lw.covariance_;q=np.sqrt(np.diag(c));c=c/np.outer(q,q);np.fill_diagonal(c,1)
  covs[d]=c
  rawcovs[d]=(z.T@z)/len(z)
  qc.append(dict(donor=d,n_cells=len(x),n_nonconstant_genes=int((sd>0).sum()),shrinkage=float(lw.shrinkage_)))
 return covs,rawcovs,qc

def make_graph(C,edge,k):
 W=np.maximum(C,0) if edge=='positive' else np.abs(C)
 np.fill_diagonal(W,0)
 keep=np.zeros(W.shape,bool)
 ix=np.argpartition(W,-k,axis=1)[:,-k:]
 keep[np.arange(len(W))[:,None],ix]=True
 W=np.where(keep|keep.T,W,0)
 deg=W.sum(1);isolated=deg<=0
 W[isolated,isolated]=1;deg=W.sum(1)
 return W,deg

def bridge_scores(W,deg,restart,genes,profiles):
 valid=[g for g in CANDIDATES if g in genes]; ix=[genes.index(g) for g in valid]
 # f = r * (I-(1-r)P)^-1 seed, P column-stochastic.
 P=W/deg[None,:]
 A=np.eye(len(W))-(1-restart)*P
 rhs=np.eye(len(W))[:,ix]
 solution=lu_solve(lu_factor(A.T,check_finite=False),rhs,check_finite=False)
 residual=float(np.max(np.abs(A.T@solution-rhs)))
 assert residual<1e-10 and np.max(np.abs(P.sum(0)-1))<1e-12
 reach=restart*solution.T
 res=[]
 for profile,(v,null) in profiles.items():
  obs=reach@v;nv=reach@null
  mu=nv.mean(1);sd=nv.std(1,ddof=1)
  z=(obs-mu)/np.where(sd>0,sd,np.nan)
  pp=(1+(nv>=obs[:,None]).sum(1))/(1+N_NULL)
  for j,g in enumerate(valid):res.append(dict(gene_symbol=g,weight_profile=profile,diffusion_score=obs[j],null_mean=mu[j],null_sd=sd[j],matched_null_z=z[j],empirical_p_enrichment=pp[j],linear_solve_max_abs_residual=residual,weighted_degree=deg[ix[j]],n_nonzero_edges=int((W[ix[j]]>0).sum())))
 # primary clinical anchor contributions, additive to observed score
 primary=profiles['discovery_DO_vs_HV'][0]
 contrib=[]
 for j,g in enumerate(valid):
  for a in ANCHORS:
   if a in genes:
    ai=genes.index(a);contrib.append(dict(gene_symbol=g,anchor_gene=a,seed_weight=primary[ai],reachability_coefficient=reach[j,ai],weighted_contribution=reach[j,ai]*primary[ai]))
 return res,contrib

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--package',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
 package=Path(args.package).resolve();out=Path(args.output).resolve();out.mkdir(exist_ok=True,parents=True)
 t=time.time();src=package/'Source_Data'
 effects=pd.read_csv(src/'Human_DILI_Proteomics/effect_estimates.csv')
 weights={}
 for cohort,contrast in itertools.product(['discovery','confirmatory'],['DO_vs_HV','DO_vs_NDO']):
  v=effects[effects.cohort.eq(cohort)&effects.contrast.eq(contrast)].set_index('protein').cliffs_delta_group1_minus_group2.abs()
  weights[f'{cohort}_{contrast}']=v.to_dict()
 meta=pd.read_csv(src/'Liver_Single_Cell/cell_metadata.csv').set_index('CellName')
 all_rows=[];lodo_rows=[];nullpool=[];qcs=[];wr=[];contributions=[];availability=[];inputs=[];ablations=[];matching_sensitivity=[]
 primary_profile='discovery_DO_vs_HV'
 for comp in COMPARTMENTS:
  path=src/f'Liver_Single_Cell/{comp.lower()}_network_matrix.csv.gz';d=pd.read_csv(path,index_col=0);genes=d.index.tolist();X=d.T.to_numpy(float)
  donors=meta.loc[d.columns,'Sample'].to_numpy()
  assert meta.loc[d.columns,'compartment'].eq(comp).all()
  assert np.isfinite(X).all() and len(np.unique(donors))==5
  inputs.append(dict(file=str(path.relative_to(package)),sha256=sha(path),n_genes=X.shape[1],n_cells=X.shape[0],n_donors=len(np.unique(donors))))
  profiles,mean,detect,pools,wt=weighted_nulls(X,donors,genes,weights)
  nullpool.extend(dict(compartment=comp,**x) for x in pools);wr.extend(dict(compartment=comp,**x) for x in wt)
  covs,rawcovs,qc=donor_covariances(X,donors);qcs.extend(dict(compartment=comp,**x) for x in qc)
  for g in CANDIDATES:availability.append(dict(compartment=comp,gene_symbol=g,network_available=g in genes,donor_mean_expression=mean[genes.index(g)] if g in genes else np.nan,donor_mean_detection=detect[genes.index(g)] if g in genes else np.nan))
  C=np.mean(list(covs.values()),axis=0)
  Wprimary,dprimary=make_graph(C,PRIMARY[0],PRIMARY[1])
  for pool_size in [25,50,100]:
   degree_profiles,_,_,_,_=weighted_nulls(X,donors,genes,weights,degree=dprimary,pool_size=pool_size)
   rr,_=bridge_scores(Wprimary,dprimary,PRIMARY[2],genes,degree_profiles)
   matching_sensitivity.extend(dict(compartment=comp,null_match='expression_detection_degree',match_pool_size=pool_size,**r) for r in rr)
  Wraw,draw=make_graph(np.mean(list(rawcovs.values()),axis=0),PRIMARY[0],PRIMARY[1])
  rr,_=bridge_scores(Wraw,draw,PRIMARY[2],genes,profiles)
  ablations.extend(dict(compartment=comp,covariance_estimator='unshrunk_Pearson',**r) for r in rr)
  for edge,k,restart in GRID:
   W,degree=make_graph(C,edge,k)
   rows,cont=bridge_scores(W,degree,restart,genes,profiles)
   spec=f'{edge}_k{k}_r{restart}'
   all_rows.extend(dict(compartment=comp,specification=spec,edge_mode=edge,k_neighbors=k,restart_probability=restart,is_primary_spec=(edge,k,restart)==PRIMARY,**r) for r in rows)
   if (edge,k,restart)==PRIMARY:contributions.extend(dict(compartment=comp,**r) for r in cont)
  for omit in sorted(covs):
   C=np.mean([v for key,v in covs.items() if key!=omit],axis=0)
   # The primary network setting is fixed before analysis and subjected to LODO.
   keep=donors!=omit
   lodo_profiles,_,_,_,_=weighted_nulls(X[keep],donors[keep],genes,weights)
   W,degree=make_graph(C,PRIMARY[0],PRIMARY[1]);rows,_=bridge_scores(W,degree,PRIMARY[2],genes,{primary_profile:lodo_profiles[primary_profile]})
   lodo_rows.extend(dict(compartment=comp,omitted_donor=omit,**r) for r in rows)
  print(comp,'complete',round(time.time()-t,1),'sec',flush=True)
 scores=pd.DataFrame(all_rows);lodo=pd.DataFrame(lodo_rows)
 # Each specification/profile is a sensitivity analysis; within it all eligible target-compartment hypotheses form one family.
 scores['empirical_q_bh']=scores.groupby(['specification','weight_profile']).empirical_p_enrichment.transform(bh)
 scores['candidate_rank']=scores.groupby(['compartment','specification','weight_profile']).matched_null_z.rank(ascending=False,method='average')
 scores['candidate_rank_fraction']=scores.groupby(['compartment','specification','weight_profile']).matched_null_z.transform(lambda x:x.rank(ascending=False,method='average')/x.notna().sum())
 lodo['empirical_q_bh']=lodo.groupby('omitted_donor').empirical_p_enrichment.transform(bh)
 lodo['candidate_rank']=lodo.groupby(['compartment','omitted_donor']).matched_null_z.rank(ascending=False,method='average')
 lodo['candidate_rank_fraction']=lodo.groupby(['compartment','omitted_donor']).matched_null_z.transform(lambda x:x.rank(ascending=False,method='average')/x.notna().sum())
 prim=scores[scores.is_primary_spec&scores.weight_profile.eq(primary_profile)].copy()
 stab=scores[scores.weight_profile.eq(primary_profile)].groupby(['compartment','gene_symbol']).agg(n_specs=('matched_null_z','size'),median_spec_z=('matched_null_z','median'),min_spec_z=('matched_null_z','min'),max_spec_z=('matched_null_z','max'),spec_top_quartile_frequency=('candidate_rank_fraction',lambda x:float((x<=.25).mean())),spec_q_lt_0_05_frequency=('empirical_q_bh',lambda x:float((x<.05).mean()))).reset_index()
 ds=lodo.groupby(['compartment','gene_symbol']).agg(n_lodo=('matched_null_z','size'),median_lodo_z=('matched_null_z','median'),min_lodo_z=('matched_null_z','min'),max_lodo_z=('matched_null_z','max'),lodo_top_quartile_count=('candidate_rank_fraction',lambda x:int((x<=.25).sum())),lodo_q_lt_0_05_count=('empirical_q_bh',lambda x:int((x<.05).sum())),median_lodo_rank_fraction=('candidate_rank_fraction','median')).reset_index()
 result=pd.DataFrame(availability).merge(prim,on=['compartment','gene_symbol'],how='left').merge(stab,on=['compartment','gene_symbol'],how='left').merge(ds,on=['compartment','gene_symbol'],how='left')
 result['prespecified_stability_gate']=result.network_available & result.empirical_q_bh.lt(.05)&result.spec_top_quartile_frequency.ge(.75)&result.lodo_top_quartile_count.ge(4)
 result['inherited_expression_support_gate']=result.donor_mean_detection.ge(.10)
 result['expression_supported_stability_gate']=result.prespecified_stability_gate & result.inherited_expression_support_gate
 result['prespecified_rank_robust_gate']=result.network_available & result.spec_top_quartile_frequency.ge(.75)&result.lodo_top_quartile_count.ge(4)
 old=pd.read_csv(src/'Virtual_Knockout/virtual_knockout_seed_summary.csv');old=old[old.ko_class.eq('drug_target')&old.ko_label.isin(CANDIDATES)]
 oldsum=old.groupby(['compartment','ko_label']).agg(old_scTenifold_median_anchor_percentile=('dili_anchor_mean_percentile','median'),old_scTenifold_min_anchor_percentile=('dili_anchor_mean_percentile','min'),old_scTenifold_max_anchor_percentile=('dili_anchor_mean_percentile','max'),old_n_seeds=('seed','nunique')).reset_index().rename(columns={'ko_label':'gene_symbol'})
 compare=result.merge(oldsum,on=['compartment','gene_symbol'],how='left')
 oldnull=pd.read_csv(src/'Integration/virtual_ko_expression_matched_null_stability.csv').rename(columns={'ko_label':'gene_symbol'})
 compare=compare.merge(oldnull[['compartment','gene_symbol','median_matched_z']],on=['compartment','gene_symbol'],how='left').rename(columns={'median_matched_z':'old_scTenifold_matched_null_z'})
 metrics=[]
 for comp in COMPARTMENTS:
  a=compare[compare.compartment.eq(comp)]
  for field in ['old_scTenifold_median_anchor_percentile','old_scTenifold_matched_null_z']:
   b=a[['matched_null_z',field]].dropna();r,p=spearmanr(b.matched_null_z,b[field]);metrics.append(dict(compartment=comp,comparison=field,n_matched_candidates=len(b),spearman_rho=r,spearman_p_descriptive=p))
 # LODO rank correlation and source-weight agreement are descriptive, not predictive validation.
 for comp in COMPARTMENTS:
  for omit in sorted(lodo.omitted_donor.unique()):
   b=prim[prim.compartment.eq(comp)][['gene_symbol','matched_null_z']].merge(lodo[lodo.compartment.eq(comp)&lodo.omitted_donor.eq(omit)][['gene_symbol','matched_null_z']],on='gene_symbol',suffixes=('_full','_lodo'))
   b=b.dropna();r,p=spearmanr(b.matched_null_z_full,b.matched_null_z_lodo);metrics.append(dict(compartment=comp,comparison=f'LODO_{omit}',n_matched_candidates=len(b),spearman_rho=r,spearman_p_descriptive=p))
 ablation=pd.DataFrame(ablations)
 ablation['empirical_q_bh']=ablation.groupby('weight_profile').empirical_p_enrichment.transform(bh)
 ablation['candidate_rank']=ablation.groupby(['compartment','weight_profile']).matched_null_z.rank(ascending=False)
 ablation=ablation.merge(scores[scores.is_primary_spec][['compartment','gene_symbol','weight_profile','matched_null_z','candidate_rank']],on=['compartment','gene_symbol','weight_profile'],suffixes=('_unshrunk','_shrinkage'))
 match=pd.DataFrame(matching_sensitivity)
 match['empirical_q_bh']=match.groupby(['null_match','match_pool_size','weight_profile']).empirical_p_enrichment.transform(bh)
 pharm=pd.read_csv(src/'Pharmacology/chembl_target_summary.csv')[['gene_symbol','max_pchembl','median_pchembl','n_documents','n_assays','pharmacology_tier']]
 result=result.merge(pharm,on='gene_symbol',how='left')
 outputs={'network_null_matching_sensitivity.csv':match,'network_shrinkage_ablation.csv':ablation,'candidate_network_summary.csv':result,'network_specification_scores.csv':scores,'network_donor_lodo_scores.csv':lodo,'network_old_algorithm_comparison.csv':compare,'network_algorithm_comparison_metrics.csv':pd.DataFrame(metrics),'network_expression_matched_null_pools.csv':pd.DataFrame(nullpool),'network_clinical_seed_weights.csv':pd.DataFrame(wr),'network_donor_covariance_qc.csv':pd.DataFrame(qcs),'network_target_anchor_contributions.csv':pd.DataFrame(contributions)}
 for name,table in outputs.items():table.to_csv(out/name,index=False)
 summary=dict(method='donor_equal_LedoitWolf_clinical_seed_RWR',n_null=N_NULL,random_seed=SEED,candidate_list=CANDIDATES,clinical_seeds=ANCHORS,primary_clinical_weight='absolute discovery DO_vs_HV Cliff delta normalized to sum one',primary_spec=dict(edge='positive',k_neighbors=20,restart=.5),specification_count=18,inputs=inputs,n_eligible_target_compartment=int(result.network_available.sum()),n_primary_q_lt_0_05=int(result.empirical_q_bh.lt(.05).sum()),n_expression_supported_stability_passes=int(result.expression_supported_stability_gate.sum()),max_linear_solve_residual=float(scores.linear_solve_max_abs_residual.max()),stability_passes=result.loc[result.prespecified_stability_gate,['compartment','gene_symbol','matched_null_z','empirical_q_bh','spec_top_quartile_frequency','lodo_top_quartile_count']].to_dict('records'),rank_robust_candidates=result.loc[result.prespecified_rank_robust_gate,['compartment','gene_symbol','matched_null_z','empirical_q_bh','spec_top_quartile_frequency','lodo_top_quartile_count']].to_dict('records'),execution_seconds=time.time()-t,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,interpretation='Expression-adjusted proximity to a human DILI protein signature; undirected and noncausal. LODO is within-atlas donor robustness, not independent validation. Fixed prefiltered gene universe.')
 (out/'network_diffusion_summary.json').write_text(json.dumps(summary,indent=2))
 print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=1):main()
