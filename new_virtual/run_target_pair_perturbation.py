#!/usr/bin/env python3
"""Exact graph interventions on inherited clinical-weighted donor-balanced networks.
These are structural network sensitivity experiments, not biological causal effects.
All parameters are frozen before primary-output computation. No target-selection tuning.
"""
import os
for v in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
 os.environ[v]='1'
from pathlib import Path
import json, hashlib, time, platform, itertools
import numpy as np
import pandas as pd
from scipy.linalg import solve
from scipy.stats import false_discovery_control
from sklearn.covariance import LedoitWolf
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent
SC=ROOT/'data/Additional_file_5_Liver_and_Lung_Single_Cell_Data/Source_Data/Liver_Single_Cell'
NET=ROOT/'data/Additional_file_4_Clinical_Evidence_and_Analysis_Code/Analysis/network_diffusion'
TARGS=['PDGFRA','FLT4']
COMPS=['Hepatocyte','Endothelial','Macrophage']
CANDIDATES=['FLT4','PDGFRA','FGFR1','KDR','JAK1','LYN','FLT1','FGFR2','FGFR3','ABL1','CDK4','FLT3','FGFR4','TYK2','AXL','LRRK2','MET']
ANCHORS=['ACO1','ASS1','FAH','CPS1','ALDOB','HPD','OTC','DMGDH','GSTA1','FBP1','PCK2','CES1','LECT2']
CONFIG={'analysis':'target-pair structural perturbation of inherited clinical-weighted RWR graph','primary_graph':{'covariance':'within-donor LedoitWolf of standardized genes; equal donor weighting','edges':'positive correlation; union of each node top20 neighbors','restart':0.5},'targets':TARGS,'comparators':['KDR','FGFR1'],'compartments':COMPS,'clinical_weights':'inherited discovery DO_vs_HV absolute Cliffs delta, normalized','primary_metric':'clinical-weighted absolute propagation-score change divided by clinical-weighted baseline score','intervention':'remove all incident edges of selected targets, add isolated-node self-loop, re-normalize columns; retain clinical restart seeds','nonadditivity':'f_AB-f_A-f_B+f_0; clinical weighted L1 relative to baseline','null_match_features':['log1p equal-donor mean expression','logit equal-donor detection with 0.005 pseudocount','log1p weighted degree'],'null_pool_size':50,'pair_null_draws':1999,'null_pool_exclusions':ANCHORS+CANDIDATES,'seed':20261002,'multiplicity':'BH across 3 interventions x 3 compartments for magnitude; separately 3 paired nonadditivity tests','donor_sensitivity':'all five leave-one-donor-out graph reconstructions; fixed inherited gene universe','attenuation_grid':[0,0.25,0.5,0.75,1],'attenuation_rule':'each edge incident to one or both target nodes is scaled once by 1-fraction','interpretation':'Undirected transcript covariance intervention is not receptor inhibition, State prediction, clinical causal identification, synergy, or a biological rescue.'}
(OUT/'analysis_specification.json').write_text(json.dumps(CONFIG,indent=2))

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def graph(C):
 W=np.maximum(C,0);np.fill_diagonal(W,0)
 ix=np.argpartition(W,-20,axis=1)[:,-20:]
 keep=np.zeros_like(W,dtype=bool);keep[np.arange(len(W))[:,None],ix]=True
 W=np.where(keep|keep.T,W,0)
 W[np.diag_indices_from(W)]=(W.sum(0)<=0).astype(float)
 return W

def propagate(W,v,indices=(),fraction=1.):
 V=W.copy()
 mask=np.zeros_like(V,dtype=bool)
 for i in indices:
  mask[i,:]=True;mask[:,i]=True
 V[mask]*=1-fraction
 deg=V.sum(0);iso=deg<=0
 V[iso,iso]=1;deg=V.sum(0)
 P=V/deg[None,:];A=np.eye(len(W))-.5*P
 f=solve(A,.5*v,assume_a='gen',check_finite=False)
 err=float(np.max(np.abs(A@f-.5*v)))
 assert err<1e-10 and np.min(f)>-1e-12 and abs(f.sum()-1)<1e-10
 return f,err

def donor_graphs(X,d):
 covs={};qc=[]
 for donor in sorted(set(d)):
  x=X[d==donor];sd=x.std(0);z=(x-x.mean(0))/np.where(sd>0,sd,1)
  lw=LedoitWolf(assume_centered=True).fit(z);c=lw.covariance_;s=np.sqrt(np.diag(c));c=c/np.outer(s,s);np.fill_diagonal(c,1)
  covs[donor]=c;qc.append({'donor':donor,'n_cells':len(x),'n_nonconstant_genes':int((sd>0).sum()),'shrinkage':lw.shrinkage_})
 return covs,qc

def main():
 t=time.time();rng=np.random.default_rng(CONFIG['seed'])
 meta=pd.read_csv(SC/'cell_metadata.csv').set_index('CellName')
 weights=pd.read_csv(NET/'network_clinical_seed_weights.csv')
 rows=[];genesout=[];nullrows=[];poolrows=[];qcrows=[];inputrows=[];doseout=[];nonadd=[];local=[];donordet=[]
 for comp in COMPS:
  path=SC/f'{comp.lower()}_network_matrix.csv.gz';frame=pd.read_csv(path,index_col=0);genes=frame.index.tolist();X=frame.T.to_numpy(float);d=meta.loc[frame.columns,'Sample'].to_numpy()
  inputrows.append({'file':str(path.relative_to(ROOT)),'sha256':sha(path),'n_cells':len(X),'n_genes':len(genes),'n_donors':len(set(d))})
  aw=weights[weights.compartment.eq(comp)&weights.weight_profile.eq('discovery_DO_vs_HV')].set_index('gene_symbol').weight
  v=np.zeros(len(genes));ai=[genes.index(g) for g in ANCHORS if g in genes]
  v[ai]=[aw[genes[i]] for i in ai];v/=v.sum();clinical=v.copy()
  ti=[genes.index(g) for g in TARGS]
  covs,qc=donor_graphs(X,d);qcrows.extend(dict(compartment=comp,**r) for r in qc)
  for donor in sorted(set(d)):
   for tg,i in zip(TARGS,ti):
    donordet.append({'compartment':comp,'donor':donor,'target':tg,'n_cells':int((d==donor).sum()),'n_detected_cells':int((X[d==donor,i]>0).sum()),'detection_fraction':float((X[d==donor,i]>0).mean()),'mean_expression':float(X[d==donor,i].mean())})
  for omit in ['none']+sorted(covs):
   keepd=[key for key in covs if key!=omit];W=graph(np.mean([covs[k] for k in keepd],axis=0));f0,err0=propagate(W,v)
   den=float(clinical@f0);sample=X[np.isin(d,keepd)]
   vecs={};metrics={}
   for label,inds in [('PDGFRA',[ti[0]]),('FLT4',[ti[1]]),('PDGFRA+FLT4',ti),('KDR',[genes.index('KDR')]),('FGFR1',[genes.index('FGFR1')]),('sham',[])]:
    f,err=propagate(W,v,inds);delta=f-f0;mag=float(clinical@np.abs(delta)/den);metrics[label]=mag;vecs[label]=f
    rows.append({'compartment':comp,'omitted_donor':omit,'intervention':label,'n_donors':len(keepd),'n_cells':len(sample),'clinical_weighted_relative_change':mag,'clinical_weighted_signed_change':float(clinical@delta/den),'global_L1_change_excluding_perturbed_nodes':float(np.abs(np.delete(delta,inds)).sum()),'max_abs_gene_change':float(np.abs(delta).max()),'solve_residual':err,'clinical_baseline_score':den})
    if omit=='none':
     for g,k in zip(genes,range(len(genes))):
      genesout.append({'compartment':comp,'intervention':label,'gene':g,'is_clinical_anchor':g in ANCHORS,'is_perturbed_gene':k in inds,'baseline_score':f0[k],'perturbed_score':f[k],'score_change':delta[k],'relative_score_change':delta[k]/f0[k] if f0[k]>0 else np.nan})
   interaction=vecs['PDGFRA+FLT4']-vecs['PDGFRA']-vecs['FLT4']+f0
   nonmag=float(clinical@np.abs(interaction)/den)
   nonadd.append({'compartment':comp,'omitted_donor':omit,'clinical_nonadditivity_relative_L1':nonmag,'nonadditivity_fraction_of_pair_change':nonmag/metrics['PDGFRA+FLT4'] if metrics['PDGFRA+FLT4'] else np.nan,'signed_clinical_nonadditivity':float(clinical@interaction/den)})
   if omit=='none':
    for a,b in [('PDGFRA','FLT4'),('FLT4','PDGFRA')]:
     bi=genes.index(b);local.append({'compartment':comp,'deleted_target':a,'readout_target':b,'baseline_score':f0[bi],'perturbed_score':vecs[a][bi],'relative_change':(vecs[a][bi]-f0[bi])/f0[bi]})
    # Numerical restoration is an operator test, not a biological rescue.
    fr,err=propagate(W,v);rest=float(np.max(np.abs(fr-f0)))
    rows.append({'compartment':comp,'omitted_donor':omit,'intervention':'restored_graph','n_donors':len(keepd),'n_cells':len(sample),'clinical_weighted_relative_change':float(clinical@np.abs(fr-f0)/den),'clinical_weighted_signed_change':float(clinical@(fr-f0)/den),'global_L1_change_excluding_perturbed_nodes':float(np.abs(fr-f0).sum()),'max_abs_gene_change':rest,'solve_residual':err,'clinical_baseline_score':den})
    for fraction in CONFIG['attenuation_grid']:
     for label,inds in [('PDGFRA',[ti[0]]),('FLT4',[ti[1]]),('PDGFRA+FLT4',ti)]:
      f,_=propagate(W,v,inds,fraction)
      doseout.append({'compartment':comp,'intervention':label,'fraction_incident_edge_attenuation':fraction,'clinical_weighted_relative_change':float(clinical@np.abs(f-f0)/den)})
    mean=np.mean([X[d==donor].mean(0) for donor in sorted(set(d))],0)
    detect=np.mean([(X[d==donor]>0).mean(0) for donor in sorted(set(d))],0)
    feat=np.c_[np.log1p(mean),np.log((detect+.005)/(1-detect+.005)),np.log1p(W.sum(0))]
    feat=(feat-feat.mean(0))/np.maximum(feat.std(0),1e-8)
    bg=np.array([i for i,g in enumerate(genes) if g not in ANCHORS+CANDIDATES]);pools=[]
    for tg,i in zip(TARGS,ti):
     dist=np.linalg.norm(feat[bg]-feat[i],axis=1);inds=bg[np.argsort(dist,kind='stable')[:50]];pools.append(inds)
     for j in inds:poolrows.append({'compartment':comp,'target':tg,'matched_gene':genes[j],'match_distance':float(np.linalg.norm(feat[j]-feat[i])),'target_mean':mean[i],'matched_mean':mean[j],'target_detection':detect[i],'matched_detection':detect[j],'target_degree':W.sum(0)[i],'matched_degree':W.sum(0)[j]})
    singles={}
    for j in sorted(set(pools[0])|set(pools[1])):
     singles[j]=propagate(W,v,[j])[0]
    for tg,pool in zip(TARGS,pools):
     for ni,j in enumerate(pool):
      nullrows.append({'compartment':comp,'intervention':tg,'null_index':ni,'control_1':genes[j],'control_2':'','clinical_weighted_relative_change':float(clinical@np.abs(singles[j]-f0)/den),'clinical_nonadditivity_relative_L1':np.nan})
    pairs=sorted(set((min(i,j),max(i,j)) for i in pools[0] for j in pools[1] if i!=j))
    n=min(len(pairs),CONFIG['pair_null_draws']);sel=rng.choice(len(pairs),size=n,replace=False)
    for ni,ij in enumerate(sel):
     i,j=pairs[ij];f,_=propagate(W,v,[i,j]);inter=f-singles[i]-singles[j]+f0
     nullrows.append({'compartment':comp,'intervention':'PDGFRA+FLT4','null_index':ni,'control_1':genes[i],'control_2':genes[j],'clinical_weighted_relative_change':float(clinical@np.abs(f-f0)/den),'clinical_nonadditivity_relative_L1':float(clinical@np.abs(inter)/den)})
   print(comp,omit,'elapsed',round(time.time()-t,1),flush=True)
  # Per-compartment checkpoint supports recovery from interruption.
  pd.DataFrame(rows).to_csv(OUT/'perturbation_metrics.csv',index=False)
  pd.DataFrame(nullrows).to_csv(OUT/'matched_null_metrics.csv.gz',index=False,compression='gzip')
 f=pd.DataFrame(rows);null=pd.DataFrame(nullrows);inter=pd.DataFrame(nonadd)
 primary=f[f.omitted_donor.eq('none')&f.intervention.isin(TARGS+['PDGFRA+FLT4'])].copy()
 cal=[]
 for _,r in primary.iterrows():
  n=null[null.compartment.eq(r.compartment)&null.intervention.eq(r.intervention)].clinical_weighted_relative_change.to_numpy();obs=r.clinical_weighted_relative_change
  cal.append(dict(r,empirical_p=(1+int((n>=obs).sum()))/(len(n)+1),n_null=len(n),null_median=float(np.median(n)),null_min=float(n.min()),null_max=float(n.max()),null_percentile=float((n<obs).mean())))
 cal=pd.DataFrame(cal);cal['q_bh_9']=false_discovery_control(cal.empirical_p)
 ncal=[]
 for _,r in inter[inter.omitted_donor.eq('none')].iterrows():
  n=null[null.compartment.eq(r.compartment)&null.intervention.eq('PDGFRA+FLT4')].clinical_nonadditivity_relative_L1.to_numpy();obs=r.clinical_nonadditivity_relative_L1
  ncal.append(dict(r,empirical_p=(1+int((n>=obs).sum()))/(len(n)+1),n_null=len(n),null_median=float(np.median(n))))
 ncal=pd.DataFrame(ncal);ncal['q_bh_3']=false_discovery_control(ncal.empirical_p)
 outputs={'perturbation_metrics.csv':f,'primary_calibration.csv':cal,'nonadditivity_by_donor_omission.csv':inter,'nonadditivity_calibration.csv':ncal,'gene_level_perturbations.csv.gz':pd.DataFrame(genesout),'matched_null_metrics.csv.gz':null,'matched_control_pools.csv':pd.DataFrame(poolrows),'donor_target_expression.csv':pd.DataFrame(donordet),'donor_covariance_qc.csv':pd.DataFrame(qcrows),'target_to_target_redistribution.csv':pd.DataFrame(local),'edge_attenuation_sensitivity.csv':pd.DataFrame(doseout),'input_manifest.csv':pd.DataFrame(inputrows)}
 for name,x in outputs.items(): x.to_csv(OUT/name,index=False,compression='gzip' if name.endswith('.gz') else None)
 summary={'config_sha256':sha(OUT/'analysis_specification.json'),'n_primary_tests':len(cal),'n_primary_fdr_discoveries':int(cal.q_bh_9.lt(.05).sum()),'n_nonadditivity_fdr_discoveries':int(ncal.q_bh_3.lt(.05).sum()),'max_solver_residual':float(f.solve_residual.max()),'max_sham_change':float(f[f.intervention.eq('sham')].max_abs_gene_change.max()),'max_restoration_residual':float(f[f.intervention.eq('restored_graph')].max_abs_gene_change.max()),'wall_seconds':time.time()-t,'python':platform.python_version(),'numpy':np.__version__,'interpretation':CONFIG['interpretation']}
 (OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2));print(cal.to_string(index=False));print(ncal.to_string(index=False));print(json.dumps(summary,indent=2))
if __name__=='__main__':
 with threadpool_limits(limits=1):main()
