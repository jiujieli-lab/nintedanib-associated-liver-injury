#!/usr/bin/env python3
"""Independent arithmetic/provenance checks on frozen integration outputs."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
p=argparse.ArgumentParser();p.add_argument('--input-root',type=Path,required=True);p.add_argument('--network-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'primary');a=p.parse_args();O=a.output_dir;checks=[]
def chk(name,result,detail=''):checks.append({'check':name,'pass':bool(result),'detail':detail})
r=pd.read_csv(O/'candidate_integrated_summary.csv');raw=pd.read_csv(O/'integrated_source_features.csv').set_index('target');scores=pd.read_csv(O/'integrated_family_percentile_scores.csv').set_index('target');net=pd.read_csv(a.network_dir/'candidate_network_summary.csv');G=raw.index.tolist()
chk('17 distinct candidates retained',len(G)==17 and len(set(G))==17)
# Direct by-candidate arithmetic without calling the production functions.
for g in G:
 z=net[net.gene_symbol.eq(g)].dropna(subset=['matched_null_z','donor_mean_detection']);v=(z.matched_null_z*z.donor_mean_detection).sum()/z.donor_mean_detection.sum()
 chk('network pooling '+g,np.isclose(v,raw.loc[g,'clinical_network'],rtol=1e-12,atol=1e-12))
for c in raw.columns:
 z=raw[c];observed=z.dropna();manual=(observed.rank()-.5)/len(observed)
 chk('percentile feature '+c,np.allclose(manual.to_numpy(),scores.loc[manual.index,c].to_numpy(),rtol=1e-12))
manual=np.power(scores.fillna(.5/17).prod(axis=1),.25);ranks=manual.rank(ascending=False)
reported=r.set_index('target').loc[G,'geometric'];chk('independent geometric ranks',np.array_equal(ranks,reported))
chk('primary top five',r.sort_values('geometric').head(5).target.tolist()==['PDGFRA','FLT4','FGFR2','JAK1','FGFR3'])
old=pd.read_csv(a.input_root/'Source_Data/Integration/candidate_consensus_ranking.csv');old=old[old.candidate_class.eq('exposure_proximal_target')].set_index('gene_symbol').loc[G];domains=old[['pharmacology','liver_context','HepG2_perturbation','STRING_context']].rank(ascending=False,na_option='bottom');temp=pd.DataFrame({'m':domains.median(axis=1),'a':domains.mean(axis=1)});tu=list(zip(temp.m,temp.a));manual_legacy=pd.Series(tu,index=G).rank(method='dense');chk('exact original median consensus',np.array_equal(manual_legacy,old.consensus_rank))
held=pd.read_csv(O/'algorithm_family_omission_evaluation.csv');s=pd.read_csv(O/'algorithm_benchmark_summary.csv').set_index('algorithm')
for method,z in held.groupby('algorithm'):
 chk('omission means '+method,np.isclose(z.rank_spearman_vs_full.mean(),s.loc[method,'mean_omission_rank_rho']) and np.isclose(z.top5_jaccard_vs_full.mean(),s.loc[method,'mean_omission_top5_jaccard']))
cal=pd.read_csv(O/'algorithm_block_permutation_calibration.csv');chk('all 85 calibrated tests present',len(cal)==85 and cal[['target','algorithm']].drop_duplicates().shape[0]==85);chk('no FDR-confirmed fusion hits',bool((cal.bh_q>.05).all()))
ws=pd.read_csv(O/'algorithm_weight_sensitivity.csv');chk('all 34 weight sensitivity outputs',len(ws)==34 and ws.n_weight_draws.eq(5000).all() and ws.top5_frequency.between(0,1).all())
sens=pd.read_csv(O/'algorithm_specification_donor_missingness_sensitivity.csv');geo=sens[sens.algorithm.eq('geometric')];chk('72 network specifications',geo[geo.sensitivity_type.eq('network_specification')].specification.nunique()==72);chk('five donor omissions',geo[geo.sensitivity_type.eq('donor_omission')].specification.nunique()==5)
node=pd.read_csv(O/'unified_candidate_gene_nodes.csv');chk('30 connected biological nodes',len(node)==30 and node.gene_symbol.nunique()==30)
e=pd.read_csv(O/'target_human_protein_bridge.csv');chk('17 by 13 target protein links',len(e)==221 and np.allclose(e.groupby('gene_symbol').within_target_contribution_fraction.sum(),1))
for _,z in pd.read_csv(O/'input_manifest.csv').iterrows():
 candidates=[a.network_dir/z.filename]+list((a.input_root/'Source_Data').glob('**/'+z.filename));paths=[q for q in candidates if q.exists()];chk('input provenance '+z.filename,any(hashlib.sha256(q.read_bytes()).hexdigest()==z.sha256 for q in paths))
v=pd.read_csv(a.input_root/'Source_Data/Virtual_Knockout/virtual_knockout_gene_results.csv.gz');mx=v.groupby(['ko_label','compartment','seed','ko_class']).distance.max().reset_index();zero=mx.distance<=1e-12;chk('210 run and 132 numeric zero audit',len(mx)==210 and zero.sum()==132);d=mx[mx.ko_class.eq('drug_target')];chk('102 drug runs and 73 numeric zero',len(d)==102 and (d.distance<=1e-12).sum()==73);chk('numerical threshold separation',np.array_equal(mx.distance<=1e-14,zero) and np.array_equal(mx.distance<=1e-10,zero))
res={'decision':'PASS' if all(z['pass'] for z in checks) else 'FAIL','n_checks':len(checks),'n_passed':sum(z['pass'] for z in checks),'checks':checks,'interpretation':'Arithmetic, source identity, and output consistency checks; not certification of biological causality.'};(O/'integration_validation.json').write_text(json.dumps(res,indent=2));print(res['decision'],res['n_passed'],'/',res['n_checks']);print([z for z in checks if not z['pass']])
