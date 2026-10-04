#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd,numpy as np
p=argparse.ArgumentParser();p.add_argument('--input-root',type=Path,required=True);p.add_argument('--network-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'primary');a=p.parse_args();O=a.output_dir
old=pd.read_csv(a.input_root/'Source_Data/Integration/candidate_consensus_ranking.csv');t=old[old.candidate_class.eq('exposure_proximal_target')];pro=pd.read_csv(a.input_root/'Source_Data/Human_DILI_Proteomics/candidate_evidence_table.csv')
node=pd.DataFrame([{'gene_symbol':g,'biological_quantity':'candidate drug-binding protein','evaluated_by':'pharmacology, liver-network propagation, cell perturbation','causal_status':'hypothesis to test'} for g in t.gene_symbol]+[{'gene_symbol':g,'biological_quantity':'circulating liver-injury protein readout','evaluated_by':'human DILI clinical effect and liver-network position','causal_status':'phenotype readout'} for g in pro.protein]);node.to_csv(O/'unified_candidate_gene_nodes.csv',index=False)
e=pd.read_csv(a.network_dir/'network_target_anchor_contributions.csv');sc=pd.read_csv(a.network_dir/'candidate_network_summary.csv');e=e.merge(sc[['gene_symbol','compartment','donor_mean_detection']],on=['gene_symbol','compartment'],validate='many_to_one');e['compartment_weight']=e.donor_mean_detection/e.groupby(['gene_symbol','anchor_gene']).donor_mean_detection.transform('sum');e['pooled_contribution']=e.weighted_contribution*e.compartment_weight
q=e.groupby(['gene_symbol','anchor_gene']).agg(detection_weighted_network_contribution=('pooled_contribution','sum'),n_available_compartments=('compartment','nunique')).reset_index();q['within_target_contribution_fraction']=q.detection_weighted_network_contribution/q.groupby('gene_symbol').detection_weighted_network_contribution.transform('sum');q['edge_interpretation']='unsigned network reachability; not regulatory sign or mediation'
cols={'protein':'anchor_gene','confirmatory_DO_vs_HV_cliffs_delta_group1_minus_group2':'human_DO_vs_HV_delta','confirmatory_DO_vs_NDO_cliffs_delta_group1_minus_group2':'human_DO_vs_NDO_delta','confirmatory_DO_vs_NDO_mann_whitney_q_bh':'human_DO_vs_NDO_q'}
q=q.merge(pro[list(cols)].rename(columns=cols),on='anchor_gene');q.to_csv(O/'target_human_protein_bridge.csv',index=False)
q.sort_values(['gene_symbol','within_target_contribution_fraction'],ascending=[True,False]).groupby('gene_symbol').head(3).to_csv(O/'top_network_contributors_per_target.csv',index=False)

# Extend provenance with clinical and target-anchor matrices used by this bridge map.
import hashlib
mp=O/'input_manifest.csv';manifest=pd.read_csv(mp)
for fp in [a.network_dir/'network_target_anchor_contributions.csv',a.input_root/'Source_Data/Human_DILI_Proteomics/candidate_evidence_table.csv']:
 record={'filename':fp.name,'sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'bytes':fp.stat().st_size}
 if not manifest.filename.eq(fp.name).any():manifest=pd.concat([manifest,pd.DataFrame([record])],ignore_index=True)
manifest.to_csv(mp,index=False)
