#!/usr/bin/env python3
"""Annotate raw archive comparisons and compare only numerically evaluable vKO runs."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr

def main():
 p=argparse.ArgumentParser();p.add_argument('--package',required=True);p.add_argument('--audit',required=True);p.add_argument('--output',required=True);a=p.parse_args();out=Path(a.output);src=Path(a.package)/'Source_Data'
 qc=pd.read_csv(a.audit)
 old=pd.read_csv(src/'Virtual_Knockout/virtual_knockout_seed_summary.csv').merge(qc[['compartment','seed','ko_label','numerically_nonzero']],on=['compartment','seed','ko_label'],how='left',validate='one_to_one')
 old=old[old.ko_class.eq('drug_target')]
 ns=old.groupby(['compartment','ko_label']).agg(n_archived_runs=('seed','size'),n_evaluable_runs=('numerically_nonzero','sum')).reset_index().rename(columns={'ko_label':'gene_symbol'})
 valid=old[old.numerically_nonzero]
 ss=valid.groupby(['compartment','ko_label']).agg(evaluable_vko_mean_anchor_percentile=('dili_anchor_mean_percentile','median')).reset_index().rename(columns={'ko_label':'gene_symbol'})
 data=pd.read_csv(out/'candidate_network_summary.csv').merge(ns,on=['compartment','gene_symbol'],how='left').merge(ss,on=['compartment','gene_symbol'],how='left')
 data['old_three_seed_evaluable']=data.n_evaluable_runs.eq(3)
 data['benchmark_interpretation']='RWR proximity vs valid vKO unsigned displacement; same reference atlas; no independent validation'
 data.to_csv(out/'network_evaluable_legacy_comparison.csv',index=False)
 rows=[]
 for comp in data.compartment.unique():
  for min_valid in [1,2,3]:
   q=data[data.compartment.eq(comp)&data.n_evaluable_runs.ge(min_valid)][['matched_null_z','evaluable_vko_mean_anchor_percentile']].dropna()
   r,pv=spearmanr(q.iloc[:,0],q.iloc[:,1]) if len(q)>=3 else (np.nan,np.nan)
   rows.append(dict(compartment=comp,minimum_numerically_evaluable_vko_runs=min_valid,n_common_candidates=len(q),spearman_rho=r,spearman_p_descriptive=pv,interpretation='conditional exploratory comparison; different network estimands; not predictive validation'))
 pd.DataFrame(rows).to_csv(out/'network_evaluable_legacy_metrics.csv',index=False)
 raw=pd.read_csv(out/'network_old_algorithm_comparison.csv');raw['legacy_rank_status']='UNFILTERED_ARCHIVE_INCLUDES_NUMERICAL_ZERO_INVALID_FOR_PRIMARY_COMPARISON';raw.to_csv(out/'network_old_algorithm_comparison.csv',index=False)
 metrics=pd.read_csv(out/'network_algorithm_comparison_metrics.csv');metrics['comparison_status']=np.where(metrics.comparison.str.startswith('old_'),'UNFILTERED_ARCHIVE_INCLUDES_NUMERICAL_ZERO_INVALID_FOR_PRIMARY_COMPARISON','VALID_DONOR_DELETION');metrics.to_csv(out/'network_algorithm_comparison_metrics.csv',index=False)
 print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
