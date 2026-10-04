#!/usr/bin/env python3
"""Diagnose finite-precision zero displacements before interpreting percentile ranks."""
import argparse,json
from pathlib import Path
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--input-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent);a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
x=pd.read_csv(a.input_root/'Source_Data/Virtual_Knockout/virtual_knockout_gene_results.csv.gz')
z=x.groupby(['ko_label','compartment','seed','ko_class']).distance.agg(['min','median','max']).reset_index();z['numerically_nonzero']=z['max']>1e-12
z.to_csv(a.output_dir/'vko_numerical_run_audit.csv',index=False)
thresholds=[]
for t in [1e-14,1e-12,1e-10]:
 for c,q in z.groupby('ko_class'):
  thresholds.append({'numeric_zero_threshold':t,'ko_class':c,'total_runs':len(q),'numerical_zero_runs':int((q['max']<=t).sum()),'numerically_nonzero_runs':int((q['max']>t).sum()),'same_classification_as_primary':bool(((q['max']>t)==q.numerically_nonzero).all())})
pd.DataFrame(thresholds).to_csv(a.output_dir/'vko_numeric_threshold_sensitivity.csv',index=False)
q=z.loc[z.ko_class.eq('drug_target')];s=q.groupby(['ko_label','compartment']).agg(n_runs=('seed','size'),n_numerically_nonzero=('numerically_nonzero','sum'),minimum_max_distance=('max','min'),maximum_max_distance=('max','max')).reset_index();s.to_csv(a.output_dir/'vko_target_compartment_numeric_validity.csv',index=False)
res={'total_gene_rows':len(x),'total_runs':len(z),'numeric_zero_runs':int((~z.numerically_nonzero).sum()),'nonzero_runs':int(z.numerically_nonzero.sum()),'primary_threshold':1e-12,'maximum_displacement_in_zero_cluster':float(z.loc[~z.numerically_nonzero,'max'].max()),'minimum_max_displacement_in_nonzero_cluster':float(z.loc[z.numerically_nonzero,'max'].min()),'drug_target_total_runs':len(q),'drug_target_zero_runs':int((~q.numerically_nonzero).sum()),'drug_targets_total':q.ko_label.nunique(),'drug_targets_any_valid_run':q.loc[q.numerically_nonzero].ko_label.nunique(),'drug_targets_with_valid_three_seed_compartment':s.loc[s.n_numerically_nonzero.eq(3)].ko_label.unique().tolist(),'interpretation':'A percentile rank from a run with only floating-point-scale displacement cannot establish a biological network perturbation.'}
(a.output_dir/'vko_numerical_audit_summary.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
