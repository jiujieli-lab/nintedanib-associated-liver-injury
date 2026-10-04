from pathlib import Path
import numpy as np,pandas as pd,itertools,json
from scipy import stats,special,optimize
P=Path(__file__).resolve().parent
meta=pd.read_csv(P/'CEL_file_manifest.csv', dtype={'inferred_CEL_run_date': str}).set_index('gsm')
assert set(meta.inferred_CEL_run_date)=={'180216','190216'}
meta.to_csv(P/'CEL_file_manifest.csv');y=pd.read_csv(P/'expression_RMA_symbols.csv.gz',index_col=0);meta=meta.loc[y.columns];g=meta.group.to_numpy();batch=(meta.inferred_CEL_run_date=='190216').to_numpy();X=np.stack([g=='CONTROL',g=='DISEASE',g=='NINT',batch],axis=1).astype(float);ix=np.linalg.inv(X.T@X);Y=y.T.to_numpy();b=ix@X.T@Y;df=len(g)-X.shape[1];s2=((Y-X@b)**2).sum(0)/df;z=np.log(s2);v=z.var(ddof=1)-special.polygamma(1,df/2);d0=2*optimize.brentq(lambda a:special.polygamma(1,a)-v,1e-5,1e8);s0=np.exp(z.mean()-(special.digamma(df/2)-np.log(df/2)-special.digamma(d0/2)+np.log(d0/2)));post=(d0*s0+df*s2)/(d0+df);c=np.array([0,-1,1,0]);ef=c@b;se=np.sqrt(post*(c@ix@c));pv=2*stats.t.sf(abs(ef/se),df+d0);t=pd.DataFrame({'gene':y.index,'log2FC_batch_adjusted':ef,'SE_moderated_batch':se,'p_batch_moderated':pv,'q_batch_genome':stats.false_discovery_control(pv)}).set_index('gene');(P/'batch_adjusted_gene_sensitivity.csv').write_text(t.to_csv())
# Conditional randomization preserves observed treatment counts in both inferred CEL dates.
DT=np.flatnonzero(g!='CONTROL');nD=sum(g=='DISEASE');nT=sum(g=='NINT');groups=[np.flatnonzero(batch[DT]==z) for z in [0,1]];counts=[sum(g[DT][i]=='NINT') for i in groups];alloc=[]
for a in itertools.combinations(groups[0],counts[0]):
 for b in itertools.combinations(groups[1],counts[1]):alloc.append(a+b)
a=np.array(alloc,dtype=np.int16);W=np.full((len(a),len(DT)),-1/nD);W[np.arange(len(a))[:,None],a]=1/nT
panel=pd.read_csv(P/'fixed_targets_and_clinical_proteins.csv');panel=panel[panel.contrast=='NINT_vs_DISEASE'].copy()
for i,r in panel.iterrows():
 if r.measurement_available:
  panel.loc[i,'p_exact_CEL_date_stratified']=np.mean(abs(W@y.loc[r.gene].to_numpy()[DT])>=abs(r.log2FC)-1e-12)
  for col in t:panel.loc[i,col]=t.loc[r.gene,col]
for role in panel.role.unique():
 ix=panel.role==role;panel.loc[ix,'q_exact_CEL_date_stratified_full_family']=stats.false_discovery_control(panel.loc[ix,'p_exact_CEL_date_stratified'].fillna(1));panel.loc[ix&~panel.measurement_available,'q_exact_CEL_date_stratified_full_family']=np.nan
panel.to_csv(P/'CEL_run_date_fixed_panel_sensitivity.csv',index=False)
summary={'inferred_CEL_run_date_group_crosstab':pd.crosstab(meta.inferred_CEL_run_date,meta.group).to_dict(),'conditional_allocations':len(W),'interpretation':'Date inferred solely from CEL filenames; sensitivity variable, not an independently documented batch covariate.'};(P/'batch_sensitivity_summary.json').write_text(json.dumps(summary,indent=2));print(panel.loc[panel.gene.isin(['KDR','FGFR1','FLT1']),['gene','log2FC_batch_adjusted','q_batch_genome','p_exact_CEL_date_stratified','q_exact_CEL_date_stratified_full_family']].to_string(index=False));print(summary)
