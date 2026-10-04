"""Fixed-candidate reanalysis of GSE125975; all archived samples retained.
Empirical-Bayes moderated linear model with scaled inverse-chi-square variance prior.
Independent-animal label permutations provide candidate/module sensitivity.
"""
from pathlib import Path
import gzip,io,json,re,ast,hashlib,itertools
import numpy as np,pandas as pd
from scipy import stats,special,optimize
from scipy.stats import false_discovery_control
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
threadpool_limits(1)
P=Path(__file__).resolve().parent
ROOT=P.parents[1]
CANDIDATES=['FLT4','PDGFRA','FGFR1','KDR','JAK1','LYN','FLT1','FGFR2','FGFR3','ABL1','CDK4','FLT3','FGFR4','TYK2','AXL','LRRK2','MET']
ANCHORS=['ACO1','ASS1','FAH','CPS1','ALDOB','HPD','OTC','DMGDH','GSTA1','FBP1','PCK2','CES1','LECT2']
MODULES={
'ECM/fibrosis':['COL1A1','COL1A2','COL3A1','FN1','POSTN','CTGF','TIMP1','LOX','SERPINE1','ACTA2','TGFBI','MMP2'],
'Macrophage–lipid':['CD68','ADGRE1','CSF1R','LPL','APOE','TREM2','LGALS3','CTSD','ABCA1','FABP5','SPP1','LIPA','PLIN2'],
'Urea/arginine':['ARG1','OTC','CPS1','ASS1','ASL','OAT'],
'Nintedanib targets':['FGFR1','FGFR3','FLT1'],
'Hepatic anchors':['OTC','GSTA1','FBP1','CES1A','CES1D','CES1F','CPS1','ASS1'],
'Oxidative stress':['NQO1','HMOX1','SOD1','SOD2','GPX1','GPX3','GSTA1','CAT','PRDX1','PRDX2'],
'Drug efflux':['ABCB1A','ABCB1B','ABCG2','ABCC1','ABCC2']}
# Resolve one established alias transparently without altering member identity.
ALIASES={'CTGF':'CCN2'}
def bh(x):
 x=np.array(x,float);o=np.full(x.shape,np.nan);ok=np.isfinite(x);o[ok]=false_discovery_control(x[ok]);return o
s=gzip.open(P/'matrix.txt.gz','rt').read();block=s.split('!series_matrix_table_begin\n')[1].split('!series_matrix_table_end')[0];raw=pd.read_csv(io.StringIO(block),sep='\t',index_col=0)
info=pd.read_csv(P/'gene_info.gz',sep='\t',dtype=str,compression='gzip',usecols=['GeneID','Symbol','description'])
info=info.drop_duplicates('GeneID').set_index('GeneID');ent=raw.index.str.replace('_at','',regex=False)
anno=pd.DataFrame({'probe_id':raw.index,'entrez_id':ent});anno['mouse_symbol']=[info.loc[x,'Symbol'] if x in info.index else '' for x in ent];anno['symbol']=anno.mouse_symbol.str.upper();anno['mapped']=anno.symbol!='';anno.loc[~anno.mapped,'symbol']='ENTREZ_'+anno.loc[~anno.mapped,'entrez_id'];anno.to_csv(P/'feature_annotation.csv',index=False)
y=raw.copy();y.index=anno.symbol;y=y.groupby(level=0,sort=True).mean();y.to_csv(P/'expression_RMA_symbols.csv.gz',compression='gzip')
meta=pd.read_csv(P/'sample_metadata.csv').set_index('gsm').loc[y.columns];g=meta.group.to_numpy();levels=['CONTROL','DISEASE','NINT'];X=np.stack([g==z for z in levels],axis=1).astype(float);Y=y.T.to_numpy(float);ix=np.linalg.inv(X.T@X);beta=ix@X.T@Y;resid=Y-X@beta;df=len(g)-len(levels);s2=np.sum(resid**2,axis=0)/df;assert (s2>0).all();z=np.log(s2);excess=np.var(z,ddof=1)-special.polygamma(1,df/2)
if excess>0:d0=2*optimize.brentq(lambda a:special.polygamma(1,a)-excess,1e-5,1e8)
else:d0=1e8
elogF=special.digamma(df/2)-np.log(df/2)-special.digamma(d0/2)+np.log(d0/2);s0=np.exp(z.mean()-elogF);post=(d0*s0+df*s2)/(d0+df);totaldf=df+d0
contrasts={'NINT_vs_DISEASE':('NINT','DISEASE'),'DISEASE_vs_CONTROL':('DISEASE','CONTROL'),'NINT_vs_CONTROL':('NINT','CONTROL')}
tabs={}
for name,(hi,lo) in contrasts.items():
 c=np.array([int(t==hi)-int(t==lo) for t in levels]);ef=c@beta;se=np.sqrt(post*(c@ix@c));t=ef/se;pval=2*stats.t.sf(abs(t),totaldf);crit=stats.t.ppf(.975,totaldf)
 a=Y[g==hi];b=Y[g==lo];ws=np.sqrt(a.var(0,ddof=1)/len(a)+b.var(0,ddof=1)/len(b));wdf=(a.var(0,ddof=1)/len(a)+b.var(0,ddof=1)/len(b))**2/((a.var(0,ddof=1)/len(a))**2/(len(a)-1)+(b.var(0,ddof=1)/len(b))**2/(len(b)-1));wp=2*stats.t.sf(abs(ef/ws),wdf)
 tab=pd.DataFrame({'gene':y.index,'contrast':name,'log2FC':ef,'SE_moderated':se,'CI95_low':ef-crit*se,'CI95_high':ef+crit*se,'t_moderated':t,'p_value':pval,'q_value_BH_genome':bh(pval),'SE_Welch':ws,'p_Welch':wp,'q_Welch_genome':bh(wp),'mean_log2_RMA':Y.mean(0),'n_numerator':len(a),'n_denominator':len(b)}).set_index('gene');
 tabs[name]=tab;(P/f'{name}_all_genes.csv').write_text(tab.to_csv())
# Exact treatment allocations for independent animal-level fixed panel and modules.
D=np.flatnonzero(g=='DISEASE');T=np.flatnonzero(g=='NINT');C=np.flatnonzero(g=='CONTROL');DT=np.r_[D,T];alloc=np.array(list(itertools.combinations(range(len(DT)),len(T))),dtype=np.int16);W=np.full((len(alloc),len(DT)),-1/len(D));W[np.arange(len(alloc))[:,None],alloc]=1/len(T)
tracks=[];L=[]
for role,genes in [('pharmacology_candidate',CANDIDATES),('clinical_protein',ANCHORS),('mouse_CES1_family_sensitivity',['CES1A','CES1B','CES1C','CES1D','CES1E','CES1F','CES1G'])]:
 for gene in genes:
  for name,tab in tabs.items():
   row={'gene':gene,'role':role,'contrast':name,'measurement_available':gene in tab.index}
   if gene in tab.index:
    row.update(tab.loc[gene].to_dict());vals=y.loc[gene].to_numpy();row.update({'mean_CONTROL':float(vals[C].mean()),'mean_DISEASE':float(vals[D].mean()),'mean_NINT':float(vals[T].mean())})
    if name=='NINT_vs_DISEASE':
     null=W@vals[DT];row['p_exact_animal_permutation']=np.mean(abs(null)>=abs(row['log2FC'])-1e-12)
     leave=[np.mean(np.delete(vals[T],i))-vals[D].mean() for i in range(len(T))]+[vals[T].mean()-np.mean(np.delete(vals[D],i)) for i in range(len(D))]
     row['leave_one_animal_min']=min(leave);row['leave_one_animal_max']=max(leave)
   tracks.append(row)
panel=pd.DataFrame(tracks)
for role in panel.role.unique():
 q=(panel.role==role)&(panel.contrast=='NINT_vs_DISEASE');panel.loc[q,'q_exact_panel']=bh(panel.loc[q,'p_exact_animal_permutation']);panel.loc[q,'q_exact_full_fixed_family']=false_discovery_control(panel.loc[q,'p_exact_animal_permutation'].fillna(1).to_numpy());panel.loc[q & ~panel.measurement_available,'q_exact_full_fixed_family']=np.nan
panel.to_csv(P/'fixed_targets_and_clinical_proteins.csv',index=False)
# Module score, using fixed genes already specified in prior public-preclinical workflow.
members=[];scores={}
for m,genes in MODULES.items():
 present=[]
 for gene in genes:
  lookup=ALIASES.get(gene,gene);found=lookup in y.index;members.append({'module':m,'requested_gene':gene,'resolved_gene':lookup,'available':found})
  if found:present.append(lookup)
 values=y.loc[present].to_numpy();values=(values-values.mean(1,keepdims=True))/values.std(1,ddof=1,keepdims=True);scores[m]=values.mean(0)
pd.DataFrame(members).to_csv(P/'module_gene_membership.csv',index=False);score=pd.DataFrame(scores,index=y.columns);score.to_csv(P/'module_scores_by_animal.csv');mr=[]
for name,(hi,lo) in contrasts.items():
 for m in MODULES:
  v=score[m].to_numpy();a=v[g==hi];b=v[g==lo];effect=a.mean()-b.mean();ss=a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b);se=np.sqrt(ss);dd=ss**2/((a.var(ddof=1)/len(a))**2/(len(a)-1)+(b.var(ddof=1)/len(b))**2/(len(b)-1));crit=stats.t.ppf(.975,dd);r={'module':m,'contrast':name,'effect_standardized_module':effect,'SE_Welch':se,'CI95_low':effect-crit*se,'CI95_high':effect+crit*se,'p_Welch':float(2*stats.t.sf(abs(effect/se),dd)),'n_genes':sum(t['available'] for t in members if t['module']==m)}
  if name=='NINT_vs_DISEASE':r['p_exact']=float(np.mean(abs(W@v[DT])>=abs(effect)-1e-12))
  mr.append(r)
modules=pd.DataFrame(mr)
for name in contrasts:
 ixm=modules.contrast==name;modules.loc[ixm,'q_Welch']=bh(modules.loc[ixm,'p_Welch']);modules.loc[ixm,'q_exact']=bh(modules.loc[ixm,'p_exact'])
modules.to_csv(P/'module_contrasts.csv',index=False)
# Transcriptome geometry and shared-reference sensitivity.
d=tabs['DISEASE_vs_CONTROL'];t=tabs['NINT_vs_DISEASE'];sel=(d.q_value_BH_genome<.05)&(abs(d.log2FC)>=.5);x=d.log2FC[sel].to_numpy();v=t.loc[sel,'log2FC'].to_numpy();ge={'disease_q05_fc05_n':int(sel.sum()),'opposite_direction_fraction':float(np.mean(x*v<0)),'closer_to_control_fraction':float(np.mean(abs(x+v)<abs(x))),'spearman_effects':float(stats.spearmanr(x,v).statistic)}
# Using all genes with absolute disease effect >=0.5 in separate disease-reference sets avoids duplicated denominator noise.
rows=[];ym=y.to_numpy();cm=ym[:,C].mean(1)
for a in itertools.combinations(D,4):
 A=np.array(a);B=np.array([i for i in D if i not in a]);xx=ym[:,A].mean(1)-cm;vv=ym[:,T].mean(1)-ym[:,B].mean(1);keep=abs(xx)>=.5
 rows.append({'disease_reference_animals':','.join(meta.iloc[A].animal_id),'treatment_reference_animals':','.join(meta.iloc[B].animal_id),'n_selected':int(keep.sum()),'opposite_direction_fraction':float(np.mean(xx[keep]*vv[keep]<0)),'spearman':float(stats.spearmanr(xx[keep],vv[keep]).statistic)})
pd.DataFrame(rows).to_csv(P/'disjoint_disease_reference_sensitivity.csv',index=False)
# PCA and sample QC, all samples retained.
var=y.var(axis=1);top=var.nlargest(min(2000,len(var))).index;pca=PCA(n_components=4,svd_solver="full");pc=pca.fit_transform(y.loc[top].T);pcdf=pd.DataFrame(pc,index=y.columns,columns=['PC1','PC2','PC3','PC4']);pcdf=pcdf.join(meta[['animal_id','group']]);pcdf.to_csv(P/'sample_PCA.csv');qc=pd.DataFrame({'sample':y.columns,'group':g,'median_log2_RMA':np.median(Y,axis=1),'IQR_log2_RMA':np.quantile(Y,.75,axis=1)-np.quantile(Y,.25,axis=1)});qc.to_csv(P/'sample_QC.csv',index=False)
summary={'analysis_version':'deposited_RMA_matrix_Python_EB','dataset':'GSE125975','group_counts':meta.group.value_counts().to_dict(),'unique_animals':meta.animal_id.nunique(),'dose_mg_kg':50,'treatment_weeks':None,'treatment_duration_note':'Specific treatment duration not independently resolved from available GEO metadata and full-text narrative; not used in analysis','age_weeks':17,'sex':'male','strain':'C57BL/6J','tissue':'non-tumor liver parenchyma','model':'DEN followed by weekly CCl4; established fibrosis before randomization','raw_matrix_features':len(raw),'analyzed_gene_features':len(y),'mapped_features':int(anno.mapped.sum()),'unmapped_features':int((~anno.mapped).sum()),'collapsed_duplicate_symbols':len(raw)-len(y),'residual_df':int(df),'prior_df':float(d0),'prior_variance':float(s0),'total_df':float(totaldf),'mean_expression_log_variance_spearman':float(stats.spearmanr(Y.mean(0),z).statistic),'exact_treatment_label_allocations':len(W),'pca_variance':pca.explained_variance_ratio_.tolist(),'DE_counts':{k:{'BH_q05':int((tb.q_value_BH_genome<.05).sum()),'BH_q05_FC05':int(((tb.q_value_BH_genome<.05)&(abs(tb.log2FC)>=.5)).sum()),'Welch_q05':int((tb.q_Welch_genome<.05).sum())} for k,tb in tabs.items()},'treatment_modules_exact_q05':modules.query('contrast=="NINT_vs_DISEASE" and q_exact<0.05').module.tolist(),'geometry':ge,'disjoint_reference':{'splits':len(rows),'median_reversal':float(pd.DataFrame(rows).opposite_direction_fraction.median()),'min_reversal':float(pd.DataFrame(rows).opposite_direction_fraction.min()),'max_reversal':float(pd.DataFrame(rows).opposite_direction_fraction.max()),'median_spearman':float(pd.DataFrame(rows).spearman.median())},'source_primary_doi':'10.1053/j.gastro.2019.07.028','interpretation':'Independent hepatic pharmacodynamic-response evaluation. Original source reported no serum ALT/AST evidence of hepatotoxicity. Not a DILI-positive model or target-causal validation.'}
(P/'analysis_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(panel.query('contrast=="NINT_vs_DISEASE"')[['gene','measurement_available','log2FC','q_value_BH_genome','q_exact_panel']].to_string(index=False));print(modules.query('contrast=="NINT_vs_DISEASE"').to_string(index=False))
