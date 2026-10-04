from pathlib import Path
import string,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
R=Path(__file__).resolve().parent;N=R/'network_diffusion';B=R/'integrated_analysis/primary'
O=R/'figures_new';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
C={'Hepatocyte':'#168C8C','Endothelial':'#3E6AAB','Macrophage':'#BC7042'}
ns=pd.read_csv(N/'candidate_network_summary.csv');bs=pd.read_csv(B/'candidate_integrated_summary.csv').sort_values('geometric');order=bs.target.tolist();alg=['legacy_consensus','median','mean','geometric','maximin','RRA'];alabel=['Baseline','Median','Mean','Geometric','Maximin','RRA']
maprows=[]
def axesgrid(title):
 f,ax=plt.subplots(3,4,figsize=(13.2,12.4),layout='constrained');f.suptitle(title,fontsize=14,fontweight='bold');
 for a,l in zip(ax.ravel(),string.ascii_uppercase):a.text(-.14,1.06,l,transform=a.transAxes,weight='bold',fontsize=13)
 return f,ax.ravel()
def hm(a,df,cmap='viridis',vmin=None,vmax=None,annot=False,fmt='.1f',bar=True):
 if annot is True and fmt=='.0f':
  annot=df.map(lambda v:'' if pd.isna(v) else str(int(v)) if float(v).is_integer() else f'{v:g}')
  fmt=''
 sns.heatmap(df,ax=a,cmap=cmap,vmin=vmin,vmax=vmax,annot=annot,fmt=fmt,cbar=bar,cbar_kws={'shrink':.65,'pad':.03},mask=df.isna(),linewidths=.25,linecolor='white');a.set_xlabel('');a.set_ylabel('');a.tick_params(length=0);a.set_yticklabels(a.get_yticklabels(),rotation=0)
def export(f,stem):
 for ext in ['pdf','svg']:f.savefig(O/f'{stem}.{ext}',bbox_inches='tight',pad_inches=.1)
 f.savefig(O/f'{stem}.png',dpi=140,bbox_inches='tight',pad_inches=.1)
 f.savefig(O/f'{stem}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'},bbox_inches='tight',pad_inches=.1)
 plt.close(f)
def source(fig,p,file,unit='candidate',analysis='descriptive'):
 maprows.append(dict(figure=fig,panel=p,source_file=file,statistical_unit=unit,analysis=analysis))
f,a=axesgrid('Figure 5 | Phenotype-guided network integration tests candidate priority and robustness')
z=ns.pivot(index='gene_symbol',columns='compartment',values='matched_null_z').reindex(order)[list(C)]
hm(a[0],z,'RdBu_r',-3,3);a[0].set_title('Clinical-protein network proximity');a[0].set_xticklabels(['Hep','Endo','Macro']);a[0].texts.append if False else None
a[0].text(.5,-.13,'Color limited to −3…3; full z values in tables',ha='center',fontsize=7.5,transform=a[0].transAxes)
for j,c in enumerate(C):
 for i,g in enumerate(order):
  q=ns.loc[ns.compartment.eq(c)&ns.gene_symbol.eq(g),'empirical_q_bh'].iloc[0]
  if q<.05:a[0].text(j+.5,i+.5,'*',ha='center',va='center',color='black',fontsize=12)
source('Figure 5','A','Network/candidate_network_summary.csv','candidate–compartment','matched-null z; BH across 50 hypotheses')
x=ns.dropna(subset=['empirical_q_bh']);
for c,col in C.items():
 t=x[x.compartment==c];a[1].scatter(t.donor_mean_detection*100,-np.log10(t.empirical_q_bh),s=24,c=col,label=c,alpha=.8)
hit_label_positions={'PDGFRA':(4,2.36),'LRRK2':(16,2.22),'TYK2':(12,1.96)}
for _,t in x[x.empirical_q_bh<.05].iterrows():
 a[1].annotate(t.gene_symbol,(t.donor_mean_detection*100,-np.log10(t.empirical_q_bh)),xytext=hit_label_positions[t.gene_symbol],textcoords='data',fontsize=8,arrowprops={'arrowstyle':'-','color':'#666','lw':.6})
a[1].set_ylim(-.08,2.52)
a[1].axvline(10,color='#666',ls='--',lw=.8);a[1].axhline(-np.log10(.05),color='#666',ls=':',lw=.8);a[1].set(xlabel='Mean donor detection (%)',ylabel='−log10(BH q)',title='Expression qualifies network hits');a[1].legend(fontsize=7,frameon=False,loc='upper right')
source('Figure 5','B','Network/candidate_network_summary.csv','candidate–compartment','detection versus empirical q')
co=pd.read_csv(N/'network_target_anchor_contributions.csv');cc=co[(co.compartment=='Macrophage')&(co.gene_symbol=='PDGFRA')].copy();cc['share']=cc.weighted_contribution/cc.weighted_contribution.sum();cc=cc.sort_values('share',ascending=False).head(7)
a[2].barh(cc.anchor_gene.iloc[::-1],cc.share.iloc[::-1]*100,color='#168C8C');a[2].set(xlabel='Contribution to diffusion score (%)',title='PDGFRA macrophage signal');source('Figure 5','C','Network/network_target_anchor_contributions.csv','protein seed','additive contribution')
lo=pd.read_csv(N/'network_donor_lodo_scores.csv');picked=[('Endothelial','LRRK2'),('Macrophage','PDGFRA'),('Endothelial','FGFR1'),('Endothelial','FLT4')]
for i,(c,g) in enumerate(picked):
 vals=lo[(lo.compartment==c)&(lo.gene_symbol==g)].matched_null_z.dropna().values;a[3].scatter(vals,np.full(len(vals),i)+np.linspace(-.12,.12,len(vals)),s=20,color=C[c]);a[3].plot([vals.min(),vals.max()],[i,i],color=C[c],lw=1)
a[3].set_yticks(range(len(picked)),[g+' ('+c[:4]+')' for c,g in picked]);a[3].set(xlabel='Matched-null z after donor omission',title='Five biological donor omissions');a[3].axvline(0,ls=':',c='gray',lw=.7)
source('Figure 5','D','Network/network_donor_lodo_scores.csv','omitted donor','fixed primary network, donor omission')
rr=pd.read_csv(B/'algorithm_candidate_rankings.csv').pivot(index='target',columns='algorithm',values='rank').reindex(order)[alg];rr.columns=alabel;hm(a[4],rr,'YlGnBu_r',1,17,True,'.0f',False);a[4].set_xticklabels(alabel,rotation=45,ha='right');a[4].set_title('Six algorithms on a fixed target set');source('Figure 5','E','Integration/algorithm_candidate_rankings.csv','candidate','rank; same 17 targets')
a[5].scatter(bs.legacy_consensus_rank,bs.geometric,color='#168C8C',s=28)
a[5].plot([1,17],[1,17],ls='--',c='#999',lw=.8)
for g in ['FGFR1','PDGFRA','FLT4','FGFR2']:
 t=bs[bs.target==g].iloc[0];a[5].annotate(g,(t.legacy_consensus_rank,t.geometric),xytext=(5,5),textcoords='offset points',fontsize=8)
a[5].set(xlabel='Baseline consensus rank',ylabel='New geometric rank',title='Reassessment changes priority',xlim=(0,19),ylim=(19,0));source('Figure 5','F','Integration/candidate_integrated_summary.csv')
be=pd.read_csv(B/'algorithm_benchmark_summary.csv').set_index('algorithm').reindex(alg)
a[6].barh(alabel[::-1],be.mean_omission_rank_rho.values[::-1],color=['#168C8C' if n=='Geometric' else '#869DAB' for n in alabel[::-1]]);a[6].set(xlim=(0,1),xlabel='Mean rank correlation after omission',title='Component-omission stability');source('Figure 5','G','Integration/algorithm_benchmark_summary.csv','algorithm × omitted family','rank stability, not prediction accuracy')
he=pd.read_csv(B/'algorithm_family_omission_evaluation.csv').pivot(index='algorithm',columns='omitted_family',values='held_feature_spearman').reindex(alg);he.index=alabel;he.columns=['Clinical\nnetwork' if c=='clinical_network' else 'Drug\nresponse' if c=='HepG2_perturbation' else 'Pharm.' if c=='pharmacology' else 'STRING' for c in he.columns];hm(a[7],he,'RdBu_r',-.5,.5,True,'.2f',False);a[7].set_title('Agreement with withheld features');source('Figure 5','H','Integration/algorithm_family_omission_evaluation.csv','candidate rankings','withheld component concordance')
a[8].barh(order[::-1],bs.top5_frequency.values[::-1],color=['#168C8C' if g in ['PDGFRA','FLT4'] else '#869DAB' for g in order[::-1]]);a[8].set(xlim=(0,1.03),xlabel='Top-five frequency in 5,000 draws',title='Sensitivity to integration weights');source('Figure 5','I','Integration/algorithm_weight_sensitivity.csv','Dirichlet weight draw','top-five frequency')
ss=pd.read_csv(B/'algorithm_specification_donor_missingness_sensitivity.csv');d=ss[(ss.algorithm=='geometric')&(ss.sensitivity_type=='donor_omission')]
if d.empty:
 print('SENS_TYPES',ss.sensitivity_type.unique());d=ss[(ss.algorithm=='geometric')&ss.sensitivity_type.str.contains('donor')]
mat=d.pivot(index='target',columns='specification',values='rank').reindex(order);mat.columns=['D'+str(i+1) for i in range(len(mat.columns))];hm(a[9],mat,'YlGnBu_r',1,17,True,'.0f',False);a[9].set_title('Candidate ranks after donor omission');source('Figure 5','J','Integration/algorithm_specification_donor_missingness_sensitivity.csv','omitted donor','geometric fusion with retained-donor detection weights')
ff=pd.read_csv(B/'integrated_family_percentile_scores.csv').set_index('target').reindex(order);ff.columns=['Pharm.','Clinical\nnetwork','Drug\nresponse','STRING'];hm(a[10],ff,'YlGnBu',0,1,False,bar=False);a[10].set_title('Inputs to the integrated score');source('Figure 5','K','Integration/integrated_family_percentile_scores.csv','candidate × source family','rank-scaled score')
qc=ns[['empirical_q_bh','prespecified_stability_gate','expression_supported_stability_gate']]
v=[50,int((qc.empirical_q_bh<.05).sum()),int(qc.prespecified_stability_gate.fillna(False).sum()),int(qc.expression_supported_stability_gate.fillna(False).sum())]
a[11].bar(['Evaluable','q < .05','Stable','Expression\nsupported'],v,color=['#869DAB','#D49A43','#168C8C','#3E6AAB']);a[11].set(ylabel='Candidate–compartment pairs',title='Joint criteria for network evidence',ylim=(0,56));
for i,n in enumerate(v):a[11].text(i,n+.8,str(n),ha='center',fontsize=10)
a[11].tick_params(axis='x',labelsize=7.5);source('Figure 5','L','Network/candidate_network_summary.csv','candidate–compartment','joint statistical, specification, donor and expression criteria')
export(f,'Figure_5_Phenotype_Guided_Network_Integration')
# Extended sensitivity figure.
f,a=axesgrid('Figure S8 | Network and integration sensitivities quantify target-selection uncertainty')
sp=ns.pivot(index='gene_symbol',columns='compartment',values='spec_top_quartile_frequency').reindex(order)[list(C)];hm(a[0],sp,'YlGnBu',0,1,bar=False);a[0].set_title('Top-quartile frequency in 18 networks');a[0].set_xticklabels(['Hep','Endo','Macro']);source('Figure S8','A','Network/candidate_network_summary.csv')
lm=ns.pivot(index='gene_symbol',columns='compartment',values='lodo_top_quartile_count').reindex(order)[list(C)];hm(a[1],lm,'YlGnBu',0,5,True,'.0f',False);a[1].set_title('Top-quartile donor-omission counts');a[1].set_xticklabels(['Hep','Endo','Macro']);source('Figure S8','B','Network/candidate_network_summary.csv')
m=pd.read_csv(N/'network_null_matching_sensitivity.csv');m=m[m.weight_profile=='discovery_DO_vs_HV'];
for c,g in [('Endothelial','LRRK2'),('Endothelial','TYK2'),('Macrophage','PDGFRA')]:
 t=m[(m.compartment==c)&(m.gene_symbol==g)].sort_values('match_pool_size');a[2].plot(t.match_pool_size,t.empirical_q_bh,marker='o',label=g+' '+c[:4])
a[2].axhline(.05,ls='--',color='gray',lw=.8);a[2].set(xlabel='Matched background pool size',ylabel='BH q',title='Expression plus degree matched null');a[2].legend(fontsize=7,frameon=False);source('Figure S8','C','Network/network_null_matching_sensitivity.csv')
ab=pd.read_csv(N/'network_shrinkage_ablation.csv');ab=ab[ab.weight_profile=='discovery_DO_vs_HV']
for c,col in C.items():
 t=ab[ab.compartment==c];a[3].scatter(t.candidate_rank_unshrunk,t.candidate_rank_shrinkage,c=col,label=c,s=20)
a[3].plot([1,17],[1,17],ls='--',c='gray');a[3].set(xlabel='Unshrunk Pearson rank',ylabel='Shrinkage rank',title='Covariance-estimator ablation');a[3].legend(fontsize=7,frameon=False);source('Figure S8','D','Network/network_shrinkage_ablation.csv')
sscore=pd.read_csv(N/'network_specification_scores.csv');pr=sscore[sscore.is_primary_spec]
for c,col in C.items():
 t=pr[pr.compartment==c].pivot(index='gene_symbol',columns='weight_profile',values='matched_null_z');a[4].scatter(t.discovery_DO_vs_HV,t.confirmatory_DO_vs_HV,s=20,color=col)
a[4].set(xlabel='Discovery-weighted z',ylabel='Confirmation-weighted z',title='Clinical-weight sensitivity');source('Figure S8','E','Network/network_specification_scores.csv')
mat=pr[(pr.gene_symbol=='FGFR1')].pivot(index='compartment',columns='weight_profile',values='matched_null_z');mat.columns=['Confirm\ninjury' if x=='confirmatory_DO_vs_HV' else 'Confirm\netiology' if x=='confirmatory_DO_vs_NDO' else 'Discovery\ninjury' if x=='discovery_DO_vs_HV' else 'Discovery\netiology' for x in mat.columns];hm(a[5],mat,'RdBu_r',-3,3,True,'.1f',False);a[5].set_title('FGFR1 under all clinical weights');source('Figure S8','F','Network/network_specification_scores.csv')
for i,g in enumerate(order[:8]):
 t=bs[bs.target==g].iloc[0];a[6].plot([t.network_specification_rank_min,t.network_specification_rank_max],[i,i],c='#168C8C');a[6].scatter(t.geometric,i,s=22,c='#168C8C')
a[6].set_yticks(range(8),order[:8]);a[6].invert_yaxis();a[6].set(xlabel='Rank across network specifications',title='Integrated specification sensitivity');source('Figure S8','G','Integration/candidate_integrated_summary.csv')
for i,g in enumerate(order[:8]):
 t=bs[bs.target==g].iloc[0];a[7].plot([t.rank_q025,t.rank_q975],[i,i],c='#168C8C');a[7].scatter(t.rank_median,i,s=22,c='#168C8C')
a[7].set_yticks(range(8),order[:8]);a[7].invert_yaxis();a[7].set(xlabel='Central 95% weight-perturbation ranks',title='Weight sensitivity ranges');source('Figure S8','H','Integration/algorithm_weight_sensitivity.csv')
mi=ss[(ss.algorithm=='geometric')&ss.sensitivity_type.str.contains('missing')];mat=mi.pivot(index='target',columns='specification',values='rank').reindex(order);hm(a[8],mat,'YlGnBu_r',1,17,True,'.0f',False);a[8].set_xticklabels(a[8].get_xticklabels(),rotation=45,ha='right',fontsize=7);a[8].set_title('Explicit missingness conventions');source('Figure S8','I','Integration/algorithm_specification_donor_missingness_sensitivity.csv')
mt=ss[(ss.algorithm=='geometric')&ss.sensitivity_type.str.contains('measurement')];mat=mt.pivot(index='target',columns='specification',values='rank').reindex(order);hm(a[9],mat,'YlGnBu_r',1,17,True,'.0f',False);a[9].set_xticklabels(a[9].get_xticklabels(),rotation=45,ha='right',fontsize=7);a[9].set_title('L1000 measurement restrictions');source('Figure S8','J','Integration/algorithm_specification_donor_missingness_sensitivity.csv')
bp=pd.read_csv(B/'algorithm_block_permutation_calibration.csv');mat=bp.pivot(index='target',columns='algorithm',values='bh_q').reindex(order);hm(a[10],mat,'YlGnBu_r',0,1,True,'.2f',False);a[10].set_xticklabels(a[10].get_xticklabels(),rotation=45,ha='right');a[10].set_title('Block-permuted aggregation q');source('Figure S8','K','Integration/algorithm_block_permutation_calibration.csv')
q=pd.read_csv(N/'network_donor_covariance_qc.csv');sns.barplot(q,x='donor',y='shrinkage',hue='compartment',palette=C,ax=a[11]);a[11].set(xlabel='Donor',ylabel='Shrinkage coefficient',title='Donor-specific regularization');a[11].set_xticklabels(['D'+str(i+1) for i in range(5)]);a[11].legend(fontsize=7,frameon=False);source('Figure S8','L','Network/network_donor_covariance_qc.csv','donor','Ledoit–Wolf coefficient')
export(f,'Supplementary_Figure_S8_Network_and_Algorithm_Robustness')
pd.DataFrame(maprows).to_csv(O/'new_figure_panel_source_map.csv',index=False)
print('Saved Figure 5, S8 and source map',flush=True)
