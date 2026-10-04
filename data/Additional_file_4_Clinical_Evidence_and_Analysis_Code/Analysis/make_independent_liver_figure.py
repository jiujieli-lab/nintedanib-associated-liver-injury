from pathlib import Path
import os,string,json,hashlib,subprocess
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from matplotlib.ticker import PercentFormatter
from PIL import Image
R=Path(__file__).resolve().parent;ROOT=R/'independent_liver'
S=Path(os.environ.get('INDEPENDENT_LIVER_FIGURE_SOURCE',str(ROOT)))
O=R/'figures_new';O.mkdir(exist_ok=True)
summary=json.loads((S/'analysis_summary.json').read_text())
de=pd.read_csv(S/'NINT_vs_DISEASE_all_genes.csv').set_index('gene')
disease=pd.read_csv(S/'DISEASE_vs_CONTROL_all_genes.csv').set_index('gene')
fixed=pd.read_csv(S/'fixed_targets_and_clinical_proteins.csv')
cand=fixed[fixed.role.eq('pharmacology_candidate')&fixed.contrast.eq('NINT_vs_DISEASE')].set_index('gene')
rank=pd.read_csv(R/'integrated_analysis/primary/candidate_integrated_summary.csv').set_index('target').sort_values('geometric')
cand=cand.reindex(rank.index)
exprpath=S/'expression_RMA_symbols.csv.gz'
if not exprpath.exists():exprpath=ROOT/'expression_RMA_symbols.csv.gz'
expr=pd.read_csv(exprpath,index_col=0)
meta=pd.read_csv(ROOT/'sample_metadata.csv').set_index('gsm').loc[expr.columns]
pca=pd.read_csv(S/'sample_PCA.csv',index_col=0)
modules=pd.read_csv(S/'module_contrasts.csv').query('contrast == "NINT_vs_DISEASE"').set_index('module')
members=pd.read_csv(S/'module_gene_membership.csv')
splits=pd.read_csv(S/'disjoint_disease_reference_sensitivity.csv')
for g in ['KDR','FLT1','FGFR1']:
 for group in ['CONTROL','DISEASE','NINT']:
  assert np.isclose(expr.loc[g,meta.index[meta.group.eq(group)]].mean(),cand.loc[g,'mean_'+group],atol=1e-9)
assert len(meta)==22 and meta.group.value_counts().to_dict()=={'NINT':10,'DISEASE':9,'CONTROL':3}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.titlesize':10,
 'axes.labelsize':8.5,'xtick.labelsize':7.4,'ytick.labelsize':7.4,'pdf.fonttype':42,
 'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig,aa=plt.subplots(3,4,figsize=(13.5,12.3),layout='constrained');a=aa.ravel()
fig.get_layout_engine().set(w_pad=.08,h_pad=.13,wspace=.07,hspace=.13)
fig.suptitle('Figure 6 | An independent liver experiment tests the ranked nintedanib target hypotheses',fontweight='bold',fontsize=14)
for ax,letter in zip(a,string.ascii_uppercase):ax.text(-.16,1.07,letter,transform=ax.transAxes,fontweight='bold',fontsize=13)
groupcolors={'CONTROL':'#90959c','DISEASE':'#cd8b41','NINT':'#168c8c'}
grouplabels={'CONTROL':'Healthy (n=3)','DISEASE':'Disease (n=9)','NINT':'Nintedanib (n=10)'}
maprows=[]
def source(panel,file,unit,analysis):maprows.append({'figure':'Figure 6','panel':panel,'source_file':'Independent_liver/'+file,'statistical_unit':unit,'analysis':analysis})

# A
for group in groupcolors:
 t=pca[pca.group.eq(group)];a[0].scatter(t.PC1,t.PC2,s=33,color=groupcolors[group],label=grouplabels[group],edgecolors='white',linewidths=.5)
a[0].set(xlabel=f'PC1 ({summary["pca_variance"][0]*100:.1f}%)',ylabel=f'PC2 ({summary["pca_variance"][1]*100:.1f}%)',title='Independent liver transcriptomes')
a[0].legend(frameon=False,fontsize=7,loc='best');source('A','sample_PCA.csv','animal; all 22','PCA of 2,000 variable genes')
# B
sig=de.q_value_BH_genome.lt(.05)
a[1].scatter(de.log2FC,-np.log10(de.q_value_BH_genome.clip(lower=1e-300)),s=5,color='#d1d5d8',alpha=.6,rasterized=True)
a[1].scatter(de.loc[sig,'log2FC'],-np.log10(de.loc[sig,'q_value_BH_genome']),s=9,c=np.where(de.loc[sig,'log2FC']<0,'#cd8b41','#168c8c'),alpha=.8,rasterized=True)
a[1].axhline(-np.log10(.05),ls='--',lw=.7,color='#777');a[1].axvline(0,lw=.6,color='#bbb')
a[1].set(xlabel='Nintedanib − disease (log2 RMA)',ylabel='−log10(genome-wide BH q)',title='Hepatic drug-response transcriptome')
a[1].text(.02,.98,f'{sig.sum():,} features at q < .05\nNintedanib n=10; disease n=9',va='top',transform=a[1].transAxes,fontsize=7.3)
source('B','NINT_vs_DISEASE_all_genes.csv','19 treated/disease animals; variance model all 22','moderated-test BH across 13,392 features')
# C
for j,(g,t) in enumerate(cand.iterrows()):
 if t.measurement_available:
  a[2].plot([t.CI95_low,t.CI95_high],[j,j],c='#9cabb5',lw=1)
  a[2].scatter(t.log2FC,j,s=22,color='#cd8b41' if t.q_exact_full_fixed_family<.05 else '#168c8c',zorder=3)
 else:a[2].text(.55,j,'Unavailable',color='#999',fontsize=6.4,ha='right',va='center')
a[2].axvline(0,lw=.7,color='#aaa');a[2].set_yticks(range(len(cand)),[g+' *' if bool(cand.loc[g,'measurement_available']) and cand.loc[g,'q_exact_full_fixed_family']<.05 else g for g in cand.index]);a[2].set_ylim(len(cand)-.2,-.8);a[2].set(xlim=(-.75,.60),xlabel='Treatment effect (moderated 95% CI)',title='Seventeen frozen drug-target candidates');a[2].tick_params(axis='y',labelsize=7)
source('C','fixed_targets_and_clinical_proteins.csv','19 treated/disease animals; variance model all 22','moderated CI; star exact BH across full fixed 17-candidate family')
# D-F: individual animals, frozen genes, no inferred pseudoreplicates.
animalrows=[]
for ax,g,letter in zip(a[3:6],['KDR','FLT1','FGFR1'],['D','E','F']):
 for j,group in enumerate(groupcolors):
  ids=meta.index[meta.group.eq(group)];vals=expr.loc[g,ids].values
  ax.scatter(j+np.linspace(-.18,.18,len(vals)),vals,s=26,color=groupcolors[group],edgecolors='white',linewidths=.45,zorder=3)
  ax.plot([j-.23,j+.23],[vals.mean(),vals.mean()],c='#333',lw=1.2)
  for gsm,val in zip(ids,vals):animalrows.append({'panel':letter,'gene':g,'sample':gsm,'animal_id':meta.loc[gsm,'animal_id'],'group':group,'log2_RMA':val})
 row=cand.loc[g];ax.set_xticks([0,1,2],['Healthy\nn=3','Disease\nn=9','NINT\nn=10']);ax.set(xlim=(-.5,2.5),ylabel='Deposited log2 RMA signal',title=g.title()+' expression by animal');ax.margins(y=.32)
 ax.text(.02,.98,f'Exact P={row.p_exact_animal_permutation:.4f}\nFixed-family q={row.q_exact_full_fixed_family:.4f}',transform=ax.transAxes,va='top',fontsize=7.1)
 source(letter,'Figure6_animal_expression_source.csv','22 individual animals; exact treatment test 19','points individual animals, bars group means; exact BH family 17')
# G
available=cand[cand.measurement_available]
for j,(g,t) in enumerate(available.iterrows()):
 a[6].plot([t.leave_one_animal_min,t.leave_one_animal_max],[j,j],color='#9cabb5',lw=1.5);a[6].scatter(t.log2FC,j,s=23,color='#cd8b41' if g=='KDR' else '#168c8c')
a[6].set_yticks(range(len(available)),available.index);a[6].invert_yaxis();a[6].axvline(0,lw=.7,c='#999');a[6].set(xlabel='Log2 effect and leave-one-animal range',title='Stability across nineteen animal omissions')
source('G','fixed_targets_and_clinical_proteins.csv','each omitted disease/treated animal','full-sample effect and 19 omission estimate min/max; not confidence intervals')
# H
totals=members.groupby('module').size()
for j,(m,t) in enumerate(modules.iterrows()):
 a[7].plot([t.CI95_low,t.CI95_high],[j,j],c='#9cabb5',lw=1);a[7].scatter(t.effect_standardized_module,j,s=24,c='#168c8c')
labels=[('DILI proteins' if m=='Hepatic anchors' else m)+f' ({int(t.n_genes)}/{int(totals[m])})' for m,t in modules.iterrows()]
a[7].set_yticks(range(len(modules)),labels);a[7].set_ylim(len(modules),-.5);a[7].axvline(0,c='#aaa',lw=.7);a[7].set(xlabel='Module-score effect (Welch 95% CI)',title='Programs evaluated at observed coverage');a[7].text(.97,.03,'Exact q across 7 modules: all ≥ '+f'{modules.q_exact.min():.3f}',ha='right',transform=a[7].transAxes,fontsize=6.8)
source('H','module_contrasts.csv; module_gene_membership.csv','19 disease/treated animals; standardized over all 22','Welch CI; exact BH across seven modules; measured/requested member coverage')
# I
names=['NINT_vs_DISEASE','DISEASE_vs_CONTROL','NINT_vs_CONTROL'];keys=['BH_q05','BH_q05_FC05','Welch_q05'];labs=['Moderated q < .05','+ |log2 effect| ≥ .5','Welch q < .05'];cs=['#168c8c','#617f91','#cd8b41']
for k,key in enumerate(keys):
 vals=[summary['DE_counts'][n][key] for n in names];a[8].bar(np.arange(3)+(k-1)*.25,vals,width=.23,color=cs[k],label=labs[k])
 for j,v in enumerate(vals):a[8].text(j+(k-1)*.25,v+8,str(v),ha='center',fontsize=6)
a[8].set_xticks(range(3),['NINT−\ndisease','Disease−\nhealthy','NINT−\nhealthy']);a[8].set(ylabel='Genome-wide significant features',title='Genome-wide contrast comparison');a[8].margins(y=.25);a[8].legend(frameon=False,fontsize=6.8,loc='upper left')
source('I','analysis_summary.json','22 independent animals','three contrast-specific genome-wide multiplicity families')
# J
jdata=disease[['log2FC','q_value_BH_genome']].rename(columns={'log2FC':'disease_log2FC'}).join(de.log2FC.rename('treatment_log2FC'));jdata=jdata[jdata.q_value_BH_genome.lt(.05)&jdata.disease_log2FC.abs().ge(.5)]
a[9].scatter(jdata.disease_log2FC,jdata.treatment_log2FC,s=10,c='#739bab',alpha=.6,rasterized=True);a[9].axhline(0,lw=.6,c='#bbb');a[9].axvline(0,lw=.6,c='#bbb');lim=max(jdata.disease_log2FC.abs().max(),jdata.treatment_log2FC.abs().max());a[9].plot([-lim,lim],[lim,-lim],ls='--',c='#999',lw=.7)
opp=(jdata.disease_log2FC*jdata.treatment_log2FC<0).mean();j_rho=spearmanr(jdata.disease_log2FC,jdata.treatment_log2FC).statistic
a[9].set(xlabel='Disease − healthy (log2 RMA)',ylabel='Nintedanib − disease (log2 RMA)',title='Descriptive disease–treatment geometry');a[9].text(.02,.98,f'{len(jdata)} disease-selected features\nOpposite direction {opp:.1%}; ρ={j_rho:.3f}',va='top',transform=a[9].transAxes,fontsize=7,bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':1})
source('J','DISEASE_vs_CONTROL_all_genes.csv; NINT_vs_DISEASE_all_genes.csv','genes descriptive; shared disease animals','disease q<.05 and abs log2 effect≥.5; no independent-feature inference')
# K
v=np.sort(splits.opposite_direction_fraction);a[10].step(v,np.arange(1,len(v)+1)/len(v),where='post',lw=1.7,color='#168c8c');a[10].xaxis.set_major_formatter(PercentFormatter(1,decimals=0));a[10].set(xlabel='Opposite-direction fraction',ylabel='Empirical cumulative fraction',title='Disjoint disease-reference sensitivity',ylim=(0,1.06));a[10].text(.02,.98,f'{len(v)} splits; median {np.median(v):.1%}\nRange {v.min():.1%}–{v.max():.1%}',transform=a[10].transAxes,va='top',fontsize=7.2)
source('K','disjoint_disease_reference_sensitivity.csv','oriented 4/5 disease-animal reference partition','126 overlapping splits, effect-only feature selection; descriptive sensitivity')
# L: independent, exact animal-randomization assessment of all six frozen algorithms.
heldout_path=R/'heldout_rank_association/deposited_matrix/external_rank_association_summary.csv'
alg=['legacy_consensus','median','mean','geometric','maximin','RRA'];alabel=['Baseline','Median','Mean','Geometric','Maximin','RRA']
l=pd.read_csv(heldout_path).set_index('algorithm').reindex(alg)
for j,(name,t) in enumerate(l.iterrows()):
 rho=t.spearman_priority_vs_abs_drug_effect
 a[11].plot([0,rho],[j,j],c='#aebec6',lw=1)
 a[11].scatter(rho,j,s=27,c='#168c8c')
 a[11].text(.48,j,f'q={t.q_BH_six_algorithms_two_sided:.3f}',ha='left',va='center',fontsize=7)
a[11].set_yticks(range(6),alabel);a[11].set_ylim(5.7,-.7);a[11].axvline(0,c='#bbb',lw=.7)
a[11].set(xlim=(-.55,.88),xlabel='Priority–|liver effect| Spearman ρ',title='Independent hepatic rank association')
a[11].text(.98,.02,'12/17 targets; 92,378 animal allocations',ha='right',transform=a[11].transAxes,fontsize=6.8)
source('L','Heldout_rank/external_rank_association_summary.csv','12 measured candidates; 19 disease/treated animals','frozen negative-rank priority versus absolute effect; two-sided exact animal-label P, BH across six algorithms')
stem='Figure_6_Independent_Hepatic_Pharmacodynamics'
for ext in ['svg','pdf']:fig.savefig(O/(stem+'.'+ext),bbox_inches='tight',pad_inches=.1)
fig.savefig(O/(stem+'.png'),dpi=145,bbox_inches='tight',pad_inches=.1);plt.close(fig)
temp=R/'figure_qa'/f'{stem}_600dpi.png'
subprocess.run(['inkscape',str(O/(stem+'.svg')),'--export-type=png','--export-area-page','--export-dpi=600','--export-background=white','--export-background-opacity=1',f'--export-filename={temp}'],check=True,capture_output=True)
subprocess.run(['convert',str(temp),'-background','white','-alpha','remove','-alpha','off','-units','PixelsPerInch','-density','600','-compress','LZW',str(O/(stem+'.tiff'))],check=True,capture_output=True)
with Image.open(O/(stem+'.tiff')) as im:im.load();assert tuple(map(float,im.info['dpi']))==(600.,600.)
temp.unlink()
pd.DataFrame(maprows).to_csv(O/'Figure6_panel_source_map.csv',index=False)
pd.DataFrame(animalrows).to_csv(O/'Figure6_animal_expression_source.csv',index=False)
jdata.to_csv(O/'Figure6_J_descriptive_geometry_source.csv');l.to_csv(O/'Figure6_L_independent_rank_association_source.csv')
obsolete=O/'Figure6_L_frozen_rank_source.csv'
if obsolete.exists():obsolete.unlink()
qa={'source_directory':str(S),'total_animals':len(meta),'treatment_test_animals':19,'candidate_family':17,
    'available_candidates':len(available),'module_family':7,'PCA_variance':summary['pca_variance'][:2],
    'genome_wide_treatment_discoveries':int(sig.sum()),'J_opposite_fraction':opp,'J_descriptive_rho':j_rho,
    'L_available_targets':12,'L_algorithms':6,'L_exact_allocations':92378,
    'L_minimum_q_six_algorithms':float(l.q_BH_six_algorithms_two_sided.min()),
    'expression_source_sha256':hashlib.sha256(exprpath.read_bytes()).hexdigest()}
(R/'figure_qa/Figure6_numerical_source_qa.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa,indent=2),flush=True)
