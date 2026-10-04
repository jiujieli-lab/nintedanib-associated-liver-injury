#!/usr/bin/env python3
"""Assay-stratified reanalysis of released nintedanib evidence; no docking or MD.

Chart contract: static publication figure; raw assay points and endpoint-separated
summaries address whether both targets have experimental engagement evidence.
Two target colors, endpoint marker shapes, hollow sensitivity points. 26 raw
records, 11 frozen eligible records, 3 additional BAO-format sensitivity records.
No null-result removal or inferential pooling across heterogeneous endpoints.
"""
from pathlib import Path
import json, hashlib, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
raw=pd.read_csv(ROOT/'chembl_all_nintedanib_activities.csv')
frozen=pd.read_csv(ROOT/'chembl_human_single_protein_activities.csv')
targetmap={'CHEMBL2007':'PDGFRA','CHEMBL1955':'FLT4'}
x=raw[raw.target_chembl_id.isin(targetmap)].copy()
x['gene_symbol']=x.target_chembl_id.map(targetmap)
x['in_frozen_analysis']=x.activity_id.isin(frozen.activity_id)
valid=x.target_organism.eq('Homo sapiens') & x.standard_relation.eq('=') & x.pchembl_value.notna() & x.data_validity_comment.isna() & x.potential_duplicate.ne(1)
x['included_format_sensitivity']=valid & x.standard_type.isin(['IC50','Kd','Ki'])
x['audit_category']=np.select([x.in_frozen_analysis,x.included_format_sensitivity,x.potential_duplicate.eq(1),x.standard_type.isin(['kon','k_off'])],['Frozen eligible','Additional assay format','Duplicate flag','Kinetics without value'],default='Single-concentration inhibition')
x['p_value_note']='pChEMBL is a potency transform; not a statistical P value.'
x.to_csv(OUT/'Binding_All_Target_Records_Audited.csv',index=False)
primary=x[x.in_frozen_analysis].copy();expanded=x[x.included_format_sensitivity].copy()
assert len(x)==26 and len(primary)==11 and len(expanded)==14
rows=[];loo=[]
for scope,df in [('Frozen primary',primary),('BAO-format sensitivity',expanded)]:
 for (g,e),v in df.groupby(['gene_symbol','standard_type']):
  d=v.groupby('document_chembl_id').pchembl_value.median()
  rows.append(dict(scope=scope,gene_symbol=g,endpoint=e,n_records=len(v),n_assays=v.assay_chembl_id.nunique(),n_documents=len(d),min_nM=v.standard_value.min(),median_nM=v.standard_value.median(),max_nM=v.standard_value.max(),median_pchembl=v.pchembl_value.median(),min_pchembl=v.pchembl_value.min(),max_pchembl=v.pchembl_value.max(),document_balanced_median_pchembl=d.median(),all_below_1000nM=bool(v.standard_value.lt(1000).all())))
  if len(d)>=2:
   for drop in d.index:
    vv=d.drop(drop);loo.append(dict(scope=scope,gene_symbol=g,endpoint=e,omitted_document=drop,remaining_documents=len(vv),median_pchembl=vv.median(),median_equivalent_nM=10**(9-vv.median())))
summary=pd.DataFrame(rows);summary.to_csv(OUT/'Binding_Assay_Stratified_Summary.csv',index=False)
loo=pd.DataFrame(loo);loo.to_csv(OUT/'Binding_Leave_One_Document_Out.csv',index=False)
pairs=[]
for (doc,endpoint),d in expanded.groupby(['document_chembl_id','standard_type']):
 if d.gene_symbol.nunique()!=2:continue
 med=d.groupby('gene_symbol').standard_value.median()
 pairs.append(dict(document_chembl_id=doc,document_year=int(d.document_year.iloc[0]),endpoint=endpoint,PDGFRA_nM=med['PDGFRA'],FLT4_nM=med['FLT4'],log10_PDGFRA_over_FLT4=np.log10(med['PDGFRA']/med['FLT4']),fold_PDGFRA_over_FLT4=med['PDGFRA']/med['FLT4'],both_in_frozen=bool(d.in_frozen_analysis.all())))
pairs=pd.DataFrame(pairs);pairs.to_csv(OUT/'Binding_Within_Document_Contrasts.csv',index=False)
audit=pd.crosstab(x.gene_symbol,x.audit_category);audit.to_csv(OUT/'Binding_Record_Audit_Counts.csv')
raw_pred=json.loads((ROOT/'Inductive_Bio_Raw_Predictions.json').read_text());smiles=raw.canonical_smiles.dropna().iloc[0]
fresh=json.loads((ROOT/'Nintedanib_ChEMBL.json').read_text())
assert smiles==fresh['molecule_structures']['canonical_smiles']
pred=[]
for b in raw_pred['content']:
 a=json.loads(b['text']);q=a['predictions'][0];assert q['smiles']==smiles
 pred.append(dict(compound='Nintedanib neutral parent, Z geometry',exact_smiles=smiles,model_id=a['model_id'],model_version=a['model_version'],config_version=a['config_version'],predicted_value=q['continuous_prediction'],returned_lower_bound=q['continuous_prediction_low'],returned_upper_bound=q['continuous_prediction_high'],interval_coverage='not stated by connector; not labeled as 95% CI',units='dimensionless',out_of_domain_flag=q['out_of_domain_flag'],low_confidence_flag=q['low_confidence_flag'],status=q['status'],retrieved_UTC='2026-10-02'))
pd.DataFrame(pred).to_csv(OUT/'Inductive_Physicochemical_Predictions.csv',index=False)
(OUT/'SMILES_Verification.json').write_text(json.dumps({'all_source_and_prediction_strings_identical':True,'canonical_smiles_length':len(smiles),'stereobond_backslash_count':smiles.count(chr(92)),'source_standard_inchi_key':fresh['molecule_structures']['standard_inchi_key'],'source_heavy_atoms':fresh['molecule_properties']['heavy_atoms'],'interpretation':'Byte-identical molecular strings; JSON backslash escaping does not change the submitted stereobond.'},indent=2)+'\n')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':10,'axes.titleweight':'bold','axes.linewidth':.7,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
color={'PDGFRA':'#286B8C','FLT4':'#C57E34'}
fig,axes=plt.subplots(3,2,figsize=(8.4,8.4),layout='constrained')
def style(ax,letter,title):
 ax.set_title(title,loc='left',pad=9)
 ax.text(-.14,1.05,letter,transform=ax.transAxes,fontsize=13,weight='bold')
 ax.spines[['top','right']].set_visible(False)
 ax.grid(axis='x',color='#E6E8EA',linewidth=.55,zorder=0)
def assay_panel(ax,gene,letter):
 d=expanded[expanded.gene_symbol.eq(gene)].sort_values(['document_year','document_chembl_id','standard_type']).reset_index(drop=True)
 labels=[]
 for i,r in d.iterrows():
  mark={'IC50':'o','Kd':'s','Ki':'^'}[r.standard_type]
  ax.scatter(r.standard_value,i,s=42,marker=mark,facecolor=color[gene] if r.in_frozen_analysis else 'white',edgecolor=color[gene],linewidth=1.1,zorder=3)
  ax.annotate(f'{r.standard_value:g}',(r.standard_value,i),xytext=(6,0),textcoords='offset points',va='center',fontsize=7.5)
  labels.append(f'{int(r.document_year)}  {r.standard_type}')
 ax.set_yticks(range(len(d)),labels);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlim(.9,1800);ax.set_xticks([1,10,100,1000],[1,10,100,1000]);ax.set_xlabel('Deposited assay value (nM; log scale)')
 style(ax,letter,f'{gene} reported measurements')
assay_panel(axes[0,0],'PDGFRA','A');assay_panel(axes[0,1],'FLT4','B')
ax=axes[1,0];d=pairs[pairs.endpoint.eq('IC50')].sort_values(['document_year','document_chembl_id']).reset_index(drop=True)
for i,r in d.iterrows():
 val=r.log10_PDGFRA_over_FLT4;ax.plot([0,val],[i,i],color='#A4A9AF',lw=1.2)
 ax.scatter(val,i,s=45,facecolor='#286B8C' if r.both_in_frozen else 'white',edgecolor='#286B8C',linewidth=1.1,zorder=3)
 ax.annotate(f'{r.fold_PDGFRA_over_FLT4:.2f}×',(val,i),xytext=(5,7),textcoords='offset points',fontsize=7)
ax.axvline(0,color='#777',lw=.7);ax.set_xlim(-.5,1);ax.set_yticks(range(len(d)),[f'{int(r.document_year)} | {r.document_chembl_id.replace("CHEMBL","")}' for _,r in d.iterrows()]);ax.invert_yaxis();ax.set_xlabel('log₁₀ (PDGFRA IC₅₀ / FLT4 IC₅₀)')
style(ax,'C','Matched-document IC₅₀ contrasts')

ax=axes[1,1];idx=0;labels=[]
for gene,ep in [('PDGFRA','IC50'),('PDGFRA','Kd'),('PDGFRA','Ki'),('FLT4','IC50'),('FLT4','Kd')]:
 for scope,offset,filled in [('Frozen primary',-.12,True),('BAO-format sensitivity',.12,False)]:
  r=summary[(summary.gene_symbol==gene)&(summary.endpoint==ep)&(summary.scope==scope)].iloc[0]
  ax.plot([r.min_pchembl,r.max_pchembl],[idx+offset]*2,color=color[gene],lw=1.1,alpha=.7)
  ax.scatter(r.median_pchembl,idx+offset,s=27,marker='o',facecolor=color[gene] if filled else 'white',edgecolor=color[gene],zorder=3)
 labels.append(f'{gene} {ep}');idx+=1
ax.set_yticks(range(idx),labels);ax.invert_yaxis();ax.set_xlim(6,9);ax.set_xlabel('pChEMBL; median and observed range')
style(ax,'D','Endpoint-stratified sensitivity')

ax=axes[2,0];cats=['Frozen eligible','Additional assay format','Duplicate flag','Kinetics without value','Single-concentration inhibition'];cols=['#286B8C','#9FB7C3','#C57E34','#C9CDD1','#E4E7E9'];left=np.zeros(2)
for cat,c in zip(cats,cols):
 vals=audit.reindex(['PDGFRA','FLT4'])[cat].to_numpy();ax.barh([0,1],vals,left=left,color=c,edgecolor='white',height=.5,label=cat)
 for i,v in enumerate(vals):
  if v:ax.text(left[i]+v/2,i,str(v),ha='center',va='center',fontsize=8,color='white' if cat=='Frozen eligible' else '#333')
 left+=vals
ax.set_yticks([0,1],['PDGFRA','FLT4']);ax.invert_yaxis();ax.set_xlim(0,20);ax.set_xticks([0,5,10,15,20]);ax.set_xlabel('Raw ChEMBL activity records (count)')
style(ax,'E','Complete record disposition')
ax.legend(loc='upper left',bbox_to_anchor=(-.03,-.21),frameon=False,ncol=2,fontsize=7,handlelength=1.2,labelspacing=.3)

ax=axes[2,1];labels=[];ys=[];i=0
for gene in ['PDGFRA','FLT4']:
 for scope in ['Frozen primary','BAO-format sensitivity']:
  d=loo[(loo.gene_symbol==gene)&(loo.endpoint=='IC50')&(loo.scope==scope)]
  y=i;ys.append(y);labels.append(f'{gene}\n'+('Frozen' if scope=='Frozen primary' else 'Sensitivity'))
  jitter=np.linspace(-.09,.09,len(d));ax.scatter(d.median_pchembl,y+jitter,s=25,facecolor=color[gene] if scope=='Frozen primary' else 'white',edgecolor=color[gene],linewidth=.85,zorder=3)
  med=summary[(summary.gene_symbol==gene)&(summary.endpoint=='IC50')&(summary.scope==scope)].median_pchembl.iloc[0];ax.plot([med,med],[y-.2,y+.2],color='#222',lw=1.2)
  i+=1
ax.set_yticks(ys,labels);ax.invert_yaxis();ax.set_xlim(6.8,8.8);ax.set_xlabel('Median pIC₅₀ after document omission')
style(ax,'F','Leave-one-document-out IC₅₀')

fig.canvas.draw()
panels=OUT/'Individual_Panels';panels.mkdir(exist_ok=True)
renderer=fig.canvas.get_renderer()
for letter,ax in zip('ABCDEF',axes.flat):
 extent=ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted()).expanded(1.04,1.06)
 fig.savefig(panels/f'Figure_4_Panel_{letter}.png',dpi=600,bbox_inches=extent)
 fig.savefig(panels/f'Figure_4_Panel_{letter}.tiff',dpi=600,bbox_inches=extent,pil_kwargs={'compression':'tiff_lzw'})
 fig.savefig(panels/f'Figure_4_Panel_{letter}.pdf',bbox_inches=extent)
fig.savefig(OUT/'Figure_Binding_Evidence.png',dpi=250)
fig.savefig(OUT/'Figure_Binding_Evidence.pdf')
fig.savefig(OUT/'Figure_Binding_Evidence.svg')
fig.savefig(OUT/'Figure_Binding_Evidence_1200dpi.tiff',dpi=1200,pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)

manifest={'analysis_date_UTC':'2026-10-02','frozen_rows':len(primary),'sensitivity_rows':len(expanded),'raw_target_rows':len(x),'methods':'Exact relations, human target, valid pChEMBL, no duplicate flags; original BAO filter preserved as primary and relaxed only in labeled sensitivity; no inferential cross-assay pooling.','not_performed':['Conventional docking','Boltz structure/binding prediction: authentication pending','Molecular dynamics: no trajectories or MD engine available'],'source_hashes':{q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in [ROOT/'chembl_all_nintedanib_activities.csv',ROOT/'chembl_human_single_protein_activities.csv',ROOT/'Inductive_Bio_Raw_Predictions.json']}}
(OUT/'Binding_Analysis_Manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(summary.to_string(index=False));print(pairs.to_string(index=False));print('Outputs:',OUT)
