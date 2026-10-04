from pathlib import Path
import json,hashlib,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lxml import etree as E

BASE=Path(__file__).resolve().parent
SRC=BASE/'integrated_analysis'
OUT=BASE/'figures_updated/Figures'
QA=BASE/'figure_qa'
OUT.mkdir(parents=True,exist_ok=True)
QA.mkdir(parents=True,exist_ok=True)
audit=pd.read_csv(SRC/'vko_numerical_run_audit.csv')
projection=pd.read_csv(SRC/'target_protein_projection_by_seed.csv')
assert audit.shape[0]==210
assert (audit['numerically_nonzero']==(audit['max']>1e-12)).all()
for threshold in [1e-14,1e-12,1e-10]:assert ((audit['max']>threshold)==audit['numerically_nonzero']).all()
drug=audit[audit.ko_class.eq('drug_target')].copy()
valid=projection[projection.cohort.eq('confirmatory')&projection.contrast.eq('DO_vs_NDO')&projection.numerically_nonzero].copy()
assert len(drug)==102 and drug.numerically_nonzero.sum()==29 and len(valid)==29
fgfr=audit[audit.ko_label.eq('FGFR1')&audit.compartment.eq('Endothelial')].sort_values('seed')
assert len(fgfr)==3 and fgfr.numerically_nonzero.sum()==1

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':7.7,
 'axes.titleweight':'bold','axes.labelsize':7,'xtick.labelsize':6.3,'ytick.labelsize':6.3,
 'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.55,
 'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2.5,'ytick.major.size':2.5,
 'axes.labelcolor':'#262626','text.color':'#262626','axes.edgecolor':'#262626'})
W,H=585.764912,294
fig=plt.figure(figsize=(W/72,H/72))
axes=[fig.add_axes([.074,.57,.385,.31]),fig.add_axes([.61,.57,.34,.31]),
      fig.add_axes([.074,.105,.385,.31]),fig.add_axes([.61,.105,.34,.31])]
colors={'Hepatocyte':'#d28a78','Endothelial':'#70b0b5','Macrophage':'#7786b5'}
comps=['Hepatocyte','Endothelial','Macrophage']
def title(ax,letter,text):
 ax.set_title(text,loc='left',pad=8)
 ax.text(-.11,1.09,letter,transform=ax.transAxes,fontweight='bold',fontsize=9.6)

# G: Unranked original distances reveal the finite-precision cluster.
ax=axes[0]
for j,c in enumerate(comps):
 d=audit[audit.compartment.eq(c)].sort_values(['ko_class','ko_label','seed'])
 jitter=np.linspace(-.19,.19,len(d))
 ax.scatter(j+jitter,d['max'],s=9,c=np.where(d.numerically_nonzero,colors[c],'#c5c8ca'),alpha=.85,linewidths=0)
ax.axhline(1e-12,color='#ad655e',ls='--',lw=.75)
ax.set_yscale('log');ax.set_ylim(2e-16,.03);ax.set_yticks([1e-15,1e-12,1e-9,1e-6,1e-3])
ax.set_xticks(range(3),['Hepatocyte','Endothelial','Macrophage']);ax.set_ylabel('Maximum unsigned distance')
ax.text(.01,.98,'132 / 210 runs at numerical zero',va='top',transform=ax.transAxes,fontsize=6.1)
ax.text(2.43,1e-12,'10⁻¹²',va='center',ha='right',fontsize=5.8,color='#ad655e',bbox={'facecolor':'white','edgecolor':'none','pad':.3})
title(ax,'G','Raw distances reveal a numerical-zero cluster')

# H: Run denominators are preserved; no zero-distance rank is retained.
ax=axes[1]
counts=drug.groupby('compartment').numerically_nonzero.agg(['sum','count']).reindex(comps)
for j,(c,row) in enumerate(counts.iterrows()):
 ax.barh(j,row['count'],height=.55,color='#dedfe1')
 ax.barh(j,row['sum'],height=.55,color=colors[c])
 ax.text(row['count']+.5,j,f"{int(row['sum'])} / {int(row['count'])}",va='center',fontsize=6.2)
ax.set_yticks(range(3),comps);ax.invert_yaxis();ax.set_xlim(0,43);ax.set_xticks([0,10,20,30,40]);ax.set_xlabel('Drug-target runs (nonzero / supplied)')
ax.text(.99,.98,'29 / 102 retained',ha='right',va='top',transform=ax.transAxes,fontsize=6.1)
title(ax,'H','Numerical calibration reduces eligible runs')

# I: FGFR1 has only one nonzero endothelial seed.
ax=axes[2]
ax.scatter(range(3),fgfr['max'],s=25,c=np.where(fgfr.numerically_nonzero,colors['Endothelial'],'#a6a9ac'),edgecolors='white',linewidths=.4,zorder=3)
ax.axhline(1e-12,color='#ad655e',ls='--',lw=.75);ax.set_yscale('log');ax.set_ylim(2e-16,.001);ax.set_yticks([1e-15,1e-12,1e-9,1e-6,1e-3]);ax.set_xlim(-.4,2.4)
ax.set_xticks(range(3),[str(x) for x in fgfr.seed]);ax.set_xlabel('Endothelial network seed');ax.set_ylabel('Maximum FGFR1-KO distance')
for j,row in enumerate(fgfr.itertuples()):
 ax.annotate(f'{row.max:.2e}',(j,row.max),xytext=(0,7),textcoords='offset points',ha='center',fontsize=6)
ax.text(.99,.96,'1 / 3 seeds nonzero',ha='right',va='top',transform=ax.transAxes,fontsize=6.1)
title(ax,'I','FGFR1 lacks cross-seed numerical support')

# J: Clinical-effect-weighted injury-protein projection against 4,000 matched sets.
ax=axes[3]
for j,c in enumerate(comps):
 d=valid[valid.compartment.eq(c)].sort_values(['target','seed'])
 xx=j+np.linspace(-.17,.17,len(d)) if len(d)>1 else np.array([j])
 ax.scatter(xx,d.matched_z,s=13,color=colors[c],edgecolors='white',linewidths=.35,zorder=3)
 fg=d.target.eq('FGFR1')
 if fg.any():
  x=xx[np.where(fg)[0][0]];y=d.loc[fg,'matched_z'].iloc[0]
  ax.annotate('FGFR1',(x,y),xytext=(6,-2),textcoords='offset points',fontsize=6.1)
ax.axhline(0,color='#a6a9ac',lw=.65);ax.set_xlim(-.4,2.4);ax.set_ylim(-2.65,.7)
ax.set_xticks(range(3),comps);ax.set_ylabel('Matched-null z score');ax.set_xlabel('Numerically nonzero drug-target runs')
ax.text(.02,.98,f"29 runs; minimum P = {valid.matched_empirical_p.min():.3f}",va='top',transform=ax.transAxes,fontsize=6.1)
title(ax,'J','No enriched weighted injury-protein response')
fragment=QA/'Figure4_QC_lower_panels.svg';fig.savefig(fragment,format='svg');plt.close(fig)

# Retain original A-F and their color bars, replacing only obsolete G-J and title.
old=OUT/'Figure_4_Liver_Localization_and_Virtual_Knockout.svg'
original=BASE/'Templates/Figure_4_Original.svg'
tree=E.parse(str(original));ns={'s':'http://www.w3.org/2000/svg'}
for n in tree.xpath('//s:text',namespaces=ns):
 if n.text=='Direct-target localization across liver compartments':n.text='Drug-target localization across liver compartments'
 if n.text=='DILI-anchor localization across liver compartments':n.text='Injury-protein localization across liver compartments'
keep={}
for a in tree.xpath('//s:g[starts-with(@id,"axes_")]',namespaces=ns):
 if a.get('id') in ['axes_7','axes_8','axes_9','axes_10','axes_13','axes_14']:
  a.getparent().remove(a)
 else:keep[a.get('id')]=hashlib.sha256(E.tostring(a)).hexdigest()
for n in tree.xpath('//s:text',namespaces=ns):
 if (n.text or '').startswith('Figure 4 |'):
  n.text='Figure 4 | Liver-cell localization and numerical calibration of virtual perturbation'
root=tree.getroot();root.set('height','840pt');root.set('viewBox','0 0 585.764912 840')
sub=E.parse(str(fragment)).getroot()
# Prefix every generated ID and local reference to avoid clip/marker collisions.
ids={n.get('id'):'qc4_'+n.get('id') for n in sub.iter() if n.get('id')}
for n in sub.iter():
 for k,v in list(n.attrib.items()):
  if k=='id':n.set(k,ids[v])
  else:
   for source,target in ids.items():v=v.replace('url(#'+source+')','url(#'+target+')')
   if v.startswith('#') and v[1:] in ids:v='#'+ids[v[1:]]
   n.set(k,v)
sub.set('x','0');sub.set('y','540');sub.set('width',str(W));sub.set('height',str(H));root.append(sub)
for a in tree.xpath('//s:g[starts-with(@id,"axes_")]',namespaces=ns):
 assert hashlib.sha256(E.tostring(a)).hexdigest()==keep[a.get('id')]
new=OUT/'Figure_4_Liver_Localization_and_Perturbation_QC.svg'
tree.write(str(new),xml_declaration=True,encoding='utf-8')
for ext in ['svg','pdf','tiff']:
 stale=old.with_suffix('.'+ext)
 if stale.exists():stale.unlink()
report={'A_F_and_colorbar_geometry_identical':True,'all_runs':len(audit),
 'all_numerical_zero':int((~audit.numerically_nonzero).sum()),'drug_target_runs':len(drug),
 'drug_target_nonzero':int(drug.numerically_nonzero.sum()),'FGFR1_endothelial_nonzero':int(fgfr.numerically_nonzero.sum()),
 'threshold':1e-12,'identical_threshold_classification':[1e-14,1e-12,1e-10],
 'minimum_weighted_null_P':float(valid.matched_empirical_p.min()),
 'sources':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [SRC/'vko_numerical_run_audit.csv',SRC/'target_protein_projection_by_seed.csv']}}
(QA/'figure4_numeric_qc_provenance.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
