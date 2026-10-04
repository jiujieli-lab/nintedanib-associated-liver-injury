#!/usr/bin/env python3
"""Publication supplement, exact quantitative output only; no generated data."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/virtual-perturbation-mpl')
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import LogLocator
OUT=Path(__file__).resolve().parent
COMPS=['Hepatocyte','Endothelial','Macrophage'];CLAB={'Hepatocyte':'Hep','Endothelial':'Endo','Macrophage':'Mac'}
TARGS=['PDGFRA','FLT4','PDGFRA+FLT4'];TLAB={'PDGFRA':'PDGFRA','FLT4':'FLT4','PDGFRA+FLT4':'Dual'}
COL={'PDGFRA':'#27678F','FLT4':'#C47D31','PDGFRA+FLT4':'#333C45'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':.6,'savefig.facecolor':'white','pdf.fonttype':42,'svg.fonttype':'none'})
contract={'question':'Do target-specific graph deletions produce clinical-signature propagation changes beyond matched graph controls?','takeaway':'All nine magnitude calibrations and all three paired nonadditivity calibrations are nonsignificant; sparse target detection limits biological applicability.','surface':'standalone publication supplementary TIFF, PDF and preview PNG','figure':'Supplementary Figure S11','palette':'blue and gold for single targets; charcoal for joint; blue/gold signed scale; neutral calibration','source_tables':['primary_calibration.csv','perturbation_metrics.csv','donor_target_expression.csv','gene_level_perturbations.csv.gz','edge_attenuation_sensitivity.csv','nonadditivity_by_donor_omission.csv','nonadditivity_calibration.csv','target_to_target_redistribution.csv'],'panels':{'A':'5-donor detection heatmap, 6 target-context rows','B':'9 primary magnitudes plus five-donor omission ranges','C':'9 BH q values; magnitude tests','D':'13 clinical gene score changes x 3 interventions, endothelial','E':'five prespecified graph attenuation settings, endothelial','F':'4 target/comparator deletions x3 compartments','G':'3 paired nonadditivity magnitudes with five donor omission ranges','H':'three paired nonadditivity nominal p and q values','I':'6 reciprocal target-score shifts; descriptive, undirected'},'export_dimensions':'8.4 x 8.4 inches, 1200 dpi TIFF (10080 x 10080 pixels); 220 dpi preview; vector PDF'}
(OUT/'figure_s11_chart_contract.json').write_text(json.dumps(contract,indent=2))
cal=pd.read_csv(OUT/'primary_calibration.csv');met=pd.read_csv(OUT/'perturbation_metrics.csv');det=pd.read_csv(OUT/'donor_target_expression.csv');g=pd.read_csv(OUT/'gene_level_perturbations.csv.gz');atten=pd.read_csv(OUT/'edge_attenuation_sensitivity.csv');inter=pd.read_csv(OUT/'nonadditivity_by_donor_omission.csv');ncal=pd.read_csv(OUT/'nonadditivity_calibration.csv');rec=pd.read_csv(OUT/'target_to_target_redistribution.csv')
fig,axs=plt.subplots(3,3,figsize=(8.4,8.4));fig.subplots_adjust(left=.09,right=.98,bottom=.065,top=.965,wspace=.56,hspace=.69)
def panel(ax,label,title):
 ax.set_title(title,loc='left',pad=8,fontweight='normal');ax.text(-.26,1.12,label,transform=ax.transAxes,fontweight='bold',fontsize=11,va='top');ax.spines[['top','right']].set_visible(False)
def heat(ax,z,rows,cols,cmap,fmt='.3g',vmin=None,vmax=None):
 im=ax.imshow(z,aspect='auto',cmap=cmap,vmin=vmin,vmax=vmax)
 ax.set_xticks(range(len(cols)),cols);ax.set_yticks(range(len(rows)),rows)
 for i in range(len(rows)):
  for j in range(len(cols)):
   val=z[i,j];ax.text(j,i,format(val,fmt),ha='center',va='center',fontsize=5.8,color='white' if vmin is not None and vmax is not None and ((val<vmin+.32*(vmax-vmin)) if str(cmap).endswith('_r') else (val>vmin+.68*(vmax-vmin))) else '#18202A')
 ax.tick_params(length=0);return im
# A: every donor retained
ax=axs[0,0];panel(ax,'A','Target detection across donors (%)')
rows=[(c,t) for c in COMPS for t in ['PDGFRA','FLT4']];donors=sorted(det.donor.unique());z=np.array([[det[(det.compartment==c)&(det.target==t)&(det.donor==d)].detection_fraction.iloc[0]*100 for d in donors] for c,t in rows]);heat(ax,z,[CLAB[c]+' / '+t for c,t in rows],[d.replace('TLH','') for d in donors],'Blues','.1f',0,35)
# B full and donor omission
ax=axs[0,1];panel(ax,'B','Clinical-score perturbation magnitude')
for y,(c,t) in enumerate((c,t) for c in COMPS for t in TARGS):
 v=cal[(cal.compartment==c)&(cal.intervention==t)].clinical_weighted_relative_change.iloc[0]*100
 ds=met[(met.compartment==c)&(met.intervention==t)&(met.omitted_donor!='none')].clinical_weighted_relative_change.to_numpy()*100
 ax.plot([ds.min(),ds.max()],[y,y],color=COL[t],lw=1.4);ax.scatter([v],[y],s=17,color=COL[t],zorder=3);ax.scatter(ds,[y]*len(ds),s=5,color=COL[t],alpha=.5)
ax.set_yticks(range(9),[CLAB[c]+' / '+TLAB[t] for c in COMPS for t in TARGS]);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlabel('Weighted score change (%)');ax.set_xlim(2e-5,.08);ax.xaxis.set_major_locator(LogLocator(base=10,numticks=5));ax.grid(axis='x',alpha=.2)
# C statistical calibration
ax=axs[0,2];panel(ax,'C','Matched-control calibration (BH q)')
z=np.array([[cal[(cal.compartment==c)&(cal.intervention==t)].q_bh_9.iloc[0] for t in TARGS] for c in COMPS]);heat(ax,z,[CLAB[c] for c in COMPS],['PDGFRA','FLT4','Dual'],'Greys_r','.3f',0,1)
ax.text(.5,-.27,'9-test family; no q < 0.05',transform=ax.transAxes,ha='center',fontsize=6.5)
# D endothelial clinical anchor response
ax=axs[1,0];panel(ax,'D','Endothelial clinical-gene scores')
sub=g[g.compartment.eq('Endothelial')&g.is_clinical_anchor&g.intervention.isin(TARGS)];anchors=['ACO1','ASS1','FAH','CPS1','ALDOB','HPD','OTC','DMGDH','GSTA1','FBP1','PCK2','CES1','LECT2'];z=sub.pivot(index='gene',columns='intervention',values='relative_score_change').reindex(index=anchors,columns=TARGS).to_numpy()*100
lim=float(np.abs(z).max());cmap=LinearSegmentedColormap.from_list('blue_white_gold',['#27678F','#FFFFFF','#C47D31']);im=ax.imshow(z,aspect='auto',cmap=cmap,vmin=-lim,vmax=lim);ax.set_yticks(range(13),anchors);ax.set_xticks(range(3),['PDGFRA','FLT4','Dual']);ax.tick_params(length=0);cb=fig.colorbar(im,ax=ax,fraction=.05,pad=.03);cb.ax.tick_params(labelsize=5);cb.set_label('Score change (%)',fontsize=6)
# E prescribed operator attenuation
ax=axs[1,1];panel(ax,'E','Endothelial edge attenuation')
for t,marker in zip(TARGS,['o','s','^']):
 d=atten[(atten.compartment=='Endothelial')&(atten.intervention==t)].sort_values('fraction_incident_edge_attenuation');ax.plot(d.fraction_incident_edge_attenuation,d.clinical_weighted_relative_change*100,label=TLAB[t],color=COL[t],marker=marker,ms=3,lw=1)
ax.set_xlabel('Incident-edge attenuation');ax.set_ylabel('Weighted score change (%)');ax.set_xticks([0,.25,.5,.75,1]);ax.legend(frameon=False,fontsize=6,loc='upper left');ax.set_ylim(bottom=0)
# F comparator context
ax=axs[1,2];panel(ax,'F','Comparator perturbations (%)')
labels=['PDGFRA','FLT4','KDR','FGFR1'];z=np.array([[met[(met.compartment==c)&(met.intervention==t)&(met.omitted_donor=='none')].clinical_weighted_relative_change.iloc[0]*100 for t in labels] for c in COMPS]);heat(ax,z,[CLAB[c] for c in COMPS],labels,'Blues','.3g',0,float(z.max()));ax.tick_params(axis='x',labelrotation=45)
# G nonadditivity magnitude with ranges
ax=axs[2,0];panel(ax,'G','Joint-response nonadditivity')
for y,c in enumerate(COMPS):
 obs=inter[(inter.compartment==c)&(inter.omitted_donor=='none')].clinical_nonadditivity_relative_L1.iloc[0];ds=inter[(inter.compartment==c)&(inter.omitted_donor!='none')].clinical_nonadditivity_relative_L1.to_numpy();ax.plot([ds.min(),ds.max()],[y,y],color='#333C45',lw=1.4);ax.scatter([obs],[y],s=20,color='#333C45');ax.scatter(ds,[y]*len(ds),s=6,color='#333C45',alpha=.5)
ax.set_yticks(range(3),[CLAB[c] for c in COMPS]);ax.set_xscale('log');ax.set_xlim(1e-12,1e-4);ax.invert_yaxis();ax.set_xlabel('Weighted absolute interaction / baseline');ax.xaxis.set_major_locator(LogLocator(base=10,numticks=4));ax.grid(axis='x',alpha=.2)
# H nonadditivity calibration
ax=axs[2,1];panel(ax,'H','Nonadditivity calibration')
y=np.arange(3)
for shift,field,col,label in [(-.17,'empirical_p','#B8BEC4','Nominal P'),(.17,'q_bh_3','#333C45','BH q')]:
 vals=[ncal[ncal.compartment==c][field].iloc[0] for c in COMPS];ax.barh(y+shift,vals,height=.3,color=col,label=label);[ax.text(min(v+.02,.95),yy,f'{v:.3f}',ha='left' if v<.9 else 'right',va='center',fontsize=5.8,color='white' if (v>.9 and field=='q_bh_3') else '#18202A') for v,yy in zip(vals,y+shift)]
ax.set_yticks(y,[CLAB[c] for c in COMPS]);ax.invert_yaxis();ax.set_xlim(0,1.04);ax.set_xlabel('Matched-control probability');ax.axvline(.05,color='#666666',ls=':',lw=.7);ax.legend(frameon=False,fontsize=6,loc='lower right')
# I reciprocal redistribution
ax=axs[2,2];panel(ax,'I','Reciprocal target-score shifts')
for shift,(deleted,readout),col,lab in [(-.17,('PDGFRA','FLT4'),COL['PDGFRA'],'PDGFRA deleted'),(.17,('FLT4','PDGFRA'),COL['FLT4'],'FLT4 deleted')]:
 vals=[rec[(rec.compartment==c)&(rec.deleted_target==deleted)&(rec.readout_target==readout)].relative_change.iloc[0]*100 for c in COMPS];ax.barh(y+shift,vals,height=.3,color=col,label=lab)
ax.set_yticks(y,[CLAB[c] for c in COMPS]);ax.invert_yaxis();ax.set_xlabel('Propagation-score change (%)');ax.axvline(0,color='#555555',lw=.6);ax.legend(frameon=False,fontsize=5.5,loc='upper right')
from matplotlib.transforms import Bbox
fig.canvas.draw();renderer=fig.canvas.get_renderer();paneldir=OUT/'Individual_Panels';paneldir.mkdir(exist_ok=True)
for letter,ax in zip('ABCDEFGHI',axs.ravel()):
 boxes=[ax.get_tightbbox(renderer)]
 if letter=='D': boxes.append(cb.ax.get_tightbbox(renderer))
 bbox=Bbox.union(boxes).transformed(fig.dpi_scale_trans.inverted()).padded(.04)
 fig.savefig(paneldir/f'Figure_S11_Panel_{letter}_600dpi.png',dpi=600,bbox_inches=bbox)
fig.savefig(OUT/'Figure_S11_Virtual_Perturbation.pdf');fig.savefig(OUT/'Figure_S11_Virtual_Perturbation_Preview.png',dpi=220)
fig.savefig(OUT/'Figure_S11_Virtual_Perturbation_1200dpi.tiff',dpi=1200,pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
print('Exported Figure S11')
