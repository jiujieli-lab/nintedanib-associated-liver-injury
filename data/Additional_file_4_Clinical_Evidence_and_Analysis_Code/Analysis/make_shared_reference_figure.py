from pathlib import Path
import json,string,hashlib,subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from PIL import Image

R=Path(__file__).resolve().parent
S=R/'pulmonary_sensitivity';O=R/'figures_new';O.mkdir(exist_ok=True)
summary=pd.read_csv(S/'shared_reference_summary.csv').set_index('dataset')
split=pd.read_csv(S/'disjoint_reference_splits.csv')
null=pd.read_csv(S/'animal_label_permutations.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.titlesize':10,
 'axes.titleweight':'bold','axes.labelsize':8.5,'xtick.labelsize':7.5,'ytick.labelsize':7.5,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,
 'axes.spines.right':False,'axes.linewidth':.7,'xtick.major.width':.6,'ytick.major.width':.6})
datasets=[('GSE278200','Rat'),('GSE308578','Mouse')]
metrics=[('reversal_fraction','Reversal fraction'),('rho','Disease–treatment Spearman ρ'),
 ('median_restoration','Median restoration (log2 CPM)')]
f,axes=plt.subplots(4,3,figsize=(12.4,12.8),layout='constrained')
f.get_layout_engine().set(h_pad=.14,w_pad=.10,hspace=.12,wspace=.09)
f.suptitle('Figure S10 | Shared-reference sensitivity of pulmonary nintedanib responses',fontweight='bold',fontsize=14)
panelrows=[];datarows=[]
def display(value,metric):return f'{100*value:.1f}%' if metric=='reversal_fraction' else f'{value:.3f}'
def label(a,letter):a.text(-.13,1.12,letter,transform=a.transAxes,fontweight='bold',fontsize=13)
for row,(dataset,animal) in enumerate(datasets):
 obs=summary.loc[dataset];ds=split[split.dataset.eq(dataset)]
 assert len(ds)==obs.n_splits
 for col,(metric,xlabel) in enumerate(metrics):
  a=axes[row,col];panel=string.ascii_uppercase[row*3+col];v=np.sort(ds[metric].to_numpy())
  yy=np.arange(1,len(v)+1)/len(v)
  a.step(v,yy,where='post',color='#168c8c',lw=1.6,label='Disjoint disease references')
  if len(v)<=10:a.scatter(v,yy,s=16,color='#168c8c',zorder=3)
  a.axvline(obs[metric],color='#c27236',lw=1.2,ls='--',label='Shared disease reference')
  full=np.r_[v,obs[metric]];span=max(np.ptp(full),.1)
  a.set_xlim(full.min()-.09*span,full.max()+.10*span);a.set_ylim(0,1.06)
  if metric=='reversal_fraction':a.xaxis.set_major_formatter(PercentFormatter(1,decimals=0))
  a.set(xlabel=xlabel,ylabel='Empirical cumulative fraction' if col==0 else '',
        title=f'{animal}: disjoint-reference {"reversal" if col==0 else "correlation" if col==1 else "restoration"}')
  a.text(.03,.97,f'{len(v)} splits\nMedian {display(np.median(v),metric)}\nRange {display(v.min(),metric)} to {display(v.max(),metric)}',
         transform=a.transAxes,va='top',fontsize=7.7,bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':2})
  if row==0 and col==0:a.legend(frameon=False,fontsize=6.8,loc='lower right')
  label(a,panel)
  panelrows.append({'figure':'Figure S10','panel':panel,'source':'Pulmonary/disjoint_reference_splits.csv',
                    'dataset':dataset,'metric':metric,'unit':'oriented disease-animal reference split',
                    'n':len(ds),'summary_source':'Pulmonary/shared_reference_summary.csv'})
  for x in ds.itertuples():datarows.append({'figure':'Figure S10','panel':panel,'dataset':dataset,'metric':metric,
      'record_type':'disjoint_reference_split','record_id':x.split,'value':getattr(x,metric)})
  datarows.append({'figure':'Figure S10','panel':panel,'dataset':dataset,'metric':metric,
      'record_type':'shared_reference_observed','record_id':-1,'value':obs[metric]})
for row,(dataset,animal) in enumerate(datasets,start=2):
 obs=summary.loc[dataset];dn=null[null.dataset.eq(dataset)]
 for col,(metric,xlabel) in enumerate(metrics):
  a=axes[row,col];panel=string.ascii_uppercase[row*3+col];v=dn[metric].to_numpy()
  a.hist(v,bins=14 if len(v)<100 else 30,color='#a7b5c3',edgecolor='white',linewidth=.5)
  a.axvline(obs[metric],color='#c27236',lw=1.4,label='Observed allocation')
  a.axvline(np.median(v),color='#617383',lw=1,ls=':',label='Null median')
  p=obs[metric+'_animal_permutation_p'];q=obs[metric+'_animal_permutation_q']
  a.text(.98,.97,f'Observed {display(obs[metric],metric)}\nP = {p:.3f}; q = {q:.3f}',
         transform=a.transAxes,ha='right',va='top',fontsize=7.7,
         bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':2})
  a.margins(y=.24)
  if metric=='reversal_fraction':a.xaxis.set_major_formatter(PercentFormatter(1,decimals=0))
  a.set(xlabel=xlabel,ylabel='Whole-animal allocations' if col==0 else '',
        title=f'{animal}: {len(v):,} animal-label allocations')
  if row==2 and col==0:a.legend(frameon=False,fontsize=6.8,loc='upper left')
  label(a,panel)
  panelrows.append({'figure':'Figure S10','panel':panel,'source':'Pulmonary/animal_label_permutations.csv',
                    'dataset':dataset,'metric':metric,'unit':'whole-animal treatment allocation',
                    'n':len(dn),'summary_source':'Pulmonary/shared_reference_summary.csv'})
  for x in dn.itertuples():datarows.append({'figure':'Figure S10','panel':panel,'dataset':dataset,'metric':metric,
      'record_type':'animal_label_allocation','record_id':x.allocation,'value':getattr(x,metric)})
  datarows.append({'figure':'Figure S10','panel':panel,'dataset':dataset,'metric':metric,
      'record_type':'shared_reference_observed','record_id':-1,'value':obs[metric]})
stem='Supplementary_Figure_S10_Shared_Reference_Sensitivity'
for ext in ['svg','pdf']:f.savefig(O/(stem+'.'+ext),bbox_inches='tight',pad_inches=.1)
f.savefig(O/(stem+'.png'),dpi=145,bbox_inches='tight',pad_inches=.1)
plt.close(f)
tmp=R/'figure_qa'/f'{stem}_600dpi.png'
subprocess.run(['inkscape',str(O/(stem+'.svg')),'--export-type=png','--export-area-page','--export-dpi=600',
                '--export-background=white','--export-background-opacity=1',f'--export-filename={tmp}'],check=True,capture_output=True)
subprocess.run(['convert',str(tmp),'-background','white','-alpha','remove','-alpha','off','-units','PixelsPerInch',
                '-density','600','-compress','LZW',str(O/(stem+'.tiff'))],check=True,capture_output=True)
with Image.open(O/(stem+'.tiff')) as im:im.load();assert tuple(map(float,im.info['dpi']))==(600.,600.)
tmp.unlink()
pd.DataFrame(panelrows).to_csv(O/'S10_panel_source_map.csv',index=False)
pd.DataFrame(datarows).to_csv(O/'S10_panel_source_data.csv',index=False)
(R/'figure_qa/S10_source_checksums.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [S/'shared_reference_summary.csv',S/'disjoint_reference_splits.csv',S/'animal_label_permutations.csv']},indent=2))
print(f'Created Figure S10: 12 panels, {len(datarows)} source-value rows',flush=True)
