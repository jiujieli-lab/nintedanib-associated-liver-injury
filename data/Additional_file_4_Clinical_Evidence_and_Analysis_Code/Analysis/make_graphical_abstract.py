from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle, FancyArrowPatch, Polygon
from PIL import Image
O=Path(__file__).resolve().parent/'figures_new';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
f,ax=plt.subplots(figsize=(13.2,7.1));ax.set(xlim=(0,13.2),ylim=(0,7.1));ax.axis('off');f.subplots_adjust(0,0,1,1)
navy='#18374A';blue='#3886A3';red='#C86D62';gold='#C49C42';gray='#5D6973'
ax.text(.35,6.68,'Clinical proteins guide the test of a hepatic injury mechanism',fontsize=18,fontweight='bold',color=navy)
ax.text(.35,6.25,'Nintedanib for fibrosing interstitial lung disease',fontsize=12,color=gray)
# pulmonary setting
for x,sgn in [(1.0,-1),(1.7,1)]:
 ax.add_patch(Ellipse((x,4.55),.8,1.65,angle=-sgn*10,facecolor='#DCECF0',edgecolor=blue,lw=1.5))
 ax.plot([1.35,x,x-.08*sgn],[5.58,4.8,4.1],color=blue,lw=2)
ax.text(1.35,3.13,'Preserve pulmonary\nantifibrotic activity',ha='center',fontsize=11,color=navy)
ax.add_patch(FancyBboxPatch((.5,2.2),1.7,.55,boxstyle='round,pad=.09,rounding_size=.22',fc=navy,ec=navy));ax.text(1.35,2.47,'NINTEDANIB',color='white',ha='center',va='center',weight='bold')
# liver microenvironment background
ax.add_patch(FancyBboxPatch((3,2.0),6.2,3.9,boxstyle='round,pad=.12',fc='#F4F7F8',ec='#D5DFE3'))
ax.text(6.1,5.64,'Human liver target-to-phenotype model',ha='center',fontsize=12,weight='bold',color=navy)
# nonparenchymal cell
ax.add_patch(Ellipse((4.3,4.4),1.65,.65,fc='#DDECF2',ec=blue,lw=1.5));ax.add_patch(Ellipse((4.3,4.4),.35,.24,fc=blue,alpha=.5,ec='none'))
for x in [3.95,4.25,4.55]:
 ax.plot([x,x],[4.65,4.92],color=navy,lw=2);ax.plot([x-.08,x+.08],[4.94,4.94],color=navy,lw=2)
ax.text(4.3,3.76,'PDGFRA · FLT4',ha='center',weight='bold',color=navy,fontsize=12)
ax.text(4.3,3.35,'Stable experimental priorities\nCell-specific expression to verify',ha='center',fontsize=9.6,color=gray)
# hepatocyte hexagon
import numpy as np
angles=np.arange(6)*np.pi/3
pts=np.column_stack([7.5+.82*np.cos(angles),4.35+.75*np.sin(angles)])
ax.add_patch(Polygon(pts,closed=True,fc='#F4DDD6',ec=red,lw=1.6));ax.add_patch(Circle((7.5,4.35),.23,fc=red,alpha=.55,ec='none'))
ax.text(7.5,3.31,'Hepatocyte metabolism\nand injury-protein release',ha='center',fontsize=10,color=navy)
# inferred unresolved link
ax.add_patch(FancyArrowPatch((5.18,4.4),(6.58,4.4),arrowstyle='->',mutation_scale=15,lw=1.6,ls='--',color=gold))
ax.text(5.89,4.92,'Intercellular\nmechanism?',ha='center',fontsize=10,color='#927123')
ax.add_patch(FancyArrowPatch((2.31,2.48),(3.8,4.04),arrowstyle='->',mutation_scale=14,lw=1.6,color=navy,connectionstyle='arc3,rad=-.15'))
ax.text(6.1,2.49,'Protein-weighted liver networks + numerical calibration',ha='center',fontsize=10,color=navy)
# clinical endpoint
for k,(lab,col) in enumerate([('FBP1',red),('GSTA1',blue),('OTC',gold)]):
 x=10.0+.9*k
 ax.add_patch(Circle((x,4.4),.34,fc=col,alpha=.15,ec=col,lw=1.8));ax.text(x,4.4,lab,ha='center',va='center',fontsize=9,color=navy,weight='bold')
ax.add_patch(FancyArrowPatch((8.39,4.35),(9.54,4.35),arrowstyle='->',mutation_scale=15,lw=1.6,ls='--',color=gold))
ax.text(10.9,3.44,'Measured human serum\ninjury phenotype',ha='center',fontsize=11,color=navy)
ax.text(10.9,2.6,'13 clinical proteins\nMulti-drug DILI cohort',ha='center',fontsize=9.7,color=gray)
# lower integrative test across organ
ax.add_patch(FancyBboxPatch((.45,.53),12.25,.85,boxstyle='round,pad=.08',fc='#FFF9E9',ec='#DFC887'))
ax.text(6.6,1.07,'CAUSAL TEST: target engagement → exposure-matched perturbation → reciprocal rescue',ha='center',fontsize=11,color=navy,weight='bold')
ax.text(6.6,.72,'Require lower hepatic injury and retained pulmonary response. No DILI-causal target is confirmed by the public-data analysis.',ha='center',fontsize=9.3,color=gray)
ax.text(.45,.15,'Dashed arrows: proposed biological links, not demonstrated causal effects.',fontsize=8.5,color=gray)
name='Graphical_Abstract_Clinical_Protein_Guided_Target_Prioritization'
f.savefig(O/(name+'.pdf'));f.savefig(O/(name+'.svg'));f.savefig(O/(name+'.png'),dpi=160)
f.savefig(O/(name+'.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
plt.close(f)
(O/'Graphical_Abstract_caption.md').write_text('Graphical abstract. Clinical proteins guide the test of a hepatic injury mechanism. Nintedanib pharmacology, the human serum injury phenotype, and donor-balanced liver networks form one candidate-to-phenotype model. PDGFRA and FLT4 are stable relative experimental priorities, while the responsible intercellular signal and intervention direction remain unestablished. Dashed arrows denote hypotheses. The proposed confirmation experiment must demonstrate exposure-matched hepatic rescue and preserve pulmonary antifibrotic activity. This schematic is separate from the six data-only main figures and presents no completed perturbation-and-rescue experiment.\n')
print(name)
