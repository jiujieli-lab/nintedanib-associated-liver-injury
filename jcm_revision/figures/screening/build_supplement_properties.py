"""Properties, explicit nintedanib poses, and separately completed receptor sensitivity."""
from pathlib import Path
import os,io,json,hashlib
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'.mplconfig'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
from plot_fixed_evidence import HERE,ROOT,TEAL,CORAL,NAVY,GRAY,property_profile,style,atomic

OUT=HERE/'final';OUT.mkdir(exist_ok=True)
SEEDS=[104729,130363,155921]

def build():
    source=ROOT/'docking/PDGFRA_Nintedanib_Receptor_Sensitivity.csv'
    d=pd.read_csv(source)
    assert len(d)==6 and d.receptor.nunique()==2, 'Separate receptor sensitivity is incomplete'
    primary='Primary_6JOL';other=next(x for x in d.receptor.unique() if x!=primary)
    assert '8PQJ' in other
    for receptor in [primary,other]:assert sorted(d[d.receptor==receptor].seed.tolist())==SEEDS
    fig=plt.figure(figsize=(7.2,7.35))
    axes=[]
    for x,w in [(.18,.24),(.46,.235),(.74,.235)]:axes.append(fig.add_axes([x,.55,w,.38]))
    for i,(ax,model,title) in enumerate(zip(axes,['mcp_public_logd','mcp_public_apka','mcp_public_bpka'],['Lipophilicity','Acidity','Basicity'])):
        property_profile(ax,model,title,showlabels=i==0)
        fig.text([.025,.43,.71][i],.977,'ABC'[i],fontsize=12,fontweight='bold',va='top')
    handles=[Line2D([],[],marker='o',ls='',ms=4,color=GRAY,label='Within domain'),
             Line2D([],[],marker='o',ls='',ms=4,mfc='white',mec=GRAY,label='Open: out of domain'),
             Line2D([],[],marker='D',ls='',ms=4,color=GRAY,label='Low confidence')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.56,.477),ncol=3,frameon=False,fontsize=7)
    fig.text(.18,.463,'Segments: provider model bounds; coverage unspecified (not confidence intervals).',fontsize=6.8)
    for i,name in enumerate(['Panel_F_PDGFRA_pose_ambiguity.png','Panel_G_FLT4_predicted_pocket.png']):
        x=.012+i*.333
        ax=fig.add_axes([x+.008,.035,.3125,2.4/7.35])
        im=Image.open(HERE/'molecular'/name);assert im.size==(1350,1440)
        ax.imshow(im);ax.set_axis_off()
    titles=['Nintedanib: PDGFRA poses','Nintedanib: FLT4 predicted pose','PDGFRA receptor sensitivity']
    for i,title in enumerate(titles):
        x=.012+i*.333
        fig.text(x,.409,'DEF'[i],fontsize=12,fontweight='bold',va='top')
        fig.text(x+.025,.404,title,fontsize=7.6,fontweight='bold',va='top')
    seedhandles=[Line2D([],[],marker=m,ls='',ms=3,color=c,label=str(k)) for k,c,m in zip(SEEDS,[TEAL,CORAL,NAVY],['o','s','^'])]
    fig.legend(handles=seedhandles,loc='lower center',bbox_to_anchor=(.876,.366),ncol=3,frameon=False,fontsize=5.2,handletextpad=.1,columnspacing=.4)
    axscore=fig.add_axes([.775,.244,.205,.117])
    axdist=fig.add_axes([.775,.076,.205,.117])
    for j,(seed,color,marker) in enumerate(zip(SEEDS,[TEAL,CORAL,NAVY],['o','s','^'])):
        for i,receptor in enumerate([primary,other]):
            r=d[(d.receptor==receptor)&(d.seed==seed)].iloc[0]
            axscore.scatter(i+(j-1)*.12,r.vina_score_kcal_mol,color=color,marker=marker,s=20,zorder=3)
            axdist.scatter(i+(j-1)*.12,r.minimum_Cys677_distance_A,color=color,marker=marker,s=20,zorder=3)
    for ax in [axscore,axdist]:
        ax.set_xlim(-.4,1.4);ax.set_xticks([0,1],['6JOL','8PQJ']);ax.tick_params(labelsize=7)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        ax.grid(axis='y',color='#E6EBEE',lw=.55,zorder=0)
    axscore.set_ylabel('Vina score\n(kcal/mol)',fontsize=7);axscore.tick_params(axis='x',labelbottom=False)
    axdist.set_ylabel('Minimum Cys677\ndistance (Å)',fontsize=7)
    axdist.set_xlabel('PDGFRA receptor structure',fontsize=7)
    fig.text(.775,.221,'Three search seeds per structure',fontsize=6.4,color=GRAY)
    fig.text(.775,.013,'8PQJ: separate apo-state sensitivity',fontsize=6.2,color=GRAY)
    stem=OUT/'Figure_S14_Properties_and_Receptor_Sensitivity'
    b=io.BytesIO();fig.savefig(b,format='png',dpi=600);atomic(stem.with_suffix('.png'),b.getvalue())
    im=Image.open(io.BytesIO(b.getvalue())).convert('RGB');t=io.BytesIO();im.save(t,format='TIFF',compression='tiff_lzw',dpi=(600,600));atomic(stem.with_suffix('.tiff'),t.getvalue())
    b=io.BytesIO();fig.savefig(b,format='pdf',dpi=600);atomic(stem.with_suffix('.pdf'),b.getvalue());plt.close(fig)
    meta={'sensitivity_rows':len(d),'receptors':[primary,other],'seeds':SEEDS,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
          'not_part_of_primary_96_rankings':True}
    (HERE/'Properties_Sensitivity_Source_Manifest.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2))

if __name__=='__main__':build()
