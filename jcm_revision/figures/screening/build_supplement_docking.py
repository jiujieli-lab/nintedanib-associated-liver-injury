"""Complete-campaign-only supplemental charts. Never render partial results."""
from pathlib import Path
import os, io, json, hashlib
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'.mplconfig'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from plot_fixed_evidence import ORDER, ROOT, HERE, TEAL, CORAL, NAVY, GRAY, style, endpoint_coverage, atomic

OUT=HERE/'final';OUT.mkdir(exist_ok=True)
SEEDS=[104729,130363,155921]

def complete_data():
    d=pd.read_csv(ROOT/'docking/Docking_All_Completed_Runs.csv')
    s=pd.read_csv(ROOT/'docking/Docking_Summary.csv')
    expected={(c,t,k) for c in ORDER for t in ['PDGFRA','FLT4'] for k in SEEDS}
    got=list(zip(d.name,d.target,d.seed.astype(int)))
    assert len(got)==len(set(got)), 'Duplicate docking identity keys'
    assert set(got)==expected, f'Campaign incomplete: {len(set(got)&expected)}/96 expected searches'
    assert len(s)==32 and s.completed_searches.eq(3).all(), 'Summary incomplete'
    assert d.vina_score_kcal_mol.notna().all()
    for _,r in s.iterrows():
        vals=d.loc[(d.name==r['name'])&(d.target==r.target),'vina_score_kcal_mol'].values
        assert np.isclose(np.median(vals),r.median_vina_score_kcal_mol)
        assert np.isclose(min(vals),r.min_vina_score_kcal_mol)
        assert np.isclose(max(vals),r.max_vina_score_kcal_mol)
    return d,s

def seed_heatmap(ax,d,target):
    p=d[d.target==target].pivot(index='name',columns='seed',values='vina_score_kcal_mol').loc[ORDER,SEEDS]
    values=p.to_numpy();color=TEAL if target=='PDGFRA' else CORAL
    cm=LinearSegmentedColormap.from_list(target,[color,'#F8FAFB'])
    vmin=np.floor(values.min());vmax=np.ceil(values.max())
    im=ax.imshow(values,aspect='auto',cmap=cm,vmin=vmin,vmax=vmax)
    for i in range(16):
        for j in range(3):
            v=values[i,j];f=(v-vmin)/(vmax-vmin)
            ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=7,color='white' if target=='PDGFRA' and f<.38 else '#253A43')
    ax.set_yticks(range(16),[x.capitalize() if x!='MAZ51' else x for x in ORDER]);ax.tick_params(length=0)
    ax.set_xticks(range(3),['104729','130363','155921']);ax.set_xlabel('Search seed')
    ax.set_title(f'{target}: three seeded searches',pad=10)
    ax.spines[:].set_visible(False)
    return im

def consistency(ax,s,target):
    d=s[s.target==target].set_index('name').loc[ORDER]
    vals=d.max_pairwise_top_pose_RMSD_A
    assert vals.notna().all()
    col=TEAL if target=='PDGFRA' else CORAL
    ax.axhspan(-.45,.45,color='#EAF1F6',zorder=0)
    ax.scatter(vals,range(16),s=17,color=col,zorder=3)
    ax.set_yticks(range(16),[x.capitalize() if x!='MAZ51' else x for x in ORDER])
    ax.set_ylim(15.6,-.7);ax.set_xlim(0,max(2,float(vals.max())*1.10))
    ax.set_xlabel('Maximum top-pose RMSD (Å)');ax.set_title(f'{target}: pose consistency',pad=10);style(ax)

def pose_sensitivity(ax):
    qc=pd.read_csv(ROOT/'docking/Nintedanib_Pose_QC_All_Searches.csv')
    qc=qc[qc.target=='PDGFRA'].copy()
    assert len(qc)==30 and set(qc.seed)==set(SEEDS)
    for seed,color,marker in zip(SEEDS,[TEAL,CORAL,NAVY],['o','s','^']):
        q=qc[qc.seed==seed]
        assert sorted(q.pose_rank.tolist())==list(range(1,11))
        ax.scatter(q.vina_score_kcal_mol,q.minimum_gatekeeper_distance_A,s=15,color=color,marker=marker,alpha=.78,label=str(seed),zorder=2)
        top=q[q.pose_rank==1].iloc[0]
        ax.scatter(top.vina_score_kcal_mol,top.minimum_gatekeeper_distance_A,s=39,facecolor=color,edgecolor='#142D39',lw=.8,marker=marker,zorder=4)
    ax.axhline(4,color='#BDC7CC',lw=.8,ls='--')
    ax.text(.97,4.15,'4 Å proximity',transform=ax.get_yaxis_transform(),ha='right',fontsize=6.5,color=GRAY)
    ax.set_xlabel('Vina score (kcal/mol)');ax.set_ylabel('Minimum Thr674 distance (Å)')
    ax.set_title('PDGFRA nintedanib poses',pad=30)
    ax.legend(frameon=False,fontsize=5.8,loc='lower left',bbox_to_anchor=(0,1.01),ncol=3,handletextpad=.15,columnspacing=.4,borderaxespad=0)
    style(ax)

def export_composite(fig,stem):
    buf=io.BytesIO();fig.savefig(buf,format='png',dpi=600)
    atomic(OUT/(stem+'.png'),buf.getvalue())
    with Image.open(io.BytesIO(buf.getvalue())) as src:
        rgb=src.convert('RGB');b=io.BytesIO();rgb.save(b,format='TIFF',compression='tiff_lzw',dpi=(600,600));atomic(OUT/(stem+'.tiff'),b.getvalue())
    b=io.BytesIO();fig.savefig(b,format='pdf',dpi=600);atomic(OUT/(stem+'.pdf'),b.getvalue())

def build():
    d,s=complete_data()
    fig=plt.figure(figsize=(7.2,10.15))
    grid=fig.add_gridspec(3,2,left=.17,right=.91,top=.956,bottom=.11,hspace=.47,wspace=.68,
                         height_ratios=[1.08,1.05,1])
    axs=[fig.add_subplot(grid[i,j]) for i in range(3) for j in range(2)]
    for i,target in enumerate(['PDGFRA','FLT4']):
        im=seed_heatmap(axs[i],d,target)
        cb=fig.colorbar(im,ax=axs[i],fraction=.04,pad=.04)
        cb.set_label('Vina score (kcal/mol)',fontsize=7);cb.ax.tick_params(labelsize=6.5,length=2)
        consistency(axs[i+2],s,target)
    endpoint_coverage(axs[4]);axs[4].set_title('Biochemical evidence coverage',pad=30)
    pose_sensitivity(axs[5])
    for i,ax in enumerate(axs):ax.text(-.43 if i in [0,1,2,3,4] else -.25,1.21 if i in [4,5] else 1.075,'ABCDEF'[i],transform=ax.transAxes,fontsize=12,fontweight='bold')
    fig.text(.17,.036,'Docking scores are interpreted within each receptor; seeds quantify search consistency, not biological replication.',fontsize=7)
    fig.text(.17,.016,'PDGFRA: experimental 6JOL; FLT4: AlphaFold prediction. E: eligible record counts; –: unavailable.',fontsize=6.65)
    export_composite(fig,'Figure_S13_Docking_and_Evidence');plt.close(fig)
    files=[ROOT/'docking/Docking_All_Completed_Runs.csv',ROOT/'docking/Docking_Summary.csv',ROOT/'docking/PDGFRA_nintedanib_pose_hinge_QC.json',ROOT/'docking/Nintedanib_Pose_QC_All_Searches.csv']
    meta={'run_count':len(d),'compound_count':d.name.nunique(),'targets':sorted(d.target.unique()),'seeds':SEEDS,
          'all_groups_have_three_searches':bool(s.completed_searches.eq(3).all()),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    (HERE/'Supplement_Docking_Source_Manifest.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2))

if __name__=='__main__':build()
