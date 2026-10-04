"""Publication Figure 7; requires a complete input-disposition BoltzMol table."""
from pathlib import Path
import os, argparse, io, json, hashlib
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'.mplconfig'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from plot_fixed_evidence import ORDER, ROOT, HERE, TEAL, CORAL, NAVY, GRAY, style, matched_kd, atomic

OUT=HERE/'final';OUT.mkdir(exist_ok=True)

def get_data(path):
    d=pd.read_csv(path)
    if 'compound' not in d and 'name' in d:d=d.rename(columns={'name':'compound'})
    assert {'compound','target','status','binding_confidence'}<=set(d.columns),d.columns.tolist()
    expected={(c,t) for c in ORDER for t in ['PDGFRA','FLT4']}
    keys=list(zip(d.compound,d.target))
    assert len(keys)==32 and len(set(keys))==32 and set(keys)==expected
    d['binding_confidence']=pd.to_numeric(d.binding_confidence,errors='coerce')
    scored=d.status.str.lower().eq('scored')
    filtered=d.status.str.lower().eq('filtered')
    assert (scored|filtered).all(),'Unresolved model outputs cannot be finalized'
    assert d.loc[scored,'binding_confidence'].notna().all()
    assert d.loc[filtered,'binding_confidence'].isna().all(),'Filtered inputs must not be assigned a numerical score'
    assert d.loc[scored,'binding_confidence'].between(0,1).all()
    return d

def binding_profile(ax,d,target):
    t=d[d.target==target].copy()
    score=t[t.status.eq('scored')].sort_values(['binding_confidence','compound'],ascending=[False,True])
    filtered=t[t.status.eq('filtered')].set_index('compound').loc[[c for c in ORDER if c in set(t[t.status.eq('filtered')].compound)]].reset_index()
    order=score.compound.tolist()+filtered.compound.tolist()
    idx=order.index('nintedanib');ax.axhspan(idx-.44,idx+.44,color='#EAF1F6',zorder=0)
    color=TEAL if target=='PDGFRA' else CORAL
    for y,r in enumerate(score.itertuples()):
        ax.scatter(r.binding_confidence,y,s=18,color=color,zorder=3)
    for y in range(len(score),len(order)):
        ax.text(.52,y,'Filtered',transform=ax.get_yaxis_transform(),ha='center',va='center',fontsize=7,color='#6D7A83')
    if len(filtered):ax.axhline(len(score)-.5,color='#BFCBD1',lw=.65,ls=(0,(2,3)))
    ax.set_yticks(range(16),[x.capitalize() if x!='MAZ51' else x for x in order])
    ax.set_ylim(15.6,-.7);ax.set_xlim(0,1.03);ax.set_xticks([0,.5,1],['0','0.5','1.0'])
    ax.set_xlabel('Binding confidence',fontsize=8);style(ax)
    ax.tick_params(axis='y',labelsize=7)
    return {'target':target,'scored':len(score),'filtered':len(filtered),'plot_order':order}

def molecular_image(fig,path,box):
    ax=fig.add_axes(box);im=Image.open(path);assert im.size==(1350,1440),im.size
    ax.imshow(im);ax.set_axis_off()

def build(path):
    csv_path=Path(path)
    d=get_data(csv_path)
    assets={
        'D':HERE/'molecular/Panel_E_redocking.png',
        'E':HERE/'molecular/Panel_E_BoltzMol_PDGFRA_imatinib.png',
        'F':HERE/'molecular/Panel_F_BoltzMol_FLT4_tivozanib.png'}
    assert all(p.exists() for p in assets.values()),'Verified molecular panels not yet available'
    fig=plt.figure(figsize=(7.2,6.75))
    # Each column is ~2.3 in, matching molecular-panel native 600-dpi widths.
    cols=[(.012,.318),(.345,.318),(.676,.312)]
    top_bottom=.51;top_height=.403
    axA=fig.add_axes([cols[0][0]+.096,top_bottom,cols[0][1]-.111,top_height])
    matched_kd(axA,footnote=False)
    axA.tick_params(axis='y',labelsize=7.2);axA.tick_params(axis='x',labelsize=6.8)
    axA.set_xlabel(r'Experimental $K_d$ (nM)',fontsize=8)
    axA.set_xticks([1,10,100,1000,10000],['1','10','100',r'$10^3$',r'$10^4$'])
    # The same-study annotation remains visible but compact; detailed flags are in legend.
    fig.text(cols[0][0]+.095,.445,'Davis et al., 2011\nOpen: duplicate-flagged records',fontsize=6.3,va='top')
    summaries=[]
    for i,target in enumerate(['PDGFRA','FLT4'],1):
        x,w=cols[i]
        ax=fig.add_axes([x+.101,top_bottom,w-.108,top_height])
        summaries.append(binding_profile(ax,d,target))
        fig.text(x+.055,.42,f'{summaries[-1]["scored"]} scored; {summaries[-1]["filtered"]} filtered',fontsize=6.8,color=GRAY)
    titles=['Same-study binding','PDGFRA: BoltzMol prediction','FLT4: BoltzMol prediction']
    for i,(x,w) in enumerate(cols):
        fig.text(x,.969,'ABC'[i],fontsize=12,fontweight='bold',va='top')
        fig.text(x+.025,.962,titles[i],fontsize=8.2,fontweight='bold',va='top')
    for i,(letter,path) in enumerate(assets.items()):
        x,w=cols[i]
        fig.text(x,.395,letter,fontsize=12,fontweight='bold',va='top')
        # Native model illustrations occupy ~2.25 in; no invented geometry.
        molecular_image(fig,path,[x+.013,.014,.3125,.35556])
    for i,title in enumerate(['Imatinib redocking control','Imatinib–PDGFRA predicted pocket','Tivozanib–FLT4 predicted pocket']):
        x,w=cols[i];fig.text(x+.025,.392,title,fontsize=7.65,fontweight='bold',va='top')
    # Final exports: raster 600 dpi, vector plot text in PDF plus native image panels.
    b=io.BytesIO();fig.savefig(b,format='png',dpi=600);atomic(OUT/'Figure_7_Compound_Comparison.png',b.getvalue())
    im=Image.open(io.BytesIO(b.getvalue())).convert('RGB');tb=io.BytesIO();im.save(tb,format='TIFF',compression='tiff_lzw',dpi=(600,600));atomic(OUT/'Figure_7_Compound_Comparison.tiff',tb.getvalue())
    pb=io.BytesIO();fig.savefig(pb,format='pdf',dpi=600);atomic(OUT/'Figure_7_Compound_Comparison.pdf',pb.getvalue())
    plt.close(fig)
    meta={'model_profiles':summaries,'score_field':'binding_confidence','score_interpretation':'dimensionless provider score, not a calibrated probability or measured affinity',
          'source_sha256':hashlib.sha256(csv_path.read_bytes()).hexdigest(),'molecular_panel_sha256':{k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in assets.items()}}
    (HERE/'Main_Figure_Source_Manifest.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('boltz_summary');args=a.parse_args();build(args.boltz_summary)
