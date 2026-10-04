from pathlib import Path
import io, os, json, hashlib
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parent/'.mplconfig'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE/'fixed_panels'
OUT.mkdir(exist_ok=True)
TEAL='#177E89'; CORAL='#D66A52'; NAVY='#244F70'; GRAY='#6A7880'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
    'axes.titlesize':9,'axes.titleweight':'bold','xtick.labelsize':7.5,
    'ytick.labelsize':7.5,'axes.spines.top':False,'axes.spines.right':False,
    'axes.linewidth':0.65,'xtick.major.width':0.65,'ytick.major.width':0.65,
    'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})

LIB = pd.read_csv(ROOT/'candidate_library/Candidate_Library_Evidence.csv')
ORDER=LIB.compound.tolist()
KD = pd.read_csv(ROOT/'candidate_library/Davis2011_Matched_Kd_Benchmark.csv')
PROP=pd.read_csv(ROOT/'screening/Inductive_Properties_All16.csv')
ENDPOINT=pd.read_csv(ROOT/'candidate_library/Endpoint_Separated_Experimental_Summary.csv')

def atomic(path,data):
    tmp=path.with_name(path.name+'.partial')
    with tmp.open('wb') as fh:
        fh.write(data);fh.flush();os.fsync(fh.fileno())
    os.replace(tmp,path)

def export(fig,stem):
    # Native plot at 600 dpi; no screenshot resizing.
    buf=io.BytesIO();fig.savefig(buf,format='png',dpi=600)
    atomic(OUT/(stem+'.png'),buf.getvalue())
    with Image.open(io.BytesIO(buf.getvalue())) as im:
        rgb=im.convert('RGB');t=io.BytesIO();rgb.save(t,format='TIFF',compression='tiff_lzw',dpi=(600,600))
        atomic(OUT/(stem+'.tiff'),t.getvalue())
    svg=io.BytesIO();fig.savefig(svg,format='svg');atomic(OUT/(stem+'.svg'),svg.getvalue())

def style(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='x',color='#E6EBEE',lw=.55)
    ax.tick_params(axis='y',length=0)
    ax.spines['left'].set_visible(False)

def matched_kd(ax,footnote=True):
    order=[x for x in ORDER if x in set(KD.compound)]
    assert len(KD)==12 and len(order)==6
    ax.axhspan(-.42,.42,color='#EAF1F6',zorder=0)
    for y,c in enumerate(order):
        rows=KD[KD.compound==c]
        assert len(rows)==2
        vals=rows.set_index('gene_symbol').standard_value
        ax.plot([vals['PDGFRA'],vals['FLT4']],[y-.115,y+.115],lw=.75,color='#BDC9CD',zorder=1)
        for _,r in rows.iterrows():
            target=r.gene_symbol;color=TEAL if target=='PDGFRA' else CORAL
            marker='o' if target=='PDGFRA' else 's';yy=y+(-.115 if target=='PDGFRA' else .115)
            v=float(r.standard_value);flag=bool(r.potential_duplicate)
            ax.scatter(v,yy,s=23,marker=marker,edgecolors=color,facecolors='white' if flag else color,linewidth=.85,zorder=3)
            if r.standard_relation=='>':
                ax.annotate('',xy=(21000,yy),xytext=(10500,yy),arrowprops={'arrowstyle':'->','lw':.9,'color':color})
                label='>10,000'
            else: label=f'{v:g}'
            # Offset the two target labels vertically, retaining numerical values.
            ax.annotate(label,(v,yy),xytext=(0,-5 if target=='FLT4' else 3),textcoords='offset points',
                fontsize=6.8,color=color,ha='center',va='top' if target=='FLT4' else 'bottom')
    ax.set_xscale('log');ax.set_xlim(.2,38000);ax.set_xticks([1,10,100,1000,10000],['1','10','100','1,000','10,000'])
    ax.set_yticks(range(len(order)),[x.capitalize() if x!='MAZ51' else x for x in order]);ax.set_ylim(len(order)+.00,-.68)
    ax.set_xlabel(r'Experimental $K_d$ (nM)');style(ax)
    handles=[Line2D([],[],marker='o',ls='',mfc=TEAL,mec=TEAL,label='PDGFRA',ms=4),
             Line2D([],[],marker='s',ls='',mfc=CORAL,mec=CORAL,label='FLT4',ms=4)]
    ax.legend(handles=handles,loc='lower left',bbox_to_anchor=(0,1.00),ncol=2,frameon=False,fontsize=7.3,
              borderaxespad=0,handletextpad=.4,columnspacing=1.5)
    if footnote:
        ax.text(0,-.24,'Davis et al., 2011\nOpen symbols: duplicate-flagged records',transform=ax.transAxes,fontsize=6.5,ha='left',va='top')
    return order

def property_profile(ax,model,title,showlabels=True):
    d=PROP[PROP.model_id==model].set_index('compound').loc[ORDER]
    assert len(d)==16 and d.status.eq('success').all() and d.identity_match.all()
    ax.axhspan(-.45,.45,color='#EAF1F6',zorder=0)
    for i,(c,r) in enumerate(d.iterrows()):
        col=NAVY if c=='nintedanib' else GRAY
        ood=bool(r.out_of_domain_flag);low=bool(r.low_confidence_flag)
        marker='D' if low else 'o'
        ax.plot([r.model_lower_bound,r.model_upper_bound],[i,i],color=col,alpha=.60,lw=1.05,zorder=1)
        ax.scatter(r.value,i,s=17 if c=='nintedanib' else 12,marker=marker,
                   facecolor='white' if ood else col,edgecolor=col,lw=.7,zorder=2)
    if showlabels:ax.set_yticks(range(16),[x.capitalize() if x!='MAZ51' else x for x in ORDER])
    else:ax.set_yticks(range(16),[])
    ax.set_ylim(15.6,-.7);ax.set_title(title,pad=9);style(ax)
    if model=='mcp_public_logd':ax.set_xlim(0,5.55);ax.set_xticks([0,1,2,3,4,5]);ax.set_xlabel('Predicted logD (pH 7.4)')
    elif model=='mcp_public_apka':ax.set_xlim(5.5,12.45);ax.set_xticks([6,8,10,12]);ax.set_xlabel('Predicted acidic pKa')
    else:ax.set_xlim(1.5,11);ax.set_xticks([2,4,6,8,10]);ax.set_xlabel('Predicted basic pKa')

def endpoint_coverage(ax):
    cols=[(t,e) for t in ['PDGFRA','FLT4'] for e in ['Kd','Ki','IC50']]
    p=ENDPOINT.pivot(index='compound',columns=['gene_symbol','endpoint'],values='biochemical_records').loc[ORDER,cols]
    values=p.to_numpy(float);masked=np.ma.masked_equal(values,0)
    im=ax.imshow(masked,aspect='auto',cmap='Blues',vmin=0,vmax=max(values.max(),1))
    ax.set_facecolor('#F3F5F6')
    for i in range(16):
        for j in range(6):
            v=int(values[i,j]);ax.text(j,i,str(v) if v else '–',ha='center',va='center',fontsize=7,
                color='white' if v>=4 else '#36454E')
    ax.set_yticks(range(16),[x.capitalize() if x!='MAZ51' else x for x in ORDER])
    ax.set_xticks(range(6),['$K_d$','$K_i$','IC$_{50}$']*2)
    ax.tick_params(length=0);ax.spines[:].set_visible(False)
    ax.axvline(2.5,color='white',lw=3)
    ax.text(1, -1.25,'PDGFRA',ha='center',va='center',fontsize=8.5,fontweight='bold',color=TEAL)
    ax.text(4, -1.25,'FLT4',ha='center',va='center',fontsize=8.5,fontweight='bold',color=CORAL)
    ax.set_ylim(15.5,-.5)
    return im

if __name__=='__main__':
    fig,ax=plt.subplots(figsize=(3.45,2.95));fig.subplots_adjust(left=.27,right=.965,top=.88,bottom=.22)
    matched_kd(ax);export(fig,'Matched_Davis_Kd');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(7.2,4.15));fig.subplots_adjust(left=.18,right=.985,top=.89,bottom=.21,wspace=.18)
    for ax,mod,title in zip(axes,['mcp_public_logd','mcp_public_apka','mcp_public_bpka'],['Lipophilicity','Acidity','Basicity']):
        property_profile(ax,mod,title,showlabels=ax is axes[0])
    for i,ax in enumerate(axes):ax.text(-.10 if i else -.58,1.08,'ABC'[i],transform=ax.transAxes,fontsize=12,fontweight='bold')
    handles=[Line2D([],[],marker='o',ls='',ms=4,color=GRAY,label='Within domain'),
             Line2D([],[],marker='o',ls='',ms=4,mfc='white',mec=GRAY,label='Open: out of domain'),
             Line2D([],[],marker='D',ls='',ms=4,color=GRAY,label='Low confidence')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.56,.055),ncol=3,frameon=False,fontsize=7)
    fig.text(.18,.019,'Points: model estimates; segments: provider model bounds, coverage unspecified (not confidence intervals).',fontsize=6.7)
    export(fig,'All16_Inductive_Properties');plt.close(fig)
    fig,ax=plt.subplots(figsize=(3.35,3.75));fig.subplots_adjust(left=.31,right=.97,top=.90,bottom=.16)
    endpoint_coverage(ax);fig.text(.31,.065,'Eligible biochemical records',fontsize=8)
    fig.text(.31,.02,'–: no eligible record; not evidence of inactivity',fontsize=6.5)
    export(fig,'Endpoint_Evidence_Coverage');plt.close(fig)
    metadata={'matched_kd_records':len(KD),'matched_kd_compounds':KD.compound.nunique(),
              'potential_duplicate_records':int(KD.potential_duplicate.sum()),'property_rows':len(PROP),
              'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
                  ROOT/'candidate_library/Davis2011_Matched_Kd_Benchmark.csv',
                  ROOT/'candidate_library/Endpoint_Separated_Experimental_Summary.csv',
                  ROOT/'screening/Inductive_Properties_All16.csv']}}
    (HERE/'Fixed_Evidence_Source_Manifest.json').write_text(json.dumps(metadata,indent=2))
    print(json.dumps(metadata,indent=2))
