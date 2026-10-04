#!/usr/bin/env python3
"""Supplementary visualization of GSE299128 whole-blood sensitivity analysis."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "blood_longitudinal"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","Nimbus Sans","DejaVu Sans"],
    "font.size":6.4, "axes.titlesize":7.3, "axes.labelsize":6.5,
    "xtick.labelsize":5.6, "ytick.labelsize":5.6, "axes.linewidth":0.65,
    "pdf.fonttype":42, "svg.fonttype":"none"
})

meta=pd.read_csv(DATA/"sample_metadata.csv")
expr=pd.read_csv(DATA/"target_log2_cpm.csv").merge(meta,on="sample",how="left",validate="one_to_one")
slopes=pd.read_csv(DATA/"participant_slopes.csv")
summary=pd.read_csv(DATA/"target_longitudinal_summary.csv")

blue="#2D6F9F"; orange="#E76F51"; teal="#2A9D8F"; gray="#5B6573"; light="#D6DDE3"

def clean(ax):
    ax.spines[['top','right']].set_visible(False); ax.tick_params(length=2.2,pad=1.4)
def letter(ax,s):
    ax.text(-.12,1.06,s,transform=ax.transAxes,fontsize=9,fontweight="bold",va="top",ha="left")

fig,axs=plt.subplots(3,2,figsize=(7.2,8.7))
plt.subplots_adjust(left=.10,right=.985,top=.94,bottom=.07,hspace=.62,wspace=.42)
fig.suptitle("Supplementary Figure S5 | Whole-blood longitudinal sensitivity analysis (GSE299128)",
             x=.10,ha="left",fontsize=10,fontweight="bold",color="#173B57")

# A visit availability.
ax=axs[0,0]
parts=sorted(meta.participant.unique()); visits=sorted(meta.visit.unique())
mat=np.zeros((len(parts),len(visits)))
for i,p in enumerate(parts):
    present=set(meta.loc[meta.participant==p,'visit'])
    mat[i,:]=[1 if v in present else np.nan for v in visits]
ax.imshow(mat,aspect='auto',cmap=mpl.colors.ListedColormap([teal]),vmin=0,vmax=1)
ax.set_yticks(range(len(parts)),parts); ax.set_xticks(range(len(visits)),visits)
ax.set_xlabel("Study visit"); ax.set_ylabel("Participant")
ax.set_title("A 57-sample, 7-participant series",loc="left",fontweight="bold"); letter(ax,"A")
for s in ax.spines.values(): s.set_visible(False)

# B mean slopes with confidence intervals.
ax=axs[0,1]
d=summary.sort_values("mean_slope_log2_cpm_per_visit")
y=np.arange(len(d))
for yi,(_,r) in zip(y,d.iterrows()):
    c=orange if r.gene_symbol=="FBP1" else gray
    ax.plot([r.slope_ci_low,r.slope_ci_high],[yi,yi],color=c,lw=.8)
    ax.scatter(r.mean_slope_log2_cpm_per_visit,yi,s=14,facecolor="white",edgecolor=c,lw=.7)
ax.axvline(0,color="#999",lw=.6,ls=(0,(2,2)))
ax.set_yticks(y,d.gene_symbol)
for tick in ax.get_yticklabels():
    if tick.get_text()=="FBP1": tick.set_color(orange); tick.set_fontweight("bold")
ax.set_xlabel("Mean log2 CPM slope per visit (95% CI)")
ax.set_title("No prespecified transcript passed FDR",loc="left",fontweight="bold")
ax.text(.98,.02,"All BH q=1.00",transform=ax.transAxes,ha="right",va="bottom",fontsize=5.5)
clean(ax); letter(ax,"B")

# C participant-specific slopes.
ax=axs[1,0]
piv=slopes.pivot(index="gene_symbol",columns="participant",values="slope_per_visit")
order=summary.sort_values("mean_slope_log2_cpm_per_visit",ascending=False).gene_symbol
piv=piv.loc[order]
lim=np.nanquantile(np.abs(piv.to_numpy()),.97)
im=ax.imshow(piv,aspect='auto',cmap="RdBu_r",vmin=-lim,vmax=lim)
ax.set_yticks(range(len(piv)),piv.index); ax.set_xticks(range(len(piv.columns)),piv.columns,rotation=45,ha="right")
ax.set_title("Participant-level slope heterogeneity",loc="left",fontweight="bold")
cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.02); cb.set_label("log2 CPM / visit")
for s in ax.spines.values(): s.set_visible(False)
letter(ax,"C")

# D mean slope versus last-baseline change.
ax=axs[1,1]
ax.scatter(summary.mean_slope_log2_cpm_per_visit,summary.median_last_minus_baseline_log2_cpm,
           s=20,facecolor="white",edgecolor=gray,lw=.8)
for _,r in summary.iterrows():
    if r.gene_symbol in ["FBP1","KDR","CES1","GSTA1","NFE2L2"]:
        c=orange if r.gene_symbol=="FBP1" else gray
        ax.annotate(r.gene_symbol,(r.mean_slope_log2_cpm_per_visit,r.median_last_minus_baseline_log2_cpm),
                    xytext=(3,2),textcoords="offset points",fontsize=5.2,color=c)
ax.axhline(0,color="#AAA",lw=.6); ax.axvline(0,color="#AAA",lw=.6)
ax.set_xlabel("Mean slope per visit"); ax.set_ylabel("Median last minus baseline (log2 CPM)")
ax.set_title("Two longitudinal summaries give no coherent shift",loc="left",fontweight="bold")
clean(ax); letter(ax,"D")

# E FBP1 participant trajectories.
ax=axs[2,0]
for p,dd in expr.sort_values('visit').groupby('participant'):
    ax.plot(dd.visit,dd.FBP1,marker='o',ms=2.5,lw=.8,alpha=.75,label=p)
ax.set_xlabel("Study visit"); ax.set_ylabel("FBP1 log2 CPM")
ax.set_title("FBP1 trajectories vary between participants",loc="left",fontweight="bold")
ax.legend(ncol=2,frameon=False,loc="best",handlelength=1.2,columnspacing=.8)
clean(ax); letter(ax,"E")

# F group-average standardized trajectories.
ax=axs[2,1]
groups={
    "DILI anchors":["ACO1","ASS1","FAH","CPS1","ALDOB","HPD","OTC","DMGDH","GSTA1","FBP1","PCK2","CES1","LECT2"],
    "Direct targets":["KDR","FLT1","FLT4","FGFR1","FGFR2","FGFR3","PDGFRA","PDGFRB"],
    "Stress module":["NFE2L2","HMOX1","TXNIP","TXNRD1","GPX4","SLC7A11"]
}
palette={"DILI anchors":orange,"Direct targets":blue,"Stress module":teal}
for name,genes in groups.items():
    available=[g for g in genes if g in expr.columns]
    standardized=expr[available].apply(lambda x:(x-x.mean())/(x.std(ddof=0) if x.std(ddof=0)>0 else 1))
    temp=pd.DataFrame({"visit":expr.visit,"module":standardized.mean(axis=1)}).groupby('visit').agg(mean=('module','mean'),se=('module','sem')).reset_index()
    ax.plot(temp.visit,temp['mean'],marker='o',ms=2.5,lw=1.1,color=palette[name],label=name)
    ax.fill_between(temp.visit,temp['mean']-temp.se,temp['mean']+temp.se,color=palette[name],alpha=.12,lw=0)
ax.axhline(0,color="#AAA",lw=.6); ax.set_xlabel("Study visit"); ax.set_ylabel("Mean standardized log2 CPM")
ax.set_title("Module summaries remain near baseline",loc="left",fontweight="bold")
ax.legend(frameon=False); clean(ax); letter(ax,"F")

note="Single-arm whole-blood RNA-seq; no adjudicated DILI endpoint. Participant is the inference unit; this analysis tests systemic detectability only."
fig.text(.10,.025,note,ha="left",va="bottom",fontsize=5.7,color=gray)

stem=OUT/"Supplementary_Figure_S5_Whole_Blood_Sensitivity"
fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight')
fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight')
fig.savefig(stem.with_suffix('.png'),dpi=600,bbox_inches='tight',facecolor='white')
fig.savefig(stem.with_suffix('.tiff'),dpi=600,bbox_inches='tight',facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)

caption="""Supplementary Figure S5 | Whole-blood longitudinal sensitivity analysis in GSE299128. (A) Availability of 57 samples across nine scheduled visits in seven participants receiving nintedanib. (B) Mean participant-specific log2-CPM slope with 95% bootstrap confidence intervals for 30 prespecified genes; exact sign-flip P values were Benjamini–Hochberg adjusted and no gene passed FDR (all q=1.00). (C) Participant-specific slopes. (D) Mean slopes compared with median last-minus-baseline change. (E) FBP1 trajectories. (F) Visit-level averages of within-gene standardized expression for DILI-anchor, direct-target and stress-response modules; shading is the standard error across available samples and is descriptive. GSE299128 is a small, single-arm whole-blood study without adjudicated DILI. It evaluates systemic longitudinal detectability and neither validates nor refutes a liver-specific toxicity mechanism."""
(OUT/"Supplementary_Figure_S5_caption.md").write_text(caption+"\n",encoding='utf-8')

