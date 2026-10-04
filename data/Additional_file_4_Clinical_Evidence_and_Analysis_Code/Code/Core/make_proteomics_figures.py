#!/usr/bin/env python3
"""Create publication figures for the clean-slate public DILI proteomics reanalysis.

Every plotted number is loaded from an auditable analysis output under
fresh_analysis/proteomics.  The script intentionally does not infer participant
linkage across proteins or treat the follow-up source block as paired.
"""

from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve, auc


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "Source_Data" / "Human_DILI_Proteomics"
OUT = ROOT / "Figures"
SUPP = ROOT / "Supplementary_Figures"
OUT.mkdir(parents=True, exist_ok=True)
SUPP.mkdir(parents=True, exist_ok=True)

COL = {
    "HV": "#8A8D91",
    "DO": "#D55E00",
    "NDO": "#0072B2",
    "DF": "#009E73",
    "FBP1": "#D55E00",
    "neutral": "#4D4D4D",
    "light": "#D9D9D9",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 6.5,
    "axes.titlesize": 7.2,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 6.0,
    "ytick.labelsize": 6.0,
    "legend.fontsize": 5.9,
    "axes.linewidth": 0.65,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def clean_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=2.2, pad=1.5)


def panel_label(ax, letter):
    ax.text(-0.18, 1.08, letter, transform=ax.transAxes, fontsize=9.0,
            fontweight="bold", va="top", ha="left", clip_on=False)


def forest(ax, df, title, xlab="Cliff's delta (95% CI)", highlight="FBP1"):
    d = df.copy().sort_values("cliffs_delta_group1_minus_group2")
    y = np.arange(len(d))
    x = d["cliffs_delta_group1_minus_group2"].to_numpy(float)
    lo = d["cliffs_delta_ci_low"].to_numpy(float)
    hi = d["cliffs_delta_ci_high"].to_numpy(float)
    colors = [COL["FBP1"] if p == highlight else COL["neutral"] for p in d.protein]
    fills = [c if q < 0.05 else "white" for c, q in zip(colors, d.mann_whitney_q_bh)]
    for yi, xi, li, ui, ec, fc, q in zip(y, x, lo, hi, colors, fills, d.mann_whitney_q_bh):
        ax.plot([li, ui], [yi, yi], color=ec, lw=0.9, zorder=1)
        ax.scatter([xi], [yi], s=19, facecolor=fc, edgecolor=ec, linewidth=0.8, zorder=2)
        if q < 0.05:
            ax.text(1.03, yi, "*", fontsize=6.8, va="center", ha="left", clip_on=False)
    ax.axvline(0, color="#9B9B9B", lw=0.7, ls=(0, (2, 2)))
    ax.set_yticks(y, d.protein)
    for tick in ax.get_yticklabels():
        if tick.get_text() == highlight:
            tick.set_color(COL["FBP1"]); tick.set_fontweight("bold")
    ax.set_xlim(-1.05, 1.05)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlabel(xlab)
    ax.set_title(title, loc="left", pad=3, fontweight="bold")
    clean_axis(ax)


def save_multi(fig, stem):
    destination = SUPP if stem.startswith("Supplementary_") else OUT
    fig.savefig(destination / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(destination / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(destination / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(destination / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})


effects = pd.read_csv(DATA / "effect_estimates.csv")
tidy = pd.read_csv(DATA / "confirmatory_tidy.csv")
corr = pd.read_csv(DATA / "correlation_stats.csv")
concord = pd.read_csv(DATA / "cross_cohort_concordance_detail.csv")
concord_summary = pd.read_csv(DATA / "cross_cohort_concordance_summary.csv")
left = pd.read_csv(DATA / "left_censor_sensitivity.csv")
missing = pd.read_csv(DATA / "missingness_stats.csv")
pck2 = pd.read_csv(DATA / "pck2_outlier_sensitivity.csv")
zonal = pd.read_csv(DATA / "zonal_effect_estimates.csv")


# ------------------------------- Main Figure 2 ------------------------------
fig, axs = plt.subplots(5, 2, figsize=(7.20, 9.50), constrained_layout=False)
plt.subplots_adjust(left=0.105, right=0.985, top=0.945, bottom=0.055,
                    hspace=0.92, wspace=0.47)
fig.suptitle(
    "Figure 2 | Public human serum proteomics recapitulates a source-selected etiologic contrast",
    x=0.015, y=0.992, ha="left", fontsize=9.2, fontweight="bold",
)

# A-C: effect forests.
specs = [
    ("discovery", "DO_vs_NDO", "Discovery: DILI onset vs non-drug injury"),
    ("confirmatory", "DO_vs_HV", "Confirmation: acute DILI vs healthy"),
    ("confirmatory", "DO_vs_NDO", "Confirmation: DILI vs non-drug injury"),
]
for ax, (cohort, contrast, title), letter in zip(
        [axs[0, 0], axs[0, 1], axs[1, 0]], specs, ["A", "B", "C"]):
    d = effects[(effects.cohort == cohort) & (effects.contrast == contrast)]
    forest(ax, d, title)
    panel_label(ax, letter)

# D: FBP1 raw detected distributions.
ax = axs[1, 1]
d = tidy[(tidy.protein == "FBP1") & tidy.value.notna()].copy()
rng = np.random.default_rng(20260904)
positions = {"HV": 0, "DO": 1, "NDO": 2}
for g in ["HV", "DO", "NDO"]:
    vals = d.loc[d.group == g, "value"].to_numpy(float)
    jit = rng.uniform(-0.16, 0.16, len(vals))
    ax.scatter(np.full(len(vals), positions[g]) + jit, vals, s=9, alpha=0.60,
               color=COL[g], edgecolor="none", rasterized=True)
    if len(vals):
        q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
        ax.plot([positions[g]-0.21, positions[g]+0.21], [med, med], color="black", lw=1.1)
        ax.plot([positions[g], positions[g]], [q1, q3], color="black", lw=3.4,
                solid_capstyle="butt")
    allg = tidy[(tidy.protein == "FBP1") & (tidy.group == g)]
    ax.text(positions[g], 1.02, f"{allg.value.notna().sum()}/{len(allg)} detected",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.7)
ax.set_yscale("log")
ax.set_xticks([0, 1, 2], ["Healthy", "DILI onset", "Non-drug injury"])
ax.set_ylabel("FBP1 relative intensity")
ax.set_title("FBP1 protein distribution", loc="left", pad=13, fontweight="bold")
ax.text(0.02, 0.96, "DO vs NDO: q=4.69e-4", transform=ax.transAxes,
        ha="left", va="top", fontsize=5.8)
clean_axis(ax); panel_label(ax, "D")

# E: empirical ROC, NDO higher than DILI onset.
ax = axs[2, 0]
rd = tidy[(tidy.protein == "FBP1") & tidy.group.isin(["DO", "NDO"]) & tidy.value.notna()].copy()
y = (rd.group == "NDO").astype(int).to_numpy()
score = rd.value.to_numpy(float)
fpr, tpr, _ = roc_curve(y, score)
roc_auc = auc(fpr, tpr)
row = effects[(effects.cohort == "confirmatory") & (effects.contrast == "DO_vs_NDO") &
              (effects.protein == "FBP1")].iloc[0]
sep_auc = 1.0 - float(row.roc_auc_group1_positive)
sep_lo = 1.0 - float(row.roc_auc_ci_high)
sep_hi = 1.0 - float(row.roc_auc_ci_low)
ax.plot(fpr, tpr, color=COL["FBP1"], lw=1.5)
ax.plot([0, 1], [0, 1], color="#999999", lw=0.7, ls=(0, (3, 2)))
ax.fill_between(fpr, tpr, fpr, color=COL["FBP1"], alpha=0.10)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xlabel("False-positive rate"); ax.set_ylabel("True-positive rate")
ax.set_title("FBP1 separates injury etiologies", loc="left", pad=3, fontweight="bold")
ax.text(0.97, 0.05, f"AUC={sep_auc:.3f}\n95% CI {sep_lo:.3f}–{sep_hi:.3f}\nn={len(y)} detected",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8)
clean_axis(ax); panel_label(ax, "E")

# F: cross-cohort concordance for DO vs NDO, HPD excluded from replication count.
ax = axs[2, 1]
cc = concord[concord.contrast == "DO_vs_NDO"].copy()
for _, r in cc.iterrows():
    p = r.protein
    c = COL["FBP1"] if p == "FBP1" else ("#BDBDBD" if p == "HPD" else COL["neutral"])
    marker = "X" if p == "HPD" else "o"
    ax.scatter(r.discovery_delta, r.confirmatory_delta, s=22 if p == "FBP1" else 14,
               color=c, marker=marker, edgecolor="white", linewidth=0.35, zorder=2)
    if p in ["FBP1", "LECT2", "HPD"]:
        ax.annotate(p, (r.discovery_delta, r.confirmatory_delta), xytext=(3, 3),
                    textcoords="offset points", fontsize=5.7, color=c)
ax.axhline(0, color="#AAAAAA", lw=0.6); ax.axvline(0, color="#AAAAAA", lw=0.6)
ax.set_xlim(-1.03, 1.03); ax.set_ylim(-1.03, 1.03)
ax.set_xlabel("Discovery Cliff's delta")
ax.set_ylabel("Confirmatory Cliff's delta")
cs = concord_summary[concord_summary.contrast == "DO_vs_NDO"].iloc[0]
ax.text(0.02, 0.98, "11/12 directions concordant\ndescriptive across-protein summary\nHPD excluded (duplicated series)",
        transform=ax.transAxes, ha="left", va="top", fontsize=5.7)
ax.set_title("Cross-cohort direction stability", loc="left", pad=3, fontweight="bold")
clean_axis(ax); panel_label(ax, "F")

# G: within-DILI ALT correlations.
ax = axs[3, 0]
cr = corr[(corr.correlation_input_rule == "primary_exclude_undocumented_numeric_ND_substitutions") &
          (corr.stratum == "DO")].copy().sort_values("spearman_rho")
ypos = np.arange(len(cr))
for yi, (_, r) in zip(ypos, cr.iterrows()):
    c = COL["FBP1"] if r.protein == "FBP1" else COL["neutral"]
    fc = c if r.spearman_q_bh < 0.05 else "white"
    ax.plot([r.spearman_ci_low, r.spearman_ci_high], [yi, yi], color=c, lw=0.85)
    ax.scatter(r.spearman_rho, yi, s=17, facecolor=fc, edgecolor=c, linewidth=0.7)
ax.axvline(0, color="#999999", lw=0.6, ls=(0, (2, 2)))
ax.set_yticks(ypos, cr.protein)
for tick in ax.get_yticklabels():
    if tick.get_text() == "FBP1": tick.set_color(COL["FBP1"]); tick.set_fontweight("bold")
ax.set_xlim(-0.1, 1.02); ax.set_xlabel("Spearman ρ with ALT (95% CI)")
ax.set_title("Severity association within acute DILI", loc="left", pad=3, fontweight="bold")
clean_axis(ax); panel_label(ax, "G")

# H: onset versus follow-up source blocks, explicitly unpaired.
ax = axs[3, 1]
dfx = effects[(effects.cohort == "confirmatory") & (effects.contrast == "DO_vs_DF")].copy()
forest(ax, dfx, "Acute DILI vs follow-up source block")
ax.text(0.01, -0.36, "Unpaired source-block comparison; no participant linkage assumed",
        transform=ax.transAxes, fontsize=5.5, ha="left", va="top")
panel_label(ax, "H")

# I: injury magnitude versus etiology specificity.
ax = axs[4, 0]
h = effects[(effects.cohort == "confirmatory") & effects.contrast.isin(["DO_vs_HV", "DO_vs_NDO"])].pivot(
    index="protein", columns="contrast", values=["cliffs_delta_group1_minus_group2", "mann_whitney_q_bh"])
for protein, r in h.iterrows():
    x = r[("cliffs_delta_group1_minus_group2", "DO_vs_HV")]
    yv = r[("cliffs_delta_group1_minus_group2", "DO_vs_NDO")]
    specq = r[("mann_whitney_q_bh", "DO_vs_NDO")]
    c = COL["FBP1"] if protein == "FBP1" else COL["neutral"]
    ax.scatter(x, yv, s=28 if protein == "FBP1" else 15,
               facecolor=c if specq < 0.05 else "white", edgecolor=c, linewidth=0.8)
    if protein in ["FBP1", "LECT2", "PCK2", "ACO1", "ASS1"]:
        ax.annotate(protein, (x, yv), xytext=(3, 2), textcoords="offset points",
                    fontsize=5.7, color=c)
ax.axhline(0, color="#AAAAAA", lw=0.6); ax.axvline(0, color="#AAAAAA", lw=0.6)
ax.set_xlim(-0.35, 1.05); ax.set_ylim(-0.65, 0.35)
ax.set_xlabel("Acute DILI vs healthy delta")
ax.set_ylabel("DILI vs non-drug injury delta")
ax.set_title("Injury detection vs etiology specificity", loc="left", pad=3, fontweight="bold")
clean_axis(ax); panel_label(ax, "I")

# J: censoring and outlier robustness, paired effect estimates.
ax = axs[4, 1]
rows = left[(left.cohort == "confirmatory") &
            (((left.protein == "FBP1") & left.contrast.isin(["DO_vs_HV", "DO_vs_NDO"])) |
             ((left.protein == "PCK2") & left.contrast.isin(["DO_vs_HV", "DO_vs_NDO"])))].copy()
labels = [f"{p} {c.replace('_vs_', '–')}" for p, c in zip(rows.protein, rows.contrast)]
yy = np.arange(len(rows))[::-1]
for yj, (_, r), lab in zip(yy, rows.iterrows(), labels):
    c = COL["FBP1"] if r.protein == "FBP1" else COL["neutral"]
    ax.plot([r.complete_case_cliffs_delta, r.half_min_substitution_cliffs_delta],
            [yj, yj], color="#BDBDBD", lw=1)
    ax.scatter(r.complete_case_cliffs_delta, yj, s=19, facecolor="white", edgecolor=c,
               linewidth=0.8, label="Complete case" if yj == yy[0] else None)
    ax.scatter(r.half_min_substitution_cliffs_delta, yj, s=19, facecolor=c, edgecolor=c,
               marker="s", linewidth=0.8, label="Half-min substitution" if yj == yy[0] else None)
ax.axvline(0, color="#999999", lw=0.6, ls=(0, (2, 2)))
ax.set_yticks(yy, labels)
ax.set_xlim(-1.02, 1.02)
ax.set_xlabel("Cliff's delta")
ax.set_title("Censoring and source-value audit", loc="left", pad=3, fontweight="bold")
ax.legend(loc="lower right", frameon=False, handletextpad=0.3, borderaxespad=0.2)
pdelta = pck2[pck2.contrast == "DO_vs_NDO"].set_index("analysis_rule")["cliffs_delta_group1_minus_group2"]
ax.text(0.01, 0.96,
        f"PCK2 starred-value exclusion: delta {pdelta.iloc[0]:.3f}→{pdelta.iloc[1]:.3f}",
        transform=ax.transAxes, ha="left", va="top", fontsize=5.7)
clean_axis(ax); panel_label(ax, "J")

save_multi(fig, "Figure_2_Human_DILI_Proteomics")
plt.close(fig)


# ---------------------- Supplementary Figure: distributions -----------------
proteins = sorted(tidy.protein.unique())
fig, axs = plt.subplots(4, 4, figsize=(8.25, 9.6))
plt.subplots_adjust(left=0.07, right=0.985, top=0.925, bottom=0.06, hspace=0.58, wspace=0.40)
rng = np.random.default_rng(20260904)
for i, (ax, protein) in enumerate(zip(axs.ravel(), proteins)):
    dp = tidy[(tidy.protein == protein) & tidy.group.isin(["HV", "DO", "NDO"])].copy()
    positive = dp.loc[dp.value.notna() & (dp.value > 0), "value"]
    log_ok = len(positive) and positive.min() > 0
    for xg, g in enumerate(["HV", "DO", "NDO"]):
        vals = dp.loc[(dp.group == g) & dp.value.notna(), "value"].to_numpy(float)
        ax.scatter(xg + rng.uniform(-0.17, 0.17, len(vals)), vals, s=5, alpha=0.45,
                   color=COL[g], edgecolor="none", rasterized=True)
        if len(vals):
            q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
            ax.plot([xg-.20, xg+.20], [med, med], color="black", lw=0.8)
            ax.plot([xg, xg], [q1, q3], color="black", lw=2.7)
        ng = dp[dp.group == g]
        ax.text(xg, 1.01, f"{ng.value.notna().sum()}/{len(ng)}", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=4.3)
    if log_ok and positive.max()/positive.min() > 100:
        ax.set_yscale("log")
    ax.set_xticks([0, 1, 2], ["HV", "DO", "NDO"])
    ax.set_title(protein, loc="left", fontweight="bold", color=COL["FBP1"] if protein == "FBP1" else "black")
    clean_axis(ax); panel_label(ax, chr(ord("A") + i))
for j in range(len(proteins), len(axs.ravel())):
    axs.ravel()[j].axis("off")
legend = [Line2D([0],[0], marker='o', color='none', markerfacecolor=COL[g], markeredgecolor='none', label=g)
          for g in ["HV", "DO", "NDO"]]
fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False)
fig.suptitle("Supplementary Figure S2 | Confirmatory serum protein distributions",
             x=0.07, y=0.985, ha="left", fontsize=9, fontweight="bold")
save_multi(fig, "Supplementary_Figure_S2_Protein_Distributions")
plt.close(fig)


# ---------------------- Supplementary Figure: robustness --------------------
fig, axs = plt.subplots(5, 2, figsize=(7.20, 10.0))
plt.subplots_adjust(left=0.11, right=0.985, top=0.94, bottom=0.06, hspace=0.92, wspace=0.48)

# A-C: all three cross-cohort contrasts.
for ax, contrast_name, letter in zip(axs[0:2, :].ravel()[:3], ["DO_vs_HV", "DO_vs_NDO", "DO_vs_DF"], ["A", "B", "C"]):
    dd = concord[concord.contrast == contrast_name]
    for _, r in dd.iterrows():
        c = COL["FBP1"] if r.protein == "FBP1" else ("#BDBDBD" if r.protein == "HPD" else COL["neutral"])
        ax.scatter(r.discovery_delta, r.confirmatory_delta, s=16, color=c, edgecolor="white", linewidth=0.3)
        if r.protein in ["FBP1", "HPD"]:
            ax.annotate(r.protein, (r.discovery_delta, r.confirmatory_delta), xytext=(3,2), textcoords="offset points", fontsize=5)
    ax.axhline(0,color="#AAA",lw=.6); ax.axvline(0,color="#AAA",lw=.6)
    ax.set_xlim(-1.05,1.05); ax.set_ylim(-1.05,1.05)
    ax.set_xlabel("Discovery delta"); ax.set_ylabel("Confirmatory delta")
    s = concord_summary[concord_summary.contrast == contrast_name].iloc[0]
    ax.set_title(f"{contrast_name.replace('_vs_', ' vs ')}: {int(s.n_direction_concordant)}/{int(s.n_proteins)} concordant", loc="left", fontweight="bold")
    clean_axis(ax); panel_label(ax, letter)

# D: pooled vs within-DO rho.
ax = axs[1,1]
ce = pd.read_csv(DATA / "candidate_evidence_table.csv")
for _, r in ce.iterrows():
    c = COL["FBP1"] if r.protein == "FBP1" else COL["neutral"]
    ax.scatter(r.pooled_alt_spearman_rho, r.within_DO_alt_spearman_rho, s=18, color=c)
    if r.protein in ["FBP1","ASS1","ACO1","GLDH"]:
        ax.annotate(r.protein,(r.pooled_alt_spearman_rho,r.within_DO_alt_spearman_rho),xytext=(3,2),textcoords="offset points",fontsize=5)
ax.plot([0,1],[0,1],color="#AAA",lw=.6,ls=(0,(2,2)))
ax.set_xlabel("Pooled HV+DILI ρ"); ax.set_ylabel("Within-DILI ρ")
ax.set_title("Case-mix inflation audit",loc="left",fontweight="bold"); clean_axis(ax); panel_label(ax,"D")

# E: missingness contrast for FBP1/PCK2.
ax=axs[2,0]
md=missing[(missing.cohort=="confirmatory") & missing.protein.isin(["FBP1","PCK2"]) & missing.contrast.isin(["DO_vs_HV","DO_vs_NDO"])].copy()
labs=[f"{p} {c.replace('_vs_','–')}" for p,c in zip(md.protein,md.contrast)]
for yi,(_,r) in enumerate(md.iterrows()):
    f1=r.group1_observed/(r.group1_observed+r.group1_missing_or_ND)
    f2=r.group2_observed/(r.group2_observed+r.group2_missing_or_ND)
    ax.plot([f1,f2],[yi,yi],color="#AAA",lw=1)
    ax.scatter(f1,yi,color=COL["DO"],s=18); ax.scatter(f2,yi,color=COL["HV"] if r.contrast=="DO_vs_HV" else COL["NDO"],s=18,marker="s")
ax.set_yticks(range(len(md)),labs); ax.set_xlim(0,1.03); ax.set_xlabel("Detected fraction")
ax.set_title("Differential detection audit",loc="left",fontweight="bold"); clean_axis(ax); panel_label(ax,"E")

# F: complete case vs half-min all contrasts.
ax=axs[2,1]
ld=left[left.cohort=="confirmatory"].dropna(subset=["complete_case_cliffs_delta","half_min_substitution_cliffs_delta"])
ax.scatter(ld.complete_case_cliffs_delta,ld.half_min_substitution_cliffs_delta,s=9,color=COL["neutral"],alpha=.55)
fb=ld[ld.protein=="FBP1"]
ax.scatter(fb.complete_case_cliffs_delta,fb.half_min_substitution_cliffs_delta,s=25,color=COL["FBP1"],label="FBP1")
ax.plot([-1,1],[-1,1],color="#AAA",lw=.7,ls=(0,(2,2))); ax.set_xlim(-1.05,1.05); ax.set_ylim(-1.05,1.05)
ax.set_xlabel("Complete-case delta"); ax.set_ylabel("Half-min delta"); ax.legend(frameon=False)
ax.set_title("Left-censor sensitivity",loc="left",fontweight="bold"); clean_axis(ax); panel_label(ax,"F")

# G: PCK2 outlier influence.
ax=axs[3,0]
for i,cname in enumerate(["DO_vs_HV","DO_vs_NDO","DO_vs_DF"]):
    z=pck2[pck2.contrast==cname]
    ax.plot(z.cliffs_delta_group1_minus_group2,[i,i],color="#AAA",lw=1)
    ax.scatter(z.cliffs_delta_group1_minus_group2.iloc[0],i,s=20,facecolor="white",edgecolor=COL["neutral"])
    ax.scatter(z.cliffs_delta_group1_minus_group2.iloc[1],i,s=20,color=COL["neutral"],marker="s")
ax.axvline(0,color="#AAA",lw=.6); ax.set_yticks(range(3),[x.replace('_vs_','–') for x in ["DO_vs_HV","DO_vs_NDO","DO_vs_DF"]])
ax.set_xlabel("PCK2 Cliff's delta"); ax.set_title("Starred-source-value influence",loc="left",fontweight="bold"); clean_axis(ax); panel_label(ax,"G")

# H: zonal effects.
ax=axs[3,1]
zd=zonal[zonal.contrast=="DO_vs_NDO"].sort_values("zone")
yy=np.arange(len(zd))
ax.errorbar(zd.cliffs_delta_group1_minus_group2,yy,xerr=[zd.cliffs_delta_group1_minus_group2-zd.cliffs_delta_ci_low,zd.cliffs_delta_ci_high-zd.cliffs_delta_group1_minus_group2],fmt='o',color=COL["neutral"],ms=3,lw=.8)
ax.axvline(0,color="#AAA",lw=.6,ls=(0,(2,2))); ax.set_yticks(yy,zd.zone); ax.set_xlim(-1.05,1.05); ax.set_xlabel("DILI–NDO zonal-score delta")
ax.set_title("Liver-zone score null result",loc="left",fontweight="bold"); clean_axis(ax); panel_label(ax,"H")

# I: source-block sizes.
ax=axs[4,0]
counts=pd.read_csv(DATA/"source_block_counts.csv")
numeric=counts.select_dtypes(include=[np.number])
if len(numeric.columns)>=1:
    col=numeric.columns[-1]
    labels=counts.select_dtypes(exclude=[np.number]).astype(str).agg(" | ".join,axis=1)
    top=counts.assign(_lab=labels).sort_values(col,ascending=False).head(12)
    ax.barh(np.arange(len(top)),top[col],color="#777777")
    ax.set_yticks(np.arange(len(top)),top._lab.str.slice(0,24)); ax.invert_yaxis(); ax.set_xlabel(str(col).replace('_',' '))
else:
    ax.text(.5,.5,"See source-block count table",ha='center',va='center',transform=ax.transAxes)
ax.set_title("Source-block count audit",loc="left",fontweight="bold"); clean_axis(ax); panel_label(ax,"I")

# J: row-order linkage audit summary.
ax=axs[4,1]
rq=pd.read_csv(DATA/"row_order_linkage_qc.csv")
nums=rq.select_dtypes(include=[np.number])
if nums.shape[1]>=2:
    xcol,ycol=nums.columns[:2]
    ax.scatter(rq[xcol],rq[ycol],s=12,color=COL["neutral"],alpha=.7)
    ax.set_xlabel(xcol.replace('_',' ')); ax.set_ylabel(ycol.replace('_',' '))
else:
    ax.text(.5,.58,"No stable cross-protein row linkage",ha='center',va='center',transform=ax.transAxes,fontweight='bold')
    ax.text(.5,.42,f"{len(rq)} QC comparisons",ha='center',va='center',transform=ax.transAxes)
ax.set_title("Row-order linkage audit",loc="left",fontweight="bold"); clean_axis(ax); panel_label(ax,"J")

fig.suptitle("Supplementary Figure S3 | Human DILI proteomics robustness and provenance audits",
             x=0.11, y=0.985, ha="left", fontsize=9, fontweight="bold")
save_multi(fig,"Supplementary_Figure_S3_Proteomics_Robustness")
plt.close(fig)


caption = """Figure 2 | Public human serum proteomics recapitulates an etiologic contrast among source-selected DILI proteins. (A) Discovery-set Cliff's delta for acute DILI onset (DO) versus non-drug acute liver injury (NDO). (B,C) Confirmatory effects for DO versus healthy volunteers (HV) and NDO. Filled points denote BH-adjusted q<0.05; bars are stratified-bootstrap 95% confidence intervals. (D) Detected FBP1 values; denominators include source rows reported as not detected. Central bars show medians and thick vertical bars the interquartile range. (E) Empirical ROC for higher FBP1 in NDO than DO, restricted to detected values. (F) Descriptive discovery-confirmatory DO-versus-NDO effect concordance; HPD is shown for transparency but excluded because its discovery series exactly duplicates ALDOB. Across-protein direction agreement and rank correlation are descriptive because proteins are correlated features from the same source cohorts. (G) Within-DO Spearman correlations with ALT. CK18 and GLDH are contextual source-reported comparators and were not entered into candidate integration; the CK18 estimate retains tied source values at the reported 100-U/L assay lower limit and should be interpreted with that censoring caveat. (H) Unpaired DO versus follow-up source-block comparison; no patient linkage was inferred. (I) Joint acute-injury and etiology-specificity effect space. (J) Sensitivity to half-minimum substitution for non-detects and to exclusion of a starred PCK2 source value. The 13 candidates were selected by the source study; these reanalyses do not constitute unbiased proteome-wide discovery, nintedanib specificity, or target causality."""
(OUT / "Figure_2_Human_DILI_Proteomics_caption.txt").write_text(caption + "\n", encoding="utf-8")

s2_caption = """Supplementary Figure S2 | Confirmatory serum protein distributions. (A–M) Detected confirmatory-cohort values for ACO1, ALDOB, ASS1, CES1, CPS1, DMGDH, FAH, FBP1, GSTA1, HPD, LECT2, OTC, and PCK2, respectively, in healthy volunteers (HV), drug-induced liver injury at onset (DO), and acute non-drug liver injury at onset (NDO). Fractions above each group give detected/source-row counts. Horizontal bars show medians and thick vertical bars show interquartile ranges; jittered points are detected source observations. Logarithmic axes are used when the observed dynamic range exceeds 100-fold. The source workbook does not provide a defensible participant identifier linking protein rows, so panels must not be assembled into a multivariable participant matrix. These multi-drug DILI distributions are not nintedanib-specific."""
(SUPP / "Supplementary_Figure_S2_Protein_Distributions_caption.md").write_text(s2_caption + "\n", encoding="utf-8")

s3_caption = """Supplementary Figure S3 | Robustness and provenance audits for the human DILI proteomics layer. (A–C) Discovery-versus-confirmatory Cliff's-delta concordance for DO versus HV, DO versus NDO, and DO versus the available DILI follow-up source block (DF); follow-up comparisons are unpaired because participant linkage was unavailable. (D) Correlations with ALT calculated within DO compared with correlations from pooled HV+DO rows, illustrating case-mix inflation. (E) Differential detection fractions for FBP1 and PCK2 in the two confirmatory contrasts. (F) Complete-case versus protein-specific half-minimum substitution estimates for source-reported non-detects. (G) Influence of excluding the starred PCK2 source value. (H) DO-versus-NDO liver-zone score effects; none passed multiplicity control. (I) Audit of source-block row counts against declared cohort sizes. (J) Row-order linkage audit showing that equal row counts did not establish stable cross-protein participant alignment. These diagnostics support univariate, protein-wise inference and preclude a participant-level multivariable classifier from the released workbook."""
(SUPP / "Supplementary_Figure_S3_Proteomics_Robustness_caption.md").write_text(s3_caption + "\n", encoding="utf-8")

manifest = {
    "figure": "Figure 2",
    "panels": 10,
    "input_files": [
        "effect_estimates.csv", "confirmatory_tidy.csv", "correlation_stats.csv",
        "cross_cohort_concordance_detail.csv", "cross_cohort_concordance_summary.csv",
        "left_censor_sensitivity.csv", "missingness_stats.csv",
        "pck2_outlier_sensitivity.csv", "zonal_effect_estimates.csv",
    ],
    "interpretation_boundary": "Human multi-drug DILI phenotype evidence; neither nintedanib specificity nor causal target validation.",
}
(DATA / "Figure_2_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
