#!/usr/bin/env python3
"""Cell Press square graphical abstract; separate from numbered result figures."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT
OUT.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

# Four inches at 300 dpi gives the exact 1200 × 1200 px Cell Press raster canvas.
fig, ax = plt.subplots(figsize=(4, 4))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")
fig.patch.set_facecolor("white")

navy = "#173B57"
blue = "#2D6F9F"
teal = "#2A9D8F"
orange = "#E76F51"
gray = "#5B6573"
line = "#D8E0E6"
pale_blue = "#EAF2F7"
pale_teal = "#EFF8F6"
pale_orange = "#FFF3EE"
light = "#F6F8FA"


def card(x, y, w, h, title, edge, fill):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.7,rounding_size=1.6",
        facecolor=fill, edgecolor=edge, linewidth=1.15,
    )
    ax.add_patch(box)
    ax.text(x + 2.2, y + h - 2.0, title, fontsize=7.0, fontweight="bold", color=edge, va="top")


def arrow(x1, y1, x2, y2, color=gray):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=10,
        linewidth=1.0, color=color, connectionstyle="arc3,rad=0",
    ))


ax.text(
    50, 97.3,
    "Role-separated map of nintedanib-associated liver injury",
    fontsize=9.2, fontweight="bold", color=navy, ha="center", va="top",
)
ax.text(
    50, 92.7,
    "Fibrotic interstitial lung disease  •  role-separated evidence  •  unit-aware public-data reanalysis",
    fontsize=4.9, color=gray, ha="center", va="top",
)

# 1. Clinical problem.
card(4, 74.5, 92, 14.0, "1  Clinical problem", blue, pale_blue)
ax.plot([12.5, 12.5], [84.2, 79.5], color=navy, lw=1.2)
ax.plot([12.5, 10.7], [82.4, 81.0], color=navy, lw=0.8)
ax.plot([12.5, 14.3], [82.4, 81.0], color=navy, lw=0.8)
ax.add_patch(Ellipse((9.9, 79.1), 4.0, 6.8, angle=-9, facecolor="#B9D8E9", edgecolor=blue, lw=0.8))
ax.add_patch(Ellipse((15.1, 79.1), 4.0, 6.8, angle=9, facecolor="#B9D8E9", edgecolor=blue, lw=0.8))
ax.add_patch(FancyBboxPatch((19.0, 80.0), 8.5, 3.0, boxstyle="round,pad=0.15,rounding_size=1.5", facecolor="white", edgecolor=orange, lw=0.8))
ax.add_patch(Polygon([[23.25, 80.0], [27.5, 80.0], [27.5, 83.0], [23.25, 83.0]], closed=True, facecolor=orange, edgecolor="none"))
ax.text(31, 83.4, "Nintedanib slows lung-function loss across fibrotic ILD", fontsize=5.6, color=navy, fontweight="bold", va="top")
ax.text(31, 79.8, "Treatment-limiting DILI lacks a molecular susceptibility map", fontsize=5.6, color=orange, fontweight="bold", va="top")
ax.text(31, 76.8, "Clinical question: distinguish drug-engaged targets from downstream injury", fontsize=4.6, color=gray, va="top")

arrow(50, 74.0, 50, 70.7)

# 2. Observed public evidence.
card(4, 51.0, 92, 19.0, "2  Observed public evidence", teal, pale_teal)
evidence = [
    (7.0, 61.7, 27.0, "FAERS active comparator", "Heterogeneous signal; ILD sensitivity"),
    (36.5, 61.7, 27.0, "Human DILI proteomics", "Source-reanalysis phenotype anchors"),
    (66.0, 61.7, 27.0, "ChEMBL + HepG2 LINCS", "Engagement + six-dose response"),
    (21.5, 54.3, 27.0, "Healthy-liver single cells", "5 donors + unsigned virtual KO"),
    (51.5, 54.3, 27.0, "Rodent lung/tissue data", "Animal-, pool-, and slice-aware"),
]
for x, y, w, head, sub in evidence:
    ax.add_patch(FancyBboxPatch((x, y), w, 5.7, boxstyle="round,pad=0.25,rounding_size=0.9", facecolor="white", edgecolor=line, lw=0.65))
    ax.add_patch(Circle((x + 2.1, y + 2.85), 0.65, facecolor=teal, edgecolor="white", lw=0.4))
    ax.text(x + 3.5, y + 3.65, head, fontsize=4.6, fontweight="bold", color=navy, va="center")
    ax.text(x + 3.5, y + 1.85, sub, fontsize=4.0, color=gray, va="center")

arrow(50, 50.5, 50, 47.2)

# 3. Role-separated evidence map.
card(4, 27.0, 92, 19.5, "3  Role-separated evidence map", orange, pale_orange)
ax.plot([50, 50], [29.0, 42.7], color="#E5C9BF", lw=0.7)
ax.text(8, 41.8, "T track  |  exposure-proximal", fontsize=5.5, fontweight="bold", color=navy, va="top")
ax.text(8, 38.1, "FGFR1 — exploratory T2-gate candidate", fontsize=5.7, fontweight="bold", color=orange, va="top")
ax.text(8, 34.7, "Direct pharmacology + endothelial detection +", fontsize=4.5, color=gray, va="top")
ax.text(8, 32.2, "Cross-seed-median virtual-KO gate", fontsize=4.3, color=gray, va="top")
ax.text(8, 30.5, "Strict all-seed T2 set empty", fontsize=4.3, color=gray, va="top")
ax.text(54, 41.8, "P track  |  downstream phenotype", fontsize=5.5, fontweight="bold", color=navy, va="top")
ax.text(54, 38.1, "FBP1 — P1 human DILI anchor", fontsize=5.7, fontweight="bold", color=orange, va="top")
ax.text(54, 34.7, "GSTA1 / OTC — P2 injury–recovery readouts", fontsize=4.5, color=gray, va="top")
ax.text(54, 32.2, "No phenotype anchor is a direct nintedanib target", fontsize=4.5, color=gray, va="top")
ax.text(50, 28.0, "No T1 confirmation  •  no protective/harmful direction inferred", fontsize=4.5, color=blue, ha="center", fontweight="bold")

arrow(50, 26.5, 50, 23.2)

# 4. Prospective functional confirmation.
card(4, 5.0, 92, 17.5, "4  Stage-gated prospective confirmation", navy, light)
ax.add_patch(Ellipse((13.0, 13.5), 10.0, 6.5, angle=-8, facecolor="#D99077", edgecolor=orange, lw=0.8))
ax.add_patch(Ellipse((16.0, 13.0), 3.8, 4.2, angle=15, facecolor="#C97860", edgecolor=orange, lw=0.65))
ax.text(23, 18.1, "Bidirectional FGFR1 / FBP1 perturbation + reciprocal rescue", fontsize=5.0, color=navy, fontweight="bold", va="top")
ax.text(23, 14.8, "Donor-replicated hepatocyte–LSEC/organoid systems; CES1 and LC–MS/MS exposure control", fontsize=4.4, color=gray, va="top")
ax.text(23, 11.7, "GO: reproducible liver rescue at matched exposure + retained pulmonary response", fontsize=4.5, color=teal, fontweight="bold", va="top")
ax.text(23, 8.5, "NO-GO: absent rescue, exposure artifact, unstable direction, or weakened lung response", fontsize=4.4, color=orange, va="top")

ax.text(50, 1.9, "Output: a falsifiable target–phenotype map — not a claimed causal or protective target", fontsize=4.8, color=navy, ha="center", fontweight="bold")

fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
base = OUT / "Graphical_Abstract_Public_Data_Target_Map"
fig.savefig(base.with_suffix(".pdf"), facecolor="white")
fig.savefig(base.with_suffix(".svg"), facecolor="white")
fig.savefig(base.with_suffix(".png"), dpi=300, facecolor="white")
fig.savefig(base.with_suffix(".tiff"), dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
