#!/usr/bin/env python3
"""Reproducible pooled-library analysis of GSE151374.

The statistical unit is the independently prepared pooled 10x library
(three mice pooled per library; n=3 libraries per condition/time point).
Cells are never treated as independent replicates.  The deliberately broad,
marker-rule annotation is used only to aggregate counts and cell fractions to
the library level; this is not a de-novo cell-atlas analysis.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import sys
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from matplotlib.lines import Line2D
from scipy import sparse, stats
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
RESULTS = ROOT / "results"
VALIDATION = ROOT / "validation"
FIGURES = ROOT.parents[1] / "figures"
for directory in (RESULTS, VALIDATION, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)

SEED = 151374
np.random.seed(SEED)

# GEO sample order and experimental mapping reconstructed from the Series Matrix/
# SOFT record.  replicate is a pooled-library replicate, not an individual mouse.
SAMPLE_ROWS = [
    ("GSM4576588", 14, "Control", 1),
    ("GSM4576589", 14, "BLM", 1),
    ("GSM4576590", 14, "BLM+NINT", 1),
    ("GSM4576591", 7, "Control", 1),
    ("GSM4576592", 7, "BLM", 1),
    ("GSM4576593", 7, "BLM+NINT", 1),
    ("GSM4576594", 14, "Control", 2),
    ("GSM4576595", 14, "BLM", 2),
    ("GSM4576596", 14, "BLM+NINT", 2),
    ("GSM4576597", 7, "Control", 2),
    ("GSM4576598", 7, "BLM", 2),
    ("GSM4576599", 7, "BLM+NINT", 2),
    ("GSM4576600", 14, "Control", 3),
    ("GSM4576601", 14, "BLM", 3),
    ("GSM4576602", 14, "BLM+NINT", 3),
    ("GSM4576603", 7, "Control", 3),
    ("GSM4576604", 7, "BLM", 3),
    ("GSM4576605", 7, "BLM+NINT", 3),
]

LINEAGE_MARKERS = {
    "Macrophage": ["Adgre1", "C1qa", "C1qb", "C1qc", "Cd68", "Csf1r", "Marco", "Mertk", "Apoe", "Fcgr1", "Lgals3"],
    "Monocyte": ["Ly6c2", "Ccr2", "Sell", "Plac8", "Fcgr3", "Lilrb4a", "Lst1", "Ctss"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Retnlg", "Csf3r", "Mmp8", "Camp", "Lcn2"],
    "Dendritic": ["Flt3", "Zbtb46", "Xcr1", "Clec10a", "Cd209a", "Ccr7", "Cd74", "H2-Ab1"],
    "T": ["Cd3d", "Cd3e", "Trac", "Cd247", "Lck", "Lat"],
    "B": ["Cd79a", "Cd79b", "Ms4a1", "Cd19", "Cd22", "Cd37", "Cd74", "H2-Aa"],
    "NK": ["Nkg7", "Klrd1", "Ncr1", "Prf1", "Gzmb", "Ccl5", "Xcl1", "Xcl2"],
    "Epithelial": ["Epcam", "Krt8", "Krt18", "Krt19", "Scgb1a1", "Sftpc"],
    "Endothelial": ["Pecam1", "Kdr", "Eng", "Emcn", "Cdh5", "Klf2"],
}

MODULES = {
    "Repair macrophage": ["Mrc1", "Retnla", "Chil3", "Il10", "Mertk", "Marco", "Igf1", "Trem2", "Gpnmb", "Fabp5"],
    "Inflammatory": ["Il1b", "Tnf", "Nfkbia", "Cxcl2", "Ccl2", "Ccr2", "S100a8", "S100a9", "Nlrp3"],
    "ECM remodeling": ["Spp1", "Tgfb1", "Mmp9", "Mmp12", "Timp1", "Lgals3", "Fn1", "Col1a1", "Col3a1"],
    "Alveolar homeostasis": ["Pparg", "Car4", "Marco", "Siglecf", "Fabp4", "Ear1", "Itgax"],
    "RTK-response": ["Fgfr1", "Fgfr2", "Fgfr3", "Pdgfra", "Pdgfrb", "Kdr", "Flt1", "Mapk1", "Mapk3", "Akt1", "Stat3"],
}

CANDIDATES = [
    "Fgfr1", "Fbp1", "Gsta1", "Otc",  # manuscript target/phenotype tracks
    "Spp1", "Mrc1", "Il1b", "Tgfb1", "Pparg", "Marco",
]


def track_role(gene: str) -> str:
    if gene == "Fgfr1":
        return "exposure_proximal_target"
    if gene == "Fbp1":
        return "DILI_phenotype_anchor_P1"
    if gene in {"Gsta1", "Otc"}:
        return "DILI_phenotype_anchor_P2"
    return "lung_phenotype_context"

COLORS = {
    "Control": "#7A7A7A",
    "BLM": "#D55E00",
    "BLM+NINT": "#0072B2",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def bh(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if pv.size == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    back = np.empty_like(adjusted)
    back[order] = adjusted
    q[valid] = back
    return q


def exact_permutation_p(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sided exact difference-in-means test for n=3 versus n=3."""
    x = np.r_[np.asarray(a, float), np.asarray(b, float)]
    na = len(a)
    observed = abs(float(np.mean(a) - np.mean(b)))
    diffs = []
    for idx in itertools.combinations(range(len(x)), na):
        idx = np.asarray(idx)
        mask = np.zeros(len(x), dtype=bool)
        mask[idx] = True
        diffs.append(abs(float(np.mean(x[mask]) - np.mean(x[~mask]))))
    return float(np.mean(np.asarray(diffs) >= observed - 1e-12))


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    """Hedges g for a minus b; returns NA when pooled variance is zero."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    df = len(a) + len(b) - 2
    sp2 = ((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / df
    if not np.isfinite(sp2) or sp2 <= 0:
        return np.nan
    d = (np.mean(a) - np.mean(b)) / math.sqrt(sp2)
    correction = 1 - 3 / (4 * df - 1)
    return float(d * correction)


def mean_diff_ci(a: np.ndarray, b: np.ndarray, alpha: float = 0.05) -> tuple[float, float, float]:
    """Welch t confidence interval for mean(a)-mean(b)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    diff = float(np.mean(a) - np.mean(b))
    va, vb = np.var(a, ddof=1) / len(a), np.var(b, ddof=1) / len(b)
    se = math.sqrt(max(va + vb, 0))
    if se == 0:
        return diff, diff, diff
    df = (va + vb) ** 2 / ((va**2) / (len(a) - 1) + (vb**2) / (len(b) - 1))
    crit = stats.t.ppf(1 - alpha / 2, df)
    return diff, float(diff - crit * se), float(diff + crit * se)


def gene_indices(names: np.ndarray, genes: list[str]) -> list[int]:
    lookup: dict[str, list[int]] = {}
    for i, name in enumerate(names):
        lookup.setdefault(str(name), []).append(i)
    return [i for gene in genes for i in lookup.get(gene, [])]


def classify_cells(x: sparse.csc_matrix, totals: np.ndarray, names: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Broad marker-rule labels; ambiguous marker ties are retained as Other."""
    score_rows = []
    for genes in LINEAGE_MARKERS.values():
        idx = gene_indices(names, genes)
        values = x[idx, :].toarray().astype(np.float64, copy=False)
        values *= (1e4 / np.maximum(totals, 1))[None, :]
        np.log1p(values, out=values)
        score_rows.append(values.mean(axis=0))
    scores = np.vstack(score_rows)
    order = np.argsort(scores, axis=0)
    top_idx = order[-1]
    top = scores[top_idx, np.arange(scores.shape[1])]
    second = scores[order[-2], np.arange(scores.shape[1])]
    labels = np.asarray(list(LINEAGE_MARKERS), dtype=object)[top_idx]
    labels[(top < 0.10) | ((top - second) < 0.05)] = "Other"
    return labels, scores


def logcpm(counts: np.ndarray, prior_count: float = 0.5) -> np.ndarray:
    libs = counts.sum(axis=0).astype(float)
    return np.log2((counts + prior_count) / (libs[None, :] + prior_count * counts.shape[0]) * 1e6)


def pca_table(expression: np.ndarray, meta: pd.DataFrame, prefix: str) -> pd.DataFrame:
    variances = np.nanvar(expression, axis=1)
    keep = np.isfinite(variances) & (variances > 0)
    idx = np.where(keep)[0]
    idx = idx[np.argsort(variances[idx])[-min(2000, len(idx)):]]
    model = PCA(n_components=2, random_state=SEED)
    coords = model.fit_transform(expression[idx, :].T)
    out = meta.copy()
    out[f"{prefix}_PC1"] = coords[:, 0]
    out[f"{prefix}_PC2"] = coords[:, 1]
    out[f"{prefix}_PC1_variance_pct"] = model.explained_variance_ratio_[0] * 100
    out[f"{prefix}_PC2_variance_pct"] = model.explained_variance_ratio_[1] * 100
    return out


def de_contrast(counts: np.ndarray, gene_ids: np.ndarray, symbols: np.ndarray,
                meta: pd.DataFrame, day: int) -> pd.DataFrame:
    sel = (meta["day"] == day) & meta["group"].isin(["BLM", "BLM+NINT"])
    sub_meta = meta.loc[sel].reset_index(drop=True)
    y = logcpm(counts[:, sel.to_numpy()])
    cpm = 2**y
    keep = (cpm >= 1).sum(axis=1) >= 3
    nint = sub_meta["group"].eq("BLM+NINT").to_numpy()
    blm = sub_meta["group"].eq("BLM").to_numpy()
    a, b = y[:, nint], y[:, blm]
    logfc = np.mean(a, axis=1) - np.mean(b, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        test = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    p = np.asarray(test.pvalue, float)
    p[~np.isfinite(p)] = 1.0
    p[~keep] = np.nan
    q = bh(p)
    effect = np.asarray([hedges_g(aa, bb) for aa, bb in zip(a, b)])
    return pd.DataFrame({
        "gene_id": gene_ids,
        "symbol": symbols,
        "day": day,
        "log2FC_NINT_vs_BLM": logfc,
        "hedges_g": effect,
        "mean_log2CPM_NINT": np.mean(a, axis=1),
        "mean_log2CPM_BLM": np.mean(b, axis=1),
        "p_welch": p,
        "q_BH": q,
        "filter_CPM1_in_3": keep,
        "n_pool_NINT": 3,
        "n_pool_BLM": 3,
    }).sort_values(["q_BH", "p_welch"], na_position="last")


def module_analysis(mac_log: np.ndarray, names: np.ndarray, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Collapse duplicate gene symbols by their mean logCPM before standardization.
    expr = pd.DataFrame(mac_log, index=names, columns=meta["sample_id"])
    expr = expr.groupby(level=0).mean()
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    score_rows = []
    membership_rows = []
    for module, genes in MODULES.items():
        available = [g for g in genes if g in z.index and np.isfinite(z.loc[g]).any()]
        for g in genes:
            membership_rows.append({"module": module, "gene": g, "present_in_matrix": g in expr.index,
                                    "used_in_score": g in available})
        scores = z.loc[available].mean(axis=0) if available else pd.Series(np.nan, index=expr.columns)
        for sample_id, value in scores.items():
            score_rows.append({"sample_id": sample_id, "module": module, "score": value,
                               "n_genes_used": len(available), "genes_used": ";".join(available)})
    scores = pd.DataFrame(score_rows).merge(meta, on="sample_id", validate="many_to_one")
    pd.DataFrame(membership_rows).to_csv(RESULTS / "module_gene_membership.csv", index=False)

    contrasts = []
    for (module, day), frame in scores[scores["group"].isin(["BLM", "BLM+NINT"])].groupby(["module", "day"], sort=False):
        a = frame.loc[frame.group == "BLM+NINT", "score"].to_numpy(float)
        b = frame.loc[frame.group == "BLM", "score"].to_numpy(float)
        diff, lo, hi = mean_diff_ci(a, b)
        welch = stats.ttest_ind(a, b, equal_var=False).pvalue
        contrasts.append({
            "module": module, "day": int(day), "mean_NINT": np.mean(a), "mean_BLM": np.mean(b),
            "difference_NINT_minus_BLM": diff, "CI95_low": lo, "CI95_high": hi,
            "hedges_g": hedges_g(a, b), "p_welch": welch,
            "p_exact_permutation": exact_permutation_p(a, b), "n_pool_each_group": 3,
        })
    contrasts = pd.DataFrame(contrasts)
    contrasts["q_BH_welch_10_tests"] = bh(contrasts["p_welch"].to_numpy())
    contrasts["q_BH_exact_10_tests"] = bh(contrasts["p_exact_permutation"].to_numpy())
    return scores, contrasts


def composition_analysis(props: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cell_type, day), frame in props[props.group.isin(["BLM", "BLM+NINT"])].groupby(["cell_type", "day"], sort=False):
        a = frame.loc[frame.group == "BLM+NINT", "proportion"].to_numpy(float)
        b = frame.loc[frame.group == "BLM", "proportion"].to_numpy(float)
        diff, lo, hi = mean_diff_ci(a, b)
        rows.append({"cell_type": cell_type, "day": int(day),
                     "mean_NINT": np.mean(a), "mean_BLM": np.mean(b),
                     "difference_NINT_minus_BLM": diff, "CI95_low": lo, "CI95_high": hi,
                     "hedges_g": hedges_g(a, b),
                     "p_welch": stats.ttest_ind(a, b, equal_var=False).pvalue,
                     "p_exact_permutation": exact_permutation_p(a, b),
                     "n_pool_each_group": 3})
    out = pd.DataFrame(rows)
    out["q_BH_welch_all_celltype_day_tests"] = bh(out.p_welch.to_numpy())
    out["q_BH_exact_all_celltype_day_tests"] = bh(out.p_exact_permutation.to_numpy())
    return out


def candidate_table(all_counts: np.ndarray, mac_counts: np.ndarray, names: np.ndarray,
                    meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for compartment, counts in [("All retained cells", all_counts), ("Macrophage stratum", mac_counts)]:
        lib = counts.sum(axis=0)
        for gene in CANDIDATES:
            idx = np.where(names == gene)[0]
            summed = counts[idx].sum(axis=0) if len(idx) else np.zeros(counts.shape[1])
            cpm = summed / np.maximum(lib, 1) * 1e6
            for j, sample in meta.iterrows():
                rows.append({**sample.to_dict(), "compartment": compartment, "gene": gene,
                             "raw_count": int(summed[j]), "CPM": float(cpm[j]),
                             "log2CPM_plus_0.5": float(np.log2(cpm[j] + 0.5)),
                             "detected_in_pool": bool(summed[j] > 0)})
    return pd.DataFrame(rows)


def candidate_contrasts(candidates: pd.DataFrame) -> pd.DataFrame:
    """Pool-level contrasts for prespecified targets and phenotype tracks."""
    rows = []
    for (compartment, gene, day), frame in candidates[
        candidates.group.isin(["BLM", "BLM+NINT"])
    ].groupby(["compartment", "gene", "day"], sort=False):
        a = frame.loc[frame.group == "BLM+NINT", "log2CPM_plus_0.5"].to_numpy(float)
        b = frame.loc[frame.group == "BLM", "log2CPM_plus_0.5"].to_numpy(float)
        diff, lo, hi = mean_diff_ci(a, b)
        # Constant all-zero genes are an explicit null/detectability result.
        p_welch = 1.0 if np.var(np.r_[a, b]) == 0 else float(stats.ttest_ind(a, b, equal_var=False).pvalue)
        rows.append({
            "compartment": compartment, "track_type": track_role(gene),
            "gene": gene, "day": int(day), "mean_log2CPM_NINT": float(np.mean(a)),
            "mean_log2CPM_BLM": float(np.mean(b)), "log2CPM_difference_NINT_minus_BLM": diff,
            "CI95_low": lo, "CI95_high": hi, "hedges_g": hedges_g(a, b),
            "p_welch": p_welch, "p_exact_permutation": exact_permutation_p(a, b),
            "detected_pools_NINT": int(frame.loc[frame.group == "BLM+NINT", "detected_in_pool"].sum()),
            "detected_pools_BLM": int(frame.loc[frame.group == "BLM", "detected_in_pool"].sum()),
            "n_pool_each_group": 3,
        })
    out = pd.DataFrame(rows)
    out["q_BH_welch_all_candidate_tests"] = bh(out.p_welch.to_numpy())
    out["q_BH_exact_all_candidate_tests"] = bh(out.p_exact_permutation.to_numpy())
    return out


def plot_s6(meta: pd.DataFrame, qc: pd.DataFrame, props: pd.DataFrame,
            all_pca: pd.DataFrame, mac_pca: pd.DataFrame, de7: pd.DataFrame,
            de14: pd.DataFrame, module_scores: pd.DataFrame,
            module_contrasts: pd.DataFrame, candidates: pd.DataFrame) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 8,
        "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "axes.linewidth": 0.6, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(14.2, 12.0), constrained_layout=False)
    gs = fig.add_gridspec(4, 6, left=0.055, right=0.955, bottom=0.055, top=0.965,
                          wspace=0.85, hspace=0.95)
    axes = [
        fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]),
        fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6]),
        fig.add_subplot(gs[2, 0:2]), fig.add_subplot(gs[2, 2:4]), fig.add_subplot(gs[2, 4:6]),
        fig.add_subplot(gs[3, 0:3]), fig.add_subplot(gs[3, 3:6]),
    ]
    for letter, ax in zip("ABCDEFGHIJ", axes):
        ax.text(-0.08, 1.10, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

    # A: retained cell yield at the actual inferential unit.
    ax = axes[0]
    q = qc.sort_values(["day", "group", "replicate"])
    xlabels = []
    x = np.arange(len(q))
    for _, row in q.iterrows():
        xlabels.append(f"d{row.day}\n{row.group.replace('BLM+','+')}\nr{row.replicate}")
    ax.bar(x, q.cells_retained, color=[COLORS[g] for g in q.group], width=0.78)
    ax.set_xticks(x, xlabels, rotation=70, ha="right")
    ax.set_ylabel("Retained cells / pooled library")
    ax.set_title("Library-level cell yield")
    ax.spines[["top", "right"]].set_visible(False)

    # B: QC summaries; each point is one pooled library.
    ax = axes[1]
    for group, frame in qc.groupby("group"):
        ax.scatter(frame.median_UMI_retained, frame.median_genes_retained,
                   s=28, color=COLORS[group], label=group, alpha=0.85, edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Median UMI / retained cell")
    ax.set_ylabel("Median genes / retained cell")
    ax.set_title("Per-library QC summaries")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="best")

    # C: all-cell pseudobulk PCA.
    ax = axes[2]
    for day, marker in [(7, "o"), (14, "s")]:
        for group, frame in all_pca[all_pca.day == day].groupby("group"):
            ax.scatter(frame.all_PC1, frame.all_PC2, s=42, marker=marker, color=COLORS[group],
                       edgecolor="white", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({all_pca.all_PC1_variance_pct.iloc[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({all_pca.all_PC2_variance_pct.iloc[0]:.1f}%)")
    ax.set_title("All-cell pseudobulk PCA")
    ax.spines[["top", "right"]].set_visible(False)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS[g], label=g) for g in COLORS]
    handles += [Line2D([0], [0], marker=m, color="#333333", linestyle="none", label=f"day {d}") for d, m in [(7, "o"), (14, "s")]]
    ax.legend(handles=handles, frameon=False, ncol=2, loc="best")

    # D: broad cell composition per pooled library.
    ax = axes[3]
    order = ["Macrophage", "Monocyte", "Neutrophil", "Dendritic", "T", "B", "NK", "Endothelial", "Epithelial", "Other"]
    palette = dict(zip(order, ["#3B6FB6", "#67A9CF", "#EF8A62", "#A6761D", "#6A51A3", "#DD3497", "#1B9E77", "#66C2A5", "#FC8D62", "#BDBDBD"]))
    wide = props.pivot(index="sample_id", columns="cell_type", values="proportion").fillna(0).reindex(meta.sample_id)
    bottom = np.zeros(len(wide))
    for cell_type in order:
        val = wide[cell_type].to_numpy() if cell_type in wide else np.zeros(len(wide))
        ax.bar(np.arange(len(wide)), val, bottom=bottom, width=0.82, color=palette[cell_type], label=cell_type)
        bottom += val
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of retained cells")
    ax.set_xticks(np.arange(len(wide)), [f"d{r.day}-{r.group.replace('BLM+','+N')}-r{r.replicate}" for _, r in meta.iterrows()], rotation=70, ha="right")
    # Keep the title below the two-row legend; a large title pad causes the
    # two text layers to collide in print-sized supplementary exports.
    ax.set_title("Broad-lineage composition", pad=6)
    ax.legend(frameon=False, ncol=5, bbox_to_anchor=(0.5, 1.13), loc="lower center")
    ax.spines[["top", "right"]].set_visible(False)

    # E: macrophage fraction, preserving the three pool values.
    ax = axes[4]
    pm = props[props.cell_type == "Macrophage"]
    positions = {(7, "BLM"): 0, (7, "BLM+NINT"): 1, (14, "BLM"): 3, (14, "BLM+NINT"): 4}
    for (day, group), xpos in positions.items():
        vals = pm[(pm.day == day) & (pm.group == group)].proportion.to_numpy()
        jitter = np.linspace(-0.08, 0.08, len(vals))
        ax.scatter(xpos + jitter, vals, color=COLORS[group], s=35, zorder=3, edgecolor="white", linewidth=0.5)
        mean = vals.mean()
        sem = stats.sem(vals)
        ax.errorbar(xpos, mean, yerr=sem, fmt="_", markersize=12, color="#222222", capsize=3, zorder=4)
    ax.set_xticks([0, 1, 3, 4], ["BLM", "+NINT", "BLM", "+NINT"])
    ax.text(0.5, -0.18, "day 7", transform=ax.get_xaxis_transform(), ha="center")
    ax.text(3.5, -0.18, "day 14", transform=ax.get_xaxis_transform(), ha="center")
    ax.set_ylabel("Macrophage fraction")
    ax.set_title("Macrophage abundance (pool-level)")
    ax.spines[["top", "right"]].set_visible(False)

    # F: macrophage pseudobulk PCA.
    ax = axes[5]
    for day, marker in [(7, "o"), (14, "s")]:
        for group, frame in mac_pca[mac_pca.day == day].groupby("group"):
            ax.scatter(frame.mac_PC1, frame.mac_PC2, s=42, marker=marker, color=COLORS[group], edgecolor="white", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({mac_pca.mac_PC1_variance_pct.iloc[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({mac_pca.mac_PC2_variance_pct.iloc[0]:.1f}%)")
    ax.set_title("Macrophage pseudobulk PCA")
    ax.spines[["top", "right"]].set_visible(False)

    # G/H: exploratory gene-level pseudobulk contrasts.
    for ax, de, day in [(axes[6], de7, 7), (axes[7], de14, 14)]:
        plot = de[de.filter_CPM1_in_3].copy()
        yv = -np.log10(plot.p_welch.clip(lower=1e-300))
        sig = plot.q_BH < 0.10
        ax.scatter(plot.log2FC_NINT_vs_BLM[~sig], yv[~sig], s=6, color="#BDBDBD", alpha=0.6, rasterized=True)
        ax.scatter(plot.log2FC_NINT_vs_BLM[sig], yv[sig], s=9, color="#7B3294", alpha=0.8, rasterized=True)
        ax.axvline(0, color="#555555", lw=0.6)
        ax.axhline(-np.log10(0.05), color="#888888", lw=0.5, ls="--")
        # No gene passed the prespecified multiplicity threshold.  Therefore
        # nominal top-P labels are intentionally omitted: they add no
        # inferential information and collide at journal print scale.
        ax.set_xlabel("log2FC, NINT vs BLM")
        ax.set_ylabel("−log10 Welch P")
        ax.set_title(f"Macrophage pseudobulk, day {day}")
        ax.spines[["top", "right"]].set_visible(False)

    # I: module scores show all pool values (no cell-level tests).
    ax = axes[8]
    modules = list(MODULES)
    xbase = np.arange(len(modules))
    for day, offset, marker in [(7, -0.16, "o"), (14, 0.16, "s")]:
        frame = module_contrasts[module_contrasts.day == day].set_index("module").reindex(modules)
        ax.errorbar(xbase + offset, frame.difference_NINT_minus_BLM,
                    yerr=[frame.difference_NINT_minus_BLM - frame.CI95_low,
                          frame.CI95_high - frame.difference_NINT_minus_BLM],
                    fmt=marker, color="#2166AC" if day == 7 else "#B2182B", label=f"day {day}", capsize=2.5, ms=4)
    ax.axhline(0, color="#555555", lw=0.6)
    ax.set_xticks(xbase, [m.replace(" ", "\n") for m in modules], rotation=25, ha="right")
    ax.set_ylabel("Module score difference\nNINT − BLM (95% CI)")
    ax.set_title("Prespecified macrophage programs")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    # J: absolute candidate/phenotype-track expression; zero-count pools stay explicit.
    ax = axes[9]
    cand = candidates[candidates.compartment == "Macrophage stratum"].copy()
    mat = cand.pivot(index="gene", columns="sample_id", values="log2CPM_plus_0.5").reindex(CANDIDATES).reindex(columns=meta.sample_id)
    im = ax.imshow(mat, aspect="auto", cmap="viridis", interpolation="nearest")
    ax.set_yticks(np.arange(len(mat)), mat.index)
    ax.set_xticks(np.arange(len(meta)), [f"d{r.day}-{r.group.replace('BLM+','+N')}-r{r.replicate}" for _, r in meta.iterrows()], rotation=70, ha="right")
    ax.set_title("FGFR1 target and phenotype anchors in macrophage pools")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("log2(CPM + 0.5)")

    fig.suptitle("GSE151374: pooled-library single-cell reanalysis of nintedanib response in bleomycin-exposed mouse lung", fontsize=11, fontweight="bold", y=0.993)
    outbase = FIGURES / "Supplementary_Figure_S6_GSE151374"
    fig.savefig(outbase.with_suffix(".pdf"), dpi=600)
    fig.savefig(outbase.with_suffix(".svg"), dpi=600)
    fig.savefig(outbase.with_suffix(".png"), dpi=600)
    fig.savefig(outbase.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main() -> None:
    files = {p.name.split("_")[0]: p for p in sorted(RAW.glob("GSM*_raw_feature_bc_matrix.h5"))}
    meta = pd.DataFrame(SAMPLE_ROWS, columns=["sample_id", "day", "group", "replicate"])
    meta["condition"] = "d" + meta.day.astype(str) + "_" + meta.group.str.replace("+", "plus", regex=False)
    meta["n_mice_pooled"] = 3
    meta["statistical_unit"] = "pooled 10x library"
    meta["raw_file"] = meta.sample_id.map(lambda x: files[x].name)
    meta["raw_file_sha256"] = meta.sample_id.map(lambda x: sha256(files[x]))
    meta.to_csv(RESULTS / "library_manifest.csv", index=False)

    all_counts = None
    mac_counts = None
    class_counts_total = None
    gene_ids = names = None
    qc_rows, count_rows, validation_rows = [], [], []

    for j, sample in meta.iterrows():
        path = files[sample.sample_id]
        with h5py.File(path, "r") as handle:
            group = handle["matrix"]
            shape = tuple(int(v) for v in group["shape"][:])
            data = group["data"][:]
            indices = group["indices"][:]
            indptr = group["indptr"][:]
            this_names = np.asarray([x.decode() for x in group["features/name"][:]], dtype=object)
            this_ids = np.asarray([x.decode() for x in group["features/id"][:]], dtype=object)
        if names is None:
            names, gene_ids = this_names, this_ids
            all_counts = np.zeros((shape[0], len(meta)), dtype=np.int64)
            mac_counts = np.zeros((shape[0], len(meta)), dtype=np.int64)
            class_counts_total = {k: np.zeros(shape[0], dtype=np.int64) for k in [*LINEAGE_MARKERS, "Other"]}
        else:
            if not (np.array_equal(names, this_names) and np.array_equal(gene_ids, this_ids)):
                raise RuntimeError(f"Feature ordering differs in {path}")
        x = sparse.csc_matrix((data, indices, indptr), shape=shape)
        totals = np.asarray(x.sum(axis=0)).ravel()
        features = np.diff(indptr)
        mito_idx = np.where(np.char.startswith(names.astype(str), "mt-"))[0]
        mito_counts = np.asarray(x[mito_idx, :].sum(axis=0)).ravel()
        mito = mito_counts / np.maximum(totals, 1)
        keep = (totals >= 500) & (features >= 200) & (features <= 6000) & (mito < 0.20)
        retained = x[:, keep]
        retained_totals = totals[keep]
        labels, scores = classify_cells(retained, retained_totals, names)

        all_counts[:, j] = np.asarray(retained.sum(axis=1)).ravel()
        macrophage = labels == "Macrophage"
        mac_counts[:, j] = np.asarray(retained[:, macrophage].sum(axis=1)).ravel()
        for cell_type in [*LINEAGE_MARKERS, "Other"]:
            mask = labels == cell_type
            n = int(mask.sum())
            count_rows.append({**sample.to_dict(), "cell_type": cell_type, "cell_count": n,
                               "proportion": n / max(len(labels), 1)})
            if n:
                class_counts_total[cell_type] += np.asarray(retained[:, mask].sum(axis=1)).ravel()

        margin = np.sort(scores, axis=0)[-1] - np.sort(scores, axis=0)[-2]
        qc_rows.append({**sample.to_dict(), "raw_barcodes": shape[1], "cells_retained": int(keep.sum()),
                        "retained_fraction_of_barcodes": float(keep.mean()),
                        "median_UMI_retained": float(np.median(totals[keep])),
                        "median_genes_retained": float(np.median(features[keep])),
                        "median_mito_fraction_retained": float(np.median(mito[keep])),
                        "cells_marker_ambiguous_other": int((labels == "Other").sum()),
                        "median_top_marker_margin": float(np.median(margin))})
        for label_idx, label in enumerate(LINEAGE_MARKERS):
            assigned = labels == label
            validation_rows.append({**sample.to_dict(), "assigned_cell_type": label,
                                    "n_cells": int(assigned.sum()),
                                    "mean_own_marker_score": float(np.mean(scores[label_idx, assigned])) if assigned.any() else np.nan,
                                    "mean_highest_competing_score": float(np.mean(np.max(np.delete(scores, label_idx, axis=0)[:, assigned], axis=0))) if assigned.any() else np.nan})
        print(f"{sample.sample_id}: retained {keep.sum():,}; macrophage {macrophage.sum():,}", flush=True)

    assert all_counts is not None and mac_counts is not None and names is not None and gene_ids is not None
    qc = pd.DataFrame(qc_rows)
    counts = pd.DataFrame(count_rows)
    props = counts.copy()
    validation = pd.DataFrame(validation_rows)
    qc.to_csv(RESULTS / "qc_by_library.csv", index=False)
    counts.to_csv(RESULTS / "celltype_counts_and_proportions_by_library.csv", index=False)
    validation.to_csv(VALIDATION / "marker_score_separation_by_library.csv", index=False)

    # Compact matrix exports preserve the exact pseudobulk inputs used downstream.
    all_frame = pd.DataFrame(all_counts, columns=meta.sample_id)
    all_frame.insert(0, "symbol", names)
    all_frame.insert(0, "gene_id", gene_ids)
    all_frame.to_csv(RESULTS / "pseudobulk_all_retained_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    mac_frame = pd.DataFrame(mac_counts, columns=meta.sample_id)
    mac_frame.insert(0, "symbol", names)
    mac_frame.insert(0, "gene_id", gene_ids)
    mac_frame.to_csv(RESULTS / "pseudobulk_macrophage_counts.tsv.gz", sep="\t", index=False, compression="gzip")

    all_log = logcpm(all_counts)
    mac_log = logcpm(mac_counts)
    all_pca = pca_table(all_log, meta, "all")
    mac_pca = pca_table(mac_log, meta, "mac")
    all_pca.to_csv(RESULTS / "pca_all_retained_pseudobulk.csv", index=False)
    mac_pca.to_csv(RESULTS / "pca_macrophage_pseudobulk.csv", index=False)

    de7 = de_contrast(mac_counts, gene_ids, names, meta, 7)
    de14 = de_contrast(mac_counts, gene_ids, names, meta, 14)
    de7.to_csv(RESULTS / "macrophage_DE_NINT_vs_BLM_day7.csv", index=False)
    de14.to_csv(RESULTS / "macrophage_DE_NINT_vs_BLM_day14.csv", index=False)

    scores, module_contrasts = module_analysis(mac_log, names, meta)
    scores.to_csv(RESULTS / "macrophage_module_scores_by_library.csv", index=False)
    module_contrasts.to_csv(RESULTS / "macrophage_module_contrasts.csv", index=False)
    comp_contrasts = composition_analysis(props)
    comp_contrasts.to_csv(RESULTS / "celltype_proportion_contrasts.csv", index=False)
    candidates = candidate_table(all_counts, mac_counts, names, meta)
    candidates.to_csv(RESULTS / "candidate_and_phenotype_tracks.csv", index=False)
    candidate_stats = candidate_contrasts(candidates)
    candidate_stats.to_csv(RESULTS / "candidate_gene_contrasts.csv", index=False)

    # Aggregated marker expression provides a transparent annotation sanity check.
    marker_rows = []
    for cell_type, vec in class_counts_total.items():
        lib = max(int(vec.sum()), 1)
        for marker_class, genes in LINEAGE_MARKERS.items():
            idx = gene_indices(names, genes)
            cpm = float(vec[idx].sum() / lib * 1e6) if idx else 0.0
            marker_rows.append({"assigned_cell_type": cell_type, "marker_set": marker_class,
                                "summed_marker_CPM": cpm, "n_marker_features": len(idx)})
    pd.DataFrame(marker_rows).to_csv(VALIDATION / "aggregate_marker_expression_by_celltype.csv", index=False)

    plot_s6(meta, qc, props, all_pca, mac_pca, de7, de14, scores, module_contrasts, candidates)

    # Key checks and a machine-readable manifest.
    expected = {(d, g): 3 for d in (7, 14) for g in ("Control", "BLM", "BLM+NINT")}
    observed = meta.groupby(["day", "group"]).size().to_dict()
    checks = {
        "sample_count": int(len(meta)),
        "all_six_condition_cells_have_three_libraries": observed == expected,
        "all_libraries_retain_at_least_1000_cells": bool((qc.cells_retained >= 1000).all()),
        "all_libraries_retain_at_least_100_macrophages": bool((counts.query("cell_type == 'Macrophage'").cell_count >= 100).all()),
        "no_cell_level_inferential_tests": True,
        "exact_permutation_minimum_two_sided_p_for_3v3": 0.1,
        "gene_level_BH_present_day7": bool(de7.q_BH.notna().any()),
        "gene_level_BH_present_day14": bool(de14.q_BH.notna().any()),
        "module_BH_present": bool(module_contrasts.q_BH_welch_10_tests.notna().all()),
    }
    (VALIDATION / "validation_checks.json").write_text(json.dumps(checks, indent=2) + "\n")

    # analysis_manifest.json and its detached validator are intentionally excluded
    # from the manifest's output list to avoid circular/stale self-checksums.
    output_paths = sorted([
        *[p for p in RESULTS.glob("*") if p.name != "analysis_manifest.json"],
        *[p for p in VALIDATION.glob("*") if p.name != "manifest_hash_validation.json"],
        *FIGURES.glob("Supplementary_Figure_S6_GSE151374.*"),
    ])
    input_paths = sorted([
        *files.values(),
        RAW / "GSE151374_RAW.tar",
        ROOT / "metadata" / "GSE151374_family.soft.gz",
        Path(__file__).resolve(),
    ])
    manifest = {
        "dataset": "GSE151374 / PRJNA635636",
        "analysis_timestamp_utc": pd.Timestamp.utcnow().isoformat(),
        "random_seed": SEED,
        "statistical_unit": "pooled 10x library (three mice pooled per library; n=3 libraries per condition/time point)",
        "qc_rule": "UMI >= 500; detected genes >= 200 and <= 6000; mitochondrial fraction < 0.20",
        "annotation": "broad marker-rule label; max mean log-normalized marker score; Other if top score <0.10 or margin <0.05",
        "gene_level_test": "macrophage pseudobulk log2(CPM+prior), two-sided Welch t; BH within each day; exploratory because n=3 vs n=3",
        "module_and_composition_test": "pooled-library Welch t plus exhaustive two-sided 3-vs-3 permutation; BH across displayed families",
        "software": {
            "python": sys.version.split()[0], "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
            "sklearn": __import__("sklearn").__version__, "matplotlib": mpl.__version__, "h5py": h5py.__version__,
        },
        "inputs": [{"path": str(p.relative_to(ROOT.parents[2])), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in input_paths],
        "outputs": [{"path": str(p.relative_to(ROOT.parents[2])), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in output_paths if p.is_file()],
    }
    manifest_path = RESULTS / "analysis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    entries = manifest["inputs"] + manifest["outputs"]
    verification = []
    for item in entries:
        path = ROOT.parents[2] / item["path"]
        verification.append({"path": item["path"], "bytes_match": path.stat().st_size == item["bytes"],
                             "sha256_match": sha256(path) == item["sha256"]})
    detached = {"manifest_path": str(manifest_path.relative_to(ROOT.parents[2])),
                "manifest_sha256": sha256(manifest_path), "entries_checked": len(verification),
                "all_recorded_hashes_match": all(x["bytes_match"] and x["sha256_match"] for x in verification),
                "entries": verification}
    (VALIDATION / "manifest_hash_validation.json").write_text(json.dumps(detached, indent=2) + "\n")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
