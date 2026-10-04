#!/usr/bin/env python3
"""Unconditioned animal-level reanalysis of PXD052594 processed proteomics.

All nintedanib-treated animals are compared with all vehicle animals present in
the public processed matrix.  Post-treatment radiomic clusters are deliberately
not used for selection, stratification, adjustment, or exclusion.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]
RAW, RESULTS, VALIDATION = ROOT / "raw", ROOT / "results", ROOT / "validation"
FIGURES = ROOT.parents[1] / "figures"
for directory in (RESULTS, VALIDATION, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)
SEED = 52594
RNG = np.random.default_rng(SEED)

MODULES = {
    "ECM remodeling": ["Col1a1", "Col1a2", "Col3a1", "Col5a1", "Fn1", "Postn", "Tnc", "Sparc", "Lox", "Ltbp2", "Ltbp4", "Mmp2", "Mmp9", "Timp1", "Tgfb1"],
    "Macrophage repair": ["Mrc1", "Cd68", "Adgre1", "Csf1r", "Marco", "Mertk", "Trem2", "Gpnmb", "Fabp5", "Lgals3"],
    "Inflammatory": ["Spp1", "Lgals3", "S100a8", "S100a9", "Ccl2", "Ccr2", "Nfkb1", "Nfkbia", "Stat3", "Tnf", "Il1b", "Nlrp3"],
    "Alveolar homeostasis": ["Pparg", "Car4", "Marco", "Siglecf", "Fabp4", "Itgax"],
    "RTK-response": ["Fgfr1", "Fgfr2", "Fgfr3", "Pdgfra", "Pdgfrb", "Kdr", "Flt1", "Flt4", "Mapk1", "Mapk3", "Akt1", "Stat3"],
}
CANDIDATES = ["Fgfr1", "Fbp1", "Gsta1", "Otc", "Spp1", "Mrc1", "Tgfb1", "Marco", "Col1a1", "Fn1"]
COLORS = {"vehicle": "#D55E00", "nintedanib": "#0072B2"}


def track_role(gene: str) -> str:
    if gene == "Fgfr1":
        return "exposure_proximal_target"
    if gene == "Fbp1":
        return "DILI_phenotype_anchor_P1"
    if gene in {"Gsta1", "Otc"}:
        return "DILI_phenotype_anchor_P2"
    return "lung_phenotype_context"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def bh(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    ranked = p[order]
    qrank = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    qrank = np.minimum(qrank, 1.0)
    q = np.empty_like(qrank)
    q[order] = qrank
    return q


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    df = len(a) + len(b) - 2
    sp2 = ((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / df
    if sp2 <= 0 or not np.isfinite(sp2):
        return np.nan
    return float((np.mean(a) - np.mean(b)) / np.sqrt(sp2) * (1 - 3 / (4 * df - 1)))


def mean_diff_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    a, b = np.asarray(a, float), np.asarray(b, float)
    diff = float(np.mean(a) - np.mean(b))
    va, vb = np.var(a, ddof=1) / len(a), np.var(b, ddof=1) / len(b)
    se = np.sqrt(va + vb)
    if se == 0:
        return diff, diff, diff
    df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
    crit = stats.t.ppf(0.975, df)
    return diff, float(diff - crit * se), float(diff + crit * se)


def monte_carlo_permutation(a: np.ndarray, b: np.ndarray, n_perm: int = 100_000) -> float:
    values = np.r_[a, b].astype(float)
    observed = abs(float(np.mean(a) - np.mean(b)))
    exceed = 0
    for _ in range(n_perm):
        perm = RNG.permutation(len(values))
        diff = abs(float(np.mean(values[perm[:len(a)]]) - np.mean(values[perm[len(a):]])))
        exceed += diff >= observed - 1e-12
    return (exceed + 1) / (n_perm + 1)


def module_scores(prot: pd.DataFrame, sample_cols: list[str], sample_meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = prot.set_index("Symbol")[sample_cols]
    expression = expression.groupby(level=0).mean()
    z = expression.sub(expression.mean(axis=1), axis=0).div(expression.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    score_rows, member_rows = [], []
    for module, genes in MODULES.items():
        used = [g for g in genes if g in z.index and np.isfinite(z.loc[g]).any()]
        for gene in genes:
            member_rows.append({"module": module, "gene": gene, "detected": gene in expression.index,
                                "used": gene in used})
        score = z.loc[used].mean(axis=0)
        for sample, value in score.items():
            score_rows.append({"sample_id": sample, "module": module, "score": float(value),
                               "n_proteins_used": len(used), "proteins_used": ";".join(used)})
    pd.DataFrame(member_rows).to_csv(RESULTS / "module_protein_membership.csv", index=False)
    scores = pd.DataFrame(score_rows).merge(sample_meta, on="sample_id", validate="many_to_one")
    contrasts = []
    for module, frame in scores.groupby("module", sort=False):
        a = frame.loc[frame.treatment == "nintedanib", "score"].to_numpy()
        b = frame.loc[frame.treatment == "vehicle", "score"].to_numpy()
        diff, lo, hi = mean_diff_ci(a, b)
        contrasts.append({"module": module, "difference_NINT_minus_vehicle": diff,
                          "CI95_low": lo, "CI95_high": hi, "hedges_g": hedges_g(a, b),
                          "p_welch": float(stats.ttest_ind(a, b, equal_var=False).pvalue),
                          "p_mannwhitney": float(stats.mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue),
                          "p_permutation_100k": monte_carlo_permutation(a, b),
                          "n_NINT": len(a), "n_vehicle": len(b)})
    contrasts = pd.DataFrame(contrasts)
    for pcol in ["p_welch", "p_mannwhitney", "p_permutation_100k"]:
        contrasts[f"q_BH_{pcol}"] = bh(contrasts[pcol].to_numpy())
    return scores, contrasts


def candidate_results(prot: pd.DataFrame, sample_cols: list[str], sample_meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_rows, contrast_rows = [], []
    for gene in CANDIDATES:
        match = prot.loc[prot.Symbol == gene]
        detected = not match.empty
        if detected:
            values = match[sample_cols].mean(axis=0)
            for sample, value in values.items():
                meta = sample_meta.set_index("sample_id").loc[sample]
                expr_rows.append({"gene": gene, "track_type": track_role(gene),
                                  "sample_id": sample, "treatment": meta.treatment,
                                  "animal_id": meta.animal_id, "log2_abundance": float(value), "detected": True})
            a = values[[s for s in sample_cols if s.startswith("nintedanib_")]].to_numpy(float)
            b = values[[s for s in sample_cols if s.startswith("vehicle_")]].to_numpy(float)
            diff, lo, hi = mean_diff_ci(a, b)
            contrast_rows.append({"gene": gene, "track_type": track_role(gene),
                                  "detected": True, "log2FC_NINT_vs_vehicle": diff,
                                  "CI95_low": lo, "CI95_high": hi, "hedges_g": hedges_g(a, b),
                                  "p_welch": float(stats.ttest_ind(a, b, equal_var=False).pvalue),
                                  "p_mannwhitney": float(stats.mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue),
                                  "n_NINT": len(a), "n_vehicle": len(b)})
        else:
            contrast_rows.append({"gene": gene, "track_type": track_role(gene),
                                  "detected": False, "log2FC_NINT_vs_vehicle": np.nan,
                                  "CI95_low": np.nan, "CI95_high": np.nan, "hedges_g": np.nan,
                                  "p_welch": np.nan, "p_mannwhitney": np.nan,
                                  "n_NINT": 10, "n_vehicle": 13})
    contrasts = pd.DataFrame(contrast_rows)
    valid = contrasts.p_welch.notna()
    contrasts.loc[valid, "q_BH_welch_detected_tracks"] = bh(contrasts.loc[valid, "p_welch"].to_numpy())
    contrasts.loc[valid, "q_BH_mannwhitney_detected_tracks"] = bh(contrasts.loc[valid, "p_mannwhitney"].to_numpy())
    return pd.DataFrame(expr_rows), contrasts


def plot_s7(prot: pd.DataFrame, sample_cols: list[str], sample_meta: pd.DataFrame, de: pd.DataFrame,
            pca: pd.DataFrame, scores: pd.DataFrame, mod: pd.DataFrame,
            cand_expr: pd.DataFrame, cand_stats: pd.DataFrame) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 8,
                         "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
                         "legend.fontsize": 6, "axes.linewidth": 0.6, "pdf.fonttype": 42,
                         "ps.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(14.2, 11.8))
    gs = fig.add_gridspec(4, 6, left=0.06, right=0.955, bottom=0.06, top=0.96, wspace=0.95, hspace=0.95)
    axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]),
            fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6]),
            fig.add_subplot(gs[2, 0:2]), fig.add_subplot(gs[2, 2:4]), fig.add_subplot(gs[2, 4:6]),
            fig.add_subplot(gs[3, 0:3]), fig.add_subplot(gs[3, 3:6])]
    for letter, ax in zip("ABCDEFGHIJ", axes):
        ax.text(-0.08, 1.10, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

    # A, public processed matrix completeness at animal level.
    ax = axes[0]
    observed = prot[sample_cols].notna().sum(axis=0)
    ax.bar(np.arange(len(sample_cols)), observed, color=[COLORS[t] for t in sample_meta.treatment])
    ax.set_xticks(np.arange(len(sample_cols)), sample_meta.animal_id, rotation=70, ha="right")
    ax.set_ylabel("Proteins quantified")
    ax.set_title("Animal-level matrix completeness")
    ax.spines[["top", "right"]].set_visible(False)

    # B, PCA.
    ax = axes[1]
    for treatment, frame in pca.groupby("treatment"):
        ax.scatter(frame.PC1, frame.PC2, s=36, color=COLORS[treatment], label=treatment, edgecolor="white", linewidth=0.4)
    ax.set_xlabel(f"PC1 ({pca.PC1_variance_pct.iloc[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.PC2_variance_pct.iloc[0]:.1f}%)")
    ax.set_title("Unconditioned lung-proteome PCA")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    # C, volcano.
    ax = axes[2]
    y = -np.log10(de.p_welch.clip(lower=1e-300))
    sig = de.q_BH < 0.05
    ax.scatter(de.log2FC_NINT_vs_vehicle[~sig], y[~sig], s=7, color="#BDBDBD", alpha=0.65, rasterized=True)
    ax.scatter(de.log2FC_NINT_vs_vehicle[sig], y[sig], s=11, color="#7B3294", alpha=0.85, rasterized=True)
    ax.axvline(0, lw=0.6, color="#555555")
    ax.axhline(-np.log10(0.05), lw=0.5, ls="--", color="#888888")
    # The all-protein BH screen is null.  Nominal top-P labels are omitted
    # because they are not a discovery threshold and overlap at print scale.
    ax.set_xlabel("log2FC, NINT vs vehicle")
    ax.set_ylabel("−log10 Welch P")
    ax.set_title("All-animal protein contrast")
    ax.spines[["top", "right"]].set_visible(False)

    # D, top-effect protein heatmap, standardized only for visualization.
    ax = axes[3]
    top = de.sort_values(["q_BH", "p_welch"]).head(20)
    matrix = prot.set_index("entrez_id").loc[top.entrez_id, sample_cols].to_numpy(float)
    matrix = (matrix - matrix.mean(axis=1, keepdims=True)) / matrix.std(axis=1, ddof=1, keepdims=True)
    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_yticks(np.arange(len(top)), top.Symbol)
    ax.set_xticks(np.arange(len(sample_cols)), sample_meta.animal_id, rotation=70, ha="right")
    ax.set_title("Top protein-level contrasts")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02); cb.set_label("Protein z score")

    # E, module effect forest.
    ax = axes[4]
    ordered = list(MODULES)
    frame = mod.set_index("module").reindex(ordered)
    yy = np.arange(len(frame))
    ax.errorbar(frame.difference_NINT_minus_vehicle, yy,
                xerr=[frame.difference_NINT_minus_vehicle - frame.CI95_low,
                      frame.CI95_high - frame.difference_NINT_minus_vehicle],
                fmt="o", color="#2166AC", capsize=3)
    ax.axvline(0, color="#555555", lw=0.6)
    ax.set_yticks(yy, frame.index)
    ax.set_xlabel("Score difference, NINT − vehicle (95% CI)")
    ax.set_title("Prespecified protein programs")
    ax.spines[["top", "right"]].set_visible(False)

    # F-H, selected module pool distributions, each dot is an animal.
    for ax, module in zip(axes[5:8], ["ECM remodeling", "Macrophage repair", "RTK-response"]):
        frame = scores[scores.module == module]
        for xpos, treatment in enumerate(["vehicle", "nintedanib"]):
            vals = frame.loc[frame.treatment == treatment, "score"].to_numpy()
            jitter = np.linspace(-0.10, 0.10, len(vals))
            ax.scatter(xpos + jitter, vals, color=COLORS[treatment], s=24, alpha=0.85, edgecolor="white", linewidth=0.35)
            ax.errorbar(xpos, vals.mean(), yerr=stats.sem(vals), fmt="_", color="#222222", capsize=2.5, markersize=10)
        ax.set_xticks([0, 1], ["vehicle", "NINT"])
        ax.set_ylabel("Module score")
        ax.set_title(module)
        ax.spines[["top", "right"]].set_visible(False)

    # I, candidate/phenotype protein effect estimates including explicit non-detection.
    ax = axes[8]
    frame = cand_stats.set_index("gene").reindex(CANDIDATES)
    y = np.arange(len(frame))
    detected = frame.detected.fillna(False).to_numpy(bool)
    ax.errorbar(frame.loc[detected, "log2FC_NINT_vs_vehicle"], y[detected],
                xerr=[frame.loc[detected, "log2FC_NINT_vs_vehicle"] - frame.loc[detected, "CI95_low"],
                      frame.loc[detected, "CI95_high"] - frame.loc[detected, "log2FC_NINT_vs_vehicle"]],
                fmt="o", color="#2166AC", capsize=2.5)
    ax.scatter(np.zeros((~detected).sum()), y[~detected], marker="x", color="#999999")
    for yi in y[~detected]:
        ax.text(0.02, yi, "not quantified", va="center", fontsize=5, color="#777777")
    ax.axvline(0, color="#555555", lw=0.6)
    ax.set_yticks(y, frame.index)
    ax.invert_yaxis()
    ax.set_xlabel("Protein log2FC, NINT vs vehicle (95% CI)")
    ax.set_title("FGFR1 target and phenotype anchors")
    ax.spines[["top", "right"]].set_visible(False)

    # J, absolute abundance heatmap; gray rows are absent from the proteome.
    ax = axes[9]
    mat = np.full((len(CANDIDATES), len(sample_cols)), np.nan)
    for i, gene in enumerate(CANDIDATES):
        row = prot.loc[prot.Symbol == gene, sample_cols]
        if not row.empty:
            mat[i] = row.mean(axis=0)
    cmap = mpl.colormaps["viridis"].copy(); cmap.set_bad("#D9D9D9")
    im = ax.imshow(mat, aspect="auto", cmap=cmap)
    ax.set_yticks(np.arange(len(CANDIDATES)), CANDIDATES)
    ax.set_xticks(np.arange(len(sample_cols)), sample_meta.animal_id, rotation=70, ha="right")
    ax.set_title("Absolute processed log2 protein abundance")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02); cb.set_label("log2 abundance")

    fig.suptitle("PXD052594: unconditioned animal-level lung-proteome response to nintedanib", fontsize=11, fontweight="bold", y=0.991)
    outbase = FIGURES / "Supplementary_Figure_S7_PXD052594"
    fig.savefig(outbase.with_suffix(".pdf"), dpi=600)
    fig.savefig(outbase.with_suffix(".svg"), dpi=600)
    fig.savefig(outbase.with_suffix(".png"), dpi=600)
    fig.savefig(outbase.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main() -> None:
    workbook = RAW / "Raw_data_values_Lauer_et_al_2024.xlsx"
    gene_info_path = RAW / "Mus_musculus.gene_info.gz"
    prot = pd.read_excel(workbook, sheet_name="proteomics_log2")
    prot["entrez_id"] = prot.entrez_id.astype(str)
    sample_cols = list(prot.columns[1:])
    gene_info = pd.read_csv(gene_info_path, sep="\t", dtype=str, low_memory=False)[["GeneID", "Symbol", "description"]]
    prot = prot.merge(gene_info, how="left", left_on="entrez_id", right_on="GeneID", validate="one_to_one")
    prot["Symbol"] = prot.Symbol.fillna("Entrez_" + prot.entrez_id)

    sample_meta = pd.DataFrame({"sample_id": sample_cols})
    sample_meta["treatment"] = sample_meta.sample_id.str.split("_").str[0]
    sample_meta["animal_id"] = "M" + sample_meta.sample_id.str.split("_").str[1]
    sample_meta["statistical_unit"] = "individual mouse"
    sample_meta.to_csv(RESULTS / "animal_sample_manifest.csv", index=False)

    nint_cols = [c for c in sample_cols if c.startswith("nintedanib_")]
    veh_cols = [c for c in sample_cols if c.startswith("vehicle_")]
    a, b = prot[nint_cols].to_numpy(float), prot[veh_cols].to_numpy(float)
    tests = stats.ttest_ind(a, b, axis=1, equal_var=False)
    logfc = a.mean(axis=1) - b.mean(axis=1)
    de = prot[["entrez_id", "Symbol", "description"]].copy()
    de["log2FC_NINT_vs_vehicle"] = logfc
    de["mean_log2_NINT"] = a.mean(axis=1)
    de["mean_log2_vehicle"] = b.mean(axis=1)
    de["hedges_g"] = [hedges_g(x, y) for x, y in zip(a, b)]
    de["p_welch"] = tests.pvalue
    de["q_BH"] = bh(de.p_welch.to_numpy())
    de["n_NINT"] = len(nint_cols); de["n_vehicle"] = len(veh_cols)
    de = de.sort_values(["q_BH", "p_welch"])
    de.to_csv(RESULTS / "proteome_NINT_vs_vehicle_all_animals.csv", index=False)

    # PCA on the 2,000 most variable proteins.
    matrix = prot[sample_cols].to_numpy(float)
    idx = np.argsort(np.var(matrix, axis=1))[-2000:]
    model = PCA(n_components=2, random_state=SEED)
    coords = model.fit_transform(matrix[idx].T)
    pca = sample_meta.copy()
    pca["PC1"], pca["PC2"] = coords[:, 0], coords[:, 1]
    pca["PC1_variance_pct"] = model.explained_variance_ratio_[0] * 100
    pca["PC2_variance_pct"] = model.explained_variance_ratio_[1] * 100
    pca.to_csv(RESULTS / "proteome_PCA.csv", index=False)

    scores, mod = module_scores(prot, sample_cols, sample_meta)
    scores.to_csv(RESULTS / "module_scores_by_animal.csv", index=False)
    mod.to_csv(RESULTS / "module_contrasts.csv", index=False)
    cand_expr, cand_stats = candidate_results(prot, sample_cols, sample_meta)
    cand_expr.to_csv(RESULTS / "candidate_expression_by_animal.csv", index=False)
    cand_stats.to_csv(RESULTS / "candidate_contrasts.csv", index=False)

    # The deposited phosphoproteome is a separate random 5+5 animal subset.
    # Analyze every deposited phosphosite without radiomic-cluster conditioning.
    phospho = pd.read_excel(workbook, sheet_name="phosphoproteomics")
    phospho_cols = [c for c in phospho.columns if c.startswith("iMaxLFQ_")]
    phospho_nint = [c for c in phospho_cols if "Nintedanib" in c]
    phospho_vehicle = [c for c in phospho_cols if "Vehicle" in c]
    pa = np.log2(phospho[phospho_nint].to_numpy(float))
    pb = np.log2(phospho[phospho_vehicle].to_numpy(float))
    ptest = stats.ttest_ind(pa, pb, axis=1, equal_var=False)
    phospho_de = phospho.iloc[:, :7].copy()
    phospho_de["log2FC_NINT_vs_vehicle"] = pa.mean(axis=1) - pb.mean(axis=1)
    phospho_de["mean_log2_NINT"] = pa.mean(axis=1)
    phospho_de["mean_log2_vehicle"] = pb.mean(axis=1)
    phospho_de["hedges_g"] = [hedges_g(x, y) for x, y in zip(pa, pb)]
    phospho_de["p_welch"] = ptest.pvalue
    phospho_de["q_BH"] = bh(phospho_de.p_welch.to_numpy())
    phospho_de["n_NINT"] = len(phospho_nint); phospho_de["n_vehicle"] = len(phospho_vehicle)
    phospho_de = phospho_de.sort_values(["q_BH", "p_welch"])
    phospho_de.to_csv(RESULTS / "phosphosite_NINT_vs_vehicle_all_deposited_subset_animals.csv", index=False)
    phospho_track_genes = ["Fgfr1", "Fgfr2", "Fgfr3", "Pdgfra", "Pdgfrb", "Kdr", "Flt1", "Flt4",
                           "Mapk1", "Mapk3", "Akt1", "Stat3", "Spp1", "Mrc1", "Tgfb1"]
    track_rows = []
    for gene in phospho_track_genes:
        frame = phospho_de.loc[phospho_de.Gene.fillna("").str.lower() == gene.lower()]
        track_rows.append({"gene": gene, "n_sites_quantified": len(frame),
                           "n_sites_q_lt_0_05": int((frame.q_BH < 0.05).sum()),
                           "minimum_q_BH": float(frame.q_BH.min()) if len(frame) else np.nan,
                           "site_at_minimum_q": frame.iloc[0].ProtPhosphoLoc if len(frame) else "not quantified",
                           "log2FC_at_minimum_q": float(frame.iloc[0].log2FC_NINT_vs_vehicle) if len(frame) else np.nan})
    pd.DataFrame(track_rows).to_csv(RESULTS / "phosphosite_target_pathway_summary.csv", index=False)

    pmat = np.log2(phospho[phospho_cols].to_numpy(float))
    pidx = np.argsort(np.var(pmat, axis=1))[-2000:]
    pmodel = PCA(n_components=2, random_state=SEED)
    pcoords = pmodel.fit_transform(pmat[pidx].T)
    phospho_meta = pd.DataFrame({"sample_id": phospho_cols})
    phospho_meta["treatment"] = phospho_meta.sample_id.str.extract(r"iMaxLFQ_(Nintedanib|Vehicle)", expand=False).str.lower()
    phospho_meta["animal_id"] = "M" + phospho_meta.sample_id.str.rsplit("_", n=1).str[-1]
    phospho_meta["PC1"], phospho_meta["PC2"] = pcoords[:, 0], pcoords[:, 1]
    phospho_meta["PC1_variance_pct"] = pmodel.explained_variance_ratio_[0] * 100
    phospho_meta["PC2_variance_pct"] = pmodel.explained_variance_ratio_[1] * 100
    phospho_meta.to_csv(RESULTS / "phosphoproteome_PCA_and_sample_manifest.csv", index=False)
    plot_s7(prot, sample_cols, sample_meta, de, pca, scores, mod, cand_expr, cand_stats)

    checks = {
        "proteins": int(len(prot)), "animals": len(sample_cols),
        "nintedanib_animals": len(nint_cols), "vehicle_animals": len(veh_cols),
        "matrix_has_no_missing_values": bool(prot[sample_cols].notna().all().all()),
        "animal_ids_unique": bool(sample_meta.animal_id.is_unique),
        "post_treatment_cluster_used_in_analysis": False,
        "one_vehicle_excluded_upstream_for_sample_preparation": "M29 (reported by source; not present in public proteome matrix)",
        "all_remaining_animals_included": True,
        "BH_applied_across_all_proteins": bool(de.q_BH.notna().all()),
        "phosphosites": int(len(phospho_de)),
        "phosphoproteome_nintedanib_animals": len(phospho_nint),
        "phosphoproteome_vehicle_animals": len(phospho_vehicle),
        "phosphoproteome_has_no_missing_values": bool(phospho[phospho_cols].notna().all().all()),
        "phosphosite_BH_applied": bool(phospho_de.q_BH.notna().all()),
    }
    (VALIDATION / "validation_checks.json").write_text(json.dumps(checks, indent=2) + "\n")

    paths = sorted([
        *[p for p in RESULTS.glob("*") if p.name != "analysis_manifest.json"],
        *[p for p in VALIDATION.glob("*") if p.name != "manifest_hash_validation.json"],
        *FIGURES.glob("Supplementary_Figure_S7_PXD052594.*"),
    ])
    manifest = {
        "dataset": "PXD052594; Zenodo 10.5281/zenodo.11395642",
        "analysis_timestamp_utc": pd.Timestamp.utcnow().isoformat(),
        "statistical_unit": "individual mouse",
        "contrast": "all nintedanib (n=10) versus all vehicle animals in processed proteome matrix (n=13)",
        "selection_guardrail": "No radiomic response cluster was used for conditioning, exclusion, stratification, or adjustment.",
        "protein_test": "two-sided Welch t test on deposited processed log2 abundance; BH across 7,006 proteins",
        "phosphosite_test": "two-sided Welch t test on log2-transformed deposited iMaxLFQ values in the separate random 5+5 subset; BH across 20,043 sites",
        "module_test": "Welch t, exact Mann-Whitney, and seeded 100,000-label Monte Carlo permutation; BH across five prespecified modules",
        "software": {"python": sys.version.split()[0], "platform": platform.platform(), "numpy": np.__version__,
                     "pandas": pd.__version__, "scipy": scipy.__version__,
                     "sklearn": __import__("sklearn").__version__, "matplotlib": mpl.__version__},
        "inputs": [{"path": str(p.relative_to(ROOT.parents[2])), "bytes": p.stat().st_size, "sha256": sha256(p)}
                   for p in [workbook, gene_info_path, RAW / "R_Script_Lauer_et_al_2024.R", Path(__file__).resolve()]],
        "outputs": [{"path": str(p.relative_to(ROOT.parents[2])), "bytes": p.stat().st_size, "sha256": sha256(p)}
                    for p in paths if p.is_file()],
    }
    manifest_path = RESULTS / "analysis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    verification = []
    for item in manifest["inputs"] + manifest["outputs"]:
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
