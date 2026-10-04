#!/usr/bin/env python3
"""Concentration-resolved nintedanib perturbation analysis and Figure 3.

This analysis uses only the freshly extracted LINCS L1000 HepG2 Level 5
signatures (six 24-h concentrations), the separately prepared MCF7 Sci-Plex
signatures, and high-confidence human single-protein ChEMBL activities.  It is
explicitly an in-vitro perturbation analysis.  It neither labels the HepG2
response as clinical DILI nor treats transcriptional response of a biochemical
target as evidence of toxicity causality.

Statistical guardrails
----------------------
* There is one Level 5 consensus signature at each concentration.  Spearman
  ordering is therefore tested using all 6! exact label permutations, not an
  asymptotic p value that can spuriously approach zero with n=6.
* BH correction is applied across all 12,328 genes.  Sensitivity q values are
  also reported within the measured landmark and inferred families.
* Competitive gene-set tests match the three L1000 measurement strata
  (landmark, best-inferred/BING, other inferred), because their score
  distributions are materially different.
* Sci-Plex/Harmonizome signatures contain only the reported top 500 up- and
  top 500 down-regulated genes.  Cross-cell-line correlations are conditional
  on that truncated overlap and are not genome-wide concordance estimates.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from statsmodels.stats.proportion import proportion_confint


ROOT = Path(__file__).resolve().parents[2]
LINCS = ROOT / "fresh_analysis" / "lincs" / "output"
PHARM = ROOT / "fresh_analysis" / "pharmacology"
OUT = ROOT / "fresh_analysis" / "lincs" / "analysis"
FIG = ROOT / "fresh_analysis" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SEED = 260904
N_PERM = 25_000
DOSES = np.array([0.04, 0.12, 0.37, 1.11, 3.33, 10.0], dtype=float)
LOG_DOSES = np.log10(DOSES)
DOSE_COLUMNS = [f"z_{str(x).replace('.', 'p')}uM" for x in DOSES]

DILI_ANCHORS = [
    "ACO1", "ASS1", "FAH", "CPS1", "ALDOB", "HPD", "OTC", "DMGDH",
    "GSTA1", "FBP1", "PCK2", "CES1", "LECT2",
]

# Sentinel modules were fixed before formal testing.  They are deliberately
# compact, interpretable mechanism panels rather than post-hoc pathway results.
MECHANISM_MODULES = {
    "Oxidative/redox defense": [
        "NFE2L2", "KEAP1", "HMOX1", "NQO1", "GCLC", "GCLM", "TXNIP",
        "TXNRD1", "GPX4", "SLC7A11", "SOD2", "CAT",
    ],
    "Mitochondrial energetics": [
        "PPARGC1A", "TFAM", "OPA1", "MFN1", "MFN2", "DNM1L", "SLC25A4",
        "CPT1A", "ACADM", "ACADVL", "PCK2", "ACO2", "IDH2", "SDHA",
        "NDUFS1", "ATP5F1A",
    ],
    "Bile-acid transport": [
        "ABCB11", "ABCC2", "ABCB4", "ATP8B1", "SLC10A1", "UGT1A1",
        "NR1H4", "CYP7A1", "CYP8B1", "BAAT", "SLC51A", "SLC51B",
    ],
    "Xenobiotic handling": [
        "CYP3A4", "CYP2C9", "CYP2D6", "CES1", "CES2", "GSTP1", "GSTA1",
        "UGT1A1", "UGT2B7", "NQO1", "EPHX1", "POR", "ABCB1", "ABCC3",
    ],
    "ER stress/apoptosis": [
        "BAX", "BAK1", "BCL2", "CASP3", "CASP8", "CASP9", "FAS", "DDIT3",
        "ATF4", "XBP1", "HSPA5", "EIF2AK3", "ERN1",
    ],
    "Hepatocyte integrity": [
        "KRT8", "KRT18", "KRT19", "VIM", "JUP", "DSP", "TJP1", "OCLN",
        "CLDN1", "ALB", "HNF4A",
    ],
    "RTK–MAPK–AKT signaling": [
        "KDR", "FLT1", "FLT4", "FGFR1", "FGFR2", "FGFR3", "PDGFRA",
        "PDGFRB", "PLCG1", "MAPK1", "MAPK3", "AKT1", "STAT3", "SRC",
    ],
}

CANONICAL_TARGETS = {
    "KDR", "FLT1", "FLT4", "FGFR1", "FGFR2", "FGFR3", "PDGFRA", "PDGFRB"
}

COLORS = {
    "navy": "#1F4E79",
    "blue": "#4C78A8",
    "sky": "#72B7B2",
    "orange": "#E07B39",
    "vermillion": "#D1495B",
    "purple": "#7A5195",
    "gray": "#8C8C8C",
    "lightgray": "#D9D9D9",
    "dark": "#2F2F2F",
    "green": "#3A7D44",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_adjust(p_values: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values, preserving NaNs."""
    p = np.asarray(list(p_values), dtype=float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    pf = p[finite]
    if not len(pf):
        return out
    order = np.argsort(pf)
    ranked = pf[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    out[finite] = restored
    return out


PERMUTATION_INDICES = np.asarray(list(itertools.permutations(range(6))), dtype=int)
XRANK = stats.rankdata(LOG_DOSES)
XRANK_CENTERED = XRANK - XRANK.mean()
XRANK_NORM = np.sqrt(np.sum(XRANK_CENTERED**2))


def exact_spearman(y: np.ndarray) -> tuple[float, float, bool]:
    """Return rho and an exact two-sided 6! permutation p value."""
    y = np.asarray(y, dtype=float)
    if np.allclose(y, y[0]):
        return 0.0, 1.0, True
    yrank = stats.rankdata(y)
    yc = yrank - yrank.mean()
    denom = XRANK_NORM * np.sqrt(np.sum(yc**2))
    rho = float(np.dot(XRANK_CENTERED, yc) / denom)
    permuted = yrank[PERMUTATION_INDICES]
    pc = permuted - permuted.mean(axis=1, keepdims=True)
    perm_rho = (pc @ XRANK_CENTERED) / (
        np.sqrt(np.sum(pc**2, axis=1)) * XRANK_NORM
    )
    p = float(np.mean(np.abs(perm_rho) >= abs(rho) - 1e-12))
    return rho, p, False


def load_gene_profiles() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = LINCS / "nintedanib_hepg2_level5_all_genes_tidy.csv.gz"
    long = pd.read_csv(path)
    expected = set(DOSES)
    observed = set(long["dose_um"].unique())
    if observed != expected:
        raise RuntimeError(f"Expected doses {sorted(expected)}, observed {sorted(observed)}")
    key = ["gene_id", "gene_symbol", "gene_title", "is_landmark", "is_bing"]
    wide = long.pivot(index=key, columns="dose_um", values="level5_zscore").reset_index()
    if len(wide) != 12_328 or wide[DOSES].isna().any().any():
        raise RuntimeError("Unexpected gene count or missing HepG2 Level 5 scores")

    z = wide[DOSES].to_numpy(dtype=float)
    trend = [exact_spearman(row) for row in z]
    wide["spearman_rho"] = [x[0] for x in trend]
    wide["spearman_exact_p"] = [x[1] for x in trend]
    wide["constant_profile"] = [x[2] for x in trend]
    wide["bh_q_all_genes"] = bh_adjust(wide["spearman_exact_p"])
    wide["bh_q_measurement_family"] = np.nan
    for _, idx in wide.groupby(["is_landmark", "is_bing"], sort=False).groups.items():
        wide.loc[idx, "bh_q_measurement_family"] = bh_adjust(
            wide.loc[idx, "spearman_exact_p"]
        )

    raw_auc = np.trapezoid(z, LOG_DOSES, axis=1)
    wide["auc_log10dose"] = raw_auc
    wide["auc_log10dose_normalized"] = raw_auc / (LOG_DOSES[-1] - LOG_DOSES[0])
    peak_index = np.argmax(np.abs(z), axis=1)
    wide["max_abs_z"] = np.max(np.abs(z), axis=1)
    wide["peak_dose_um"] = DOSES[peak_index]
    wide["peak_z"] = z[np.arange(len(z)), peak_index]
    dominant = np.sign(wide["auc_log10dose_normalized"].to_numpy())
    fallback = np.sign(wide["spearman_rho"].to_numpy())
    dominant = np.where(dominant == 0, fallback, dominant)
    wide["direction_consistency"] = np.mean(np.sign(z) == dominant[:, None], axis=1)
    wide["trend_magnitude"] = np.abs(wide["spearman_rho"]) * wide["max_abs_z"]
    wide["trend_direction"] = np.select(
        [wide["spearman_rho"].gt(0), wide["spearman_rho"].lt(0)],
        ["increasing", "decreasing"], default="flat",
    )

    loo = np.empty((len(wide), 6), dtype=float)
    for omit in range(6):
        keep = np.arange(6) != omit
        x = stats.rankdata(LOG_DOSES[keep])
        xc = x - x.mean()
        xnorm = np.sqrt(np.sum(xc**2))
        for i, row in enumerate(z[:, keep]):
            if np.allclose(row, row[0]):
                loo[i, omit] = 0.0
            else:
                yr = stats.rankdata(row)
                yc = yr - yr.mean()
                loo[i, omit] = np.dot(xc, yc) / (xnorm * np.sqrt(np.sum(yc**2)))
    wide["loo_rho_median"] = np.median(loo, axis=1)
    full_sign = np.sign(wide["spearman_rho"].to_numpy())
    wide["loo_sign_stability"] = np.mean(np.sign(loo) == full_sign[:, None], axis=1)
    wide["loo_min_abs_rho"] = np.min(np.abs(loo), axis=1)

    wide["measurement_stratum"] = np.select(
        [wide["is_landmark"].eq(1), wide["is_bing"].eq(1)],
        ["measured_landmark", "best_inferred"],
        default="other_inferred",
    )
    wide["stratum_response_percentile"] = (
        wide.groupby("measurement_stratum")["trend_magnitude"]
        .rank(method="average", pct=True)
    )
    wide["stratum_auc_percentile_centered"] = (
        wide.groupby("measurement_stratum")["auc_log10dose_normalized"]
        .rank(method="average", pct=True) - 0.5
    )
    for dose, name in zip(DOSES, DOSE_COLUMNS):
        wide[name] = wide[dose]
    wide = wide.drop(columns=list(DOSES))
    wide = wide.sort_values(["trend_magnitude", "max_abs_z"], ascending=False)
    wide.to_csv(OUT / "gene_dose_response_metrics.csv.gz", index=False)
    return long, wide


def competitive_test(
    metrics: pd.DataFrame,
    genes: list[str],
    label: str,
    rng: np.random.Generator,
    n_perm: int = N_PERM,
) -> dict:
    """Measurement-stratum-matched competitive random-gene-set test."""
    unique = list(dict.fromkeys(genes))
    member = metrics.loc[metrics["gene_symbol"].isin(unique)].copy()
    missing = sorted(set(unique) - set(member["gene_symbol"]))
    if not len(member):
        raise RuntimeError(f"No genes from set {label!r} were present")
    observed = float(member["stratum_response_percentile"].mean())
    observed_signed = float(member["stratum_auc_percentile_centered"].mean())
    null = np.zeros(n_perm, dtype=float)
    null_signed = np.zeros(n_perm, dtype=float)
    total = len(member)
    strata_counts = member["measurement_stratum"].value_counts()
    for stratum, n in strata_counts.items():
        pool = metrics.loc[metrics["measurement_stratum"].eq(stratum)]
        mag = pool["stratum_response_percentile"].to_numpy(dtype=float)
        signed = pool["stratum_auc_percentile_centered"].to_numpy(dtype=float)
        # Each random set is sampled without replacement within stratum.
        for b in range(n_perm):
            indices = rng.choice(len(pool), size=int(n), replace=False)
            weight = n / total
            null[b] += weight * float(mag[indices].mean())
            null_signed[b] += weight * float(signed[indices].mean())
    sd = float(null.std(ddof=1))
    competitive_z = (observed - float(null.mean())) / sd if sd > 0 else np.nan
    p_one_sided = (1 + int(np.sum(null >= observed))) / (n_perm + 1)
    signed_p = (1 + int(np.sum(np.abs(null_signed) >= abs(observed_signed)))) / (n_perm + 1)
    return {
        "set_name": label,
        "n_requested": len(unique),
        "n_present": len(member),
        "n_measured_landmark": int(member["is_landmark"].sum()),
        "n_best_inferred": int(((member["is_landmark"] == 0) & (member["is_bing"] == 1)).sum()),
        "n_other_inferred": int(((member["is_landmark"] == 0) & (member["is_bing"] == 0)).sum()),
        "missing_genes": ";".join(missing),
        "mean_stratum_response_percentile": observed,
        "null_mean": float(null.mean()),
        "null_sd": sd,
        "competitive_z": competitive_z,
        "permutation_p_one_sided": p_one_sided,
        "mean_signed_auc_percentile_centered": observed_signed,
        "signed_auc_permutation_p_two_sided": signed_p,
        "median_spearman_rho": float(member["spearman_rho"].median()),
        "median_auc_log10dose_normalized": float(member["auc_log10dose_normalized"].median()),
        "median_max_abs_z": float(member["max_abs_z"].median()),
        "permutations": n_perm,
    }


def analyze_gene_sets(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    anchor_stat = pd.DataFrame([
        competitive_test(metrics, DILI_ANCHORS, "DILI proteomic anchors", rng)
    ])
    anchor_stat["bh_q"] = anchor_stat["permutation_p_one_sided"]
    anchor_stat.to_csv(OUT / "dili_anchor_set_statistics.csv", index=False)

    anchor_metrics = pd.DataFrame({"gene_symbol": DILI_ANCHORS}).merge(
        metrics, on="gene_symbol", how="left", validate="one_to_one"
    )
    anchor_metrics["present_in_l1000"] = anchor_metrics["gene_id"].notna()
    represented = anchor_metrics["present_in_l1000"]
    anchor_metrics.loc[represented, "bh_q_within_anchor_set"] = bh_adjust(
        anchor_metrics.loc[represented, "spearman_exact_p"]
    )
    anchor_metrics.to_csv(OUT / "dili_anchor_metrics.csv", index=False)

    rows = []
    member_rows = []
    for module, genes in MECHANISM_MODULES.items():
        rows.append(competitive_test(metrics, genes, module, rng))
        for gene in genes:
            member_rows.append({
                "module": module,
                "gene_symbol": gene,
                "definition_basis": "a_priori_mechanistic_sentinel_panel",
                "present_in_l1000": gene in set(metrics["gene_symbol"]),
            })
    module = pd.DataFrame(rows)
    module["bh_q_magnitude_family"] = bh_adjust(module["permutation_p_one_sided"])
    module["bh_q_signed_auc_family"] = bh_adjust(
        module["signed_auc_permutation_p_two_sided"]
    )
    module.to_csv(OUT / "mechanism_module_statistics.csv", index=False)
    pd.DataFrame(member_rows).to_csv(OUT / "mechanism_module_members.csv", index=False)
    return anchor_stat, anchor_metrics, module


def analyze_sciplex(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sci = pd.read_csv(PHARM / "sciplex_nintedanib_signature_long.csv")
    if sci.groupby("dose").size().to_dict() != {"10uM": 1000, "1uM": 1000}:
        raise RuntimeError("Expected a top-1000 Sci-Plex signature at each dose")
    metric_cols = [
        "gene_symbol", "is_landmark", "is_bing", "spearman_rho",
        "auc_log10dose_normalized", "max_abs_z", "trend_magnitude",
    ] + DOSE_COLUMNS
    base = metrics[metric_cols]
    rows = []
    pairs = []
    pair_specs = [(1.11, "1uM"), (1.11, "10uM"), (10.0, "1uM"), (10.0, "10uM")]
    dose_to_col = {d: c for d, c in zip(DOSES, DOSE_COLUMNS)}
    for hep_dose, sci_dose in pair_specs:
        sci_part = sci.loc[sci["dose"].eq(sci_dose), ["gene_symbol", "z", "direction"]]
        merged = base.merge(sci_part, on="gene_symbol", how="inner", validate="one_to_one")
        merged = merged.rename(columns={dose_to_col[hep_dose]: "hepg2_level5_z", "z": "sciplex_z"})
        merged["hepg2_dose_um"] = hep_dose
        merged["sciplex_dose"] = sci_dose
        merged["dose_pair"] = f"HepG2 {hep_dose:g} µM vs MCF7 {sci_dose}"
        merged["sign_concordant"] = (
            np.sign(merged["hepg2_level5_z"]) == np.sign(merged["sciplex_z"])
        )
        pairs.append(merged)
        for subset, selected in [
            ("all_overlap", merged),
            ("measured_landmark_only", merged.loc[merged["is_landmark"].eq(1)]),
        ]:
            rho, rho_p = stats.spearmanr(
                selected["hepg2_level5_z"], selected["sciplex_z"]
            )
            k = int(selected["sign_concordant"].sum())
            n = len(selected)
            lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
            sign_p = stats.binomtest(k, n, p=0.5, alternative="two-sided").pvalue
            rows.append({
                "dose_pair": merged["dose_pair"].iat[0],
                "hepg2_dose_um": hep_dose,
                "sciplex_dose": sci_dose,
                "subset": subset,
                "n_overlap": n,
                "spearman_rho": float(rho),
                "spearman_p": float(rho_p),
                "sign_concordant_n": k,
                "sign_concordance": k / n,
                "sign_concordance_ci_low": float(lo),
                "sign_concordance_ci_high": float(hi),
                "sign_binomial_p": float(sign_p),
                "sciplex_signature_truncated_top1000": True,
            })
    summary = pd.DataFrame(rows)
    for subset, idx in summary.groupby("subset").groups.items():
        summary.loc[idx, "spearman_bh_q_within_subset"] = bh_adjust(
            summary.loc[idx, "spearman_p"]
        )
        summary.loc[idx, "sign_bh_q_within_subset"] = bh_adjust(
            summary.loc[idx, "sign_binomial_p"]
        )
    pair_frame = pd.concat(pairs, ignore_index=True)
    summary.to_csv(OUT / "sciplex_hepg2_concordance.csv", index=False)
    pair_frame.to_csv(OUT / "sciplex_hepg2_gene_pairs.csv.gz", index=False)
    return summary, pair_frame


def analyze_chembl(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    chembl = pd.read_csv(PHARM / "chembl_target_summary.csv")
    cross = chembl.merge(metrics, on="gene_symbol", how="left", validate="one_to_one")
    cross["present_in_l1000"] = cross["gene_id"].notna()
    cross["is_canonical_nintedanib_target"] = cross["gene_symbol"].isin(CANONICAL_TARGETS)
    cross.to_csv(OUT / "chembl_transcript_crosswalk.csv", index=False)

    association_rows = []
    groups = {
        "canonical_targets": cross["is_canonical_nintedanib_target"],
        "multi_document_potent": cross["multi_document_potent"].eq(True),
        "all_potent_pchembl_ge_6": cross["max_pchembl"].ge(6),
        "all_mapped_targets": cross["max_pchembl"].notna(),
    }
    for label, mask in groups.items():
        part = cross.loc[mask & cross["spearman_rho"].notna()]
        for response in [
            "spearman_rho", "auc_log10dose_normalized", "max_abs_z", "trend_magnitude"
        ]:
            rho, p = stats.spearmanr(part["max_pchembl"], part[response])
            association_rows.append({
                "target_set": label,
                "response_metric": response,
                "n_targets": len(part),
                "spearman_rho": float(rho),
                "p_value": float(p),
            })
    association = pd.DataFrame(association_rows)
    for target_set, idx in association.groupby("target_set").groups.items():
        association.loc[idx, "bh_q_within_target_set"] = bh_adjust(
            association.loc[idx, "p_value"]
        )
    association.to_csv(OUT / "chembl_potency_association.csv", index=False)

    rng = np.random.default_rng(SEED + 99)
    set_rows = []
    target_gene_sets = {
        "Canonical nintedanib targets": sorted(CANONICAL_TARGETS),
        "Multi-document potent targets": sorted(
            cross.loc[cross["multi_document_potent"].eq(True), "gene_symbol"].dropna().unique()
        ),
        "All pChEMBL≥6 targets": sorted(
            cross.loc[cross["max_pchembl"].ge(6), "gene_symbol"].dropna().unique()
        ),
    }
    for name, genes in target_gene_sets.items():
        set_rows.append(competitive_test(metrics, genes, name, rng))
    target_sets = pd.DataFrame(set_rows)
    target_sets["bh_q_magnitude_family"] = bh_adjust(
        target_sets["permutation_p_one_sided"]
    )
    target_sets["bh_q_signed_auc_family"] = bh_adjust(
        target_sets["signed_auc_permutation_p_two_sided"]
    )
    target_sets.to_csv(OUT / "chembl_target_set_statistics.csv", index=False)
    return cross, association, target_sets


def make_integrated_evidence(
    metrics: pd.DataFrame, sci_pairs: pd.DataFrame, chembl: pd.DataFrame
) -> pd.DataFrame:
    selected = sorted(set(DILI_ANCHORS) | CANONICAL_TARGETS)
    base = pd.DataFrame({"gene_symbol": selected}).merge(
        metrics, on="gene_symbol", how="left", validate="one_to_one"
    )
    sci = (
        sci_pairs[["gene_symbol", "sciplex_dose", "sciplex_z"]]
        .drop_duplicates()
        .pivot(index="gene_symbol", columns="sciplex_dose", values="sciplex_z")
        .rename(columns={"1uM": "sciplex_mcf7_1uM_z", "10uM": "sciplex_mcf7_10uM_z"})
        .reset_index()
    )
    chem_cols = [
        "gene_symbol", "max_pchembl", "median_pchembl", "n_documents",
        "n_assays", "pharmacology_tier", "multi_document_potent",
    ]
    out = base.merge(sci, on="gene_symbol", how="left").merge(
        chembl[chem_cols], on="gene_symbol", how="left"
    )
    out["is_dili_anchor"] = out["gene_symbol"].isin(DILI_ANCHORS)
    out["is_canonical_nintedanib_target"] = out["gene_symbol"].isin(CANONICAL_TARGETS)
    out.to_csv(OUT / "integrated_gene_evidence.csv", index=False)
    return out


def style_axis(ax: mpl.axes.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=2.5, width=0.6, labelsize=6.4)
    ax.xaxis.label.set_size(6.8)
    ax.yaxis.label.set_size(6.8)
    ax.title.set_fontsize(7.6)
    ax.title.set_fontweight("bold")
    ax.title.set_horizontalalignment("left")
    ax.title.set_x(0.0)


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.16, 1.14, label, transform=ax.transAxes, fontsize=11,
        fontweight="bold", va="top", ha="left", color="#111111"
    )


def heatmap(
    ax: mpl.axes.Axes,
    matrix: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    title: str,
    vlim: float | None = None,
    annotate: bool = False,
) -> mpl.image.AxesImage:
    arr = np.ma.masked_invalid(np.asarray(matrix, dtype=float))
    if vlim is None:
        finite = np.abs(arr.compressed())
        vlim = float(np.quantile(finite, 0.97)) if len(finite) else 1.0
    vlim = max(vlim, 0.5)
    cmap = mpl.colormaps["RdBu_r"].copy()
    cmap.set_bad("#ECECEC")
    image = ax.imshow(arr, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(0, -vlim, vlim))
    ax.set_xticks(np.arange(len(col_labels)), col_labels, rotation=45, ha="right", fontsize=6.2)
    ax.set_yticks(np.arange(len(row_labels)), row_labels, fontsize=6.2)
    ax.set_title(title, pad=5, loc="left", fontsize=7.6, fontweight="bold")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if annotate and len(row_labels) <= 16:
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                value = arr[i, j]
                if value is not np.ma.masked:
                    ax.text(j, i, f"{float(value):.1f}", ha="center", va="center", fontsize=5.9,
                            color="white" if abs(float(value)) > 0.58 * vlim else "#222222")
    return image


def add_small_colorbar(fig: mpl.figure.Figure, ax: mpl.axes.Axes, image, label: str) -> None:
    box = ax.get_position()
    cax = fig.add_axes([box.x1 + 0.006, box.y0 + box.height * 0.15, 0.008, box.height * 0.7])
    cb = fig.colorbar(image, cax=cax)
    cb.ax.tick_params(labelsize=6.0, length=2)
    cb.set_label(label, fontsize=6.2)
    cb.outline.set_linewidth(0.4)


def save_figure(fig: mpl.figure.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    png_tmp = stem.parent / f".{stem.name}.png.tmp"
    tiff_tmp = stem.parent / f".{stem.name}.tiff.tmp"
    fig.savefig(png_tmp, format="png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(
        tiff_tmp, format="tiff", dpi=600, bbox_inches="tight", facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    # Atomic replacement prevents a concurrent document-assembly process from
    # observing a partially written high-resolution raster.
    png_tmp.replace(stem.with_suffix(".png"))
    tiff_tmp.replace(stem.with_suffix(".tiff"))


def plot_main_figure(
    long: pd.DataFrame,
    metrics: pd.DataFrame,
    anchor_stat: pd.DataFrame,
    anchors: pd.DataFrame,
    modules: pd.DataFrame,
    sci_summary: pd.DataFrame,
    sci_pairs: pd.DataFrame,
    chembl: pd.DataFrame,
    associations: pd.DataFrame,
) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 6.5,
        "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.transparent": False,
    })
    fig = plt.figure(figsize=(7.25, 10.4), facecolor="white")
    gs = fig.add_gridspec(
        5, 2, left=0.105, right=0.955, bottom=0.045, top=0.945,
        hspace=0.78, wspace=0.54,
    )
    axes = [fig.add_subplot(gs[i, j]) for i in range(5) for j in range(2)]

    # A — full-data distribution summaries, split by measurement stratum.
    ax = axes[0]
    positions, values, box_colors = [], [], []
    for i, dose in enumerate(DOSES):
        for offset, landmark in [(-0.15, 1), (0.15, 0)]:
            positions.append(i + offset)
            values.append(long.loc[(long["dose_um"] == dose) & (long["is_landmark"] == landmark), "level5_zscore"])
            box_colors.append(COLORS["navy"] if landmark else COLORS["lightgray"])
    bp = ax.boxplot(
        values, positions=positions, widths=0.24, patch_artist=True, showfliers=False,
        medianprops={"color": "white", "linewidth": 0.8},
        whiskerprops={"linewidth": 0.5}, capprops={"linewidth": 0.5},
        boxprops={"linewidth": 0.5},
    )
    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_edgecolor(COLORS["dark"] if color == COLORS["lightgray"] else color)
    ax.axhline(0, color="#777777", lw=0.5, ls="--")
    ax.set_xticks(range(6), [f"{d:g}" for d in DOSES])
    ax.set_xlabel("Nintedanib (µM; 24 h)")
    ax.set_ylabel("L1000 Level 5 z score")
    ax.set_title("Global response distribution")
    ax.legend(
        handles=[
            mpl.patches.Patch(facecolor=COLORS["navy"], label="Measured landmark (n=978)"),
            mpl.patches.Patch(facecolor=COLORS["lightgray"], edgecolor=COLORS["dark"], label="Inferred (n=11,350)"),
        ], frameon=False, fontsize=6.1, loc="lower left",
    )
    style_axis(ax); panel_label(ax, "A")

    # B — exact trend scan; no claim of FDR-positive genes if none exist.
    ax = axes[1]
    inferred = metrics["is_landmark"].eq(0)
    ax.scatter(
        metrics.loc[inferred, "spearman_rho"],
        -np.log10(metrics.loc[inferred, "spearman_exact_p"]),
        s=3, color="#B8B8B8", alpha=0.35, edgecolors="none", rasterized=True,
    )
    ax.scatter(
        metrics.loc[~inferred, "spearman_rho"],
        -np.log10(metrics.loc[~inferred, "spearman_exact_p"]),
        s=6, color=COLORS["navy"], alpha=0.75, edgecolors="none", rasterized=True,
    )
    nominal_p_ref = -math.log10(0.05)
    ax.axhline(nominal_p_ref, color=COLORS["vermillion"], lw=0.7, ls="--")
    ax.text(0.98, nominal_p_ref, "nominal exact P=0.05", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=6.0, color=COLORS["vermillion"])
    ax.axvline(0, color="#777777", lw=0.5)
    q_hits = int(metrics["bh_q_all_genes"].lt(0.05).sum())
    q_min = metrics["bh_q_all_genes"].min()
    ax.text(0.03, 0.96, f"Exact 6! test; BH q<0.05: {q_hits}\nminimum q={q_min:.3f}",
            transform=ax.transAxes, va="top", fontsize=6.1,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.0})
    ax.set_xlabel("Spearman ρ across six doses")
    ax.set_ylabel("−log₁₀ exact P")
    ax.set_title("Genome-wide concentration ordering")
    style_axis(ax); panel_label(ax, "B")

    # C — strongest measured landmark trends, clearly labeled as ranked not significant.
    ax = axes[2]
    top_landmark = metrics.loc[metrics["is_landmark"].eq(1)].nlargest(12, "trend_magnitude")
    top_landmark = top_landmark.sort_values("spearman_rho")
    hm = heatmap(
        ax, top_landmark[DOSE_COLUMNS].to_numpy(), top_landmark["gene_symbol"].tolist(),
        [f"{d:g}" for d in DOSES], "Top measured dose-responsive profiles",
    )
    add_small_colorbar(fig, ax, hm, "Level 5 z")
    panel_label(ax, "C")

    # D — all prespecified public-serum DILI anchors; missing genes remain visible.
    ax = axes[3]
    anchor_plot = anchors.set_index("gene_symbol").reindex(DILI_ANCHORS)
    hm = heatmap(
        ax, anchor_plot[DOSE_COLUMNS].to_numpy(), DILI_ANCHORS,
        [f"{d:g}" for d in DOSES], "DILI proteomic-anchor profiles", vlim=2.0,
    )
    add_small_colorbar(fig, ax, hm, "Level 5 z")
    panel_label(ax, "D")

    # E — anchor-level adjusted response percentiles and set-level test.
    ax = axes[4]
    present = anchors.loc[anchors["present_in_l1000"]].sort_values("stratum_response_percentile")
    y = np.arange(len(present))
    colors = np.where(present["spearman_rho"] >= 0, COLORS["orange"], COLORS["blue"])
    ax.hlines(y, 0.5, present["stratum_response_percentile"], color="#CACACA", lw=0.8)
    ax.scatter(present["stratum_response_percentile"], y, c=colors, s=20,
               edgecolor=np.where(present["gene_symbol"].eq("FBP1"), "black", "white"), linewidth=0.6)
    ax.axvline(0.5, color="#777777", ls="--", lw=0.6)
    ax.set_yticks(y, present["gene_symbol"], fontsize=6.1)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Response-magnitude percentile\n(within L1000 measurement stratum)")
    p = anchor_stat["permutation_p_one_sided"].iat[0]
    ax.text(0.02, 0.98, f"Set n={len(present)}; matched permutation P={p:.3f}",
            transform=ax.transAxes, va="top", fontsize=6.1,
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.90})
    ax.set_title("Anchor-set competitive test")
    style_axis(ax); panel_label(ax, "E")

    # F — all prespecified mechanism modules, corrected as a family.
    ax = axes[5]
    mod = modules.sort_values("competitive_z")
    y = np.arange(len(mod))
    ax.hlines(y, 0, mod["competitive_z"], color="#C8C8C8", lw=1)
    sig = mod["bh_q_magnitude_family"].lt(0.05)
    ax.scatter(mod["competitive_z"], y, s=25,
               c=np.where(sig, COLORS["vermillion"], COLORS["navy"]), edgecolor="white", lw=0.5)
    compact = {
        "Oxidative/redox defense": "Oxidative/redox",
        "Mitochondrial energetics": "Mito energetics",
        "Bile-acid transport": "Bile-acid transport",
        "Xenobiotic handling": "Xenobiotic handling",
        "ER stress/apoptosis": "ER stress/apoptosis",
        "Hepatocyte integrity": "Hepatocyte integrity",
        "RTK–MAPK–AKT signaling": "RTK–MAPK–AKT",
    }
    labels = [f"{compact[n]} ({k})" for n, k in zip(mod["set_name"], mod["n_present"])]
    ax.set_yticks(y, labels, fontsize=6.1)
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_xlabel("Matched competitive z score")
    ax.set_title("Prespecified mechanism modules")
    style_axis(ax); panel_label(ax, "F")

    # G/H — matched nominal concentration cross-cell-line concordance.
    for panel_index, hep_dose, sci_dose, label in [(6, 1.11, "1uM", "G"), (7, 10.0, "10uM", "H")]:
        ax = axes[panel_index]
        pair = sci_pairs.loc[
            sci_pairs["hepg2_dose_um"].eq(hep_dose) & sci_pairs["sciplex_dose"].eq(sci_dose)
        ].copy()
        ax.scatter(pair["sciplex_z"], pair["hepg2_level5_z"], s=5, color="#9A9A9A",
                   alpha=0.45, edgecolors="none", rasterized=True)
        if "FBP1" in set(pair["gene_symbol"]):
            fbp = pair.loc[pair["gene_symbol"].eq("FBP1")].iloc[0]
            ax.scatter(fbp["sciplex_z"], fbp["hepg2_level5_z"], s=28,
                       color=COLORS["vermillion"], edgecolor="white", linewidth=0.6, zorder=4)
            ax.annotate("FBP1", (fbp["sciplex_z"], fbp["hepg2_level5_z"]),
                        xytext=(4, 4), textcoords="offset points", fontsize=6.1, fontweight="bold")
        coef = np.polyfit(pair["sciplex_z"], pair["hepg2_level5_z"], 1)
        xx = np.linspace(pair["sciplex_z"].quantile(0.005), pair["sciplex_z"].quantile(0.995), 100)
        ax.plot(xx, coef[0] * xx + coef[1], color=COLORS["navy"], lw=0.8)
        ax.axhline(0, color="#BBBBBB", lw=0.4); ax.axvline(0, color="#BBBBBB", lw=0.4)
        stat = sci_summary.loc[
            sci_summary["hepg2_dose_um"].eq(hep_dose)
            & sci_summary["sciplex_dose"].eq(sci_dose)
            & sci_summary["subset"].eq("all_overlap")
        ].iloc[0]
        ax.text(0.03, 0.97,
                f"n={int(stat.n_overlap)}; ρ={stat.spearman_rho:.2f}\n"
                f"sign={stat.sign_concordance:.1%}",
                transform=ax.transAxes, va="top", fontsize=6.1,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.0})
        ax.set_xlabel(f"MCF7 Sci-Plex {sci_dose} standardized value")
        ax.set_ylabel(f"HepG2 {hep_dose:g} µM Level 5 z")
        ax.set_title("Cross-context response concordance")
        style_axis(ax); panel_label(ax, label)

    # I — observed HepG2 profiles for independently supported biochemical targets.
    ax = axes[8]
    high = chembl.loc[
        chembl["multi_document_potent"].eq(True) & chembl["present_in_l1000"]
    ].copy()
    high["target_class"] = np.where(high["is_canonical_nintedanib_target"], "canonical", "other multi-document")
    high = high.sort_values(["target_class", "max_pchembl"], ascending=[True, False])
    labels = [f"{g}{'•' if c else ''}" for g, c in zip(high["gene_symbol"], high["is_canonical_nintedanib_target"])]
    hm = heatmap(
        ax, high[DOSE_COLUMNS].to_numpy(), labels, [f"{d:g}" for d in DOSES],
        "High-confidence ChEMBL targets", vlim=2.0,
    )
    add_small_colorbar(fig, ax, hm, "Level 5 z")
    panel_label(ax, "I")

    # J — potency is independent evidence, and is not assumed to predict response magnitude.
    ax = axes[9]
    canonical = high["is_canonical_nintedanib_target"]
    sizes = 10 + 7 * high["n_documents"].clip(upper=7)
    ax.scatter(high.loc[~canonical, "max_pchembl"], high.loc[~canonical, "trend_magnitude"],
               s=sizes.loc[~canonical], marker="o", color=COLORS["sky"], edgecolor="white", lw=0.5,
               label="Other multi-document")
    ax.scatter(high.loc[canonical, "max_pchembl"], high.loc[canonical, "trend_magnitude"],
               s=sizes.loc[canonical], marker="D", color=COLORS["orange"], edgecolor="white", lw=0.5,
               label="Canonical target")
    offsets = {
        "LYN": (3, 4), "FLT3": (3, -6), "KDR": (3, 3), "PDGFRB": (3, -7),
        "PDGFRA": (3, 3), "FLT4": (3, 4), "FGFR3": (3, -7), "FLT1": (-17, 4),
        "RET": (4, 4), "SRC": (-11, 4), "AXL": (3, -7), "NTRK1": (4, -2),
    }
    for _, row in high.iterrows():
        ax.annotate(row["gene_symbol"], (row["max_pchembl"], row["trend_magnitude"]),
                    xytext=offsets.get(row["gene_symbol"], (3, 2)), textcoords="offset points", fontsize=6.0)
    assoc = associations.loc[
        associations["target_set"].eq("multi_document_potent")
        & associations["response_metric"].eq("trend_magnitude")
    ].iloc[0]
    ax.text(0.03, 0.97, f"n={int(assoc.n_targets)}; descriptive Spearman ρ={assoc.spearman_rho:.2f}",
            transform=ax.transAxes, va="top", fontsize=6.1,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.0})
    ax.set_xlabel("Maximum pChEMBL (biochemical potency)")
    ax.set_ylabel("HepG2 |ρ| × max|z|")
    ax.set_title("Potency–transcription relationship")
    ax.legend(frameon=False, fontsize=6.0, loc="lower left")
    style_axis(ax); panel_label(ax, "J")

    fig.suptitle(
        "Figure 3 | Concentration-resolved nintedanib perturbation maps selected in-vitro profiles",
        x=0.5, y=0.982, fontsize=8.8, fontweight="bold",
    )
    save_figure(fig, FIG / "Figure_3_Nintedanib_HepG2_Dose_Perturbation")
    plt.close(fig)


def plot_supplement(
    long: pd.DataFrame,
    metrics: pd.DataFrame,
    modules: pd.DataFrame,
    sci_summary: pd.DataFrame,
    chembl_sets: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(4, 2, figsize=(7.25, 9.6))
    axes = axes.ravel()

    # A
    ax = axes[0]
    bins = np.linspace(-1.02, 1.02, 37)
    ax.hist(metrics.loc[metrics["is_landmark"].eq(0), "spearman_rho"], bins=bins,
            density=True, color=COLORS["lightgray"], alpha=0.9, label="Inferred")
    ax.hist(metrics.loc[metrics["is_landmark"].eq(1), "spearman_rho"], bins=bins,
            density=True, histtype="step", lw=1.2, color=COLORS["navy"], label="Measured")
    ax.set_xlabel("Spearman ρ"); ax.set_ylabel("Density"); ax.set_title("Trend distribution")
    ax.legend(frameon=False, fontsize=5.5); style_axis(ax); panel_label(ax, "A")

    # B
    ax = axes[1]
    asym_p = []
    for _, row in metrics.iterrows():
        if row["constant_profile"]:
            asym_p.append(1.0)
        else:
            asym_p.append(stats.spearmanr(np.arange(6), [row[c] for c in DOSE_COLUMNS]).pvalue)
    asym_p = np.maximum(np.asarray(asym_p), np.finfo(float).tiny)
    raw_x = -np.log10(asym_p)
    plot_x = np.minimum(raw_x, 4.0)
    ax.scatter(plot_x, -np.log10(metrics["spearman_exact_p"]), s=3,
               color=COLORS["purple"], alpha=0.35, edgecolors="none", rasterized=True)
    ax.plot([0, 2.7], [0, 2.7], color="#777777", ls="--", lw=0.6)
    ax.set_xlim(-0.05, 4.15); ax.set_ylim(-0.05, 2.7)
    ax.text(0.98, 0.94, f"{int((raw_x > 4).sum())} profiles capped",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.2)
    ax.set_xlabel("−log₁₀ asymptotic P (capped at 4)"); ax.set_ylabel("−log₁₀ exact P")
    ax.set_title("Small-n exact-test calibration")
    style_axis(ax); panel_label(ax, "B")

    # C
    ax = axes[2]
    groups = [metrics.loc[metrics["is_landmark"].eq(1), "loo_sign_stability"],
              metrics.loc[metrics["is_landmark"].eq(0), "loo_sign_stability"]]
    bp = ax.boxplot(groups, patch_artist=True, widths=0.5, showfliers=False,
                    medianprops={"color": "white"})
    for box, col in zip(bp["boxes"], [COLORS["navy"], COLORS["gray"]]): box.set_facecolor(col)
    ax.set_xticks([1, 2], ["Measured", "Inferred"])
    ax.set_ylabel("Leave-one-dose-out sign stability")
    ax.set_title("Trend-direction robustness"); style_axis(ax); panel_label(ax, "C")

    # D
    ax = axes[3]
    mod = modules.sort_values("mean_signed_auc_percentile_centered")
    y = np.arange(len(mod))
    ax.hlines(y, 0, mod["mean_signed_auc_percentile_centered"], color="#C8C8C8", lw=1)
    ax.scatter(mod["mean_signed_auc_percentile_centered"], y, s=24,
               c=np.where(mod["mean_signed_auc_percentile_centered"] >= 0, COLORS["orange"], COLORS["blue"]),
               edgecolor="white", lw=0.5)
    ax.set_yticks(y, mod["set_name"], fontsize=5.2)
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_xlabel("Mean centered signed-AUC percentile")
    ax.set_title("Module response direction"); style_axis(ax); panel_label(ax, "D")

    # E/F/G concordance matrices.
    for panel_i, value, subset, title, label, vmin, vmax, cmap in [
        (4, "spearman_rho", "all_overlap", "Sci-Plex correlation sensitivity", "E", -0.25, 0.25, "RdBu_r"),
        (5, "sign_concordance", "all_overlap", "Sci-Plex sign concordance", "F", 0.45, 0.60, "YlGnBu"),
        (6, "spearman_rho", "measured_landmark_only", "Measured-landmark sensitivity", "G", -0.25, 0.25, "RdBu_r"),
    ]:
        ax = axes[panel_i]
        sub = sci_summary.loc[sci_summary["subset"].eq(subset)]
        mat = sub.pivot(index="hepg2_dose_um", columns="sciplex_dose", values=value).reindex(
            index=[1.11, 10.0], columns=["1uM", "10uM"]
        )
        image = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_xticks([0, 1], ["MCF7 1", "MCF7 10"])
        ax.set_yticks([0, 1], ["HepG2 1.11", "HepG2 10"])
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{mat.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(mat.iloc[i, j] - (0.525 if value == 'sign_concordance' else 0)) > (0.035 if value == 'sign_concordance' else .13) else "#222")
        ax.set_title(title); ax.tick_params(length=0, labelsize=5.5)
        for spine in ax.spines.values(): spine.set_visible(False)
        panel_label(ax, label)

    # H
    ax = axes[7]
    targets = chembl_sets.sort_values("competitive_z")
    y = np.arange(len(targets))
    ax.hlines(y, 0, targets["competitive_z"], color="#C8C8C8", lw=1)
    ax.scatter(targets["competitive_z"], y, s=28, color=COLORS["green"], edgecolor="white", lw=0.5)
    short_target_names = {
        "Canonical nintedanib targets": "Canonical targets",
        "Multi-document potent targets": "Multi-doc potent",
        "All pChEMBL≥6 targets": "All pChEMBL≥6",
    }
    ax.set_yticks(y, [f"{short_target_names[x]} (n={n})" for x, n in zip(targets["set_name"], targets["n_present"])], fontsize=5.3)
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_xlabel("Matched competitive z score")
    ax.set_title("ChEMBL target-set response enrichment")
    style_axis(ax); panel_label(ax, "H")

    fig.subplots_adjust(hspace=0.70, wspace=0.52, top=0.925, bottom=0.055)
    fig.suptitle("Supplementary robustness analyses for the nintedanib perturbation layer",
                 fontsize=8.5, fontweight="bold", y=0.985)
    save_figure(fig, FIG / "Supplementary_Figure_S3_Nintedanib_Perturbation_Robustness")
    plt.close(fig)


def write_caption(
    metrics: pd.DataFrame,
    anchor_stat: pd.DataFrame,
    modules: pd.DataFrame,
    sci_summary: pd.DataFrame,
    associations: pd.DataFrame,
) -> None:
    q_hits = int(metrics["bh_q_all_genes"].lt(0.05).sum())
    matched = sci_summary.loc[
        sci_summary["subset"].eq("all_overlap")
        & (
            ((sci_summary["hepg2_dose_um"] == 1.11) & (sci_summary["sciplex_dose"] == "1uM"))
            | ((sci_summary["hepg2_dose_um"] == 10.0) & (sci_summary["sciplex_dose"] == "10uM"))
        )
    ].sort_values("hepg2_dose_um")
    assoc = associations.loc[
        associations["target_set"].eq("multi_document_potent")
        & associations["response_metric"].eq("trend_magnitude")
    ].iloc[0]
    text = f"""# Figure 3 | Concentration-resolved nintedanib perturbation maps selected in-vitro profiles without an FDR-significant gene

**(A)** Distribution of LINCS L1000 Level 5 consensus z scores after 24 h nintedanib exposure across six concentrations in HepG2 cells, separated into directly measured landmark genes (n=978) and inferred genes (n=11,350). Boxes show medians and interquartile ranges; whiskers extend to 1.5× the interquartile range. **(B)** Gene-wise Spearman concentration-ordering coefficients and exact two-sided P values obtained by enumerating all 6! dose-label permutations. The red dashed line marks nominal P=0.05 and is not an FDR threshold. Benjamini–Hochberg correction was applied across all 12,328 genes; {q_hits} genes met q<0.05 (minimum q={metrics['bh_q_all_genes'].min():.3f}). **(C)** The 12 measured landmark genes with the largest |ρ|×maximum |Level 5 z| score. These profiles are ranked responses, not FDR-significant hits. **(D)** Dose-resolved profiles of the 13 prespecified human-serum DILI proteomic anchors; DMGDH was not represented in this L1000 matrix and is shown as missing. **(E)** Anchor-specific response-magnitude percentiles calculated within each L1000 measurement stratum. The outline identifies FBP1. The 12 represented anchors were not collectively enriched for large dose responses in a 25,000-iteration measurement-stratum-matched competitive test (P={anchor_stat['permutation_p_one_sided'].iat[0]:.3f}). **(F)** Analogous competitive tests for seven prespecified mechanistic sentinel modules; multiplicity was controlled by Benjamini–Hochberg correction. **(G, H)** Descriptive, overlap-conditional cross-context concordance between the reported top-1,000 MCF7 Sci-Plex genes and HepG2 L1000 scores at approximately matched concentrations: 1 µM versus 1.11 µM (n={int(matched.iloc[0].n_overlap)}, ρ={matched.iloc[0].spearman_rho:.2f}, sign concordance={matched.iloc[0].sign_concordance:.1%}) and 10 µM versus 10 µM (n={int(matched.iloc[1].n_overlap)}, ρ={matched.iloc[1].spearman_rho:.2f}, sign concordance={matched.iloc[1].sign_concordance:.1%}). Because Sci-Plex reports a truncated top-500 up/top-500 down signature in MCF7 and genes are correlated features, these values are not genome-wide validation estimates or gene-as-independent inference. FBP1 is marked because it was the only prespecified DILI anchor reported in both Sci-Plex dose signatures. **(I)** HepG2 profiles of ChEMBL targets supported by potent activity in at least two assays from at least two documents; dots mark the canonical VEGFR/FGFR/PDGFR families. **(J)** Maximum biochemical pChEMBL versus HepG2 transcriptional dose-response magnitude among these targets (n={int(assoc.n_targets)}; descriptive Spearman ρ={assoc.spearman_rho:.2f}). Point area encodes the number of independent ChEMBL documents. Biochemical potency and transcriptional remodeling are treated as complementary evidence layers, not interchangeable measures of target engagement or toxicity.

All tests are two-sided unless explicitly labeled one-sided. Competitive tests use a one-sided alternative of greater response magnitude and sample genes without replacement while preserving the measured-landmark, best-inferred, and other-inferred composition of each set. The underlying LINCS profiles are composite Level 5 signatures, so the concentration-ordering analysis quantifies reproducible signature geometry across deposited dose summaries rather than biological-replicate variance. HepG2 and MCF7 are in-vitro model contexts; no panel demonstrates nintedanib-specific clinical DILI or causal patient-level effects.
"""
    (OUT / "Figure_3_caption.md").write_text(text, encoding="utf-8")

    supp = """# Supplementary Figure S3 | Robustness analyses for the nintedanib perturbation layer

**(A)** Distribution of six-dose Spearman coefficients in directly measured landmark and inferred L1000 genes. **(B)** Comparison of asymptotic and exact permutation P values, illustrating the finite resolution of inference with six composite dose signatures. **(C)** Leave-one-dose-out stability of the direction of the full six-dose Spearman coefficient. **(D)** Descriptive signed response of prespecified mechanism modules, expressed as mean centered within-stratum percentile of normalized log-dose AUC; formal module tests and adjusted P values are reported in the accompanying table. **(E)** Spearman concordance for all four HepG2–MCF7 dose-pair combinations among overlapping reported genes. **(F)** Corresponding sign concordance. **(G)** Correlation sensitivity analysis restricted to directly measured L1000 landmark genes. **(H)** Measurement-stratum-matched response-magnitude enrichment of canonical, multi-document potent, and all pChEMBL≥6 ChEMBL target sets. Sci-Plex analyses are conditional on its truncated reported signatures. All panels are in-vitro evidence and should not be interpreted as clinical hepatotoxicity estimates.
"""
    (OUT / "Supplementary_Figure_S3_caption.md").write_text(supp, encoding="utf-8")


def write_manuscript_ready_text(
    metrics: pd.DataFrame,
    anchor_stat: pd.DataFrame,
    anchors: pd.DataFrame,
    modules: pd.DataFrame,
    sci_summary: pd.DataFrame,
    associations: pd.DataFrame,
    chembl_sets: pd.DataFrame,
) -> None:
    """Write drop-in Methods and Results prose using computed values only."""
    top = (
        metrics.loc[metrics["is_landmark"].eq(1)]
        .nlargest(12, "trend_magnitude")["gene_symbol"].tolist()
    )
    matched_1 = sci_summary.loc[
        sci_summary["subset"].eq("all_overlap")
        & sci_summary["hepg2_dose_um"].eq(1.11)
        & sci_summary["sciplex_dose"].eq("1uM")
    ].iloc[0]
    matched_10 = sci_summary.loc[
        sci_summary["subset"].eq("all_overlap")
        & sci_summary["hepg2_dose_um"].eq(10.0)
        & sci_summary["sciplex_dose"].eq("10uM")
    ].iloc[0]
    matched_10_landmark = sci_summary.loc[
        sci_summary["subset"].eq("measured_landmark_only")
        & sci_summary["hepg2_dose_um"].eq(10.0)
        & sci_summary["sciplex_dose"].eq("10uM")
    ].iloc[0]
    fbp = anchors.loc[anchors["gene_symbol"].eq("FBP1")].iloc[0]
    anchor_inc = anchors.loc[
        anchors["present_in_l1000"]
        & anchors["spearman_rho"].ge(0.85),
        ["gene_symbol", "spearman_rho", "spearman_exact_p", "bh_q_within_anchor_set"],
    ].sort_values("spearman_rho", ascending=False)
    anchor_phrase = "; ".join(
        f"{r.gene_symbol}, ρ={r.spearman_rho:.2f}, exact P={r.spearman_exact_p:.3f}, anchor-family q={r.bh_q_within_anchor_set:.3f}"
        for r in anchor_inc.itertuples()
    )
    best_module = modules.sort_values("permutation_p_one_sided").iloc[0]
    assoc = associations.loc[
        associations["target_set"].eq("multi_document_potent")
        & associations["response_metric"].eq("trend_magnitude")
    ].iloc[0]
    canonical = chembl_sets.loc[
        chembl_sets["set_name"].eq("Canonical nintedanib targets")
    ].iloc[0]

    text = f"""# Manuscript-ready Methods and Results: concentration-resolved perturbation layer

## Methods

### LINCS L1000 nintedanib concentration series

We analyzed the nintedanib 24-h HepG2 Level 5 consensus signatures deposited in GEO series GSE70138. Six concentrations were available (0.04, 0.12, 0.37, 1.11, 3.33, and 10 µM), spanning 2.40 log10 concentration units. Each Level 5 profile is a MODZ-derived composite signature; five concentrations summarized three deposited replicate instances and the 1.11-µM signature summarized two instances. The analysis matrix comprised 12,328 genes, including 978 directly measured L1000 landmark genes and 11,350 inferred genes. Inferred genes were further separated into best-inferred (BING) and other-inferred strata. These strata were retained throughout the analysis because inference can compress or shift the score distribution relative to directly measured probes.

For each gene, we quantified monotonic concentration ordering using Spearman's rank correlation between expression z score and log10 concentration. Because only six composite dose profiles were available, two-sided P values were calculated by enumerating all 6! dose-label permutations. Benjamini–Hochberg (BH) correction was applied across all 12,328 genes; measurement-family q values were additionally calculated within the measured-landmark, BING, and other-inferred strata as sensitivity estimates. We summarized response amplitude as the maximum absolute Level 5 z score and the trapezoidal area under the z-score curve over log10 concentration, normalized by the observed log-dose range. Direction consistency was the fraction of the six z scores sharing the sign of the normalized AUC. A composite ranking metric, |Spearman ρ| × maximum |z|, combined monotonicity and amplitude without introducing an additional hypothesis test. Robustness was assessed by omitting each concentration in turn and recording the median coefficient, minimum absolute coefficient, and fraction of leave-one-dose-out coefficients with the full-series direction.

### Prespecified DILI anchors and mechanistic sentinel modules

The 13 serum-protein candidates reported in a public multi-cohort DILI proteomics study (ACO1, ASS1, FAH, CPS1, ALDOB, HPD, OTC, DMGDH, GSTA1, FBP1, PCK2, CES1, and LECT2) were fixed as a downstream human-injury anchor set before testing. Seven compact mechanistic sentinel modules were likewise prespecified: oxidative/redox defense, mitochondrial energetics, bile-acid transport, xenobiotic handling, endoplasmic-reticulum stress/apoptosis, hepatocyte integrity, and RTK–MAPK–AKT signaling. Exact membership is reported in `mechanism_module_members.csv`.

For competitive tests, each gene's |ρ| × maximum |z| score was converted to a percentile within its L1000 measurement stratum. The set statistic was the mean stratum-specific percentile. We generated 25,000 null sets by sampling without replacement while preserving the number of measured-landmark, BING, and other-inferred genes represented in the observed set. The primary one-sided P value tested whether the observed set had greater response magnitude than measurement-matched genes. A secondary two-sided test evaluated the mean centered within-stratum percentile of signed normalized AUC. BH correction was performed across the seven mechanism modules. Missing genes were retained in audit tables but excluded from the corresponding test denominator.

### Orthogonal MCF7 Sci-Plex concordance

We compared the HepG2 profiles with independently generated 24-h nintedanib Sci-Plex signatures in MCF7 cells at 1 and 10 µM. The Harmonizome representation contains the 500 highest positive and 500 lowest negative standardized values at each concentration; consequently, all cross-context estimates are conditional on the reported top-1,000 genes and are not genome-wide concordance estimates. We examined all four combinations of HepG2 1.11 or 10 µM and MCF7 1 or 10 µM, with the approximately concentration-matched comparisons designated a priori for presentation. Spearman correlation quantified rank concordance, and sign concordance was tested against 50% with an exact two-sided binomial test and Wilson 95% confidence interval. BH correction was applied across the four dose-pair comparisons within the all-overlap and measured-landmark-only analyses separately.

### ChEMBL pharmacology–transcriptome crosswalk

Nintedanib activities were restricted to exact quantitative measurements against human single-protein targets with a valid pChEMBL value, no ChEMBL data-validity warning, and no duplicate flag. We defined a multi-document potent target as one with maximum pChEMBL ≥6 supported by at least two assays from at least two independent ChEMBL documents. Canonical targets were the VEGFR1–3, FGFR1–3, and PDGFRα/β families. Maximum pChEMBL was related to the HepG2 response metrics using Spearman correlation. ChEMBL target sets were also evaluated using the same L1000-measurement-stratum-matched competitive permutation procedure. Biochemical potency and transcript response were treated as distinct evidence layers; neither was interpreted as proof of in-vivo target engagement or clinical toxicity causality.

## Results

### Exact concentration-ordering analysis identifies ranked in-vitro responses but no FDR-significant single gene

Across 12,328 HepG2 profiles, the finite exact test produced a minimum attainable two-sided P value of 0.00278. No gene met genome-wide BH q<0.05 (minimum q={metrics['bh_q_all_genes'].min():.3f}), preventing false precision from the asymptotic Spearman approximation at n=6. The strongest measured-landmark profiles by the prespecified |ρ| × maximum |z| ranking were {', '.join(top)}. These genes define the most prominent measured in-vitro response geometry in this dose series, but none constitutes a statistically confirmed clinical hepatotoxicity target.

### Human DILI anchors are not globally enriched, although selected anchors rise at higher concentrations

Twelve of 13 serum DILI anchors were represented in L1000; DMGDH was absent. As a set, the anchors had a mean measurement-stratum-adjusted response percentile of {anchor_stat['mean_stratum_response_percentile'].iat[0]:.3f}, indistinguishable from matched random gene sets (competitive z={anchor_stat['competitive_z'].iat[0]:.2f}, P={anchor_stat['permutation_p_one_sided'].iat[0]:.3f}). Three anchors showed strong positive rank ordering before multiplicity control ({anchor_phrase}); none met q<0.05 within the 12 represented anchors. FBP1 did not display a monotonic six-dose trend (ρ={fbp.spearman_rho:.2f}, exact P={fbp.spearman_exact_p:.3f}) but increased from a HepG2 Level 5 z score of {fbp.z_1p11uM:.3f} at 1.11 µM to {fbp.z_10p0uM:.3f} at 10 µM. No prespecified mechanism module survived BH correction; the smallest nominal enrichment was observed for {best_module.set_name.lower()} (n={int(best_module.n_present)}, competitive z={best_module.competitive_z:.2f}, P={best_module.permutation_p_one_sided:.3f}, q={best_module.bh_q_magnitude_family:.3f}). Thus, the perturbation layer supports selective response nodes rather than a globally coordinated DILI-anchor program.

### Cross-cell-line concordance is modest and sensitive to the L1000 measurement layer

The approximately matched 1-µM comparison included {int(matched_1.n_overlap)} reported genes and showed weak positive rank concordance (ρ={matched_1.spearman_rho:.3f}, P={matched_1.spearman_p:.3g}, BH q={matched_1.spearman_bh_q_within_subset:.3g}) without significant sign concordance ({matched_1.sign_concordance:.1%}, 95% CI {matched_1.sign_concordance_ci_low:.1%}–{matched_1.sign_concordance_ci_high:.1%}; P={matched_1.sign_binomial_p:.3f}). At 10 µM, {int(matched_10.n_overlap)} genes showed modest rank (ρ={matched_10.spearman_rho:.3f}, P={matched_10.spearman_p:.3g}, BH q={matched_10.spearman_bh_q_within_subset:.3g}) and sign concordance ({matched_10.sign_concordance:.1%}, 95% CI {matched_10.sign_concordance_ci_low:.1%}–{matched_10.sign_concordance_ci_high:.1%}; P={matched_10.sign_binomial_p:.3g}, BH q={matched_10.sign_bh_q_within_subset:.3g}). However, the 10-µM rank correlation was not retained when restricted to the {int(matched_10_landmark.n_overlap)} directly measured L1000 landmark genes (ρ={matched_10_landmark.spearman_rho:.3f}, P={matched_10_landmark.spearman_p:.3f}). FBP1 was the only prespecified DILI anchor reported in both MCF7 Sci-Plex signatures and changed in the same positive direction in both cell contexts at the presented concentrations. This is orthogonal in-vitro support for a high-concentration FBP1 response, not evidence that FBP1 mediates patient-level nintedanib DILI.

### Biochemical potency and transcriptional remodeling contribute nonredundant evidence

Sixteen L1000-represented targets met the multi-document potent ChEMBL definition. Their maximum biochemical pChEMBL values did not correlate with the |ρ| × maximum |z| HepG2 response magnitude (ρ={assoc.spearman_rho:.3f}, P={assoc.p_value:.3f}). The eight canonical VEGFR/FGFR/PDGFR targets were not enriched for response magnitude (competitive z={canonical.competitive_z:.2f}, P={canonical.permutation_p_one_sided:.3f}, q={canonical.bh_q_magnitude_family:.3f}), although their signed AUC percentiles were collectively positive in the secondary test (P={canonical.signed_auc_permutation_p_two_sided:.4f}, q={canonical.bh_q_signed_auc_family:.4f} across three ChEMBL target sets). This separation shows why biochemical affinity, transcriptional response, and downstream human injury proteins must remain distinct components of the evidence model.

### Interpretation boundary

These analyses establish concentration-resolved and cross-context in-vitro perturbation evidence. HepG2 and MCF7 are not adjudicated DILI cohorts, the L1000 observations are consensus signatures rather than independent patient samples, and Sci-Plex concordance is conditioned on a top-1,000-gene list. Accordingly, the results prioritize measurable response nodes for integration with human proteomics, liver-cell network perturbation, and pharmacovigilance; they do not by themselves demonstrate nintedanib-specific clinical DILI or causal target mediation.
"""
    (OUT / "manuscript_ready_lincs_methods_results.md").write_text(text, encoding="utf-8")


def write_summary_and_dictionary(
    metrics: pd.DataFrame,
    anchor_stat: pd.DataFrame,
    modules: pd.DataFrame,
    sci_summary: pd.DataFrame,
    associations: pd.DataFrame,
) -> None:
    matched_10 = sci_summary.loc[
        sci_summary["subset"].eq("all_overlap")
        & sci_summary["hepg2_dose_um"].eq(10.0)
        & sci_summary["sciplex_dose"].eq("10uM")
    ].iloc[0]
    assoc = associations.loc[
        associations["target_set"].eq("multi_document_potent")
        & associations["response_metric"].eq("trend_magnitude")
    ].iloc[0]
    best_module = modules.sort_values("permutation_p_one_sided").iloc[0]
    summary = {
        "analysis_population": {
            "hepg2_level5_genes": len(metrics),
            "measured_landmark_genes": int(metrics["is_landmark"].sum()),
            "inferred_genes": int(metrics["is_landmark"].eq(0).sum()),
            "doses_um": DOSES.tolist(),
            "exposure_hours": 24,
        },
        "exact_dose_trend": {
            "minimum_p": float(metrics["spearman_exact_p"].min()),
            "minimum_bh_q": float(metrics["bh_q_all_genes"].min()),
            "bh_q_lt_0_05": int(metrics["bh_q_all_genes"].lt(0.05).sum()),
        },
        "dili_anchor_set": {
            "represented": int(anchor_stat["n_present"].iat[0]),
            "competitive_z": float(anchor_stat["competitive_z"].iat[0]),
            "permutation_p": float(anchor_stat["permutation_p_one_sided"].iat[0]),
        },
        "best_nominal_mechanism_module": {
            "name": best_module["set_name"],
            "competitive_z": float(best_module["competitive_z"]),
            "permutation_p": float(best_module["permutation_p_one_sided"]),
            "bh_q": float(best_module["bh_q_magnitude_family"]),
        },
        "matched_10uM_hepg2_mcf7": {
            "n_overlap": int(matched_10["n_overlap"]),
            "spearman_rho": float(matched_10["spearman_rho"]),
            "spearman_bh_q": float(matched_10["spearman_bh_q_within_subset"]),
            "sign_concordance": float(matched_10["sign_concordance"]),
            "sign_bh_q": float(matched_10["sign_bh_q_within_subset"]),
        },
        "multi_document_chembl_potency_vs_transcript_response": {
            "n_targets": int(assoc["n_targets"]),
            "spearman_rho": float(assoc["spearman_rho"]),
            "p_value": float(assoc["p_value"]),
        },
        "main_figure_panels": 10,
        "interpretation_boundary": (
            "The perturbation results are in-vitro prioritization evidence and do not "
            "demonstrate nintedanib-specific clinical DILI or causal target mediation."
        ),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    definitions = [
        ("spearman_rho", "Spearman rank correlation between Level 5 z score and log10 dose across six composite signatures."),
        ("spearman_exact_p", "Two-sided exact P value from all 6! dose-label permutations."),
        ("bh_q_all_genes", "Benjamini-Hochberg q value across 12,328 genes."),
        ("bh_q_measurement_family", "BH q value within measured-landmark, BING, or other-inferred family."),
        ("auc_log10dose", "Trapezoidal integral of Level 5 z score over log10 concentration."),
        ("auc_log10dose_normalized", "Log-dose AUC divided by the observed log10-dose range; interpretable as log-dose-weighted mean z."),
        ("max_abs_z", "Maximum absolute Level 5 z score across six doses."),
        ("direction_consistency", "Fraction of six z scores sharing the sign of normalized log-dose AUC."),
        ("trend_magnitude", "Descriptive ranking score |Spearman rho| multiplied by maximum |z|; not a test statistic."),
        ("loo_sign_stability", "Fraction of six leave-one-dose-out correlations sharing the full-series direction."),
        ("stratum_response_percentile", "Percentile of trend_magnitude within L1000 measurement stratum."),
        ("competitive_z", "Observed gene-set mean response percentile minus matched-null mean, divided by matched-null SD."),
        ("permutation_p_one_sided", "Plus-one corrected proportion of 25,000 measurement-stratum-matched null sets at least as large as observed."),
        ("sciplex sign_concordance", "Fraction of overlapping genes with the same response sign in HepG2 and the truncated MCF7 Sci-Plex signature."),
        ("max_pchembl", "Maximum valid exact pChEMBL value for a human single-protein target."),
    ]
    pd.DataFrame(definitions, columns=["field", "definition"]).to_csv(
        OUT / "data_dictionary.csv", index=False
    )


def write_manifest() -> None:
    source_paths = [
        LINCS / "nintedanib_hepg2_level5_all_genes_tidy.csv.gz",
        LINCS / "nintedanib_hepg2_signature_metadata.csv",
        PHARM / "sciplex_nintedanib_signature_long.csv",
        PHARM / "chembl_target_summary.csv",
        PHARM / "chembl_human_single_protein_activities.csv",
    ]
    outputs = sorted([p for p in OUT.iterdir() if p.is_file()])
    figures = sorted(FIG.glob("Figure_3_Nintedanib_HepG2_Dose_Perturbation.*"))
    figures += sorted(FIG.glob("Supplementary_Figure_S3_Nintedanib_Perturbation_Robustness.*"))
    payload = {
        "analysis": "de_novo_nintedanib_hepg2_concentration_response",
        "random_seed": SEED,
        "gene_set_permutations": N_PERM,
        "doses_um": DOSES.tolist(),
        "exposure_hours": 24,
        "sources": [
            {"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in source_paths
        ],
        "outputs": [
            {"path": str(p.relative_to(ROOT)), "sha256": sha256(p)}
            for p in outputs + figures if p.name != "analysis_manifest.json"
        ],
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": mpl.__version__,
        },
        "interpretation_boundary": (
            "HepG2 L1000 and MCF7 Sci-Plex are in-vitro drug-response models. "
            "They do not identify clinical DILI cases or prove toxicity causality; "
            "ChEMBL potency does not establish in-vivo target engagement."
        ),
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    long, metrics = load_gene_profiles()
    anchor_stat, anchors, modules = analyze_gene_sets(metrics)
    sci_summary, sci_pairs = analyze_sciplex(metrics)
    chembl, associations, chembl_sets = analyze_chembl(metrics)
    make_integrated_evidence(metrics, sci_pairs, chembl)
    plot_main_figure(
        long, metrics, anchor_stat, anchors, modules,
        sci_summary, sci_pairs, chembl, associations,
    )
    plot_supplement(long, metrics, modules, sci_summary, chembl_sets)
    write_caption(metrics, anchor_stat, modules, sci_summary, associations)
    write_manuscript_ready_text(
        metrics, anchor_stat, anchors, modules, sci_summary,
        associations, chembl_sets,
    )
    write_summary_and_dictionary(metrics, anchor_stat, modules, sci_summary, associations)
    write_manifest()
    print(json.dumps({
        "genes": len(metrics),
        "bh_q_lt_0_05": int(metrics["bh_q_all_genes"].lt(0.05).sum()),
        "represented_dili_anchors": int(anchors["present_in_l1000"].sum()),
        "main_figure_panels": 10,
        "figure_stem": str((FIG / "Figure_3_Nintedanib_HepG2_Dose_Perturbation").relative_to(ROOT)),
    }, indent=2))


if __name__ == "__main__":
    main()
