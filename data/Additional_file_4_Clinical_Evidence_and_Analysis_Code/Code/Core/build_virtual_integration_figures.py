#!/usr/bin/env python3
"""Build transparent multi-evidence integration and Figures 4-5.

This script deliberately separates exposure-proximal pharmacologic targets from
downstream circulating DILI phenotype anchors. scTenifoldKnk outputs are used as
unsigned network-displacement ranks only. They are never interpreted as
protective, harmful, causal, or as expression fold changes.

Inputs are de novo public-data outputs under ``fresh_analysis``. Outputs include
tidy audit tables, leave-one-source-out rankings, Pareto fronts, and two
ten-panel result figures in vector and 600-dpi raster formats.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import rankdata, spearmanr


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "fresh_analysis"
OUT = ANALYSIS / "integration"
FIG_OUT = ANALYSIS / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG_OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260904
RNG = np.random.default_rng(SEED)

DILI_ANCHORS = [
    "ACO1", "ASS1", "FAH", "CPS1", "ALDOB", "HPD", "OTC", "DMGDH",
    "GSTA1", "FBP1", "PCK2", "CES1", "LECT2",
]
CANONICAL_TARGETS = ["KDR", "FLT1", "FLT4", "FGFR1", "FGFR2", "FGFR3", "PDGFRA", "PDGFRB"]
PRIMARY_COMPARTMENTS = ["Hepatocyte", "Endothelial", "Macrophage"]
ALL_COMPARTMENTS = ["Hepatocyte", "Endothelial", "Macrophage", "Cholangiocyte", "Stellate", "Immune_other"]

PALETTE = {
    "teal": "#007C83",
    "blue": "#3B6FB6",
    "orange": "#D97904",
    "red": "#B33A3A",
    "purple": "#7B5AA6",
    "green": "#4A8C5D",
    "gray": "#7A7A7A",
    "light_gray": "#E9ECEF",
    "dark": "#242424",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_adjust(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    x = p[ok]
    order = np.argsort(x)
    ranked = x[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    out[ok] = restored
    return out


def percentile_score(series: pd.Series, higher_better: bool = True) -> pd.Series:
    """Empirical [0,1] score with average ranks; preserves missingness."""
    out = pd.Series(np.nan, index=series.index, dtype=float)
    mask = series.notna()
    if mask.sum() == 1:
        out.loc[mask] = 1.0
    elif mask.sum() > 1:
        r = series.loc[mask].rank(method="average", ascending=higher_better)
        out.loc[mask] = (r - 1) / (mask.sum() - 1)
    return out


def signed_logit_fraction(x: pd.Series) -> pd.Series:
    return np.log((x.clip(1e-4, 1 - 1e-4)) / (1 - x.clip(1e-4, 1 - 1e-4)))


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.11, 1.08, label, transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="top", ha="left")


def clean_axis(ax: mpl.axes.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def load_inputs() -> dict[str, pd.DataFrame]:
    paths = {
        "sc_expr": ANALYSIS / "liver_singlecell" / "target_expression_by_compartment.csv",
        "sc_donor": ANALYSIS / "liver_singlecell" / "donor_level_target_expression.csv",
        "sc_type": ANALYSIS / "liver_singlecell" / "donor_celltype_target_expression.csv",
        "pca": ANALYSIS / "liver_singlecell" / "descriptive_pca_coordinates.csv",
        "pharm": ANALYSIS / "pharmacology" / "chembl_target_summary.csv",
        "sciplex": ANALYSIS / "pharmacology" / "sciplex_dose_concordance.csv",
        "prot": ANALYSIS / "proteomics" / "candidate_evidence_table.csv",
        "vko_gene": ANALYSIS / "virtual_ko" / "virtual_knockout_gene_results.csv.gz",
        "vko_seed": ANALYSIS / "virtual_ko" / "virtual_knockout_seed_summary.csv",
        "vko_stable": ANALYSIS / "virtual_ko" / "virtual_knockout_stability_summary.csv",
        "vko_down": ANALYSIS / "virtual_ko" / "stable_top_downstream_genes.csv",
        "lincs": ANALYSIS / "lincs" / "output" / "nintedanib_hepg2_level5_all_genes_tidy.csv.gz",
        "string": ANALYSIS / "inputs" / "string_candidate_network.tsv",
        "string_expanded": ANALYSIS / "inputs" / "string_candidate_network_expanded.tsv",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required inputs: {missing}")
    return {key: pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",") for key, path in paths.items()}


def derive_lincs_metrics(lincs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Use agent-generated metrics if available; otherwise derive equivalent metrics."""
    metrics_path = ANALYSIS / "lincs" / "analysis" / "gene_dose_response_metrics.csv.gz"
    if metrics_path.exists():
        metrics = pd.read_csv(metrics_path)
        # Harmonize anticipated names.
        aliases = {
            "max_abs_zscore": "max_abs_z",
            "rho": "spearman_rho",
            "spearman_exact_p": "exact_p",
            "bh_q_all_genes": "bh_q",
            "auc_log10dose": "auc_logdose",
            "auc_log10dose_normalized": "auc_logdose_normalized",
        }
        metrics = metrics.rename(columns={old: new for old, new in aliases.items() if old in metrics and new not in metrics})
        zcols = [c for c in metrics.columns if c.startswith("z_")]
        long = metrics.melt(
            id_vars=[c for c in ["gene_symbol", "gene_id", "is_landmark", "is_bing"] if c in metrics.columns],
            value_vars=zcols,
            var_name="dose_label",
            value_name="level5_zscore",
        )
        long["dose_um"] = long["dose_label"].str.extract(r"z_([0-9]+(?:p[0-9]+)?)")[0].str.replace("p", ".", regex=False).astype(float)
        return metrics, long

    rows = []
    doses = sorted(lincs["dose_um"].unique())
    for gene, group in lincs.groupby("gene_symbol", sort=False):
        group = group.sort_values("dose_um")
        if group["dose_um"].nunique() < 4:
            continue
        rho, p = spearmanr(np.log10(group["dose_um"]), group["level5_zscore"])
        z = group["level5_zscore"].to_numpy(float)
        logdose = np.log10(group["dose_um"].to_numpy(float))
        auc = float(np.trapezoid(z, logdose))
        row = {
            "gene_symbol": gene,
            "gene_id": group["gene_id"].iloc[0],
            "is_landmark": int(group["is_landmark"].iloc[0]),
            "is_bing": int(group["is_bing"].iloc[0]),
            "spearman_rho": rho,
            "exact_p": p,
            "auc_logdose": auc,
            "max_abs_z": float(np.max(np.abs(z))),
            "direction_consistency": float(max(np.mean(z >= 0), np.mean(z <= 0))),
            "peak_dose_um": float(group.iloc[np.argmax(np.abs(z))]["dose_um"]),
            "peak_z": float(group.iloc[np.argmax(np.abs(z))]["level5_zscore"]),
        }
        for dose in doses:
            vals = group.loc[group["dose_um"].eq(dose), "level5_zscore"]
            row[f"z_{str(dose).replace('.', 'p')}uM"] = float(vals.iloc[0]) if len(vals) else np.nan
        rows.append(row)
    metrics = pd.DataFrame(rows)
    metrics["bh_q"] = bh_adjust(metrics["exact_p"])
    return metrics, lincs.copy()


def matched_null_analysis(data: dict[str, pd.DataFrame], n_permutations: int = 3000) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Expression/detection-matched null for aggregate DILI-anchor displacement.

    Matching is done separately for each compartment using nearest neighbors in
    standardized log-mean-expression and logit-detection space. This does not
    turn the model into a causal test; it only calibrates whether the fixed DILI
    anchor set is unusually ranked relative to similarly expressed genes.
    """
    vko = data["vko_gene"].copy()
    drug_runs = vko.loc[vko["ko_class"].isin(["drug_target", "canonical_poly_target"])].copy()
    rows: list[dict[str, object]] = []

    for compartment in PRIMARY_COMPARTMENTS:
        matrix_path = ANALYSIS / "liver_singlecell" / f"{compartment.lower()}_network_matrix.csv.gz"
        matrix = pd.read_csv(matrix_path, index_col=0)
        expr = pd.DataFrame({
            "gene": matrix.index.astype(str),
            "mean_expression": matrix.mean(axis=1).to_numpy(float),
            "detection_fraction": matrix.gt(0).mean(axis=1).to_numpy(float),
        }).set_index("gene")
        feature = pd.DataFrame(index=expr.index)
        feature["x"] = np.log1p(expr["mean_expression"])
        feature["y"] = signed_logit_fraction(expr["detection_fraction"])
        feature = (feature - feature.mean()) / feature.std(ddof=0).replace(0, 1)

        sub = drug_runs.loc[drug_runs["compartment"].eq(compartment)]
        for (seed, ko_label, ko_class), run in sub.groupby(["seed", "ko_label", "ko_class"], sort=False):
            run = run.loc[~run["is_knocked_gene"].astype(bool)].copy()
            present = set(run["gene"]) & set(feature.index)
            anchors = [g for g in DILI_ANCHORS if g in present]
            background = sorted(present - set(DILI_ANCHORS) - set(str(ko_label).split("+")))
            if len(anchors) < 4 or len(background) < 50:
                continue
            by_gene = run.set_index("gene")["distance_percentile"]
            candidate_sets: dict[str, np.ndarray] = {}
            bg_features = feature.loc[background].to_numpy(float)
            for anchor in anchors:
                target = feature.loc[anchor].to_numpy(float)
                dist = np.sum((bg_features - target) ** 2, axis=1)
                k = min(30, len(background))
                nearest = np.asarray(background, dtype=object)[np.argpartition(dist, k - 1)[:k]]
                candidate_sets[anchor] = nearest

            local_rng = np.random.default_rng(int(seed) + sum(ord(c) for c in str(ko_label)))
            null = np.empty(n_permutations, dtype=float)
            for i in range(n_permutations):
                selected: list[str] = []
                used: set[str] = set()
                for anchor in anchors:
                    pool = [g for g in candidate_sets[anchor] if g not in used]
                    if not pool:
                        pool = candidate_sets[anchor].tolist()
                    choice = str(local_rng.choice(pool))
                    selected.append(choice)
                    used.add(choice)
                null[i] = float(by_gene.reindex(selected).mean())
            observed = float(by_gene.reindex(anchors).mean())
            null_mean = float(np.mean(null))
            null_sd = float(np.std(null, ddof=1))
            z = (observed - null_mean) / null_sd if null_sd > 0 else np.nan
            p_two = (1 + np.sum(np.abs(null - null_mean) >= abs(observed - null_mean))) / (n_permutations + 1)
            rows.append({
                "compartment": compartment,
                "seed": int(seed),
                "ko_label": ko_label,
                "ko_class": ko_class,
                "n_anchors": len(anchors),
                "n_background": len(background),
                "n_permutations": n_permutations,
                "observed_anchor_mean_percentile": observed,
                "matched_null_mean": null_mean,
                "matched_null_sd": null_sd,
                "matched_null_z": z,
                "empirical_p_two_sided": p_two,
            })

    seed_df = pd.DataFrame(rows)
    seed_df["empirical_q_bh_within_compartment"] = seed_df.groupby("compartment")["empirical_p_two_sided"].transform(bh_adjust)
    seed_df.to_csv(OUT / "virtual_ko_expression_matched_null_by_seed.csv", index=False)
    stable = (
        seed_df.groupby(["compartment", "ko_label", "ko_class"], as_index=False)
        .agg(
            n_seeds=("seed", "nunique"),
            median_matched_z=("matched_null_z", "median"),
            min_matched_z=("matched_null_z", "min"),
            max_matched_z=("matched_null_z", "max"),
            n_seed_nominal_p_lt_0_05=("empirical_p_two_sided", lambda x: int((x < 0.05).sum())),
            n_seed_bh_q_lt_0_05=("empirical_q_bh_within_compartment", lambda x: int((x < 0.05).sum())),
        )
    )
    stable.to_csv(OUT / "virtual_ko_expression_matched_null_stability.csv", index=False)
    return seed_df, stable


def top20_jaccard(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    vko = data["vko_gene"].copy()
    vko = vko.loc[~vko["is_knocked_gene"].astype(bool)]
    rows = []
    for keys, group in vko.groupby(["compartment", "ko_label", "ko_class"], sort=False):
        sets = {}
        for seed, run in group.groupby("seed"):
            sets[int(seed)] = set(run.nlargest(20, "distance")["gene"].astype(str))
        pair_values = []
        for a, b in itertools.combinations(sorted(sets), 2):
            union = sets[a] | sets[b]
            j = len(sets[a] & sets[b]) / len(union) if union else np.nan
            pair_values.append(j)
            rows.append({
                "compartment": keys[0], "ko_label": keys[1], "ko_class": keys[2],
                "seed_a": a, "seed_b": b, "jaccard_top20": j,
            })
        if len(sets) == 1:
            rows.append({
                "compartment": keys[0], "ko_label": keys[1], "ko_class": keys[2],
                "seed_a": next(iter(sets)), "seed_b": np.nan, "jaccard_top20": np.nan,
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "virtual_ko_top20_pairwise_jaccard.csv", index=False)
    return out


def string_graph_metrics(data: dict[str, pd.DataFrame], candidate_genes: list[str]) -> tuple[pd.DataFrame, nx.Graph]:
    # Union the compact candidate query and the expanded neighborhood query.
    # The compact result retains lower-score direct candidate edges that are
    # absent from the expanded high-confidence neighborhood response.
    edges = pd.concat([data["string"], data["string_expanded"]], ignore_index=True)
    edges["_edge_key"] = edges.apply(
        lambda r: "||".join(sorted([str(r["preferredName_A"]), str(r["preferredName_B"])])), axis=1
    )
    edges = edges.sort_values("score", ascending=False).drop_duplicates("_edge_key")
    edges = edges.loc[edges["score"].ge(0.4)].copy()
    graph = nx.Graph()
    for row in edges.itertuples(index=False):
        a, b, score = str(row.preferredName_A), str(row.preferredName_B), float(row.score)
        graph.add_edge(a, b, score=score, distance=-math.log(max(score, 1e-6)),
                       experimental=float(row.escore), database=float(row.dscore), text=float(row.tscore))
    rows = []
    for gene in candidate_genes:
        if gene not in graph:
            rows.append({"gene_symbol": gene, "string_present": False, "max_path_product_to_opposite_class": 0.0,
                         "weighted_degree": 0.0, "nearest_opposite_node": "", "shortest_hops": np.nan})
            continue
        opposite = DILI_ANCHORS if gene not in DILI_ANCHORS else CANONICAL_TARGETS
        best_score, best_node, best_hops = 0.0, "", np.nan
        for target in opposite:
            if target not in graph or target == gene:
                continue
            try:
                path = nx.shortest_path(graph, gene, target, weight="distance")
            except nx.NetworkXNoPath:
                continue
            product = 1.0
            for u, v in zip(path[:-1], path[1:]):
                product *= float(graph[u][v]["score"])
            if product > best_score:
                best_score, best_node, best_hops = product, target, len(path) - 1
        rows.append({
            "gene_symbol": gene,
            "string_present": True,
            "max_path_product_to_opposite_class": best_score,
            "weighted_degree": float(sum(float(d["score"]) for _, _, d in graph.edges(gene, data=True))),
            "nearest_opposite_node": best_node,
            "shortest_hops": best_hops,
        })
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "string_candidate_path_metrics.csv", index=False)
    return metrics, graph


def pareto_front_ids(values: pd.DataFrame, maximize: list[str]) -> pd.Series:
    """Iterative non-dominated sorting; lower front number is better."""
    x = values[maximize].to_numpy(float)
    fronts = np.full(len(values), np.nan)
    remaining = list(range(len(values)))
    front = 1
    while remaining:
        nondominated = []
        for i in remaining:
            dominated = False
            for j in remaining:
                if i == j:
                    continue
                if np.all(x[j] >= x[i]) and np.any(x[j] > x[i]):
                    dominated = True
                    break
            if not dominated:
                nondominated.append(i)
        fronts[nondominated] = front
        remaining = [i for i in remaining if i not in nondominated]
        front += 1
    return pd.Series(fronts.astype(int), index=values.index)


def consensus_and_loso(frame: pd.DataFrame, domains: list[str], candidate_class: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    x = frame.copy()
    rank_cols = []
    for domain in domains:
        col = f"{domain}_rank"
        x[col] = x[domain].rank(method="average", ascending=False, na_option="bottom")
        rank_cols.append(col)
    x["consensus_rank_statistic"] = x[rank_cols].median(axis=1)
    x["consensus_rank_mean_tiebreaker"] = x[rank_cols].mean(axis=1)
    x["consensus_rank"] = x[["consensus_rank_statistic", "consensus_rank_mean_tiebreaker"]].apply(tuple, axis=1).rank(method="dense").astype(int)
    loso_rows = []
    for omitted in domains:
        remaining = [f"{d}_rank" for d in domains if d != omitted]
        stat = x[remaining].median(axis=1)
        tie = x[remaining].mean(axis=1)
        temp = pd.DataFrame({"stat": stat, "tie": tie}, index=x.index).sort_values(["stat", "tie"])
        temp["rank"] = np.arange(1, len(temp) + 1)
        for idx, row in temp.iterrows():
            loso_rows.append({
                "candidate_class": candidate_class,
                "gene_symbol": x.loc[idx, "gene_symbol"],
                "omitted_source_group": omitted,
                "leave_one_source_out_rank": int(row["rank"]),
            })
    loso = pd.DataFrame(loso_rows)
    stability = loso.groupby("gene_symbol")["leave_one_source_out_rank"].agg(["min", "median", "max"]).rename(
        columns={"min": "loso_rank_best", "median": "loso_rank_median", "max": "loso_rank_worst"}
    )
    x = x.merge(stability, left_on="gene_symbol", right_index=True, how="left")
    return x, loso


def integrate_candidates(
    data: dict[str, pd.DataFrame],
    lincs_metrics: pd.DataFrame,
    string_metrics: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prot = data["prot"].copy()
    sc = data["sc_expr"].copy()
    stable = data["vko_stable"].copy()
    pharm = data["pharm"].copy()
    sciplex = data["sciplex"].copy()

    # --- Downstream phenotype anchors: four independent source groups. ---
    a = pd.DataFrame({"gene_symbol": DILI_ANCHORS}).merge(prot, left_on="gene_symbol", right_on="protein", how="left")
    clinical_parts = pd.DataFrame(index=a.index)
    clinical_parts["etiology"] = percentile_score(a["confirmatory_DO_vs_NDO_roc_auc_separation"])
    clinical_parts["acute"] = percentile_score(a["confirmatory_DO_vs_HV_roc_auc_separation"])
    clinical_parts["recovery"] = percentile_score(a["confirmatory_DO_vs_DF_roc_auc_separation"])
    clinical_parts["severity"] = percentile_score(a["within_DO_alt_spearman_rho"].abs())
    # Etiology discrimination is the prespecified primary clinical property.
    a["clinical_DILI"] = 0.45 * clinical_parts["etiology"] + 0.20 * clinical_parts["acute"] + 0.20 * clinical_parts["recovery"] + 0.15 * clinical_parts["severity"]
    missing_recovery = clinical_parts["recovery"].isna()
    if missing_recovery.any():
        a.loc[missing_recovery, "clinical_DILI"] = (
            0.45 * clinical_parts.loc[missing_recovery, "etiology"]
            + 0.25 * clinical_parts.loc[missing_recovery, "acute"]
            + 0.30 * clinical_parts.loc[missing_recovery, "severity"]
        )

    hep = sc.loc[sc["compartment"].eq("Hepatocyte"), ["gene_symbol", "mean_detection_fraction", "mean_expression"]]
    anchor_vko = stable.loc[(stable["compartment"].eq("Hepatocyte")) & (stable["ko_class"].eq("dili_anchor")),
                            ["ko_label", "median_dili_anchor_percentile", "n_seeds"]].rename(columns={"ko_label": "gene_symbol"})
    a = a.merge(hep, on="gene_symbol", how="left").merge(anchor_vko, on="gene_symbol", how="left")
    a["liver_expression_component"] = percentile_score(a["mean_detection_fraction"])
    a["liver_network_component"] = percentile_score(a["median_dili_anchor_percentile"])
    a["liver_context"] = 0.65 * a["liver_expression_component"] + 0.35 * a["liver_network_component"]

    lcols = [c for c in ["gene_symbol", "spearman_rho", "exact_p", "bh_q", "max_abs_z", "is_landmark", "is_bing", "constant_profile"] if c in lincs_metrics.columns]
    a = a.merge(lincs_metrics[lcols].drop_duplicates("gene_symbol"), on="gene_symbol", how="left")
    reliability = np.where(a["is_landmark"].fillna(0).eq(1), 1.0, np.where(a["is_bing"].fillna(0).eq(1), 0.8, 0.5))
    perturb_raw = (0.60 * percentile_score(a["spearman_rho"].abs()) + 0.40 * percentile_score(a["max_abs_z"])) * reliability
    a["HepG2_perturbation"] = perturb_raw.where(a["spearman_rho"].notna())

    a = a.merge(string_metrics, on="gene_symbol", how="left")
    a["STRING_context"] = 0.65 * percentile_score(a["max_path_product_to_opposite_class"]) + 0.35 * percentile_score(a["weighted_degree"])
    sci = sciplex[["gene_symbol", "observed_in_both_doses", "direction_concordant", "mean_z"]].copy()
    a = a.merge(sci, on="gene_symbol", how="left")
    a["SciPlex_replication_flag"] = a["observed_in_both_doses"].fillna(False).astype(bool) & a["direction_concordant"].fillna(False).astype(bool)

    anchor_domains = ["clinical_DILI", "liver_context", "HepG2_perturbation", "STRING_context"]
    eligible = a[anchor_domains].notna().all(axis=1)
    a["pareto_front"] = np.nan
    a.loc[eligible, "pareto_front"] = pareto_front_ids(a.loc[eligible], anchor_domains)
    a, anchor_loso = consensus_and_loso(a, anchor_domains, "downstream_DILI_anchor")

    # Phenotype anchors remain on a P track and are never presented as direct
    # nintedanib-binding targets.
    specific = a["dili_vs_non_dili_specificity_signal"].fillna(False).astype(bool) & a["specificity_replicates_discovery_direction"].fillna(False).astype(bool)
    liver_localized = a["mean_detection_fraction"].ge(0.25)
    a["proteomics_tier_original"] = a["evidence_tier"].astype(str)
    a["evidence_tier"] = "P3: limited/inconsistent phenotype anchor"
    a.loc[a["proteomics_tier_original"].str.startswith("C:"), "evidence_tier"] = "P2: general injury/recovery phenotype anchor"
    a.loc[specific & liver_localized, "evidence_tier"] = "P1: DILI-discriminating phenotype anchor"
    a["candidate_role"] = "downstream circulating phenotype anchor"
    a["causal_status"] = "not established"

    # --- Exposure-proximal pharmacologic targets. ---
    screened = stable.loc[stable["ko_class"].eq("drug_target"), "ko_label"].drop_duplicates().tolist()
    t = pharm.loc[pharm["gene_symbol"].isin(screened)].copy()
    t = t.sort_values(["multi_document_potent", "max_pchembl"], ascending=False).drop_duplicates("gene_symbol")
    tier_value = t["pharmacology_tier"].map({
        "canonical_target": 1.0,
        "multi_document_noncanonical": 0.75,
        "single_document_or_single_assay_screen": 0.35,
    }).fillna(0.0)
    t["pharmacology"] = 0.55 * percentile_score(t["max_pchembl"]) + 0.25 * percentile_score(t["median_pchembl"]) + 0.20 * tier_value

    expr_primary = sc.loc[sc["compartment"].isin(PRIMARY_COMPARTMENTS),
                          ["gene_symbol", "compartment", "mean_detection_fraction"]].copy()
    max_detection = expr_primary.groupby("gene_symbol", as_index=False)["mean_detection_fraction"].max().rename(
        columns={"mean_detection_fraction": "max_primary_detection"}
    )
    vko_drug = stable.loc[stable["ko_class"].eq("drug_target")].copy()
    vko_drug["vko_within_compartment_rank"] = vko_drug.groupby("compartment")["median_dili_anchor_percentile"].rank(pct=True)
    pair = vko_drug.merge(
        expr_primary,
        left_on=["ko_label", "compartment"],
        right_on=["gene_symbol", "compartment"],
        how="left",
    )
    pair["detection_within_compartment_rank"] = pair.groupby("compartment")["mean_detection_fraction"].rank(pct=True)
    # Non-compensatory maximin context: expression and virtual-KO support must
    # arise in the same compartment; a high value in one cannot offset a low
    # value in the other.
    pair["same_compartment_context_score"] = pair[["detection_within_compartment_rank", "vko_within_compartment_rank"]].min(axis=1)
    pair["same_compartment_primary_gate"] = (
        pair["mean_detection_fraction"].ge(0.10)
        & pair["n_seeds"].eq(3)
        & pair["vko_within_compartment_rank"].ge(0.75)
    )

    seed_values = data["vko_seed"].loc[data["vko_seed"]["ko_class"].eq("drug_target")].copy()
    seed_values["seed_vko_within_compartment_rank"] = seed_values.groupby(
        ["compartment", "seed"]
    )["dili_anchor_mean_percentile"].rank(pct=True)
    seed_sensitivity = (
        seed_values.groupby(["ko_label", "compartment"], as_index=False)
        .agg(
            n_seed_results=("seed", "nunique"),
            n_seeds_top_quartile=("seed_vko_within_compartment_rank", lambda x: int((x >= 0.75).sum())),
            min_seed_vko_rank=("seed_vko_within_compartment_rank", "min"),
            median_seed_vko_rank=("seed_vko_within_compartment_rank", "median"),
        )
    )
    pair = pair.merge(seed_sensitivity, on=["ko_label", "compartment"], how="left")
    pair["same_compartment_all_seed_sensitivity_gate"] = (
        pair["mean_detection_fraction"].ge(0.10)
        & pair["n_seed_results"].eq(3)
        & pair["n_seeds_top_quartile"].eq(3)
    )
    pair.to_csv(OUT / "target_same_compartment_context_all_pairs.csv", index=False)

    # Select the best joint compartment. Passing the primary gate takes
    # precedence, followed by the maximin context score; separate maxima are
    # never spliced across cell types.
    joint = (
        pair.sort_values(
            ["ko_label", "same_compartment_primary_gate", "same_compartment_context_score", "mean_detection_fraction"],
            ascending=[True, False, False, False],
        )
        .drop_duplicates("ko_label")
        [[
            "ko_label", "compartment", "mean_detection_fraction", "median_dili_anchor_percentile",
            "vko_within_compartment_rank", "detection_within_compartment_rank",
            "same_compartment_context_score", "same_compartment_primary_gate",
            "same_compartment_all_seed_sensitivity_gate", "n_seeds", "n_seed_results",
            "n_seeds_top_quartile", "min_seed_vko_rank", "median_seed_vko_rank",
        ]]
        .rename(columns={
            "ko_label": "gene_symbol",
            "compartment": "joint_context_compartment",
            "mean_detection_fraction": "joint_context_detection_fraction",
        })
    )
    t = t.merge(max_detection, on="gene_symbol", how="left").merge(joint, on="gene_symbol", how="left")
    t["liver_context"] = t["same_compartment_context_score"]
    t = t.merge(lincs_metrics[lcols].drop_duplicates("gene_symbol"), on="gene_symbol", how="left")
    reliability_t = np.where(t["is_landmark"].fillna(0).eq(1), 1.0, np.where(t["is_bing"].fillna(0).eq(1), 0.8, 0.5))
    t["HepG2_perturbation"] = ((0.60 * percentile_score(t["spearman_rho"].abs()) + 0.40 * percentile_score(t["max_abs_z"])) * reliability_t).where(t["spearman_rho"].notna())
    t = t.merge(string_metrics, on="gene_symbol", how="left")
    t["STRING_context"] = 0.65 * percentile_score(t["max_path_product_to_opposite_class"]) + 0.35 * percentile_score(t["weighted_degree"])
    t = t.merge(sci, on="gene_symbol", how="left")
    t["SciPlex_replication_flag"] = t["observed_in_both_doses"].fillna(False).astype(bool) & t["direction_concordant"].fillna(False).astype(bool)
    target_domains = ["pharmacology", "liver_context", "HepG2_perturbation", "STRING_context"]
    eligible_t = t[target_domains].notna().all(axis=1)
    t["pareto_front"] = np.nan
    t.loc[eligible_t, "pareto_front"] = pareto_front_ids(t.loc[eligible_t], target_domains)
    t, target_loso = consensus_and_loso(t, target_domains, "exposure_proximal_target")
    strong_pharm = t["pharmacology_tier"].isin(["canonical_target", "multi_document_noncanonical"])
    t["hepg2_measurement_valid"] = (
        t["is_landmark"].fillna(0).eq(1)
        & ~t.get("constant_profile", pd.Series(False, index=t.index)).fillna(True).astype(bool)
        & t["spearman_rho"].notna()
    )
    t["hepg2_alternative_gate"] = t["hepg2_measurement_valid"] & t["spearman_rho"].abs().ge(0.70)
    t["strong_pharmacology_gate"] = strong_pharm
    t["primary_T2_gate"] = strong_pharm & (
        t["same_compartment_primary_gate"].fillna(False)
        | t["hepg2_alternative_gate"].fillna(False)
    )
    t["strict_all_seed_T2_sensitivity_gate"] = strong_pharm & (
        t["same_compartment_all_seed_sensitivity_gate"].fillna(False)
        | t["hepg2_alternative_gate"].fillna(False)
    )
    t["evidence_tier"] = "T3: screen/context target hypothesis"
    t.loc[t["primary_T2_gate"], "evidence_tier"] = "T2: same-compartment target-network hypothesis"
    t["candidate_role"] = "exposure-proximal pharmacologic target"
    t["causal_status"] = "not established"
    gate_cols = [
        "gene_symbol", "pharmacology_tier", "strong_pharmacology_gate",
        "joint_context_compartment", "joint_context_detection_fraction",
        "vko_within_compartment_rank", "n_seeds", "same_compartment_primary_gate",
        "n_seeds_top_quartile", "min_seed_vko_rank",
        "same_compartment_all_seed_sensitivity_gate", "is_landmark", "constant_profile",
        "spearman_rho", "hepg2_measurement_valid", "hepg2_alternative_gate",
        "primary_T2_gate", "strict_all_seed_T2_sensitivity_gate", "evidence_tier",
    ]
    t[[c for c in gate_cols if c in t.columns]].sort_values(
        ["primary_T2_gate", "strict_all_seed_T2_sensitivity_gate", "gene_symbol"],
        ascending=[False, False, True],
    ).to_csv(OUT / "target_tier_gate_sensitivity.csv", index=False)

    loso = pd.concat([anchor_loso, target_loso], ignore_index=True)
    a["candidate_class"] = "downstream_DILI_anchor"
    t["candidate_class"] = "exposure_proximal_target"
    combined = pd.concat([a, t], ignore_index=True, sort=False)

    # Long-form domain table keeps source group boundaries explicit.
    domain_rows = []
    for row in combined.itertuples(index=False):
        gene = row.gene_symbol
        cls = row.candidate_class
        domains = anchor_domains if cls == "downstream_DILI_anchor" else target_domains
        for domain in domains:
            domain_rows.append({
                "candidate_class": cls,
                "gene_symbol": gene,
                "source_group": domain,
                "domain_score_0_to_1": getattr(row, domain),
                "independence_note": {
                    "clinical_DILI": "one human DILI proteomics publication; subcohorts/features are not counted as independent sources",
                    "pharmacology": "one ChEMBL domain; multiple documents strengthen within-domain confidence",
                    "liver_context": "GSE115469 expression and scTenifold displacement are grouped as one source and combined by a same-compartment maximin score",
                    "HepG2_perturbation": "LINCS L1000 HepG2 six-dose Level-5 signatures",
                    "STRING_context": "STRING association network; path evidence is contextual, not causal",
                }[domain],
            })
    domain_long = pd.DataFrame(domain_rows)

    keep = [
        "candidate_class", "candidate_role", "gene_symbol", "evidence_tier", "causal_status",
        "pareto_front", "consensus_rank", "loso_rank_best", "loso_rank_median", "loso_rank_worst",
        "clinical_DILI", "pharmacology", "liver_context", "HepG2_perturbation", "STRING_context",
        "SciPlex_replication_flag", "dili_vs_non_dili_specificity_signal", "specificity_replicates_discovery_direction",
        "proteomics_tier_original", "mean_detection_fraction", "max_primary_detection", "median_dili_anchor_percentile",
        "joint_context_compartment", "joint_context_detection_fraction", "vko_within_compartment_rank",
        "detection_within_compartment_rank", "same_compartment_context_score",
        "same_compartment_primary_gate", "same_compartment_all_seed_sensitivity_gate",
        "n_seed_results", "n_seeds_top_quartile", "min_seed_vko_rank", "n_seeds",
        "hepg2_measurement_valid", "hepg2_alternative_gate", "strong_pharmacology_gate",
        "primary_T2_gate", "strict_all_seed_T2_sensitivity_gate", "max_pchembl", "n_documents",
        "pharmacology_tier", "spearman_rho", "exact_p", "bh_q", "max_abs_z", "is_landmark", "is_bing", "constant_profile",
        "max_path_product_to_opposite_class", "nearest_opposite_node", "shortest_hops",
    ]
    keep = [c for c in keep if c in combined.columns]
    combined[keep].sort_values(["candidate_class", "consensus_rank", "gene_symbol"]).to_csv(OUT / "candidate_consensus_ranking.csv", index=False)
    domain_long.to_csv(OUT / "candidate_domain_scores_long.csv", index=False)
    loso.to_csv(OUT / "candidate_leave_one_source_out_ranks.csv", index=False)
    shortlist = combined.loc[combined["evidence_tier"].str.startswith("T2"), [c for c in keep if c in combined.columns]]
    shortlist.to_csv(OUT / "tier2_candidate_shortlist.csv", index=False)
    phenotype = combined.loc[combined["candidate_class"].eq("downstream_DILI_anchor"), [c for c in keep if c in combined.columns]]
    phenotype.to_csv(OUT / "phenotype_anchor_tiers.csv", index=False)
    return combined, domain_long, loso


def configure_style() -> None:
    sns.set_theme(style="ticks", context="paper")
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 7.4,
        "axes.titlesize": 7.7,
        "axes.labelsize": 7.1,
        "xtick.labelsize": 6.8,
        "ytick.labelsize": 6.8,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.2,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def dotmap(ax: mpl.axes.Axes, frame: pd.DataFrame, genes: list[str], compartments: list[str], title: str) -> None:
    sub = frame.loc[frame["gene_symbol"].isin(genes) & frame["compartment"].isin(compartments)].copy()
    sub["gene_symbol"] = pd.Categorical(sub["gene_symbol"], categories=genes[::-1], ordered=True)
    sub["compartment"] = pd.Categorical(sub["compartment"], categories=compartments, ordered=True)
    x = sub["compartment"].cat.codes
    y = sub["gene_symbol"].cat.codes
    size = 12 + 180 * sub["mean_detection_fraction"].clip(0, 1)
    color = np.log1p(sub["mean_expression"])
    sca = ax.scatter(x, y, s=size, c=color, cmap="viridis", edgecolor="white", linewidth=0.35)
    ax.set_xticks(range(len(compartments)), [c.replace("Immune_other", "Other immune") for c in compartments], rotation=35, ha="right")
    ax.set_yticks(range(len(genes)), genes[::-1])
    ax.set_xlim(-0.6, len(compartments) - 0.4)
    ax.set_ylim(-0.6, len(genes) - 0.4)
    ax.grid(color="#EEEEEE", linewidth=0.6)
    ax.set_title(title, loc="left", fontweight="bold")
    cbar = plt.colorbar(sca, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("log(1 + mean expression)", fontsize=6.8)
    handles = [plt.scatter([], [], s=12 + 180 * s, facecolor="#777777", edgecolor="white", label=f"{int(s*100)}%") for s in [0.1, 0.5, 0.9]]
    ax.legend(handles=handles, title="Detected cells", ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.40), frameon=False, handletextpad=0.2, columnspacing=0.7)


def save_figure(fig: mpl.figure.Figure, stem: str) -> dict[str, str]:
    outputs = {}
    for ext in ["svg", "pdf"]:
        path = FIG_OUT / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight")
        outputs[ext] = sha256(path)
    preview = FIG_OUT / f"{stem}_preview.png"
    fig.savefig(preview, dpi=180, bbox_inches="tight")
    outputs["preview_png"] = sha256(preview)
    png = FIG_OUT / f"{stem}_600dpi.png"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    outputs["png_600dpi"] = sha256(png)
    tif = FIG_OUT / f"{stem}_600dpi.tiff"
    fig.savefig(tif, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    outputs["tiff_600dpi"] = sha256(tif)
    return outputs


def figure4(data: dict[str, pd.DataFrame], matched_stable: pd.DataFrame, jaccard: pd.DataFrame) -> dict[str, str]:
    fig = plt.figure(figsize=(8.0, 11.2), constrained_layout=True)
    gs = fig.add_gridspec(5, 2, height_ratios=[1.0, 1.25, 1.15, 1.0, 1.15])
    axes = [fig.add_subplot(gs[i, j]) for i in range(5) for j in range(2)]

    # A-C: within-compartment descriptive PCA; embeddings were fitted separately.
    pca = data["pca"]
    comp_palette = {
        "Hepatocyte_1": "#0B6E75", "Hepatocyte_2": "#2A9D8F", "Hepatocyte_3": "#62B6A6",
        "Hepatocyte_4": "#E9C46A", "Hepatocyte_5": "#F4A261", "Hepatocyte_6": "#E76F51",
        "Central_venous_LSECs": "#3B6FB6", "Periportal_LSECs": "#76A5D5", "Portal_endothelial_Cells": "#B7D4EF",
        "Inflammatory_Macrophage": "#B33A3A", "Non-inflammatory_Macrophage": "#E89A8D",
    }
    for ax, comp, letter in zip(axes[:3], PRIMARY_COMPARTMENTS, ["A", "B", "C"]):
        sub = pca.loc[pca["compartment"].eq(comp)]
        for ct, g in sub.groupby("CellType"):
            ax.scatter(g["PC1"], g["PC2"], s=4.5, alpha=0.55, color=comp_palette.get(ct, "#777777"), label=ct.replace("_", " "), linewidth=0)
        ax.set_title(f"{comp}: descriptive PCA", loc="left", fontweight="bold")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.legend(frameon=False, markerscale=2.2, ncol=2 if comp == "Hepatocyte" else 1, loc="best")
        clean_axis(ax)
        panel_label(ax, letter)

    sc = data["sc_expr"]
    direct = CANONICAL_TARGETS + ["FLT3", "AXL", "LYN", "JAK1", "ABL1", "TYK2"]
    dotmap(axes[3], sc, direct, ALL_COMPARTMENTS, "Direct-target localization across liver compartments")
    panel_label(axes[3], "D")
    dotmap(axes[4], sc, DILI_ANCHORS, ALL_COMPARTMENTS, "DILI-anchor localization across liver compartments")
    panel_label(axes[4], "E")

    # F: donor-level detection reproducibility for compartment-relevant canonical targets.
    donor = data["sc_donor"]
    focus = [("KDR", "Endothelial"), ("FLT1", "Endothelial"), ("FLT4", "Endothelial"),
             ("FGFR1", "Endothelial"), ("FGFR2", "Hepatocyte"), ("FGFR3", "Hepatocyte"),
             ("PDGFRA", "Stellate"), ("PDGFRB", "Stellate")]
    frows = []
    for gene, comp in focus:
        z = donor.loc[(donor["gene_symbol"].eq(gene)) & (donor["compartment"].eq(comp))]
        if len(z):
            frows.append({"label": f"{gene} | {comp}", "mean": z["detection_fraction"].mean(),
                          "min": z["detection_fraction"].min(), "max": z["detection_fraction"].max(), "n_donors": z["donor"].nunique()})
    f = pd.DataFrame(frows).sort_values("mean")
    ax = axes[5]
    y = np.arange(len(f))
    ax.hlines(y, f["min"], f["max"], color="#AAB3BB", linewidth=2)
    ax.scatter(f["mean"], y, color=PALETTE["blue"], s=32, zorder=3)
    ax.set_yticks(y, f["label"])
    ax.set_xlabel("Detection fraction across donors (mean; range)")
    ax.set_xlim(-0.02, 1.0)
    ax.set_title("Target detection across liver donors", loc="left", fontweight="bold")
    ax.text(0.99, 0.02, "Stellate n=37 cells: exploratory", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8, color=PALETTE["gray"])
    clean_axis(ax); panel_label(ax, "F")

    # G: median aggregate DILI-anchor displacement percentiles.
    stable = data["vko_stable"].loc[data["vko_stable"]["ko_class"].eq("drug_target")].copy()
    genes = stable.groupby("ko_label")["median_dili_anchor_percentile"].max().nlargest(13).index.tolist()
    heat = stable.loc[stable["ko_label"].isin(genes)].pivot(index="ko_label", columns="compartment", values="median_dili_anchor_percentile").reindex(index=genes, columns=PRIMARY_COMPARTMENTS)
    sns.heatmap(heat, ax=axes[6], cmap="mako", vmin=0, vmax=0.65, annot=True, fmt=".2f", cbar_kws={"label": "Median displacement percentile"}, linewidths=0.35, linecolor="white")
    axes[6].set_xlabel("Blank = no eligible KO after compartment-specific filtering", fontsize=6.8, color=PALETTE["gray"]); axes[6].set_ylabel("")
    axes[6].set_title("Compartment-dependent unsigned virtual-KO ranking", loc="left", fontweight="bold")
    panel_label(axes[6], "G")

    # H: seed-level values plus top-20 Jaccard stability.
    seed = data["vko_seed"]
    chosen = ["JAK1", "MET", "TYK2", "ABL1", "FLT1", "KDR+FLT1+FGFR2+FGFR3"]
    sub = seed.loc[(seed["compartment"].eq("Hepatocyte")) & seed["ko_label"].isin(chosen)].copy()
    for label, g in sub.groupby("ko_label"):
        short = "canonical poly-target" if "+" in label else label
        axes[7].plot(g["seed"].astype(str), g["dili_anchor_mean_percentile"], marker="o", ms=4, label=short)
    axes[7].axhline(0.5, color="#999999", linestyle="--", linewidth=0.9, label="50th percentile")
    axes[7].set_ylabel("Mean DILI-anchor displacement percentile")
    axes[7].set_xlabel("Network seed")
    axes[7].set_ylim(0, 0.72)
    axes[7].set_title("Three-seed sensitivity in hepatocytes", loc="left", fontweight="bold")
    axes[7].legend(frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=6.2)
    clean_axis(axes[7]); panel_label(axes[7], "H")

    # I: matched-null calibration.
    m = matched_stable.loc[matched_stable["ko_class"].eq("drug_target")].copy()
    mgenes = m.groupby("ko_label")["median_matched_z"].apply(lambda x: x.abs().max()).nlargest(12).index.tolist()
    mh = m.loc[m["ko_label"].isin(mgenes)].pivot(index="ko_label", columns="compartment", values="median_matched_z").reindex(index=mgenes, columns=PRIMARY_COMPARTMENTS)
    lim = max(1.5, float(np.nanmax(np.abs(mh.to_numpy()))))
    sns.heatmap(mh, ax=axes[8], cmap="vlag", center=0, vmin=-lim, vmax=lim, annot=True, fmt=".1f", cbar_kws={"label": "Median matched-null z"}, linewidths=0.35, linecolor="white")
    axes[8].set_xlabel("Blank = no eligible KO after compartment-specific filtering", fontsize=6.8, color=PALETTE["gray"]); axes[8].set_ylabel("")
    axes[8].set_title("Exploratory matched-null calibration of anchor ranking", loc="left", fontweight="bold")
    panel_label(axes[8], "I")

    # J: pairwise top-20 overlap exposes gene-level stability limits.
    j = jaccard.loc[jaccard["ko_class"].isin(["drug_target", "canonical_poly_target"])].copy()
    sns.boxplot(data=j, x="compartment", y="jaccard_top20", hue="ko_class", ax=axes[9],
                palette={"drug_target": PALETTE["blue"], "canonical_poly_target": PALETTE["orange"]},
                showfliers=False, linewidth=0.8)
    sns.stripplot(data=j, x="compartment", y="jaccard_top20", hue="ko_class", ax=axes[9], dodge=True,
                  palette={"drug_target": PALETTE["blue"], "canonical_poly_target": PALETTE["orange"]}, alpha=0.45, size=2.5, linewidth=0)
    handles, labels = axes[9].get_legend_handles_labels()
    axes[9].legend(handles[:2], ["single target", "canonical poly-target"], frameon=False, loc="upper right")
    axes[9].set_ylim(-0.02, 1.02)
    axes[9].set_xlabel("")
    axes[9].set_ylabel("Pairwise Jaccard index of top-20 genes")
    axes[9].set_title("Top-gene stability across network seeds", loc="left", fontweight="bold")
    clean_axis(axes[9]); panel_label(axes[9], "J")

    fig.suptitle("Figure 4 | Human-liver localization and uncertainty-aware virtual perturbation", fontsize=11, fontweight="bold", x=0.01, ha="left")
    return save_figure(fig, "Figure_4_Liver_Localization_and_Virtual_Knockout")


def figure5(
    data: dict[str, pd.DataFrame],
    lincs_metrics: pd.DataFrame,
    lincs_long: pd.DataFrame,
    combined: pd.DataFrame,
    loso: pd.DataFrame,
    graph: nx.Graph,
) -> dict[str, str]:
    fig = plt.figure(figsize=(8.0, 11.2), constrained_layout=True)
    gs = fig.add_gridspec(5, 2, height_ratios=[1.05, 1.0, 1.15, 1.0, 1.15])
    axes = [fig.add_subplot(gs[i, j]) for i in range(5) for j in range(2)]
    prot = data["prot"].set_index("protein").reindex(DILI_ANCHORS)

    # A: clinical evidence dimensions (one source group; shown separately, not counted as independent datasets).
    clinical = pd.DataFrame({
        "Etiology\nDO vs NDO": prot["confirmatory_DO_vs_NDO_roc_auc_separation"],
        "Acute injury\nDO vs HV": prot["confirmatory_DO_vs_HV_roc_auc_separation"],
        "Recovery\nDO vs DF": prot["confirmatory_DO_vs_DF_roc_auc_separation"],
        "ALT severity\n|rho|": prot["within_DO_alt_spearman_rho"].abs(),
    })
    clinical = clinical.sort_values("Etiology\nDO vs NDO", ascending=False)
    sns.heatmap(clinical, ax=axes[0], cmap="crest", vmin=0.0, vmax=1.0, annot=True, fmt=".2f",
                linewidths=0.35, linecolor="white", cbar_kws={"label": "AUC / |Spearman rho|"})
    axes[0].set_xlabel(""); axes[0].set_ylabel("")
    axes[0].set_title("Human DILI evidence dimensions", loc="left", fontweight="bold")
    panel_label(axes[0], "A")

    # B: liver expression versus anchor-network displacement.
    anchors = combined.loc[combined["candidate_class"].eq("downstream_DILI_anchor")].copy()
    ax = axes[1]
    size = 35 + 170 * anchors["clinical_DILI"].fillna(0)
    ax.scatter(anchors["mean_detection_fraction"], anchors["median_dili_anchor_percentile"],
               s=size, c=anchors["clinical_DILI"], cmap="viridis", edgecolor="white", linewidth=0.6)
    for row in anchors.itertuples(index=False):
        if row.evidence_tier.startswith("P1") or row.gene_symbol in ["GSTA1", "OTC", "ALDOB", "ASS1"]:
            ax.text(row.mean_detection_fraction + 0.01, row.median_dili_anchor_percentile + 0.005, row.gene_symbol, fontsize=6.8)
    ax.axhline(0.5, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Hepatocyte detection fraction")
    ax.set_ylabel("Median anchor-displacement percentile")
    ax.set_title("Anchor localization and network rank", loc="left", fontweight="bold")
    clean_axis(ax); panel_label(ax, "B")

    # C: six-dose HepG2 profiles for DILI anchors.
    z = lincs_long.loc[lincs_long["gene_symbol"].isin(DILI_ANCHORS)].copy()
    heat = z.pivot_table(index="gene_symbol", columns="dose_um", values="level5_zscore", aggfunc="first")
    order = anchors.sort_values("consensus_rank")["gene_symbol"].tolist()
    heat = heat.reindex(index=[g for g in order if g in heat.index])
    lim = max(2.0, float(np.nanmax(np.abs(heat.to_numpy()))))
    sns.heatmap(heat, ax=axes[2], cmap="vlag", center=0, vmin=-lim, vmax=lim, annot=True, fmt=".1f",
                linewidths=0.35, linecolor="white", cbar_kws={"label": "Level-5 z"})
    axes[2].set_xlabel("Nintedanib dose (µM)"); axes[2].set_ylabel("")
    axes[2].set_title("Observed HepG2 transcriptional response across six doses", loc="left", fontweight="bold")
    panel_label(axes[2], "C")

    # D: P1 plus exposure-responsive P2 anchors; no smoothing or direction claim.
    profile_genes = anchors.loc[
        anchors["evidence_tier"].str.startswith("P1")
        | (anchors["evidence_tier"].str.startswith("P2") & anchors["spearman_rho"].abs().ge(0.80)),
        "gene_symbol",
    ].tolist()
    for gene in profile_genes:
        g = z.loc[z["gene_symbol"].eq(gene)].sort_values("dose_um")
        axes[3].plot(g["dose_um"], g["level5_zscore"], marker="o", ms=4, label=gene)
    axes[3].axhline(0, color="#999999", linewidth=0.8)
    axes[3].set_xscale("log")
    axes[3].set_xlabel("Nintedanib dose (µM; log scale)")
    axes[3].set_ylabel("LINCS Level-5 z-score")
    axes[3].set_title("Heterogeneous P1/P2 dose profiles", loc="left", fontweight="bold")
    axes[3].legend(frameon=False, ncol=2)
    clean_axis(axes[3]); panel_label(axes[3], "D")

    # E: all overlap-conditional genes establish the global cross-model boundary;
    # FBP1 is highlighted because it is the sole DILI anchor present in both
    # truncated Sci-Plex dose lists.
    sci = data["sciplex"].copy()
    pairs = [(1.11, "1uM", PALETTE["blue"]), (10.0, "10uM", PALETTE["orange"])]
    annotations = []
    for dose, sci_col, color in pairs:
        h = lincs_long.loc[lincs_long["dose_um"].eq(dose), ["gene_symbol", "level5_zscore"]]
        m = h.merge(sci[["gene_symbol", sci_col]].dropna(), on="gene_symbol")
        rho, p = spearmanr(m["level5_zscore"], m[sci_col]) if len(m) >= 3 else (np.nan, np.nan)
        axes[4].scatter(m["level5_zscore"], m[sci_col], s=7, alpha=0.25, color=color, edgecolor="none", label=f"{dose:g} vs {sci_col} (n={len(m)})")
        annotations.append(f"{dose:g} µM: rho={rho:.2f}")
        fbp = m.loc[m["gene_symbol"].eq("FBP1")]
        if len(fbp):
            axes[4].scatter(fbp["level5_zscore"], fbp[sci_col], s=55, facecolor="none", edgecolor=color, linewidth=1.4)
            axes[4].text(float(fbp["level5_zscore"].iloc[0]), float(fbp[sci_col].iloc[0]), " FBP1", fontsize=6.8)
    axes[4].axhline(0, color="#BBBBBB", linewidth=0.6); axes[4].axvline(0, color="#BBBBBB", linewidth=0.6)
    axes[4].set_xlabel("HepG2 LINCS Level-5 z-score")
    axes[4].set_ylabel("MCF7 Sci-Plex standardized value")
    axes[4].set_title("FBP1 is the sole cross-context anchor overlap", loc="left", fontweight="bold")
    axes[4].text(0.02, 0.98, "\n".join(annotations), transform=axes[4].transAxes, va="top", ha="left", fontsize=6.8)
    axes[4].legend(frameon=False, loc="lower right")
    clean_axis(axes[4]); panel_label(axes[4], "E")

    # F: explicit source-group matrix for shortlisted candidates from both roles.
    targets = combined.loc[combined["candidate_class"].eq("exposure_proximal_target")].sort_values("consensus_rank")
    top_targets = targets.head(8)["gene_symbol"].tolist()
    top_anchors = anchors.sort_values("consensus_rank").head(8)["gene_symbol"].tolist()
    display = combined.loc[combined["gene_symbol"].isin(top_anchors + top_targets)].copy()
    display["label"] = display["gene_symbol"] + np.where(display["candidate_class"].eq("downstream_DILI_anchor"), " [R]", " [T]")
    matrix = display.set_index("label")[["clinical_DILI", "pharmacology", "liver_context", "HepG2_perturbation", "STRING_context"]]
    matrix = matrix.reindex([f"{g} [R]" for g in top_anchors] + [f"{g} [T]" for g in top_targets])
    matrix = matrix.rename(columns={
        "clinical_DILI": "Clinical DILI", "pharmacology": "Pharmacology",
        "liver_context": "Liver context", "HepG2_perturbation": "HepG2",
        "STRING_context": "STRING",
    })
    cmap = LinearSegmentedColormap.from_list("evidence", ["#F4F4F4", "#BFDDE0", "#007C83"])
    axes[5].set_facecolor("#ECEFF1")
    sns.heatmap(matrix, ax=axes[5], cmap=cmap, vmin=0, vmax=1, annot=True, fmt=".2f", mask=matrix.isna(),
                linewidths=0.35, linecolor="white", cbar_kws={"label": "Within-role domain score"})
    axes[5].set_xlabel(""); axes[5].set_ylabel("")
    axes[5].set_title("Independent evidence domains preserve candidate role", loc="left", fontweight="bold")
    # Role suffixes are already encoded in row labels; a footer here crowds the
    # rotated domain labels and is therefore intentionally omitted.
    panel_label(axes[5], "F")

    # G: non-compensatory Pareto view for response anchors.
    ax = axes[6]
    eligible = anchors.dropna(subset=["pareto_front"]).copy()
    yscore = eligible[["liver_context", "HepG2_perturbation", "STRING_context"]].mean(axis=1)
    colors = eligible["pareto_front"].map({1: PALETTE["red"], 2: PALETTE["orange"]}).fillna(PALETTE["gray"])
    ax.scatter(eligible["clinical_DILI"], yscore, c=colors, s=55, edgecolor="white", linewidth=0.5)
    for gene, xval, yval, front in zip(eligible["gene_symbol"], eligible["clinical_DILI"], yscore, eligible["pareto_front"]):
        if front <= 2 or gene == "FBP1":
            ax.text(xval + 0.01, yval + 0.008, gene, fontsize=6.8)
    ax.set_xlabel("Clinical DILI domain score")
    ax.set_ylabel("Mean contextual/perturbation score")
    ax.set_title("Pareto fronts avoid compensating away weak evidence", loc="left", fontweight="bold")
    ax.legend(handles=[Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE['red'], label='Front 1', markersize=6),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE['orange'], label='Front 2', markersize=6),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE['gray'], label='Later fronts', markersize=6)],
              frameon=False, loc="best")
    clean_axis(ax); panel_label(ax, "G")

    # H-I: leave-one-source-out rank intervals.
    def rank_interval_panel(ax: mpl.axes.Axes, frame: pd.DataFrame, title: str, letter: str, role: str) -> None:
        show = frame.sort_values("consensus_rank").head(12).sort_values("consensus_rank", ascending=False)
        y = np.arange(len(show))
        ax.hlines(y, show["loso_rank_best"], show["loso_rank_worst"], color="#AAB3BB", linewidth=2)
        if role == "anchor":
            colors = np.where(
                show["evidence_tier"].str.startswith("P1"), PALETTE["red"],
                np.where(show["evidence_tier"].str.startswith("P2"), PALETTE["orange"], PALETTE["blue"]),
            )
        else:
            colors = np.where(show["evidence_tier"].str.startswith("T2"), PALETTE["red"], PALETTE["blue"])
        ax.scatter(show["consensus_rank"], y, c=colors, s=35, zorder=3, edgecolor="white", linewidth=0.4)
        ax.set_yticks(y, show["gene_symbol"])
        ax.invert_xaxis()
        ax.set_xlabel("Rank (lower is better; interval = leave-one-source-out)")
        ax.set_title(title, loc="left", fontweight="bold")
        clean_axis(ax); panel_label(ax, letter)
    rank_interval_panel(axes[7], anchors, "P-track ranks under leave-one-source-out", "H", "anchor")
    rank_interval_panel(axes[8], targets, "T-track ranks use matched liver compartments", "I", "target")

    # J: evidence network restricted to P1, exposure-responsive P2, and T2.
    p1 = anchors.loc[anchors["evidence_tier"].str.startswith("P1"), "gene_symbol"].tolist()
    p2_exposure = anchors.loc[
        anchors["evidence_tier"].str.startswith("P2") & anchors["spearman_rho"].abs().ge(0.80),
        "gene_symbol",
    ].tolist()
    tier2_t = targets.loc[targets["evidence_tier"].str.startswith("T2"), "gene_symbol"].tolist()
    selected = set(p1 + p2_exposure + tier2_t)
    # Add at most one best STRING path per prioritized response/target pair.
    for a_gene in p1 + p2_exposure:
        for t_gene in tier2_t:
            if a_gene in graph and t_gene in graph:
                try:
                    path = nx.shortest_path(graph, a_gene, t_gene, weight="distance")
                except nx.NetworkXNoPath:
                    continue
                # Include the prespecified five-hop target–anchor routes reported
                # in the source table; longer paths remain excluded for legibility.
                if len(path) <= 6:
                    selected.update(path)
    subgraph = graph.subgraph(selected).copy()
    ax = axes[9]
    if len(subgraph):
        pos = nx.spring_layout(subgraph, seed=SEED, weight="score", k=1.1 / math.sqrt(max(len(subgraph), 1)))
        node_colors = []
        node_sizes = []
        for node in subgraph.nodes:
            if node in p1:
                node_colors.append(PALETTE["red"]); node_sizes.append(260)
            elif node in p2_exposure:
                node_colors.append(PALETTE["orange"]); node_sizes.append(240)
            elif node in tier2_t:
                node_colors.append(PALETTE["blue"]); node_sizes.append(260)
            else:
                node_colors.append("#D5DADD"); node_sizes.append(150)
        widths = [0.5 + 2.6 * subgraph[u][v]["score"] for u, v in subgraph.edges]
        nx.draw_networkx_edges(subgraph, pos, ax=ax, width=widths, alpha=0.55, edge_color="#737B80")
        nx.draw_networkx_nodes(subgraph, pos, ax=ax, node_color=node_colors, node_size=node_sizes, edgecolors="white", linewidths=0.7)
        nx.draw_networkx_labels(subgraph, pos, ax=ax, font_size=6.8)
        cross_edges = [
            (u, v) for u, v in subgraph.edges
            if ((u in (p1 + p2_exposure) and v in tier2_t) or (v in (p1 + p2_exposure) and u in tier2_t))
        ]
        if not cross_edges:
            ax.text(0.99, 0.02, "No direct cross-role edge at STRING score ≥0.4",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8, color=PALETTE["gray"])
    ax.set_axis_off()
    ax.set_title("Sparse STRING associations among prioritized nodes", loc="left", fontweight="bold")
    ax.legend(handles=[Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE['red'], label='P1 phenotype anchor', markersize=7),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE['orange'], label='Exposure-responsive P2 anchor', markersize=7),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE['blue'], label='Exploratory T2-gate candidate', markersize=7),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='#D5DADD', label='STRING bridge', markersize=7)],
              frameon=False, loc="upper left", fontsize=6.0)
    panel_label(ax, "J")

    fig.suptitle("Figure 5 | Role-separated integration prioritizes one unconfirmed exposure-proximal hypothesis", fontsize=10.5, fontweight="bold", x=0.01, ha="left")
    return save_figure(fig, "Figure_5_Multi_Evidence_Target_Prioritization")


def write_captions(combined: pd.DataFrame, matched: pd.DataFrame, jaccard: pd.DataFrame) -> None:
    p1 = combined.loc[(combined["candidate_class"].eq("downstream_DILI_anchor")) & combined["evidence_tier"].str.startswith("P1"), "gene_symbol"].tolist()
    p2 = combined.loc[(combined["candidate_class"].eq("downstream_DILI_anchor")) & combined["evidence_tier"].str.startswith("P2"), "gene_symbol"].tolist()
    p3 = combined.loc[(combined["candidate_class"].eq("downstream_DILI_anchor")) & combined["evidence_tier"].str.startswith("P3"), "gene_symbol"].tolist()
    tier2_t = combined.loc[(combined["candidate_class"].eq("exposure_proximal_target")) & combined["evidence_tier"].str.startswith("T2"), "gene_symbol"].tolist()
    strict_t2 = combined.loc[
        combined["candidate_class"].eq("exposure_proximal_target")
        & combined["strict_all_seed_T2_sensitivity_gate"].fillna(False),
        "gene_symbol",
    ].tolist()
    med_j = jaccard.loc[jaccard["ko_class"].eq("drug_target"), "jaccard_top20"].median()
    matched_sig = matched.loc[matched["n_seed_bh_q_lt_0_05"].gt(0), ["compartment", "ko_label", "n_seed_bh_q_lt_0_05"]]

    caption4 = f"""# Figure 4. Human-liver localization and uncertainty-aware virtual perturbation

**(A-C)** Donor-capped descriptive principal-component embeddings for hepatocytes, endothelial cells and macrophages from GSE115469 (up to 250 cells per donor for display). This cap limited dominance by large donors but did not create equal donor weights. Embeddings were fitted separately within each compartment and are not on a shared coordinate system. **(D)** Mean expression (color) and detection fraction (area) of canonical and screened noncanonical nintedanib pharmacologic targets across six liver compartments. **(E)** Corresponding localization of 13 circulating DILI protein anchors. **(F)** Donor-level target detection means and ranges; stellate-cell estimates are exploratory because only 37 stellate cells were available. **(G)** Median percentile rank of the aggregate DILI-anchor displacement after each single-gene virtual knockout across three seeds. **(H)** Seed-level hepatocyte estimates for selected pharmacologic targets and the canonical multi-target knockout; the dashed line marks the 50th percentile, not a causal null. **(I)** Exploratory median z score relative to expression- and detection-matched gene-set permutations. {len(matched_sig)} target-compartment combinations had one seed with within-compartment BH q<0.05, but none reproduced BH q<0.05 in at least two seeds. In G and I, blank heat-map cells indicate that no eligible virtual-knockout estimate was available after compartment-specific gene filtering; blanks are not zeros. **(J)** Pairwise overlap of top-20 displaced genes across seeds (median Jaccard among single-target KOs, {med_j:.2f}), documenting limited gene-level stability. scTenifoldKnk distances are unsigned network displacement measures. They do not identify protective direction, molecular fold change, or a clinical causal effect.
"""
    caption5 = f"""# Figure 5. Multi-source prioritization distinguishes response anchors from pharmacologic targets

**(A)** Four prespecified clinical properties recalculated from the public DILI proteomics workbook: etiology discrimination (DILI onset [DO] versus non-DILI liver injury onset [NDO]), acute injury detection (DO versus healthy volunteers [HV]), onset-to-follow-up recovery (DO versus DILI follow-up [DF]), and within-DO association with ALT. These are subfeatures of one clinical source group. **(B)** Hepatocyte detection versus the three-seed median percentile of aggregate DILI-anchor displacement after virtual knockout of each anchor; point area/color encodes the composite clinical DILI score. **(C)** Six-dose, 24-hour HepG2 LINCS Level-5 profiles for available DILI anchors. **(D)** Unsmoothed dose profiles for the P1 anchor and P2 anchors with |six-dose rho|≥0.80; P labels describe phenotype evidence and are not drug-target tiers. **(E)** Overlap-conditional comparison with the truncated MCF7 Sci-Plex lists. FBP1 was the sole DILI anchor appearing at both Sci-Plex doses; across all listed overlapping genes, Spearman rho was 0.086 near 1 µM and 0.140 at 10 µM. Absence from Sci-Plex denotes absence from the published top/bottom lists rather than no response. **(F)** Explicit source-group matrix for leading response anchors [R] and exposure-proximal targets [T]. ChEMBL documents are nested within one pharmacology domain. For T-track candidates, liver context is the maximum same-compartment maximin of expression rank and virtual-knockout rank; values from different compartments are never spliced. **(G)** Non-dominated Pareto fronts for response anchors. **(H)** Consensus ranks and leave-one-source-out ranges for the 12 leading P-track candidates. **(I)** Corresponding ranges for the 12 leading T-track candidates using same-compartment liver context. The prespecified exploratory T2 gate required strong pharmacology plus either detection ≥10%, three-seed availability, and a median aggregate-anchor virtual-KO rank in the top quartile in the same compartment, or a valid directly measured HepG2 landmark with |rho|≥0.70. FGFR1 alone met this exploratory gate, with endothelial seed-specific virtual-KO ranks of 0.083, 0.917, and 1.000; the repeated computations are stability probes, not biological replicates. **(J)** Sparse STRING paths among prioritized P1/exposure-responsive P2 anchors and exploratory T2 hypotheses; STRING edges are association evidence and do not establish a toxicity pathway. P1: {', '.join(p1)}. P2: {', '.join(p2)}. P3: {', '.join(p3)}. Candidate meeting the exploratory T2 gate: {', '.join(tier2_t) if tier2_t else 'none'}. Under the stricter sensitivity requiring top-quartile rank in every seed, T2: {', '.join(strict_t2) if strict_t2 else 'none'}. No T1 target was assigned because nintedanib-specific human liver perturbation and prospective functional validation were unavailable.
"""
    (FIG_OUT / "Figure_4_caption.md").write_text(caption4)
    caption5 = caption5.replace(
        "Sparse STRING paths among",
        "Sparse STRING paths of no more than five hops among",
    )
    (FIG_OUT / "Figure_5_caption.md").write_text(caption5)


def write_methods_results(combined: pd.DataFrame, matched: pd.DataFrame, jaccard: pd.DataFrame) -> None:
    anchors = combined.loc[combined["candidate_class"].eq("downstream_DILI_anchor")].copy()
    targets = combined.loc[combined["candidate_class"].eq("exposure_proximal_target")].copy()
    p1 = anchors.loc[anchors["evidence_tier"].str.startswith("P1")].sort_values("consensus_rank")
    p2 = anchors.loc[anchors["evidence_tier"].str.startswith("P2")].sort_values("consensus_rank")
    p3 = anchors.loc[anchors["evidence_tier"].str.startswith("P3")].sort_values("consensus_rank")
    tier2_t = targets.loc[targets["evidence_tier"].str.startswith("T2")].sort_values("consensus_rank")
    tier3_t = targets.loc[targets["evidence_tier"].str.startswith("T3")].sort_values("consensus_rank")
    strict_t2 = targets.loc[targets["strict_all_seed_T2_sensitivity_gate"].fillna(False)].sort_values("consensus_rank")

    def rank_phrase(frame: pd.DataFrame, gene: str) -> str:
        hit = frame.loc[frame["gene_symbol"].eq(gene)]
        if hit.empty:
            return f"{gene}, not ranked"
        r = hit.iloc[0]
        return f"{gene}, rank {int(r.consensus_rank)} (leave-one-source-out {int(r.loso_rank_best)}–{int(r.loso_rank_worst)})"
    text = f"""# Integration Methods and Results (manuscript-ready)

## Methods

### Human-liver localization and virtual perturbation calibration

Normalized single-cell expression from five healthy donor livers (GSE115469; 8,444 cells) was summarized with the donor as the inferential unit. Cell-level principal-component embeddings were used descriptively. We grouped hepatocytes (n=3,501), endothelial cells (n=844), and macrophages (n=1,192) for virtual perturbation. For each compartment and each of three prespecified seeds, no more than 60 cells per donor were sampled without replacement before construction of three principal-component-regression networks. Single-gene and expressed canonical multi-target virtual knockouts were evaluated with scTenifoldKnk. The primary output was the percentile rank of the unsigned network-displacement distance; it was not interpreted as expression direction, protection, harm, or clinical causality.

To calibrate aggregate displacement of the 13 prespecified DILI anchors, each anchor was matched to 30 background genes in standardized log-mean-expression and logit-detection space within the same compartment. We generated 3,000 matched gene sets per target, compartment, and seed, calculated an empirical two-sided P value, and controlled the false-discovery rate within compartment. Pairwise Jaccard indices of the top 20 displaced genes quantified seed-to-seed gene-level stability.

### Source-independent candidate integration

We analyzed two candidate roles separately: exposure-proximal pharmacologic targets and downstream circulating DILI phenotype anchors. Evidence was grouped by dataset, rather than by the number of derived features. ChEMBL potency and document replication formed one pharmacology group; GSE115469 expression and virtual perturbation formed one liver-context group; the public DILI proteomics workbook formed one clinical group; six-dose HepG2 LINCS formed one observed-perturbation group; and STRING formed one contextual association group. Sci-Plex MCF7 membership at both doses was retained as an external replication flag because only the published top/bottom gene lists were available.

Within each candidate role, continuous features were converted to empirical percentile scores. The clinical DILI score assigned 45% to DO-versus-NDO discrimination, 20% each to DO-versus-HV detection and DO-versus-DF recovery, and 15% to within-DO ALT correlation; for proteins without follow-up measurements, available components were renormalized (45%, 25%, and 30%). The pharmacology score combined maximum pChEMBL (55%), median pChEMBL (25%), and a prespecified confidence level for canonical multi-document, noncanonical multi-document, or single-document evidence (20%). For each T-track target and liver compartment, liver context was the minimum of the within-compartment detection percentile and the within-compartment median virtual-knockout percentile rank. We then retained the compartment with the maximum of this maximin score. Thus, high expression in one compartment could not compensate for a virtual-knockout rank obtained from another compartment. HepG2 perturbation combined absolute dose-rank correlation (60%) and maximum absolute Level-5 z score (40%), multiplied by 1.0 for landmark, 0.8 for best-inferred, or 0.5 for other inferred transcripts. STRING context combined maximum path-product to the opposite candidate class (65%) and weighted degree (35%).

We used non-dominated Pareto fronts and the median of source-specific ranks rather than a single pooled probability. Robustness was tested by recomputing the consensus rank after omitting each source group. Phenotype and target tiers were deliberately separated. P1 required replicated DILI-versus-non-DILI discrimination and hepatocyte localization; P2 denoted general injury/recovery proteins; and P3 denoted limited or inconsistent phenotype evidence. The prespecified exploratory T2 gate required canonical or multi-document noncanonical pharmacology plus either (i) at least 10% detection, results from all three seeds, and a median aggregate-anchor virtual-knockout rank in the top quartile within that same compartment, or (ii) a nonconstant, directly measured HepG2 landmark transcript with an absolute six-dose Spearman rho of at least 0.70. The alternative HepG2 route could not be satisfied by an inferred L1000 transcript. T3 included all remaining screened target hypotheses. A stricter sensitivity required the within-compartment virtual-knockout rank to remain in the top quartile in each of the three seeds. T1 required nintedanib-specific human liver evidence or prospective functional validation and was not assigned.

## Results

### Liver localization separates exposure-proximal targets from circulating injury anchors

Canonical VEGFR-family targets were concentrated in endothelial cells, whereas most circulating DILI anchors were detected predominantly in hepatocytes. This compartmental separation argues against treating plasma injury proteins as biochemical nintedanib targets. Virtual knockout yielded compartment-dependent aggregate anchor ranks, but gene-level top-20 sets were only partially stable (median single-target Jaccard {jaccard.loc[jaccard['ko_class'].eq('drug_target'), 'jaccard_top20'].median():.2f}). The matched-null analysis identified {int((matched['n_seed_bh_q_lt_0_05'] > 0).sum())} target-compartment combinations with at least one BH-significant seed; no direction of effect was inferred.

### Split-role integration yields one exploratory T2-gate hypothesis and preserves phenotype-anchor status

FBP1 was the sole P1 phenotype anchor because it uniquely reproduced DILI-versus-non-DILI discrimination and was hepatocyte localized; notably, non-DILI onset values were higher, so this is etiologic discrimination rather than a DILI-specific increase. The P-track consensus positions most relevant to the perturbation analysis were {rank_phrase(anchors, 'GSTA1')}, {rank_phrase(anchors, 'FBP1')}, and {rank_phrase(anchors, 'OTC')}. P2 comprised {', '.join(p2['gene_symbol'])}; P3 comprised {', '.join(p3['gene_symbol'])}. These labels describe downstream human phenotypes, not biochemical nintedanib targets.

After enforcing the same-compartment gate, the exploratory T2 target-network shortlist comprised {', '.join(tier2_t['gene_symbol']) if len(tier2_t) else 'no targets'}. FLT1 and FGFR3 were downgraded to T3 because their expression and favorable virtual-knockout evidence arose from different compartments under the earlier permissive rule. The directly measured HepG2 alternative route added {', '.join(targets.loc[targets['hepg2_alternative_gate'].fillna(False) & targets['strong_pharmacology_gate'].fillna(False), 'gene_symbol']) if (targets['hepg2_alternative_gate'].fillna(False) & targets['strong_pharmacology_gate'].fillna(False)).any() else 'no target'}. Under the stricter requirement that the same-compartment virtual-knockout rank remain in the top quartile in every seed, the T2 set was {', '.join(strict_t2['gene_symbol']) if len(strict_t2) else 'empty'}. Consensus ranking and tiering were retained as separate outputs: a high median source rank could not bypass the noncompensatory biological gate. No T1 target was assigned, and no candidate is established as a mediator of nintedanib-associated liver injury.

## Interpretation guardrails

- The human DILI proteomics cohorts were not nintedanib-specific; their proteins are phenotype anchors.
- Healthy-liver single-cell expression supports cellular plausibility but does not measure DILI state.
- scTenifoldKnk is an unsigned, model-derived network displacement analysis and is not a causal or protective-effect estimator.
- HepG2 and MCF7 are transformed cell lines; cross-model transcript concordance was weak.
- ChEMBL multiple documents strengthen one pharmacology source; they are not independent omics validations.
- STRING paths are associations, often including text-mining evidence, and are not a demonstrated mechanism.
"""
    (OUT / "integration_methods_results.md").write_text(text)


def main() -> None:
    configure_style()
    data = load_inputs()
    lincs_metrics, lincs_long = derive_lincs_metrics(data["lincs"])
    lincs_metrics.to_csv(OUT / "hepg2_gene_dose_metrics_used.csv.gz", index=False, compression="gzip")
    candidate_union = sorted(set(DILI_ANCHORS) | set(data["vko_stable"].loc[data["vko_stable"]["ko_class"].eq("drug_target"), "ko_label"]))
    string_metrics, graph = string_graph_metrics(data, candidate_union)
    matched_seed, matched_stable = matched_null_analysis(data)
    jaccard = top20_jaccard(data)
    combined, domain_long, loso = integrate_candidates(data, lincs_metrics, string_metrics)
    figure_hashes = {
        "figure4": figure4(data, matched_stable, jaccard),
        "figure5": figure5(data, lincs_metrics, lincs_long, combined, loso, graph),
    }
    write_captions(combined, matched_stable, jaccard)
    write_methods_results(combined, matched_stable, jaccard)

    manifest = {
        "analysis": "transparent multi-evidence integration for nintedanib-associated DILI target prioritization",
        "analysis_script": str(Path(__file__).relative_to(ROOT)),
        "analysis_script_sha256": sha256(Path(__file__)),
        "validation_script": "fresh_analysis/scripts/validate_virtual_integration_outputs.py",
        "validation_script_sha256": sha256(ANALYSIS / "scripts" / "validate_virtual_integration_outputs.py"),
        "random_seed": SEED,
        "candidate_roles_kept_separate": True,
        "source_independence_rule": "derived features from the same dataset are grouped within one evidence domain",
        "liver_context_rule": "same-compartment maximin of expression rank and median virtual-KO rank; no cross-compartment splicing",
        "primary_T2_gate": "prespecified exploratory screen: strong pharmacology AND ((same compartment detection>=0.10, n_seeds=3, median vKO rank>=0.75) OR (directly measured nonconstant HepG2 landmark and abs(rho)>=0.70))",
        "strict_T2_sensitivity": "primary gate with same-compartment vKO rank>=0.75 in every seed",
        "phenotype_tiers": "P1/P2/P3 are distinct from target T1/T2/T3",
        "virtual_knockout_interpretation": "unsigned network displacement only; no causal, protective, harmful, or expression-direction claim",
        "pareto_rule": "iterative non-dominated sorting within candidate role",
        "leave_one_source_out": True,
        "tier1_assigned": False,
        "figure_panel_counts": {"Figure_4": 10, "Figure_5": 10},
        "figure_hashes": figure_hashes,
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                ANALYSIS / "proteomics" / "candidate_evidence_table.csv",
                ANALYSIS / "liver_singlecell" / "target_expression_by_compartment.csv",
                ANALYSIS / "pharmacology" / "chembl_target_summary.csv",
                ANALYSIS / "virtual_ko" / "virtual_knockout_stability_summary.csv",
                ANALYSIS / "virtual_ko" / "virtual_knockout_seed_summary.csv",
                ANALYSIS / "virtual_ko" / "virtual_knockout_gene_results.csv.gz",
                ANALYSIS / "lincs" / "output" / "nintedanib_hepg2_level5_all_genes_tidy.csv.gz",
                ANALYSIS / "inputs" / "string_candidate_network.tsv",
                ANALYSIS / "inputs" / "string_candidate_network_expanded.tsv",
            ]
        },
        "output_table_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                OUT / "candidate_consensus_ranking.csv",
                OUT / "candidate_domain_scores_long.csv",
                OUT / "candidate_leave_one_source_out_ranks.csv",
                OUT / "tier2_candidate_shortlist.csv",
                OUT / "phenotype_anchor_tiers.csv",
                OUT / "target_same_compartment_context_all_pairs.csv",
                OUT / "target_tier_gate_sensitivity.csv",
                OUT / "virtual_ko_expression_matched_null_by_seed.csv",
                OUT / "virtual_ko_expression_matched_null_stability.csv",
                OUT / "virtual_ko_top20_pairwise_jaccard.csv",
                OUT / "string_candidate_path_metrics.csv",
            ]
        },
        "text_artifact_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                FIG_OUT / "Figure_4_caption.md",
                FIG_OUT / "Figure_5_caption.md",
                OUT / "integration_methods_results.md",
                OUT / "input_audit.md",
            ]
        },
    }
    (OUT / "integration_analysis_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({
        "P1": combined.loc[combined["evidence_tier"].str.startswith("P1"), ["gene_symbol", "consensus_rank", "loso_rank_best", "loso_rank_worst"]].to_dict("records"),
        "P2": combined.loc[combined["evidence_tier"].str.startswith("P2"), ["gene_symbol", "consensus_rank", "loso_rank_best", "loso_rank_worst"]].to_dict("records"),
        "P3": combined.loc[combined["evidence_tier"].str.startswith("P3"), ["gene_symbol", "consensus_rank", "loso_rank_best", "loso_rank_worst"]].to_dict("records"),
        "T2": combined.loc[combined["evidence_tier"].str.startswith("T2"), ["gene_symbol", "consensus_rank", "loso_rank_best", "loso_rank_worst", "joint_context_compartment"]].to_dict("records"),
        "strict_all_seed_T2": combined.loc[combined["strict_all_seed_T2_sensitivity_gate"].fillna(False), ["gene_symbol"]].to_dict("records"),
        "outputs": manifest["figure_hashes"],
    }, indent=2))


if __name__ == "__main__":
    main()
