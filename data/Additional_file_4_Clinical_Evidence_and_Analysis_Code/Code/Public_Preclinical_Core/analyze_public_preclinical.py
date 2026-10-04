#!/usr/bin/env python3
"""Reanalyse public nintedanib preclinical perturbation datasets.

The script preserves the experimental unit encoded by each source.  GSE278200
and GSE308578 are analysed at the independent-animal level.  GSE120804 and
GSE120679 are pooled-tissue slice experiments and their slice-level statistics
are explicitly labelled exploratory.  PXD024058 is descriptive because the
public report is internally inconsistent about whether its three columns per
group are biological or technical replicates.

Outputs are deterministic and written under fresh_analysis/preclinical/results
and fresh_analysis/figures.
"""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import math
import os
import platform
import re
import sys
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
import statsmodels
from matplotlib.lines import Line2D
from scipy import linalg, optimize, special, stats
from sklearn.decomposition import PCA
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
PRE = ROOT / "preclinical"
RAW = PRE / "raw"
META = PRE / "metadata"
RES = PRE / "results"
FIG = ROOT / "figures"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260904
RNG = np.random.default_rng(SEED)

COLORS = {
    "control": "#697386",
    "disease": "#D97706",
    "nintedanib": "#0F766E",
    "blue": "#2563A5",
    "purple": "#7C3AED",
    "pink": "#C24172",
    "light": "#D6DEE8",
    "dark": "#17212B",
}

TARGET_GENES = ["FGFR1", "FGFR3", "FLT1", "OTC", "GSTA1", "FBP1", "CES1D"]
CONTEXT_GENES = ["CPS1", "ASS1", "ARG1", "COL1A1", "ACTA2", "POSTN", "FN1", "CTGF", "ABCB1B"]
DISPLAY_GENES = TARGET_GENES + CONTEXT_GENES

MODULES = {
    "ECM/fibrosis": ["COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "CTGF", "TIMP1", "LOX", "SERPINE1", "ACTA2", "TGFBI", "MMP2"],
    "Macrophage–lipid": ["CD68", "ADGRE1", "CSF1R", "LPL", "APOE", "TREM2", "LGALS3", "CTSD", "ABCA1", "FABP5", "SPP1", "LIPA", "PLIN2"],
    "Urea/arginine": ["ARG1", "OTC", "CPS1", "ASS1", "ASL", "OAT"],
    "Nintedanib targets": ["FGFR1", "FGFR3", "FLT1"],
    "Hepatic anchors": ["OTC", "GSTA1", "FBP1", "CES1A", "CES1D", "CES1F", "CPS1", "ASS1"],
    "Oxidative stress": ["NQO1", "HMOX1", "SOD1", "SOD2", "GPX1", "GPX3", "GSTA1", "CAT", "PRDX1", "PRDX2"],
    "Drug efflux": ["ABCB1A", "ABCB1B", "ABCG2", "ABCC1", "ABCC2"],
}


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def parse_gtf(path: Path) -> pd.DataFrame:
    records = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            a = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            gid = a.get("gene_id", "").split(".")[0]
            if not gid:
                continue
            records.append(
                {
                    "gene_id": gid,
                    "symbol": a.get("gene_name", gid),
                    "biotype": a.get("gene_biotype", a.get("gene_type", "unknown")),
                }
            )
    return pd.DataFrame(records).drop_duplicates("gene_id")


def parse_soft_samples(path: Path) -> pd.DataFrame:
    rows, cur = [], None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if cur:
                    rows.append(cur)
                cur = {}
            elif cur is not None:
                if line.startswith("!Sample_title = "):
                    cur["title"] = line.split(" = ", 1)[1]
                elif line.startswith("!Sample_geo_accession = "):
                    cur["geo_accession"] = line.split(" = ", 1)[1]
                elif line.startswith("!Sample_characteristics_ch1 = treatment: "):
                    cur["treatment"] = line.split("treatment: ", 1)[1]
                elif line.startswith("!Sample_characteristics_ch1 = replicate group: "):
                    cur["replicate_group"] = line.split("replicate group: ", 1)[1]
                elif line.startswith("!Sample_description = Library name: "):
                    cur["library"] = line.split("Library name: ", 1)[1]
                elif line.startswith("!Sample_relation = BioSample: "):
                    cur["biosample"] = line.rsplit("/", 1)[-1]
                elif line.startswith("!Sample_relation = SRA: "):
                    cur["sra_relation"] = line.split("term=", 1)[-1]
    if cur:
        rows.append(cur)
    return pd.DataFrame(rows)


def aggregate_by_symbol(counts: pd.DataFrame, ann: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Map Ensembl IDs and sum counts for duplicated symbols."""
    x = counts.copy()
    x.index = x.index.astype(str).str.split(".").str[0]
    x = x.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    a = ann.set_index("gene_id").reindex(x.index)
    keep = a["symbol"].notna().to_numpy()
    x = x.loc[keep]
    a = a.loc[keep]
    symbols = a["symbol"].astype(str).str.upper()
    x["__symbol"] = symbols.to_numpy()
    agg = x.groupby("__symbol", sort=False).sum(numeric_only=True)
    info = (
        a.assign(gene_id=a.index.astype(str), symbol_upper=symbols.to_numpy())
        .groupby("symbol_upper", sort=False)
        .agg(gene_id=("gene_id", lambda z: ";".join(map(str, z))), biotype=("biotype", lambda z: ";".join(sorted(set(map(str, z))))))
    )
    return agg, info


def calc_tmm_factors(counts: pd.DataFrame, logratio_trim=0.30, sum_trim=0.05) -> pd.Series:
    y = counts.to_numpy(float)
    lib = y.sum(axis=0)
    # Reference library: upper quartile nearest the mean upper quartile.
    with np.errstate(divide="ignore", invalid="ignore"):
        uq = np.nanquantile(np.where(y > 0, y / lib, np.nan), 0.75, axis=0)
    ref_idx = int(np.nanargmin(np.abs(uq - np.nanmean(uq))))
    factors = np.ones(y.shape[1])
    for j in range(y.shape[1]):
        if j == ref_idx:
            continue
        obs, ref = y[:, j], y[:, ref_idx]
        keep = (obs > 0) & (ref > 0) & np.isfinite(obs) & np.isfinite(ref)
        obs, ref = obs[keep], ref[keep]
        if len(obs) < 100:
            continue
        m = np.log2((obs / lib[j]) / (ref / lib[ref_idx]))
        a = 0.5 * np.log2((obs / lib[j]) * (ref / lib[ref_idx]))
        w = 1.0 / ((lib[j] - obs) / (lib[j] * obs) + (lib[ref_idx] - ref) / (lib[ref_idx] * ref))
        ml, mh = np.quantile(m, [logratio_trim, 1 - logratio_trim])
        al, ah = np.quantile(a, [sum_trim, 1 - sum_trim])
        use = (m >= ml) & (m <= mh) & (a >= al) & (a <= ah) & np.isfinite(w) & (w > 0)
        if use.sum() >= 50:
            factors[j] = 2 ** np.average(m[use], weights=w[use])
    factors = factors / np.exp(np.mean(np.log(factors)))
    return pd.Series(factors, index=counts.columns, name="TMM_factor")


def normalize_counts(counts: pd.DataFrame, min_group_n: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    lib = counts.sum(axis=0)
    raw_cpm = counts.div(lib, axis=1) * 1e6
    keep = (raw_cpm >= 1).sum(axis=1) >= min_group_n
    filt = counts.loc[keep]
    tmm = calc_tmm_factors(filt)
    eff = lib * tmm
    logcpm = np.log2((filt + 0.5).div(eff + 1.0, axis=1) * 1e6)
    norm = pd.DataFrame({"sample": counts.columns, "library_size": lib.values, "TMM_factor": tmm.reindex(counts.columns).values, "effective_library_size": eff.values})
    filt_table = pd.DataFrame({"gene": counts.index, "retained": counts.index.isin(filt.index), "max_raw_CPM": raw_cpm.max(axis=1).values, "n_samples_CPM_ge_1": (raw_cpm >= 1).sum(axis=1).values})
    return logcpm, norm, filt_table


def estimate_prior_df(log_s2_residual: np.ndarray, residual_df: float) -> tuple[float, float]:
    """Robust moment estimator for the scaled-F prior used by moderated t tests."""
    z = np.asarray(log_s2_residual, float)
    z = z[np.isfinite(z)]
    if len(z) < 100:
        return 10.0, np.var(z)
    lo, hi = np.quantile(z, [0.05, 0.95])
    zw = np.clip(z, lo, hi)
    observed = float(np.var(zw, ddof=1))
    sampling = float(special.polygamma(1, residual_df / 2.0))
    excess = max(observed - sampling, 1e-8)
    if excess <= 1e-7:
        return 1e6, observed
    f = lambda d: special.polygamma(1, d / 2.0) - excess
    try:
        df0 = float(optimize.brentq(f, 0.1, 1e6))
    except ValueError:
        df0 = 1e6 if f(1e6) > 0 else 0.1
    return df0, observed


def moderated_group_model(logcpm: pd.DataFrame, groups: pd.Series, contrasts: dict[str, tuple[str, str]]) -> tuple[dict[str, pd.DataFrame], dict]:
    """limma-trend-style robust empirical-Bayes model on TMM logCPM."""
    groups = groups.reindex(logcpm.columns)
    levels = list(dict.fromkeys(groups.tolist()))
    X = pd.get_dummies(groups, dtype=float)[levels].to_numpy(float)
    Y = logcpm.T.to_numpy(float)  # samples x genes
    xtxi = linalg.pinv(X.T @ X)
    beta = xtxi @ X.T @ Y
    fitted = X @ beta
    resid = Y - fitted
    df = X.shape[0] - np.linalg.matrix_rank(X)
    s2 = np.sum(resid**2, axis=0) / df
    s2 = np.maximum(s2, np.finfo(float).tiny)
    ave = np.mean(Y, axis=0)
    trend = lowess(np.log(s2), ave, frac=0.40, it=2, return_sorted=False)
    log_resid = np.log(s2) - trend
    df0, winsor_var = estimate_prior_df(log_resid, df)
    elogf = special.digamma(df / 2.0) - np.log(df / 2.0) - special.digamma(df0 / 2.0) + np.log(df0 / 2.0)
    prior_var = np.exp(trend - elogf)
    post_var = (df0 * prior_var + df * s2) / (df0 + df)
    results = {}
    for name, (numerator, denominator) in contrasts.items():
        c = np.zeros(len(levels))
        c[levels.index(numerator)] = 1
        c[levels.index(denominator)] = -1
        effect = c @ beta
        unscaled = float(np.sqrt(c @ xtxi @ c))
        se = np.sqrt(post_var) * unscaled
        tval = effect / se
        p = 2 * stats.t.sf(np.abs(tval), df + df0)
        crit = stats.t.ppf(0.975, df + df0)
        tab = pd.DataFrame(
            {
                "gene": logcpm.index,
                "contrast": name,
                "log2FC": effect,
                "SE_moderated": se,
                "CI95_low": effect - crit * se,
                "CI95_high": effect + crit * se,
                "t_moderated": tval,
                "p_value": p,
                "q_value_BH_genome": bh(p),
                "average_logCPM": ave,
            }
        ).sort_values("p_value")
        results[name] = tab
    diagnostics = {"n_samples": int(X.shape[0]), "n_genes": int(Y.shape[1]), "design_rank": int(np.linalg.matrix_rank(X)), "residual_df": int(df), "prior_df": float(df0), "winsorized_log_variance": winsor_var, "group_levels": levels, "method": "TMM logCPM plus robust limma-trend-style empirical-Bayes moderated linear model"}
    return results, diagnostics


def welch_group_model(logexpr: pd.DataFrame, groups: pd.Series, contrasts: dict[str, tuple[str, str]], scope: str) -> dict[str, pd.DataFrame]:
    groups = groups.reindex(logexpr.columns)
    out = {}
    for name, (num, den) in contrasts.items():
        a = logexpr.loc[:, groups == num].to_numpy(float)
        b = logexpr.loc[:, groups == den].to_numpy(float)
        effect = np.nanmean(a, axis=1) - np.nanmean(b, axis=1)
        tt = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
        # Hedges-like probability-of-superiority direction via all pairwise comparisons.
        ps = np.array([np.mean(x[:, None] > y[None, :]) + 0.5 * np.mean(x[:, None] == y[None, :]) for x, y in zip(a, b)])
        tab = pd.DataFrame({"gene": logexpr.index, "contrast": name, "log2FC": effect, "t_Welch": tt.statistic, "p_value_slice_exploratory": tt.pvalue, "q_value_BH_slice_exploratory": bh(tt.pvalue), "probability_superiority": ps, "inferential_scope": scope})
        out[name] = tab.sort_values("p_value_slice_exploratory")
    return out


def pca_table(logexpr: pd.DataFrame, meta: pd.DataFrame, n_features=2000) -> tuple[pd.DataFrame, np.ndarray]:
    v = logexpr.var(axis=1).nlargest(min(n_features, logexpr.shape[0])).index
    model = PCA(n_components=2, random_state=SEED)
    z = model.fit_transform(logexpr.loc[v].T)
    out = meta.set_index("sample").reindex(logexpr.columns).reset_index()
    out["PC1"] = z[:, 0]
    out["PC2"] = z[:, 1]
    return out, model.explained_variance_ratio_


def exact_permutation_p(a: np.ndarray, b: np.ndarray) -> float:
    vals = np.r_[a, b]
    n_a = len(a)
    obs = abs(np.mean(a) - np.mean(b))
    if len(vals) <= 12 and math.comb(len(vals), n_a) <= 200000:
        diffs = []
        all_idx = np.arange(len(vals))
        for idx in itertools.combinations(all_idx, n_a):
            use = np.zeros(len(vals), bool)
            use[list(idx)] = True
            diffs.append(abs(np.mean(vals[use]) - np.mean(vals[~use])))
        return float(np.mean(np.asarray(diffs) >= obs - 1e-12))
    # Reproducible Monte Carlo permutation for larger experiments.
    ge = 0
    n_perm = 50000
    for _ in range(n_perm):
        perm = RNG.permutation(vals)
        ge += abs(np.mean(perm[:n_a]) - np.mean(perm[n_a:])) >= obs - 1e-12
    return (ge + 1) / (n_perm + 1)


def module_scores(logexpr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Standardise each gene across samples before averaging to prevent abundant genes dominating.
    z = logexpr.sub(logexpr.mean(axis=1), axis=0).div(logexpr.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    scores, coverage = {}, []
    for module, genes in MODULES.items():
        present = [g for g in genes if g in z.index]
        if present:
            scores[module] = z.loc[present].mean(axis=0)
        coverage.append({"module": module, "n_present": len(present), "n_prespecified": len(genes), "genes_present": ";".join(present)})
    return pd.DataFrame(scores).T, pd.DataFrame(coverage)


def module_contrasts(scores: pd.DataFrame, groups: pd.Series, contrasts: dict[str, tuple[str, str]], unit_scope: str, use_permutation: bool) -> pd.DataFrame:
    rows = []
    groups = groups.reindex(scores.columns)
    for cname, (num, den) in contrasts.items():
        for mod in scores.index:
            a = scores.loc[mod, groups == num].dropna().to_numpy(float)
            b = scores.loc[mod, groups == den].dropna().to_numpy(float)
            effect = float(np.mean(a) - np.mean(b))
            se = math.sqrt(np.var(a, ddof=1) / len(a) + np.var(b, ddof=1) / len(b))
            dof_num = (np.var(a, ddof=1) / len(a) + np.var(b, ddof=1) / len(b)) ** 2
            dof_den = (np.var(a, ddof=1) / len(a)) ** 2 / (len(a) - 1) + (np.var(b, ddof=1) / len(b)) ** 2 / (len(b) - 1)
            dof = dof_num / dof_den if dof_den > 0 else np.inf
            p_welch = 2 * stats.t.sf(abs(effect / se), dof) if se > 0 else np.nan
            p_perm = exact_permutation_p(a, b) if use_permutation else np.nan
            rows.append({"module": mod, "contrast": cname, "effect_z_score": effect, "SE_Welch": se, "CI95_low": effect - stats.t.ppf(.975, dof) * se, "CI95_high": effect + stats.t.ppf(.975, dof) * se, "p_Welch": p_welch, "p_permutation": p_perm, "n_numerator": len(a), "n_denominator": len(b), "experimental_unit_scope": unit_scope})
    out = pd.DataFrame(rows)
    for cname in out["contrast"].unique():
        idx = out["contrast"] == cname
        pcol = "p_permutation" if use_permutation else "p_Welch"
        out.loc[idx, "q_BH_within_contrast_modules"] = bh(out.loc[idx, pcol])
    return out


def candidate_results(de: dict[str, pd.DataFrame], dataset: str, p_scope: str) -> pd.DataFrame:
    rows = []
    for cname, tab in de.items():
        x = tab.set_index("gene")
        for gene in DISPLAY_GENES:
            if gene in x.index:
                r = x.loc[gene]
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                d = {"dataset": dataset, "gene": gene, "contrast": cname, "role": "target/anchor" if gene in TARGET_GENES else "context", "detected_after_filtering": True}
                d.update(r.to_dict())
                d["p_value_scope"] = p_scope
                rows.append(d)
            else:
                rows.append({"dataset": dataset, "gene": gene, "contrast": cname, "role": "target/anchor" if gene in TARGET_GENES else "context", "detected_after_filtering": False, "p_value_scope": p_scope})
    return pd.DataFrame(rows)


def reversal_summary(disease: pd.DataFrame, treatment: pd.DataFrame, q_col: str, threshold=0.5) -> tuple[pd.DataFrame, dict]:
    d = disease.set_index("gene")[["log2FC", q_col]].rename(columns={"log2FC": "disease_log2FC", q_col: "disease_q"})
    t = treatment.set_index("gene")[["log2FC", q_col]].rename(columns={"log2FC": "treatment_log2FC", q_col: "treatment_q"})
    x = d.join(t, how="inner").dropna()
    x["directional_reversal"] = np.sign(x["disease_log2FC"]) == -np.sign(x["treatment_log2FC"])
    x["reversal_index"] = -x["treatment_log2FC"] / x["disease_log2FC"].replace(0, np.nan)
    x["distance_ratio_to_control"] = np.abs(x["disease_log2FC"] + x["treatment_log2FC"]) / np.abs(x["disease_log2FC"]).replace(0, np.nan)
    selected = x[(x["disease_q"] < 0.05) & (x["disease_log2FC"].abs() >= threshold)].copy()
    directional = selected[selected["directional_reversal"]].copy()
    rho, p = stats.spearmanr(selected["disease_log2FC"], selected["treatment_log2FC"]) if len(selected) > 2 else (np.nan, np.nan)
    summ = {"n_disease_responsive": int(len(selected)), "fraction_directionally_reversed": float(selected["directional_reversal"].mean()) if len(selected) else np.nan, "fraction_closer_to_control": float((selected["distance_ratio_to_control"] < 1).mean()) if len(selected) else np.nan, "median_reversal_index_among_directional": float(directional["reversal_index"].median()) if len(directional) else np.nan, "fraction_directional_with_reversal_index_0.5_to_1.5": float(directional["reversal_index"].between(0.5, 1.5).mean()) if len(directional) else np.nan, "median_absolute_disease_log2FC": float(directional["disease_log2FC"].abs().median()) if len(directional) else np.nan, "median_absolute_treatment_log2FC": float(directional["treatment_log2FC"].abs().median()) if len(directional) else np.nan, "spearman_disease_vs_treatment_rho": float(rho), "spearman_nominal_p_gene_as_unit_descriptive": float(p), "selection": f"disease q<0.05 and |log2FC|>={threshold}", "caution": "Gene-wise geometry is descriptive because genes are correlated; the gene-wise correlation p value is not used for biological inference."}
    return x.reset_index(), summ


def read_gse278200(rat_ann):
    raw = pd.read_csv(RAW / "GSE278200_NINT_raw_counts.txt.gz", sep="\t")
    raw = raw.loc[:, ~raw.columns.str.startswith("Unnamed")]
    counts = raw.set_index("ID")
    counts, info = aggregate_by_symbol(counts, rat_ann)
    samples = list(counts.columns)
    meta = pd.DataFrame({"sample": samples, "group": [re.match(r"(SAL|BLM|NINT)", s).group(1) for s in samples]})
    meta["experimental_unit"] = "independent animal"
    meta["species"] = "Rattus norvegicus"
    meta["tissue"] = "right lung homogenate"
    meta["dose_time"] = meta["group"].map({"SAL": "saline + vehicle; day 28", "BLM": "BLM 1 U/kg IT days 0/4 + vehicle; day 28", "NINT": "BLM + nintedanib 100 mg/kg PO daily days 7–28"})
    logcpm, norm, filt = normalize_counts(counts, min_group_n=4)
    groups = meta.set_index("sample")["group"]
    contrasts = {"BLM_vs_SAL": ("BLM", "SAL"), "NINT_vs_BLM": ("NINT", "BLM"), "NINT_vs_SAL": ("NINT", "SAL")}
    de, diag = moderated_group_model(logcpm, groups, contrasts)
    pca, pve = pca_table(logcpm, meta)
    scores, coverage = module_scores(logcpm)
    mod = module_contrasts(scores, groups, contrasts, "independent animal", use_permutation=True)
    cand = candidate_results(de, "GSE278200", "genome-wide BH; moderated animal-level model")
    rev, rev_summ = reversal_summary(de["BLM_vs_SAL"], de["NINT_vs_BLM"], "q_value_BH_genome")
    return locals()


def read_gse308578(mouse_ann):
    raw = pd.read_csv(RAW / "GSE308578_Gene_counts.csv.gz")
    raw = raw.rename(columns={raw.columns[0]: "gene_id"})
    counts = raw.set_index("gene_id")
    counts, info = aggregate_by_symbol(counts, mouse_ann)
    sm = parse_soft_samples(META / "GSE308578_family.soft.gz")
    sm["library"] = sm["library"].astype(str)
    mapping = sm.set_index("library")
    absent = [s for s in counts.columns if s not in mapping.index]
    if absent:
        raise ValueError(f"GSE308578 count libraries absent from SOFT: {absent}")
    meta = mapping.reindex(counts.columns).reset_index()
    meta = meta.rename(columns={meta.columns[0]: "sample"})
    meta["group"] = meta["treatment"].map({"CTRL Vehicle": "CTRL", "BLEO-IPF Nintedanib Vehicle": "BLM_vehicle", "BLEO-IPF Nintedanib": "NINT"})
    meta["experimental_unit"] = "independent animal"
    meta["species"] = "Mus musculus"
    meta["tissue"] = "inferior lung lobe"
    meta["dose_time"] = meta["group"].map({"CTRL": "saline + vehicle", "BLM_vehicle": "BLM 2 mg/kg IT + vehicle", "NINT": "BLM + nintedanib 60 mg/kg PO BID days 8–28"})
    logcpm, norm, filt = normalize_counts(counts, min_group_n=10)
    groups = meta.set_index("sample")["group"]
    contrasts = {"BLM_vehicle_vs_CTRL": ("BLM_vehicle", "CTRL"), "NINT_vs_BLM_vehicle": ("NINT", "BLM_vehicle"), "NINT_vs_CTRL": ("NINT", "CTRL")}
    de, diag = moderated_group_model(logcpm, groups, contrasts)
    pca, pve = pca_table(logcpm, meta)
    scores, coverage = module_scores(logcpm)
    mod = module_contrasts(scores, groups, contrasts, "independent animal", use_permutation=True)
    cand = candidate_results(de, "GSE308578", "genome-wide BH; moderated animal-level model")
    rev, rev_summ = reversal_summary(de["BLM_vehicle_vs_CTRL"], de["NINT_vs_BLM_vehicle"], "q_value_BH_genome")
    # Leave-one-animal-out effect ranges for prespecified genes and modules.
    loo_rows = []
    for sample in logcpm.columns:
        g = groups[sample]
        if g == "CTRL":
            continue
        reduced = logcpm.drop(columns=sample)
        rg = groups.drop(index=sample)
        for gene in DISPLAY_GENES:
            if gene in reduced.index:
                eff = reduced.loc[gene, rg == "NINT"].mean() - reduced.loc[gene, rg == "BLM_vehicle"].mean()
                loo_rows.append({"feature_type": "gene", "feature": gene, "omitted_sample": sample, "omitted_group": g, "log2FC_NINT_vs_vehicle": eff})
        rscores, _ = module_scores(reduced)
        for module in rscores.index:
            eff = rscores.loc[module, rg == "NINT"].mean() - rscores.loc[module, rg == "BLM_vehicle"].mean()
            loo_rows.append({"feature_type": "module", "feature": module, "omitted_sample": sample, "omitted_group": g, "log2FC_NINT_vs_vehicle": eff})
    loo = pd.DataFrame(loo_rows)
    loo_summary = loo.groupby(["feature_type", "feature"], as_index=False).agg(loo_min=("log2FC_NINT_vs_vehicle", "min"), loo_max=("log2FC_NINT_vs_vehicle", "max"), loo_median=("log2FC_NINT_vs_vehicle", "median"), sign_concordance=("log2FC_NINT_vs_vehicle", lambda z: max(np.mean(z > 0), np.mean(z < 0))), n_leave_one_out=("log2FC_NINT_vs_vehicle", "size"))
    return locals()


def read_gse120804(rat_ann89):
    raw = pd.read_csv(RAW / "GSE120804_counts.txt.gz", sep="\t")
    counts = raw.set_index("GeneID")
    counts, info = aggregate_by_symbol(counts, rat_ann89)
    ann = pd.read_csv(RAW / "GSE120804_geo_sample_annotation_edit.csv.gz")
    ann["sample"] = ann["raw.file"].str.replace(r"_Run1.*$", "", regex=True)
    amap = ann.set_index("sample")
    absent = [s for s in counts.columns if s not in amap.index]
    if absent:
        raise ValueError(f"GSE120804 columns absent from annotation: {absent}")
    meta = amap.reindex(counts.columns).reset_index()
    meta = meta.rename(columns={meta.columns[0]: "sample"})
    meta["group"] = meta["Replicate.group"]
    meta["experimental_unit"] = "slice observation from pooled tissue (3 rats pooled per ex vivo study)"
    meta["species"] = "Rattus norvegicus"
    meta["tissue"] = "liver PCLS"
    meta["dose_time"] = meta["group"].map({"Sham": "sham PCLS, 48 h", "BDL": "BDL PCLS, vehicle 48 h", "BDL_Nint": "BDL PCLS + nintedanib 1.0 µM, 48 h"})
    logcpm, norm, filt = normalize_counts(counts, min_group_n=9)
    groups = meta.set_index("sample")["group"]
    contrasts = {"BDL_vs_Sham": ("BDL", "Sham"), "BDL_Nint_vs_BDL": ("BDL_Nint", "BDL"), "BDL_Nint_vs_Sham": ("BDL_Nint", "Sham")}
    de = welch_group_model(logcpm, groups, contrasts, "exploratory slice-level; tissues from three rats were pooled before randomisation")
    scores, coverage = module_scores(logcpm)
    mod = module_contrasts(scores, groups, contrasts, "pooled-tissue slice observation; not an independent animal", use_permutation=False)
    cand = candidate_results(de, "GSE120804", "slice-level exploratory BH; no animal-population inference")
    rev, rev_summ = reversal_summary(de["BDL_vs_Sham"].rename(columns={"q_value_BH_slice_exploratory": "q"}), de["BDL_Nint_vs_BDL"].rename(columns={"q_value_BH_slice_exploratory": "q"}), "q")
    return locals()


def read_gse120679(rat_ann89):
    logcpm = pd.read_csv(RAW / "GSE120679_logCPM_expressed_protein_coding.txt.gz", sep="\t").set_index("GeneID")
    # The GEO processed matrix columns and SOFT samples are source-provided in the same order.
    sm = parse_soft_samples(META / "GSE120679_family.soft.gz")
    if len(sm) != logcpm.shape[1]:
        raise ValueError("GSE120679 SOFT/matrix sample count mismatch")
    meta = pd.DataFrame({"sample": logcpm.columns, "geo_accession": sm["geo_accession"], "title": sm["title"], "group": sm["replicate_group"]})
    meta["experimental_unit"] = "slice observation from pooled tissue (3 rats pooled per ex vivo study)"
    meta["species"] = "Rattus norvegicus"
    meta["tissue"] = "lung PCLS"
    meta["dose_time"] = meta["group"].map({"None_72_None": "vehicle PCLS, 72 h", "TGFb_72_None": "TGFβ1 PCLS, 72 h", "TGFb_72_Nint": "TGFβ1 + nintedanib 0.3 µM, 72 h"})
    ann = rat_ann89.set_index("gene_id")
    symbols = ann.reindex(logcpm.index.astype(str).str.split(".").str[0])["symbol"].astype(str).str.upper()
    logcpm["__symbol"] = symbols.to_numpy()
    logcpm = logcpm[logcpm["__symbol"] != "NAN"].groupby("__symbol").mean(numeric_only=True)
    groups = meta.set_index("sample")["group"]
    contrasts = {"TGFb_vs_vehicle": ("TGFb_72_None", "None_72_None"), "TGFb_Nint_vs_TGFb": ("TGFb_72_Nint", "TGFb_72_None"), "TGFb_Nint_vs_vehicle": ("TGFb_72_Nint", "None_72_None")}
    de = welch_group_model(logcpm, groups, contrasts, "exploratory slice-level; tissues from three rats were pooled before randomisation")
    scores, coverage = module_scores(logcpm)
    mod = module_contrasts(scores, groups, contrasts, "pooled-tissue slice observation; not an independent animal", use_permutation=False)
    cand = candidate_results(de, "GSE120679", "slice-level exploratory BH; no animal-population inference")
    return locals()


def read_pxd024058():
    raw = pd.read_csv(RAW / "PXD024058_all.protein.xls", sep="\t")
    qcols = [c for c in raw.columns if re.match(r"Control_A|Model_B|Mod\+NIB_C", c)]
    raw["gene"] = raw["Gene Name"].fillna("").astype(str).str.split(";").str[0].str.upper()
    x = raw[raw["gene"] != ""].copy()
    vals = x[qcols].apply(pd.to_numeric, errors="coerce")
    vals = np.log2(vals.where(vals > 0))
    vals["gene"] = x["gene"].to_numpy()
    expr = vals.groupby("gene").median(numeric_only=True)
    rows = []
    for gene, r in expr.iterrows():
        ctr = r[[c for c in qcols if c.startswith("Control")]].dropna().to_numpy(float)
        mod = r[[c for c in qcols if c.startswith("Model")]].dropna().to_numpy(float)
        nib = r[[c for c in qcols if c.startswith("Mod+NIB")]].dropna().to_numpy(float)
        rows.append({"gene": gene, "BLM_vs_control_log2FC_descriptive": np.mean(mod) - np.mean(ctr) if len(mod) and len(ctr) else np.nan, "NIB_vs_BLM_log2FC_descriptive": np.mean(nib) - np.mean(mod) if len(nib) and len(mod) else np.nan, "n_control_columns": len(ctr), "n_model_columns": len(mod), "n_nintedanib_columns": len(nib), "direction_consistency_NIB_minus_model": np.mean(nib[:, None] - mod[None, :] > 0) if len(nib) and len(mod) else np.nan, "inferential_scope": "descriptive only; article conflicts on biological versus technical replicate provenance"})
    effects = pd.DataFrame(rows)
    cand = effects[effects["gene"].isin(DISPLAY_GENES)].copy()
    return locals()


def write_dataset_outputs(tag: str, obj: dict, independent: bool):
    obj["meta"].to_csv(RES / f"{tag}_sample_metadata.csv", index=False)
    if "norm" in obj:
        obj["norm"].to_csv(RES / f"{tag}_normalization_factors.csv", index=False)
    if "filt" in obj:
        obj["filt"].to_csv(RES / f"{tag}_expression_filter.csv.gz", index=False, compression="gzip")
    if "pca" in obj:
        obj["pca"].to_csv(RES / f"{tag}_pca_coordinates.csv", index=False)
    de_all = pd.concat(obj["de"].values(), ignore_index=True)
    de_all.to_csv(RES / f"{tag}_differential_expression.csv.gz", index=False, compression="gzip")
    obj["cand"].to_csv(RES / f"{tag}_candidate_results.csv", index=False)
    obj["scores"].T.rename_axis("sample").reset_index().to_csv(RES / f"{tag}_module_scores.csv", index=False)
    obj["coverage"].to_csv(RES / f"{tag}_module_coverage.csv", index=False)
    obj["mod"].to_csv(RES / f"{tag}_module_contrasts.csv", index=False)
    if "rev" in obj:
        obj["rev"].to_csv(RES / f"{tag}_reversal_geometry.csv.gz", index=False, compression="gzip")


def cross_model_table(objs: dict) -> pd.DataFrame:
    specs = {
        "GSE278200 rat lung RNA (animal n=4/group)": (objs["g278"]["de"]["NINT_vs_BLM"], "log2FC"),
        "GSE308578 mouse lung RNA (animal 17 vs 16)": (objs["g308"]["de"]["NINT_vs_BLM_vehicle"], "log2FC"),
        "GSE120804 rat liver PCLS RNA (pooled slices 10 vs 10)": (objs["g804"]["de"]["BDL_Nint_vs_BDL"], "log2FC"),
        "GSE120679 rat lung PCLS RNA (pooled slices 6 vs 6)": (objs["g679"]["de"]["TGFb_Nint_vs_TGFb"], "log2FC"),
        "PXD024058 mouse lung protein (3 columns/group; provenance ambiguous)": (objs["pxd"]["effects"].rename(columns={"NIB_vs_BLM_log2FC_descriptive": "log2FC"}), "log2FC"),
    }
    rows = []
    features = DISPLAY_GENES + list(MODULES)
    for model, (tab, ecol) in specs.items():
        gt = tab.set_index("gene")[ecol]
        for feature in features:
            if feature in MODULES:
                present = [g for g in MODULES[feature] if g in gt.index]
                val = float(gt.reindex(present).mean()) if present else np.nan
                n = len(present)
                kind = "module mean gene effect"
            else:
                val = float(gt.get(feature, np.nan))
                n = int(feature in gt.index)
                kind = "gene"
            rows.append({"model": model, "feature": feature, "feature_type": kind, "n_features_observed": n, "nintedanib_vs_disease_log2_effect": val})
    return pd.DataFrame(rows)


def setup_plotting():
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.4, "axes.titlesize": 8.2, "axes.labelsize": 7.4, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.7, "pdf.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white", "axes.facecolor": "white"})


def panel_label(ax, label):
    ax.text(-0.13, 1.19, label, transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="top", ha="left", color=COLORS["dark"])


def scatter_pca(ax, pca, pve, group_order, labels, colors, title, subtitle):
    for g in group_order:
        d = pca[pca["group"] == g]
        ax.scatter(d.PC1, d.PC2, s=34, color=colors[g], edgecolor="white", linewidth=.6, label=f"{labels[g]} (n={len(d)})", alpha=.95)
    ax.axhline(0, color="#E5E7EB", lw=.6, zorder=0); ax.axvline(0, color="#E5E7EB", lw=.6, zorder=0)
    ax.set_xlabel(f"PC1 ({pve[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({pve[1]*100:.1f}%)")
    ax.set_title(title, loc="left", fontweight="bold", pad=12)
    ax.text(0, 1.01, subtitle, transform=ax.transAxes, fontsize=7.0, color="#56616F", va="bottom")
    ax.legend(frameon=False, fontsize=6.8, loc="best", handletextpad=.3)


def reversal_scatter(ax, rev, summary, title, genes_to_label):
    d = rev[(rev["disease_q"] < .05) & (rev["disease_log2FC"].abs() >= .5)].copy()
    ax.scatter(d["disease_log2FC"], d["treatment_log2FC"], s=9, alpha=.32, color=COLORS["blue"], linewidth=0)
    lim = np.nanmax(np.abs(np.r_[d["disease_log2FC"], d["treatment_log2FC"]])) if len(d) else 2
    lim = min(max(lim, 1), 7)
    ax.plot([-lim, lim], [lim, -lim], ls="--", lw=.8, color="#77808C", label="complete reversal")
    ax.axhline(0, lw=.55, color="#C8CED6"); ax.axvline(0, lw=.55, color="#C8CED6")
    for g in genes_to_label:
        z = d[d["gene"] == g]
        if len(z):
            x, y = z.iloc[0][["disease_log2FC", "treatment_log2FC"]]
            ax.scatter([x], [y], s=30, color=COLORS["pink"], edgecolor="white", lw=.5, zorder=3)
            ax.text(x, y, " " + g, fontsize=6.8, va="center")
    ax.set_xlabel("Disease effect, log2FC"); ax.set_ylabel("Nintedanib effect, log2FC")
    ax.set_title(title, loc="left", fontweight="bold", pad=12)
    ax.text(0, 1.01, f"n={summary['n_disease_responsive']}; reversed={summary['fraction_directionally_reversed']*100:.1f}%; median RI={summary['median_reversal_index_among_directional']:.2f}", transform=ax.transAxes, fontsize=7.0, color="#56616F", va="bottom")


def volcano(ax, tab, title, subtitle):
    q = tab["q_value_BH_genome"].clip(lower=1e-300)
    sig = (q < .05) & (tab["log2FC"].abs() >= .5)
    ax.scatter(tab.loc[~sig, "log2FC"], -np.log10(q[~sig]), s=5, color="#BBC4CE", alpha=.45, linewidth=0)
    ax.scatter(tab.loc[sig, "log2FC"], -np.log10(q[sig]), s=8, color=COLORS["nintedanib"], alpha=.65, linewidth=0)
    lab = tab.assign(score=-np.log10(q)).sort_values("score", ascending=False)
    force = ["COL1A1", "POSTN", "CES1D"]
    picked = list(lab.loc[sig].head(2).gene) + [g for g in force if g in set(tab.gene)]
    shown = set()
    for k, g in enumerate(picked):
        if g in shown: continue
        r = lab[lab.gene == g]
        if len(r):
            rr = r.iloc[0]; shown.add(g)
            if rr["score"] > 1:
                ax.annotate(g, (rr.log2FC, rr.score), xytext=((5 if rr.log2FC >= 0 else -5), 3 + (k % 3) * 5), textcoords="offset points", fontsize=6.8, ha="left" if rr.log2FC >= 0 else "right", va="bottom", arrowprops={"arrowstyle":"-","lw":.35,"color":"#7A8490"})
    ax.axvline(-.5, color="#9099A6", ls="--", lw=.6); ax.axvline(.5, color="#9099A6", ls="--", lw=.6)
    ax.axhline(-np.log10(.05), color="#9099A6", ls="--", lw=.6)
    ax.set_ylim(bottom=-.05, top=max(3.2, ax.get_ylim()[1] * 1.12))
    ax.set_xlabel("Nintedanib vs disease, log2FC"); ax.set_ylabel("−log10 genome-wide q")
    ax.set_title(title, loc="left", fontweight="bold", pad=12)
    ax.text(0, 1.01, subtitle + f"; q<0.05 & |FC|≥0.5: {sig.sum():,}", transform=ax.transAxes, fontsize=7.0, color="#56616F", va="bottom")


def effect_forest(ax, cand, contrast, genes, title, subtitle):
    d = cand[(cand["contrast"] == contrast) & cand["gene"].isin(genes)].copy()
    d["gene"] = pd.Categorical(d["gene"], categories=genes[::-1], ordered=True)
    d = d.sort_values("gene")
    y = np.arange(len(d))
    detected=d["log2FC"].notna().to_numpy()
    if "CI95_low" in d:
        xerr = np.vstack([(d.log2FC - d.CI95_low).to_numpy()[detected], (d.CI95_high - d.log2FC).to_numpy()[detected]])
        ax.errorbar(d.log2FC.to_numpy()[detected], y[detected], xerr=xerr, fmt="o", ms=4.2, color=COLORS["nintedanib"], ecolor="#8CA7A4", capsize=2, lw=.8)
    else:
        ax.scatter(d.log2FC.to_numpy()[detected], y[detected], s=24, color=COLORS["nintedanib"])
    ax.axvline(0, color="#8A949F", lw=.7)
    ax.set_yticks(y, d.gene.astype(str)); ax.set_xlabel("Nintedanib vs disease, log2FC")
    for yi in y[~detected]:
        ax.text(.98, yi, "below expression filter", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.8, color="#8A949F")
    ax.set_title(title, loc="left", fontweight="bold", pad=12)
    ax.text(0, 1.01, subtitle, transform=ax.transAxes, fontsize=7.0, color="#56616F", va="bottom")


def module_forest(ax, mod, contrast, title, subtitle, color=None):
    d = mod[mod.contrast == contrast].copy()
    order = list(MODULES)[::-1]
    d["module"] = pd.Categorical(d.module, categories=order, ordered=True)
    d = d.sort_values("module")
    y = np.arange(len(d))
    c = color or COLORS["nintedanib"]
    ax.errorbar(d.effect_z_score, y, xerr=np.vstack([d.effect_z_score-d.CI95_low, d.CI95_high-d.effect_z_score]), fmt="o", ms=4.2, color=c, ecolor="#9BA6B1", capsize=2, lw=.8)
    for yi, (_, r) in zip(y, d.iterrows()):
        if r.q_BH_within_contrast_modules < .05:
            ax.text(r.effect_z_score, yi, "  *", va="center", fontsize=8, fontweight="bold")
    ax.axvline(0, color="#8A949F", lw=.7)
    ax.set_yticks(y, d.module.astype(str)); ax.set_xlabel("Module effect (mean gene z-score)")
    ax.set_title(title, loc="left", fontweight="bold", pad=12)
    ax.text(0, 1.01, subtitle, transform=ax.transAxes, fontsize=7.0, color="#56616F", va="bottom")


def g308_loo_panel(ax, obj):
    genes = ["ABCB1B", "FGFR1", "FGFR3", "FLT1", "OTC", "GSTA1", "FBP1", "COL1A1", "POSTN"]
    d = obj["cand"][(obj["cand"].contrast == "NINT_vs_BLM_vehicle") & obj["cand"].gene.isin(genes)].set_index("gene")
    loo = obj["loo_summary"]
    loo = loo[loo.feature_type == "gene"].set_index("feature")
    rows=[]
    for g in genes:
        if g in d.index:
            rows.append({"gene":g,"effect":d.loc[g,"log2FC"],"lo":loo.loc[g,"loo_min"] if g in loo.index else np.nan,"hi":loo.loc[g,"loo_max"] if g in loo.index else np.nan,"q":d.loc[g,"q_value_BH_genome"]})
    z=pd.DataFrame(rows).iloc[::-1].reset_index(drop=True); y=np.arange(len(z))
    valid=z.effect.notna().to_numpy()
    ax.hlines(y[valid],z.lo.to_numpy()[valid],z.hi.to_numpy()[valid],color="#9AA5B1",lw=2,label="leave-one-animal-out range")
    ax.scatter(z.effect.to_numpy()[valid],y[valid],s=30,color=[COLORS["pink"] if g=="ABCB1B" else COLORS["nintedanib"] for g in z.gene[valid]],edgecolor="white",lw=.5,zorder=3)
    ax.axvline(0,color="#8A949F",lw=.7); ax.set_yticks(y,z.gene); ax.set_xlabel("Nintedanib vs vehicle, log2FC")
    for yi in y[~valid]:
        ax.text(.98, yi, "below expression filter", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.8, color="#8A949F")
    ax.set_title("Model-boundary targets and efflux",loc="left",fontweight="bold",pad=12)
    ax.text(0,1.01,"GSE308578; points=all animals, bars=leave-one-out range",transform=ax.transAxes,fontsize=7.0,color="#56616F",va="bottom")


def liver_pcls_panel(ax, obj):
    genes=["FGFR1","FGFR3","FLT1","OTC","GSTA1","FBP1","CES1D","CPS1","ASS1","ARG1","COL1A1","ACTA2"]
    d=obj["cand"][(obj["cand"].contrast=="BDL_Nint_vs_BDL") & obj["cand"].gene.isin(genes)].copy()
    d["gene"]=pd.Categorical(d.gene,categories=genes[::-1],ordered=True); d=d.sort_values("gene"); y=np.arange(len(d))
    q=d["q_value_BH_slice_exploratory"].to_numpy(float)
    ax.scatter(d.log2FC,y,s=34,c=[COLORS["purple"] if v<.05 else "#AAB3BE" for v in q],edgecolor="white",lw=.5)
    ax.axvline(0,color="#8A949F",lw=.7); ax.set_yticks(y,d.gene.astype(str)); ax.set_xlabel("BDL+Nintedanib vs BDL, log2FC")
    ax.set_title("Hepatic-context perturbation",loc="left",fontweight="bold",pad=12)
    ax.text(0,1.01,"GSE120804 liver PCLS; pooled-slice observations (10 vs 10), exploratory",transform=ax.transAxes,fontsize=7.0,color="#56616F",va="bottom")
    ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=COLORS["purple"],markeredgecolor="white",label="slice-level BH q<0.05"),Line2D([0],[0],marker="o",color="none",markerfacecolor="#AAB3BE",markeredgecolor="white",label="q≥0.05")],frameon=False,fontsize=6.8,loc="lower right")


def heatmap_panel(ax, cross):
    features=["FGFR1","FGFR3","FLT1","OTC","GSTA1","FBP1","CES1D","ECM/fibrosis","Macrophage–lipid","Urea/arginine","Drug efflux"]
    labels={
        "GSE278200 rat lung RNA (animal n=4/group)":"Rat lung RNA\nanimal n=4/group",
        "GSE308578 mouse lung RNA (animal 17 vs 16)":"Mouse lung RNA\nanimal 17 vs 16",
        "GSE120804 rat liver PCLS RNA (pooled slices 10 vs 10)":"Rat liver PCLS RNA\npooled slices 10 vs 10",
        "GSE120679 rat lung PCLS RNA (pooled slices 6 vs 6)":"Rat lung PCLS RNA\npooled slices 6 vs 6",
        "PXD024058 mouse lung protein (3 columns/group; provenance ambiguous)":"Mouse lung protein\n3 columns; provenance ?",
    }
    mat=cross.pivot(index="model",columns="feature",values="nintedanib_vs_disease_log2_effect").reindex(index=list(labels),columns=features)
    lim=max(.5,float(np.nanpercentile(np.abs(mat.to_numpy()),95))); lim=min(lim,2.5)
    im=ax.imshow(mat.clip(-lim,lim),aspect="auto",cmap="PuOr_r",vmin=-lim,vmax=lim,interpolation="nearest")
    ax.set_xticks(np.arange(len(features)),features,rotation=55,ha="right"); ax.set_yticks(np.arange(len(labels)),[labels[x] for x in labels])
    ax.set_title("Cross-model nintedanib response",loc="left",fontweight="bold",pad=12)
    ax.text(0,1.01,"Cell: gene log2FC or module mean; grey = not observed",transform=ax.transAxes,fontsize=7.0,color="#56616F",va="bottom")
    im.cmap.set_bad("#D8DDE3")
    cb=plt.colorbar(im,ax=ax,fraction=.025,pad=.02); cb.set_label("Nintedanib vs disease\nlog2 effect",fontsize=6.8); cb.ax.tick_params(labelsize=6.8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v=mat.iloc[i,j]
            if np.isfinite(v): ax.text(j,i,f"{v:.1f}",ha="center",va="center",fontsize=6.6,color="white" if abs(v)>.55*lim else "#1F2937")


def make_main_figure(objs, cross):
    setup_plotting()
    fig, axes=plt.subplots(5,2,figsize=(8.0,11.4),constrained_layout=False)
    plt.subplots_adjust(left=.11,right=.96,top=.925,bottom=.055,hspace=.76,wspace=.48)
    fig.suptitle("Figure 6 | Public models define pulmonary-response and tissue-context boundaries",x=.11,y=.986,ha="left",fontsize=10.2,fontweight="bold",color=COLORS["dark"])
    fig.text(.11,.966,"Independent-animal inference is separated from pooled-slice and provenance-ambiguous descriptive evidence",ha="left",fontsize=7.2,color="#56616F")
    ax=axes[0,0]; scatter_pca(ax,objs["g278"]["pca"],objs["g278"]["pve"],["SAL","BLM","NINT"],{"SAL":"SAL","BLM":"BLM","NINT":"BLM + NINT"},{"SAL":COLORS["control"],"BLM":COLORS["disease"],"NINT":COLORS["nintedanib"]},"Independent-rat lung transcriptomes","GSE278200; each point is one animal"); panel_label(ax,"A")
    ax=axes[0,1]; reversal_scatter(ax,objs["g278"]["rev"],objs["g278"]["rev_summ"],"Partial reversal of BLM effects",["COL1A1","POSTN","FBP1","FGFR3","FLT1"]); panel_label(ax,"B")
    ax=axes[1,0]; volcano(ax,objs["g278"]["de"]["NINT_vs_BLM"],"Rat lung nintedanib response","Animal-level moderated model; n=4 vs 4"); panel_label(ax,"C")
    ax=axes[1,1]; effect_forest(ax,objs["g278"]["cand"],"NINT_vs_BLM",["FGFR1","FGFR3","FLT1","OTC","GSTA1","FBP1","CES1D","COL1A1","ACTA2","POSTN"],"Prespecified targets and anchors","GSE278200; 95% moderated-t CIs, genome-wide q in table"); panel_label(ax,"D")
    ax=axes[2,0]; module_forest(ax,objs["g278"]["mod"],"NINT_vs_BLM","Rat lung program effects","GSE278200; exact label-permutation q across 7 prespecified modules"); panel_label(ax,"E")
    ax=axes[2,1]; scatter_pca(ax,objs["g308"]["pca"],objs["g308"]["pve"],["CTRL","BLM_vehicle","NINT"],{"CTRL":"CTRL","BLM_vehicle":"BLM + vehicle","NINT":"BLM + NINT"},{"CTRL":COLORS["control"],"BLM_vehicle":COLORS["disease"],"NINT":COLORS["nintedanib"]},"Independent-mouse model boundary","GSE308578; each point is one animal"); panel_label(ax,"F")
    ax=axes[3,0]; reversal_scatter(ax,objs["g308"]["rev"],objs["g308"]["rev_summ"],"Lower-magnitude transcriptomic reversal",["ABCB1B","COL1A1","POSTN","FGFR1","FBP1"]); panel_label(ax,"G")
    ax=axes[3,1]; g308_loo_panel(ax,objs["g308"]); panel_label(ax,"H")
    ax=axes[4,0]; liver_pcls_panel(ax,objs["g804"]); panel_label(ax,"I")
    ax=axes[4,1]; heatmap_panel(ax,cross); panel_label(ax,"J")
    for ext,kwargs in [("pdf",{}),("svg",{}),("png",{"dpi":320}),("tiff",{"dpi":400,"pil_kwargs":{"compression":"tiff_lzw"}})]:
        fig.savefig(FIG/f"Figure_6_Public_Preclinical_Constraints.{ext}",bbox_inches="tight",facecolor="white",**kwargs)
    plt.close(fig)


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    rat103=parse_gtf(RAW/"Rattus_norvegicus.Rnor_6.0.103.gtf.gz")
    rat89=parse_gtf(RAW/"Rattus_norvegicus.Rnor_6.0.89.gtf.gz")
    mouse=parse_gtf(RAW/"Mus_musculus.GRCm39.113.gtf.gz")
    objs={}
    objs["g278"]=read_gse278200(rat103)
    objs["g308"]=read_gse308578(mouse)
    objs["g804"]=read_gse120804(rat89)
    objs["g679"]=read_gse120679(rat89)
    objs["pxd"]=read_pxd024058()

    for tag,key in [("GSE278200","g278"),("GSE308578","g308"),("GSE120804","g804"),("GSE120679","g679")]:
        write_dataset_outputs(tag,objs[key],tag in {"GSE278200","GSE308578"})
    objs["g308"]["loo"].to_csv(RES/"GSE308578_leave_one_animal_out_full.csv.gz",index=False,compression="gzip")
    objs["g308"]["loo_summary"].to_csv(RES/"GSE308578_leave_one_animal_out_summary.csv",index=False)
    objs["pxd"]["effects"].to_csv(RES/"PXD024058_protein_effects_descriptive.csv.gz",index=False,compression="gzip")
    objs["pxd"]["cand"].to_csv(RES/"PXD024058_candidate_effects_descriptive.csv",index=False)
    cross=cross_model_table(objs)
    cross.to_csv(RES/"cross_model_candidate_module_effects.csv",index=False)
    make_main_figure(objs,cross)

    summary={
        "GSE278200": {"diagnostics":objs["g278"]["diag"],"reversal":objs["g278"]["rev_summ"],"group_n":objs["g278"]["meta"].groupby("group").size().to_dict()},
        "GSE308578": {"diagnostics":objs["g308"]["diag"],"reversal":objs["g308"]["rev_summ"],"group_n":objs["g308"]["meta"].groupby("group").size().to_dict(),"metadata_resolution":"GEO SOFT and ARRIVE Supplementary Figure S1 agree on CTRL 10, vehicle 16, nintedanib 17; main Figure 6 caption in source article reverses vehicle/nintedanib n. No sample was relabelled."},
        "GSE120804": {"group_n":objs["g804"]["meta"].groupby("group").size().to_dict(),"scope":"pooled-tissue slice observations; exploratory only"},
        "GSE120679": {"group_n":objs["g679"]["meta"].groupby("group").size().to_dict(),"scope":"pooled-tissue slice observations; exploratory only"},
        "PXD024058": {"n_proteins":int(len(objs["pxd"]["effects"])),"scope":"descriptive only because replicate provenance is internally inconsistent"},
    }
    (RES/"analysis_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    # Machine-readable frozen-output validation. Scientific boundary warnings are
    # deliberate and do not constitute computational validation failures.
    checks=[]
    def check(name, condition, detail):
        checks.append({"check":name,"passed":bool(condition),"detail":detail})
    check("GSE278200 independent-animal group sizes", objs["g278"]["meta"].groupby("group").size().to_dict()=={"BLM":4,"NINT":4,"SAL":4}, objs["g278"]["meta"].groupby("group").size().to_dict())
    check("GSE308578 GEO/ARRIVE group sizes", objs["g308"]["meta"].groupby("group").size().to_dict()=={"BLM_vehicle":16,"CTRL":10,"NINT":17}, objs["g308"]["meta"].groupby("group").size().to_dict())
    check("Unique independent-animal identifiers", objs["g278"]["meta"]["sample"].is_unique and objs["g308"]["meta"]["sample"].is_unique, "GSE278200 and GSE308578")
    check("Genome-wide q values bounded", all(x["q_value_BH_genome"].dropna().between(0,1).all() for o in [objs["g278"],objs["g308"]] for x in o["de"].values()), "all independent-animal contrasts")
    check("Independent-animal module permutations present", objs["g278"]["mod"]["p_permutation"].notna().all() and objs["g308"]["mod"]["p_permutation"].notna().all(), "exact for GSE278200; 50,000 Monte Carlo permutations for GSE308578")
    check("Seven prespecified modules per contrast", len(objs["g278"]["mod"])==21 and len(objs["g308"]["mod"])==21, {"GSE278200":len(objs["g278"]["mod"]),"GSE308578":len(objs["g308"]["mod"])})
    check("Cross-model matrix complete shape", len(cross)==5*(len(DISPLAY_GENES)+len(MODULES)), {"rows":len(cross),"expected":5*(len(DISPLAY_GENES)+len(MODULES))})
    check("Figure 6 four formats", all((FIG/f"Figure_6_Public_Preclinical_Constraints.{e}").exists() and (FIG/f"Figure_6_Public_Preclinical_Constraints.{e}").stat().st_size>10000 for e in ["pdf","svg","png","tiff"]), "PDF/SVG/PNG/TIFF")
    validation={"validation_passed":all(c["passed"] for c in checks),"n_checks":len(checks),"n_failed":sum(not c["passed"] for c in checks),"checks":checks,"scientific_boundary_warnings":["GSE120804 and GSE120679 slice observations derive from pooled rat tissues and cannot support animal-population inference.","PXD024058 replicate provenance is internally inconsistent; its effects are descriptive only.","GSE308578 source main-Figure-6 caption reverses nintedanib/vehicle n values; GEO SOFT and ARRIVE Supplementary Figure S1 agree on 17/16, which was used without relabelling.","Reversal geometry treats genes as display units; gene correlation precludes interpreting its nominal correlation p value as biological inference."]}
    validation_text=json.dumps(validation,indent=2,ensure_ascii=False)+"\n"
    (RES/"validation_report.json").write_text(validation_text)
    (RES/"preclinical_quality_audit.json").write_text(validation_text)
    # The workspace may perform asynchronous lossless PNG optimization after
    # savefig returns. Allow that operation to settle before freezing hashes;
    # otherwise the manifest can capture the pre-optimization byte stream.
    time.sleep(5.0)
    manifest=[]
    manifest_path=RES/"analysis_manifest.json"
    for p in sorted([*RAW.glob("*"),*META.glob("*"),*RES.glob("*"),*(PRE/"scripts").glob("*"),*PRE.glob("*.md"),*PRE.glob("*.csv"),*FIG.glob("Figure_6_Public_Preclinical_Constraints.*")]):
        if p.is_file() and p != manifest_path:
            manifest.append({"path":str(p.relative_to(ROOT.parent)),"bytes":p.stat().st_size,"sha256":sha256(p)})
    meta_manifest={"generated_utc":"2026-09-04","seed":SEED,"python":sys.version,"platform":platform.platform(),"packages":{"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"statsmodels":statsmodels.__version__,"scikit-learn":sklearn.__version__,"matplotlib":mpl.__version__},"files":manifest}
    manifest_path.write_text(json.dumps(meta_manifest,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
