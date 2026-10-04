#!/usr/bin/env python3
"""Targeted longitudinal re-analysis of GSE299128 whole-blood RNA-seq.

The public trial has repeated samples from seven participants and no adjudicated
nintedanib DILI cases. It is therefore used only as a tissue-specificity and
systemic-exposure sensitivity analysis, never as validation of liver toxicity.
"""

from __future__ import annotations

import itertools
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "fresh_data"
OUT = ROOT / "fresh_analysis" / "blood_longitudinal"
OUT.mkdir(parents=True, exist_ok=True)

GENES = [
    "ACO1", "ASS1", "FAH", "CPS1", "ALDOB", "HPD", "OTC", "DMGDH", "GSTA1",
    "FBP1", "PCK2", "CES1", "LECT2", "KDR", "FLT1", "FLT4", "FGFR1", "FGFR2",
    "FGFR3", "PDGFRA", "PDGFRB", "ABCB11", "ABCC2", "UGT1A1", "NFE2L2",
    "HMOX1", "TXNIP", "TXNRD1", "GPX4", "SLC7A11",
]


def parse_metadata(columns: list[str]) -> pd.DataFrame:
    rows = []
    for column in columns:
        match = re.search(r"Sample_(\d+-\d+)_v(\d+)", column)
        if match is None:
            raise ValueError(f"Unexpected sample name: {column}")
        rows.append(
            {
                "sample": column,
                "participant": match.group(1),
                "visit": int(match.group(2)),
            }
        )
    return pd.DataFrame(rows)


def exact_signflip_p(slopes: np.ndarray) -> float:
    """Two-sided exact randomization P from the complete sign-flip space.

    Because all 2**n assignments are enumerated, no Monte Carlo add-one
    correction is applied; the observed assignment is already part of the
    reference distribution.
    """
    observed = abs(float(np.mean(slopes)))
    null = []
    for signs in itertools.product([-1.0, 1.0], repeat=len(slopes)):
        null.append(abs(float(np.mean(slopes * np.asarray(signs)))))
    return sum(value >= observed - 1e-15 for value in null) / len(null)


def main() -> None:
    path = DATA / "GSE299128_count_matrix.csv.gz"
    counts = pd.read_csv(path, index_col=0)
    counts.index = counts.index.astype(str)
    counts = counts.groupby(level=0).sum()
    metadata = parse_metadata(counts.columns.tolist())
    metadata.to_csv(OUT / "sample_metadata.csv", index=False)

    library_size = counts.sum(axis=0)
    log_cpm = np.log2(counts.div(library_size, axis=1) * 1_000_000 + 0.5)
    present = [gene for gene in GENES if gene in log_cpm.index]
    log_cpm.loc[present].T.reset_index(names="sample").to_csv(
        OUT / "target_log2_cpm.csv", index=False
    )

    boot_rng = np.random.default_rng(20260904)
    summaries = []
    slopes_long = []
    for gene in present:
        participant_slopes = []
        baseline_last = []
        for participant, pm in metadata.groupby("participant", sort=True):
            pm = pm.sort_values("visit")
            y = log_cpm.loc[gene, pm["sample"]].to_numpy(dtype=float)
            x = pm["visit"].to_numpy(dtype=float)
            if len(np.unique(x)) < 2:
                continue
            slope = float(np.polyfit(x, y, 1)[0])
            participant_slopes.append(slope)
            baseline_last.append(float(y[-1] - y[0]))
            slopes_long.append(
                {"gene_symbol": gene, "participant": participant, "slope_per_visit": slope,
                 "last_minus_baseline_log2_cpm": float(y[-1] - y[0]), "n_visits": len(y)}
            )
        slopes = np.asarray(participant_slopes)
        boot = np.asarray([
            np.mean(boot_rng.choice(slopes, size=len(slopes), replace=True)) for _ in range(5000)
        ])
        summaries.append(
            {
                "gene_symbol": gene,
                "n_participants": len(slopes),
                "mean_slope_log2_cpm_per_visit": float(np.mean(slopes)),
                "slope_ci_low": float(np.quantile(boot, 0.025)),
                "slope_ci_high": float(np.quantile(boot, 0.975)),
                "exact_signflip_p": exact_signflip_p(slopes),
                "median_last_minus_baseline_log2_cpm": float(np.median(baseline_last)),
            }
        )
    summaries = pd.DataFrame(summaries)
    summaries["q_bh"] = multipletests(summaries["exact_signflip_p"], method="fdr_bh")[1]
    summaries.to_csv(OUT / "target_longitudinal_summary.csv", index=False)
    pd.DataFrame(slopes_long).to_csv(OUT / "participant_slopes.csv", index=False)

    manifest = {
        "accession": "GSE299128",
        "retrieval_date": "2026-09-04",
        "n_samples": int(counts.shape[1]),
        "n_participants": int(metadata["participant"].nunique()),
        "n_genes_tested": int(len(summaries)),
        "primary_test": (
            "two-sided exact participant-level sign-flip test of the mean "
            "within-participant slope, enumerating all 2^7 assignments"
        ),
        "multiplicity": "Benjamini-Hochberg across prespecified genes",
        "interpretation_boundary": (
            "Whole blood from a small single-arm nintedanib study without adjudicated DILI; "
            "results assess systemic longitudinal detectability, not liver-specific toxicity."
        ),
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
