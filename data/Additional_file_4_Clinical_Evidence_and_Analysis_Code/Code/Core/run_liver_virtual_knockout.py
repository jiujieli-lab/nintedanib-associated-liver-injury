#!/usr/bin/env python3
"""Run donor-balanced scTenifoldKnk virtual knockouts from scratch.

The analysis is deliberately framed as network prioritisation. scTenifoldKnk
returns unsigned network displacement after an in-silico knockout; it does not
estimate a protective direction, a molecular fold-change, or a causal clinical
effect. Multiple random seeds and donor-balanced cell sampling quantify the
stability of the ranking.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from scTenifold.core._base import scTenifoldKnk


ROOT = Path(__file__).resolve().parents[2]
SC = ROOT / "fresh_analysis" / "liver_singlecell"
OUT = ROOT / "fresh_analysis" / "virtual_ko"
OUT.mkdir(parents=True, exist_ok=True)
PARTS = OUT / "parts"
PARTS.mkdir(parents=True, exist_ok=True)

SEEDS = [20260904, 20260911, 20260918]
COMPARTMENTS = ["Hepatocyte", "Endothelial", "Macrophage"]
DILI_ANCHORS = [
    "ACO1", "ASS1", "FAH", "CPS1", "ALDOB", "HPD", "OTC", "DMGDH",
    "GSTA1", "FBP1", "PCK2", "CES1", "LECT2",
]
CANONICAL_TARGETS = ["KDR", "FLT1", "FLT4", "FGFR1", "FGFR2", "FGFR3", "PDGFRA", "PDGFRB"]
MECHANISM_GENES = [
    "ABCB11", "ABCC2", "SLC10A1", "UGT1A1", "CES1", "NFE2L2", "KEAP1",
    "HMOX1", "NQO1", "GCLC", "GCLM", "TXNIP", "TXNRD1", "GPX4", "SLC7A11",
    "SOD2", "CAT", "BAX", "BCL2", "CASP3", "CASP8", "KRT8", "KRT18",
]


def candidate_genes(compartment: str, matrix: pd.DataFrame) -> tuple[list[str], list[str]]:
    expression = pd.read_csv(SC / "target_expression_by_compartment.csv")
    expression = expression.loc[expression["compartment"].eq(compartment)].set_index("gene_symbol")
    pharm = pd.read_csv(ROOT / "fresh_analysis" / "pharmacology" / "chembl_target_summary.csv")
    replicated = pharm.loc[pharm["replicated_potent"].eq(True), "gene_symbol"].dropna().tolist()

    drug_targets = []
    for gene in list(dict.fromkeys(CANONICAL_TARGETS + replicated)):
        if gene not in matrix.index or gene not in expression.index:
            continue
        detection = float(expression.at[gene, "mean_detection_fraction"])
        if detection >= (0.01 if gene in CANONICAL_TARGETS else 0.03):
            drug_targets.append(gene)

    anchors = []
    for gene in DILI_ANCHORS:
        if gene not in matrix.index or gene not in expression.index:
            continue
        if float(expression.at[gene, "mean_detection_fraction"]) >= 0.03:
            anchors.append(gene)
    return drug_targets, anchors


def balanced_cells(meta: pd.DataFrame, available: list[str], compartment: str, seed: int) -> list[str]:
    sub = meta.loc[available]
    sub = sub.loc[sub["compartment"].eq(compartment)]
    selected = []
    for offset, (donor, dm) in enumerate(sub.groupby("Sample", sort=True)):
        n = min(60, len(dm))
        selected.extend(dm.sample(n=n, replace=False, random_state=seed + offset)["CellName"].tolist())
    return selected


def annotate_result(
    result: pd.DataFrame,
    compartment: str,
    seed: int,
    ko_label: str,
    ko_class: str,
    elapsed: float,
) -> pd.DataFrame:
    annotated = result.copy()
    annotated.columns = [c.lower().replace(" ", "_").replace("-", "_") for c in annotated.columns]
    annotated["compartment"] = compartment
    annotated["seed"] = seed
    annotated["ko_label"] = ko_label
    annotated["ko_class"] = ko_class
    annotated["ko_elapsed_seconds"] = elapsed
    annotated["distance_percentile"] = rankdata(annotated["distance"], method="average") / len(annotated)
    annotated["is_dili_anchor"] = annotated["gene"].isin(DILI_ANCHORS)
    annotated["is_mechanism_gene"] = annotated["gene"].isin(MECHANISM_GENES)
    annotated["is_canonical_drug_target"] = annotated["gene"].isin(CANONICAL_TARGETS)
    annotated["is_knocked_gene"] = annotated["gene"].isin(ko_label.split("+"))
    return annotated


def run() -> None:
    meta = pd.read_csv(SC / "cell_metadata.csv").set_index("CellName", drop=False)
    manifest_path = OUT / "virtual_knockout_run_manifest.csv"
    if manifest_path.exists():
        run_rows = pd.read_csv(manifest_path).to_dict("records")
    else:
        run_rows = []
    completed = {
        (str(row["compartment"]), int(row["seed"]), str(row["ko_label"]))
        for row in run_rows
        if (PARTS / f'{row["compartment"]}__{int(row["seed"])}__{row["ko_label"]}.csv.gz').exists()
    }

    for compartment in COMPARTMENTS:
        path = SC / f"{compartment.lower()}_network_matrix.csv.gz"
        full = pd.read_csv(path, index_col=0)
        drug_targets, anchors = candidate_genes(compartment, full)
        initial_candidates = list(dict.fromkeys(drug_targets + anchors))
        for seed in SEEDS:
            cells = balanced_cells(meta, full.columns.tolist(), compartment, seed)
            sample = full.loc[:, cells]
            variance = sample.var(axis=1)
            forced = [
                gene
                for gene in list(dict.fromkeys(initial_candidates + MECHANISM_GENES))
                if gene in sample.index and variance.get(gene, 0) > 0
            ]
            hvg = variance.loc[variance.gt(0)].nlargest(420).index.tolist()
            genes = list(dict.fromkeys(forced + hvg))[:460]
            sample = sample.loc[genes]

            sc = scTenifoldKnk(
                sample,
                qc_kws={
                    "min_lib_size": 0,
                    "remove_outlier_cells": False,
                    "min_percent": 0.01,
                    "max_mito_ratio": 1,
                    "min_exp_avg": 0,
                    "min_exp_sum": 0,
                    "plot": False,
                },
                nc_kws={
                    "n_nets": 3,
                    "n_samp_cells": min(250, len(cells)),
                    "n_comp": 3,
                    "scale_scores": True,
                    "symmetric": False,
                    "q": 0.95,
                    "random_state": seed,
                    "backend": "joblib-threading",
                    "n_jobs": 3,
                    "replace": True,
                },
                td_kws={
                    "method": "parafac",
                    "n_decimal": 2,
                    "K": 3,
                    "tol": 1e-6,
                    "max_iter": 500,
                    "random_state": seed,
                },
                ma_kws={"d": 20, "tol": 1e-8},
            )
            network_start = time.perf_counter()
            sc.run_step("qc")
            sc.run_step("nc")
            sc.run_step("td")
            network_seconds = time.perf_counter() - network_start

            available = set(sc.tensor_dict["WT"].index)
            this_drug = [g for g in drug_targets if g in available]
            this_anchor = [g for g in anchors if g in available]
            specs: list[tuple[list[str], str]] = [([g], "drug_target") for g in this_drug]
            specs += [([g], "dili_anchor") for g in this_anchor]
            canonical_group = [g for g in CANONICAL_TARGETS if g in available]
            if len(canonical_group) >= 2:
                specs.append((canonical_group, "canonical_poly_target"))

            for genes_ko, ko_class in specs:
                ko_label = "+".join(genes_ko)
                part_path = PARTS / f"{compartment}__{seed}__{ko_label}.csv.gz"
                if (compartment, seed, ko_label) in completed:
                    continue
                ko_start = time.perf_counter()
                sc.run_step("ko", ko_genes=genes_ko)
                sc.run_step("ma")
                sc.run_step("dr")
                elapsed = time.perf_counter() - ko_start
                annotated = annotate_result(
                    sc.d_regulation, compartment, seed, ko_label, ko_class, elapsed
                )
                temp_part = part_path.with_suffix(".tmp.gz")
                annotated.to_csv(temp_part, index=False, compression="gzip")
                temp_part.replace(part_path)
                run_rows.append(
                    {
                        "compartment": compartment,
                        "seed": seed,
                        "ko_label": ko_label,
                        "ko_class": ko_class,
                        "n_cells": len(cells),
                        "n_genes_input": len(genes),
                        "n_genes_after_qc": len(available),
                        "n_donors": int(meta.loc[cells, "Sample"].nunique()),
                        "network_seconds": network_seconds,
                        "ko_seconds": elapsed,
                    }
                )
                pd.DataFrame(run_rows).drop_duplicates(
                    ["compartment", "seed", "ko_label"], keep="last"
                ).to_csv(manifest_path, index=False)
                completed.add((compartment, seed, ko_label))

    part_files = sorted(PARTS.glob("*.csv.gz"))
    if not part_files:
        raise RuntimeError("No completed virtual-knockout part files were produced")
    results = pd.concat([pd.read_csv(path) for path in part_files], ignore_index=True)
    final_path = OUT / "virtual_knockout_gene_results.csv.gz"
    temp_final = OUT / "virtual_knockout_gene_results.tmp.gz"
    results.to_csv(temp_final, index=False, compression="gzip")
    temp_final.replace(final_path)

    eval_frame = results.loc[~results["is_knocked_gene"]].copy()
    summaries = []
    for keys, group in eval_frame.groupby(["compartment", "seed", "ko_label", "ko_class"]):
        dili = group.loc[group["is_dili_anchor"]]
        mechanism = group.loc[group["is_mechanism_gene"]]
        top20 = group.nlargest(20, "distance")
        summaries.append(
            {
                "compartment": keys[0],
                "seed": keys[1],
                "ko_label": keys[2],
                "ko_class": keys[3],
                "n_q_lt_0_05": int(group["adjusted_p_value"].lt(0.05).sum()),
                "dili_anchor_mean_percentile": float(dili["distance_percentile"].mean()) if len(dili) else np.nan,
                "dili_anchor_top20_count": int(top20["is_dili_anchor"].sum()),
                "mechanism_mean_percentile": float(mechanism["distance_percentile"].mean()) if len(mechanism) else np.nan,
                "mechanism_top20_count": int(top20["is_mechanism_gene"].sum()),
                "median_distance": float(group["distance"].median()),
            }
        )
    summaries = pd.DataFrame(summaries)
    summaries.to_csv(OUT / "virtual_knockout_seed_summary.csv", index=False)

    stable = (
        summaries.groupby(["compartment", "ko_label", "ko_class"], as_index=False)
        .agg(
            n_seeds=("seed", "nunique"),
            median_dili_anchor_percentile=("dili_anchor_mean_percentile", "median"),
            min_dili_anchor_percentile=("dili_anchor_mean_percentile", "min"),
            max_dili_anchor_percentile=("dili_anchor_mean_percentile", "max"),
            median_mechanism_percentile=("mechanism_mean_percentile", "median"),
            median_n_q_lt_0_05=("n_q_lt_0_05", "median"),
        )
    )
    stable.to_csv(OUT / "virtual_knockout_stability_summary.csv", index=False)

    downstream = (
        eval_frame.sort_values(
            ["compartment", "seed", "ko_label", "distance"], ascending=[True, True, True, False]
        )
        .groupby(["compartment", "seed", "ko_label"], as_index=False)
        .head(20)
        .groupby(["compartment", "ko_label", "gene"], as_index=False)
        .agg(top20_seed_frequency=("seed", "nunique"), median_distance=("distance", "median"))
        .sort_values(["compartment", "ko_label", "top20_seed_frequency", "median_distance"], ascending=[True, True, False, False])
    )
    downstream.to_csv(OUT / "stable_top_downstream_genes.csv", index=False)

    provenance = {
        "software": "scTenifoldpy 0.4.0 (import module scTenifold)",
        "seeds": SEEDS,
        "sampling": "up to 60 cells per donor without replacement before each network build",
        "network_ensemble": "3 principal-component regression networks per seed",
        "target_screening": "multi-assay pChEMBL>=6 targets expressed in the compartment; independent-document support is annotated separately, and canonical targets are allowed at >=1% donor-mean detection",
        "interpretation": "unsigned network displacement used for prioritisation only; no protective direction or causal effect inferred",
    }
    (OUT / "virtual_knockout_provenance.json").write_text(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    run()
