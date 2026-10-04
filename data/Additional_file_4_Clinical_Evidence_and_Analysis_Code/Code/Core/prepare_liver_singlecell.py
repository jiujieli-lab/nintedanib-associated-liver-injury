#!/usr/bin/env python3
"""Prepare a de novo human-liver single-cell analysis for GSE115469.

The public matrix contains normalized expression for 8,444 cells from five
healthy donors. The donor, rather than the cell, is the unit for uncertainty
summaries. Cell-level embeddings are descriptive only.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "fresh_data"
OUT = ROOT / "fresh_analysis" / "liver_singlecell"
OUT.mkdir(parents=True, exist_ok=True)

DILI_ANCHORS = [
    "ACO1", "ASS1", "FAH", "CPS1", "ALDOB", "HPD", "OTC", "DMGDH",
    "GSTA1", "FBP1", "PCK2", "CES1", "LECT2",
]

CANONICAL_NINTEDANIB_TARGETS = [
    "KDR", "FLT1", "FLT4", "FGFR1", "FGFR2", "FGFR3", "PDGFRA", "PDGFRB",
]

MECHANISM_GENES = [
    "ABCB11", "ABCC2", "SLC10A1", "UGT1A1", "CES1", "CES2", "ABCB1",
    "NFE2L2", "KEAP1", "HMOX1", "NQO1", "GCLC", "GCLM", "TXNIP", "TXNRD1",
    "GPX4", "SLC7A11", "SOD2", "CAT", "BAX", "BCL2", "CASP3", "CASP8",
    "CYP3A4", "CYP2C9", "CYP2D6", "ALB", "KRT8", "KRT18", "VIM",
    "COL1A1", "COL1A2", "ACTA2", "TGFB1", "TGFB2", "SMAD2", "SMAD3",
    "VEGFA", "PDGFA", "PDGFB", "FGF2", "KDR", "FLT1", "FLT4",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def broad_compartment(cell_type: str) -> str:
    if cell_type.startswith("Hepatocyte"):
        return "Hepatocyte"
    if "LSEC" in cell_type or cell_type == "Portal_endothelial_Cells":
        return "Endothelial"
    if "Macrophage" in cell_type:
        return "Macrophage"
    if cell_type == "Cholangiocytes":
        return "Cholangiocyte"
    if cell_type == "Hepatic_Stellate_Cells":
        return "Stellate"
    return "Immune_other"


def donor_summary(expr: pd.DataFrame, meta: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    genes = [g for g in genes if g in expr.index]
    for compartment, sub in meta.groupby("compartment", sort=False):
        for donor, dm in sub.groupby("Sample", sort=False):
            cells = dm["CellName"].tolist()
            values = expr.loc[genes, cells]
            means = values.mean(axis=1)
            detects = values.gt(0).mean(axis=1)
            for gene in genes:
                rows.append(
                    {
                        "gene_symbol": gene,
                        "compartment": compartment,
                        "donor": donor,
                        "n_cells": len(cells),
                        "mean_expression": float(means[gene]),
                        "detection_fraction": float(detects[gene]),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    meta_path = DATA / "GSE115469_CellClusterType.txt.gz"
    expr_path = DATA / "GSE115469_Data.csv.gz"
    meta = pd.read_csv(meta_path, sep="\t")
    meta["compartment"] = meta["CellType"].map(broad_compartment)
    meta.to_csv(OUT / "cell_metadata.csv", index=False)

    # Explicit float32 limits memory while preserving more numerical precision
    # than is needed for the deposited normalized matrix. A per-column mapping
    # is required because the first CSV column stores gene names.
    header = pd.read_csv(expr_path, nrows=0).columns
    numeric_dtypes = {column: np.float32 for column in header[1:]}
    expr = pd.read_csv(expr_path, index_col=0, dtype=numeric_dtypes)
    expr.index = expr.index.astype(str)
    if not set(meta["CellName"]).issubset(expr.columns):
        missing = sorted(set(meta["CellName"]) - set(expr.columns))
        raise RuntimeError(f"Metadata cells absent from expression matrix: {missing[:5]}")
    expr = expr.loc[~expr.index.duplicated(keep="first")]

    pharm = pd.read_csv(ROOT / "fresh_analysis" / "pharmacology" / "chembl_target_summary.csv")
    replicated = pharm.loc[pharm["replicated_potent"].eq(True), "gene_symbol"].dropna().tolist()
    forced = sorted(set(DILI_ANCHORS + CANONICAL_NINTEDANIB_TARGETS + MECHANISM_GENES + replicated))

    present = [g for g in forced if g in expr.index]
    missing = sorted(set(forced) - set(present))
    pd.DataFrame({"requested_gene": forced, "present": [g in expr.index for g in forced]}).to_csv(
        OUT / "forced_gene_audit.csv", index=False
    )

    ds = donor_summary(expr, meta, present)
    ds.to_csv(OUT / "donor_level_target_expression.csv", index=False)
    summary = (
        ds.groupby(["gene_symbol", "compartment"], as_index=False)
        .agg(
            n_donors=("donor", "nunique"),
            n_cells=("n_cells", "sum"),
            mean_expression=("mean_expression", "mean"),
            min_donor_expression=("mean_expression", "min"),
            max_donor_expression=("mean_expression", "max"),
            mean_detection_fraction=("detection_fraction", "mean"),
            min_donor_detection=("detection_fraction", "min"),
            max_donor_detection=("detection_fraction", "max"),
        )
    )
    summary.to_csv(OUT / "target_expression_by_compartment.csv", index=False)

    # Donor-level, cell-type-specific summaries for zonation and reproducibility.
    type_rows = []
    for (cell_type, donor), sub in meta.groupby(["CellType", "Sample"], sort=False):
        cells = sub["CellName"].tolist()
        values = expr.loc[present, cells]
        for gene, mean, detect in zip(
            present, values.mean(axis=1).to_numpy(), values.gt(0).mean(axis=1).to_numpy()
        ):
            type_rows.append(
                {
                    "gene_symbol": gene,
                    "cell_type": cell_type,
                    "donor": donor,
                    "n_cells": len(cells),
                    "mean_expression": float(mean),
                    "detection_fraction": float(detect),
                }
            )
    pd.DataFrame(type_rows).to_csv(OUT / "donor_celltype_target_expression.csv", index=False)

    selected_compartments = ["Hepatocyte", "Endothelial", "Macrophage"]
    matrix_manifest = []
    embedding_cells = []
    hvg_union: set[str] = set()
    for compartment in selected_compartments:
        cells = meta.loc[meta["compartment"].eq(compartment), "CellName"].tolist()
        sub = expr.loc[:, cells]
        detected = sub.gt(0).mean(axis=1)
        variance = sub.var(axis=1)
        eligible = detected.ge(0.01) & variance.gt(0)
        top_hvg = variance.loc[eligible].nlargest(650).index.tolist()
        genes = list(dict.fromkeys(top_hvg + present))
        # Remove zero-variance forced genes because network regression cannot use them.
        genes = [g for g in genes if variance.get(g, 0) > 0]
        hvg_union.update(genes)
        matrix = sub.loc[genes].astype(np.float32)
        out_path = OUT / f"{compartment.lower()}_network_matrix.csv.gz"
        matrix.to_csv(out_path, compression="gzip")
        matrix_manifest.append(
            {
                "compartment": compartment,
                "n_cells": len(cells),
                "n_donors": int(meta.loc[meta["compartment"].eq(compartment), "Sample"].nunique()),
                "n_genes": len(genes),
                "matrix_path": str(out_path.relative_to(ROOT)),
                "sha256": sha256(out_path),
            }
        )

        # Descriptive PCA with a per-donor cap (not equal donor weighting) for a compact embedding.
        chosen = []
        for donor, dm in meta.loc[meta["compartment"].eq(compartment)].groupby("Sample"):
            n_take = min(250, len(dm))
            chosen.extend(dm.sample(n=n_take, random_state=20260904)["CellName"].tolist())
        pca_genes = top_hvg[:300]
        x = sub.loc[pca_genes, chosen].T.to_numpy(dtype=np.float32)
        x = StandardScaler(with_mean=True, with_std=True).fit_transform(x)
        coords = PCA(n_components=2, random_state=20260904).fit_transform(x)
        e = meta.set_index("CellName").loc[chosen, ["Sample", "CellType", "compartment"]].reset_index()
        e["PC1"] = coords[:, 0]
        e["PC2"] = coords[:, 1]
        embedding_cells.append(e)

    pd.DataFrame(matrix_manifest).to_csv(OUT / "network_matrix_manifest.csv", index=False)
    pd.concat(embedding_cells, ignore_index=True).to_csv(OUT / "descriptive_pca_coordinates.csv", index=False)

    manifest = {
        "accession": "GSE115469",
        "retrieval_date": "2026-09-04",
        "expression_sha256": sha256(expr_path),
        "metadata_sha256": sha256(meta_path),
        "n_cells": int(expr.shape[1]),
        "n_genes": int(expr.shape[0]),
        "n_donors": int(meta["Sample"].nunique()),
        "cell_type_counts": meta["CellType"].value_counts().to_dict(),
        "compartment_counts": meta["compartment"].value_counts().to_dict(),
        "dili_anchors": DILI_ANCHORS,
        "canonical_nintedanib_targets": CANONICAL_NINTEDANIB_TARGETS,
        "missing_forced_genes": missing,
        "statistical_unit": "donor for uncertainty summaries; cells descriptive only",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
