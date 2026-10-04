#!/usr/bin/env python3
"""Prepare de novo nintedanib pharmacology and observed perturbation evidence.

This script intentionally consumes only the freshly downloaded ChEMBL activity
pages and Harmonizome/Sci-Plex pages. It does not read any previous manuscript,
ranked target list, figure, or analysis output.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "fresh_data"
INPUT = ROOT / "fresh_analysis" / "inputs"
OUT = ROOT / "fresh_analysis" / "pharmacology"
OUT.mkdir(parents=True, exist_ok=True)
CANONICAL_TARGETS = {"KDR", "FLT1", "FLT4", "FGFR1", "FGFR2", "FGFR3", "PDGFRA", "PDGFRB"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_activities() -> pd.DataFrame:
    paths = [
        DATA / "chembl_nintedanib_activities.json",
        INPUT / "chembl_nintedanib_activities_offset1000.json",
        INPUT / "chembl_nintedanib_activities_offset2000.json",
    ]
    rows: list[dict] = []
    manifest = []
    for path in paths:
        payload = json.loads(path.read_text())
        rows.extend(payload["activities"])
        manifest.append(
            {
                "source": "ChEMBL activity API",
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "retrieval_date": "2026-09-04",
                "n_records": len(payload["activities"]),
            }
        )
    pd.DataFrame(manifest).to_csv(OUT / "source_manifest.csv", index=False)
    frame = pd.DataFrame(rows).drop_duplicates("activity_id")
    frame["pchembl_value"] = pd.to_numeric(frame["pchembl_value"], errors="coerce")
    frame["standard_value"] = pd.to_numeric(frame["standard_value"], errors="coerce")
    return frame


def fetch_target_metadata(ids: list[str]) -> list[dict]:
    cache = INPUT / "chembl_nintedanib_target_metadata.json"
    if cache.exists():
        return json.loads(cache.read_text())["targets"]
    endpoint = "https://www.ebi.ac.uk/chembl/api/data/target.json"
    all_targets: list[dict] = []
    for start in range(0, len(ids), 40):
        chunk = ids[start : start + 40]
        response = requests.get(
            endpoint,
            params={"target_chembl_id__in": ",".join(chunk), "limit": 100},
            timeout=120,
        )
        response.raise_for_status()
        all_targets.extend(response.json()["targets"])
        time.sleep(0.2)
    cache.write_text(
        json.dumps(
            {
                "retrieval_date": "2026-09-04",
                "endpoint": endpoint,
                "targets": all_targets,
            },
            indent=2,
        )
    )
    return all_targets


def gene_symbol(target: dict) -> tuple[str | None, str | None]:
    components = target.get("target_components") or []
    if len(components) != 1:
        return None, None
    component = components[0]
    symbols = [
        x["component_synonym"]
        for x in component.get("target_component_synonyms", [])
        if x.get("syn_type") == "GENE_SYMBOL"
    ]
    return (symbols[0] if symbols else None), component.get("accession")


def prepare_pharmacology(frame: pd.DataFrame) -> pd.DataFrame:
    # High-confidence means human, single-protein, exact quantitative relation,
    # valid pChEMBL value, and no ChEMBL validity warning. This includes broad
    # profiling assays, so evidence density and replicate consistency remain
    # explicit downstream rather than converting every hit to a clinical target.
    high = frame.loc[
        frame["target_organism"].eq("Homo sapiens")
        & frame["bao_label"].eq("single protein format")
        & frame["standard_relation"].eq("=")
        & frame["pchembl_value"].notna()
        & frame["data_validity_comment"].isna()
        & frame["potential_duplicate"].ne(1)
    ].copy()
    ids = sorted(high["target_chembl_id"].dropna().unique())
    targets = fetch_target_metadata(ids)
    map_rows = []
    for target in targets:
        symbol, accession = gene_symbol(target)
        map_rows.append(
            {
                "target_chembl_id": target["target_chembl_id"],
                "gene_symbol": symbol,
                "uniprot_accession": accession,
                "target_type": target.get("target_type"),
                "target_metadata_pref_name": target.get("pref_name"),
            }
        )
    mapping = pd.DataFrame(map_rows)
    high = high.merge(mapping, on="target_chembl_id", how="left")
    high.to_csv(OUT / "chembl_human_single_protein_activities.csv", index=False)

    summary = (
        high.groupby(
            [
                "target_chembl_id",
                "gene_symbol",
                "uniprot_accession",
                "target_pref_name",
            ],
            dropna=False,
        )
        .agg(
            n_activities=("activity_id", "size"),
            n_documents=("document_chembl_id", "nunique"),
            n_assays=("assay_chembl_id", "nunique"),
            max_pchembl=("pchembl_value", "max"),
            median_pchembl=("pchembl_value", "median"),
            min_pchembl=("pchembl_value", "min"),
            q25_pchembl=("pchembl_value", lambda x: x.quantile(0.25)),
            q75_pchembl=("pchembl_value", lambda x: x.quantile(0.75)),
        )
        .reset_index()
    )
    summary["potent_pchembl_ge_6"] = summary["max_pchembl"].ge(6)
    summary["replicated_potent"] = (
        summary["max_pchembl"].ge(6)
        & summary["n_activities"].ge(2)
        & summary["n_assays"].ge(2)
    )
    # The historical column name above is retained for reproducibility, but
    # multiple assays can originate from one broad kinome-screen publication.
    # Independent-source support is therefore recorded separately and is the
    # evidence used for primary pharmacology interpretation.
    summary["multi_document_potent"] = (
        summary["max_pchembl"].ge(6)
        & summary["n_documents"].ge(2)
        & summary["n_assays"].ge(2)
    )
    summary["pharmacology_tier"] = np.select(
        [
            summary["gene_symbol"].isin(CANONICAL_TARGETS) & summary["max_pchembl"].ge(6),
            summary["multi_document_potent"],
            summary["max_pchembl"].ge(6),
        ],
        [
            "canonical_target",
            "multi_document_noncanonical",
            "single_document_or_single_assay_screen",
        ],
        default="weak_or_incomplete",
    )
    summary = summary.sort_values(
        ["replicated_potent", "max_pchembl", "n_activities"],
        ascending=[False, False, False],
    )
    summary.to_csv(OUT / "chembl_target_summary.csv", index=False)
    return summary


def prepare_sciplex() -> pd.DataFrame:
    pages = {
        "1uM": DATA / "harmonizome_nintedanib_mcf7_1uM_24h.html",
        "10uM": INPUT / "harmonizome_nintedanib_mcf7_10uM_24h.html",
    }
    rows = []
    for dose, path in pages.items():
        tables = pd.read_html(path)
        gene_tables = [
            table
            for table in tables
            if {"Symbol", "Standardized Value"}.issubset(table.columns)
        ]
        if len(gene_tables) != 2:
            raise RuntimeError(f"Expected two gene tables in {path}; got {len(gene_tables)}")
        for table in gene_tables:
            part = table.copy()
            part["direction"] = np.where(part["Standardized Value"] > 0, "up", "down")
            part["dose"] = dose
            part["cell_line"] = "MCF7"
            part["duration_h"] = 24
            part["source_sha256"] = sha256(path)
            rows.append(part)
    long = pd.concat(rows, ignore_index=True).rename(
        columns={"Symbol": "gene_symbol", "Name": "gene_name", "Standardized Value": "z"}
    )
    long["abs_z"] = long["z"].abs()
    long["rank_within_dose"] = long.groupby("dose")["abs_z"].rank(
        method="min", ascending=False
    )
    long.to_csv(OUT / "sciplex_nintedanib_signature_long.csv", index=False)

    wide = long.pivot_table(index="gene_symbol", columns="dose", values="z", aggfunc="first")
    wide["observed_in_both_doses"] = wide.notna().all(axis=1)
    wide["direction_concordant"] = np.where(
        wide["observed_in_both_doses"], np.sign(wide["1uM"]) == np.sign(wide["10uM"]), np.nan
    )
    wide["mean_z"] = wide[["1uM", "10uM"]].mean(axis=1, skipna=True)
    wide["min_abs_z"] = wide[["1uM", "10uM"]].abs().min(axis=1, skipna=True)
    wide.reset_index().to_csv(OUT / "sciplex_dose_concordance.csv", index=False)
    return long


def main() -> None:
    activities = read_activities()
    activities.to_csv(OUT / "chembl_all_nintedanib_activities.csv", index=False)
    targets = prepare_pharmacology(activities)
    signature = prepare_sciplex()
    summary = {
        "n_chembl_activities": int(len(activities)),
        "n_human_single_protein_targets": int(len(targets)),
        "n_replicated_potent_targets": int(targets["replicated_potent"].sum()),
        "n_multi_document_potent_targets": int(targets["multi_document_potent"].sum()),
        "n_sciplex_rows": int(len(signature)),
        "n_sciplex_unique_genes": int(signature["gene_symbol"].nunique()),
        "interpretation_boundary": (
            "Sci-Plex signatures are observed nintedanib responses in MCF7 cells, "
            "not hepatocyte DILI effects; ChEMBL biochemical potency does not by "
            "itself establish in-vivo target engagement or toxicity causality."
        ),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
