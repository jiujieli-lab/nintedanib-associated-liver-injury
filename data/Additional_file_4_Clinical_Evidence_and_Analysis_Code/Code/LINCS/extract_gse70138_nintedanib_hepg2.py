#!/usr/bin/env python3
"""Extract the six nintedanib/HepG2/24 h Level-5 signatures from GSE70138.

The GEO Level-5 file is a GCTX (HDF5) matrix.  This script reads only the six
required columns from the on-disk matrix, joins the official GSE70138 gene and
signature metadata, and writes both all-gene and landmark-only matrices in
wide and tidy long form.  It never materializes the complete 118,050-column
matrix in memory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd


PERT_ID = "BRD-K49075727"
EXPECTED_DOSES_UM = (0.04, 0.12, 0.37, 1.11, 3.33, 10.0)
EXPECTED_CELL = "HEPG2"
EXPECTED_TIME = "24 h"
GEO_ACCESSION = "GSE70138"
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/"
    "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
)


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block_size):
            digest.update(chunk)
    return digest.hexdigest()


def decode_h5_strings(values: Iterable[object]) -> list[str]:
    decoded: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))
    return decoded


def dose_label(dose_um: float) -> str:
    return f"z_{dose_um:g}uM_24h".replace(".", "p")


def select_signature_metadata(sig_info_path: Path) -> pd.DataFrame:
    sig = pd.read_csv(sig_info_path, sep="\t", compression="infer", dtype=str)
    chosen = sig.loc[
        (sig["pert_id"] == PERT_ID)
        & (sig["cell_id"] == EXPECTED_CELL)
        & (sig["pert_itime"] == EXPECTED_TIME)
    ].copy()
    chosen["dose_um"] = pd.to_numeric(
        chosen["pert_idose"].str.replace(" um", "", regex=False), errors="raise"
    )
    chosen["replicate_count"] = chosen["distil_id"].str.split("|").str.len()
    chosen = chosen.sort_values("dose_um", kind="stable").reset_index(drop=True)

    actual_doses = tuple(chosen["dose_um"].tolist())
    if len(chosen) != 6 or not np.allclose(actual_doses, EXPECTED_DOSES_UM):
        raise ValueError(
            "Expected exactly six HepG2/24 h signatures at doses "
            f"{EXPECTED_DOSES_UM}; found {actual_doses}."
        )
    if chosen["sig_id"].duplicated().any():
        raise ValueError("Duplicate sig_id values in selected signature metadata.")
    chosen["output_column"] = chosen["dose_um"].map(dose_label)
    return chosen


def load_gene_info(gene_info_path: Path) -> pd.DataFrame:
    genes = pd.read_csv(gene_info_path, sep="\t", compression="infer", dtype=str)
    required = {
        "pr_gene_id",
        "pr_gene_symbol",
        "pr_gene_title",
        "pr_is_lm",
        "pr_is_bing",
    }
    missing = required.difference(genes.columns)
    if missing:
        raise ValueError(f"Missing gene-info columns: {sorted(missing)}")
    if genes["pr_gene_id"].duplicated().any():
        raise ValueError("Official gene-info table contains duplicated pr_gene_id values.")
    return genes.set_index("pr_gene_id", drop=False)


def extract_six_columns(
    gctx_path: Path, signature_ids: list[str]
) -> tuple[np.ndarray, list[str], dict[str, object]]:
    with h5py.File(gctx_path, "r") as handle:
        matrix_path = "/0/DATA/0/matrix"
        row_id_path = "/0/META/ROW/id"
        col_id_path = "/0/META/COL/id"
        for required_path in (matrix_path, row_id_path, col_id_path):
            if required_path not in handle:
                raise KeyError(f"Required GCTX object is absent: {required_path}")

        matrix = handle[matrix_path]
        row_ids = decode_h5_strings(handle[row_id_path][...])
        col_ids = decode_h5_strings(handle[col_id_path][...])
        col_lookup = {value: index for index, value in enumerate(col_ids)}
        absent = [value for value in signature_ids if value not in col_lookup]
        if absent:
            raise KeyError(f"Selected sig_id values absent from GCTX columns: {absent}")

        requested_indices = [col_lookup[value] for value in signature_ids]
        sorted_order = np.argsort(requested_indices)
        sorted_indices = np.asarray(requested_indices, dtype=np.int64)[sorted_order]

        if matrix.shape == (len(row_ids), len(col_ids)):
            sorted_values = np.asarray(matrix[:, sorted_indices], dtype=np.float32)
        elif matrix.shape == (len(col_ids), len(row_ids)):
            sorted_values = np.asarray(matrix[sorted_indices, :], dtype=np.float32).T
        else:
            raise ValueError(
                "GCTX matrix dimensions do not match ROW/COL identifiers: "
                f"matrix={matrix.shape}, rows={len(row_ids)}, cols={len(col_ids)}"
            )

        restore_order = np.argsort(sorted_order)
        values = sorted_values[:, restore_order]
        structure = {
            "matrix_dataset": matrix_path,
            "matrix_shape_on_disk": list(matrix.shape),
            "matrix_dtype_on_disk": str(matrix.dtype),
            "row_id_dataset": row_id_path,
            "column_id_dataset": col_id_path,
            "row_count": len(row_ids),
            "column_count": len(col_ids),
            "selected_column_indices_zero_based": requested_indices,
        }
    return values, row_ids, structure


def attach_gene_metadata(
    values: np.ndarray,
    row_ids: list[str],
    selected: pd.DataFrame,
    genes: pd.DataFrame,
) -> pd.DataFrame:
    missing_gene_ids = [gene_id for gene_id in row_ids if gene_id not in genes.index]
    if missing_gene_ids:
        raise KeyError(
            f"{len(missing_gene_ids)} GCTX row IDs are absent from gene-info metadata; "
            f"first values: {missing_gene_ids[:10]}"
        )
    aligned = genes.loc[row_ids].reset_index(drop=True).copy()
    aligned = aligned.rename(
        columns={
            "pr_gene_id": "gene_id",
            "pr_gene_symbol": "gene_symbol",
            "pr_gene_title": "gene_title",
            "pr_is_lm": "is_landmark",
            "pr_is_bing": "is_bing",
        }
    )
    aligned["is_landmark"] = aligned["is_landmark"].astype(int)
    aligned["is_bing"] = aligned["is_bing"].astype(int)
    for position, output_column in enumerate(selected["output_column"]):
        aligned[output_column] = values[:, position]
    return aligned


def wide_to_long(wide: pd.DataFrame, selected: pd.DataFrame) -> pd.DataFrame:
    identifiers = ["gene_id", "gene_symbol", "gene_title", "is_landmark", "is_bing"]
    signature_columns = selected["output_column"].tolist()
    long = wide.melt(
        id_vars=identifiers,
        value_vars=signature_columns,
        var_name="output_column",
        value_name="level5_zscore",
    )
    signature_annotation = selected[
        [
            "output_column",
            "sig_id",
            "dose_um",
            "pert_itime",
            "cell_id",
            "pert_id",
            "pert_iname",
            "replicate_count",
            "distil_id",
        ]
    ].rename(columns={"sig_id": "signature_id", "pert_itime": "exposure_time"})
    long = long.merge(signature_annotation, on="output_column", how="left", validate="many_to_one")
    return long[
        identifiers
        + [
            "signature_id",
            "pert_id",
            "pert_iname",
            "cell_id",
            "dose_um",
            "exposure_time",
            "replicate_count",
            "distil_id",
            "level5_zscore",
        ]
    ]


def write_checksum_manifest(paths: list[Path], manifest_path: Path) -> None:
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        for path in sorted(paths, key=lambda value: value.name):
            handle.write(f"{sha256_file(path)}  {path.name}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gctx", required=True, type=Path, help="Decompressed GSE70138 Level-5 .gctx file")
    parser.add_argument("--sig-info", required=True, type=Path)
    parser.add_argument("--gene-info", required=True, type=Path)
    parser.add_argument("--compressed-source", type=Path, default=None)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected = select_signature_metadata(args.sig_info)
    genes = load_gene_info(args.gene_info)
    values, row_ids, matrix_structure = extract_six_columns(
        args.gctx, selected["sig_id"].tolist()
    )
    wide = attach_gene_metadata(values, row_ids, selected, genes)
    landmark = wide.loc[wide["is_landmark"] == 1].copy()

    if len(wide) != 12_328:
        raise ValueError(f"Expected 12,328 GCTX genes; found {len(wide):,}.")
    if len(landmark) != 978:
        raise ValueError(f"Expected all 978 landmark genes; found {len(landmark):,}.")
    if values.shape != (12_328, 6) or not np.isfinite(values).all():
        raise ValueError(f"Extracted values failed shape/finiteness checks: {values.shape}")

    signature_metadata_path = args.output_dir / "nintedanib_hepg2_signature_metadata.csv"
    all_wide_path = args.output_dir / "nintedanib_hepg2_level5_all_genes_wide.csv.gz"
    all_long_path = args.output_dir / "nintedanib_hepg2_level5_all_genes_tidy.csv.gz"
    lm_wide_path = args.output_dir / "nintedanib_hepg2_level5_landmark_wide.csv"
    lm_long_path = args.output_dir / "nintedanib_hepg2_level5_landmark_tidy.csv"
    manifest_json_path = args.output_dir / "extraction_manifest.json"

    selected.rename(columns={"sig_id": "signature_id"}).to_csv(
        signature_metadata_path, index=False, quoting=csv.QUOTE_MINIMAL
    )
    wide.to_csv(all_wide_path, index=False, compression="gzip", float_format="%.7g")
    wide_to_long(wide, selected).to_csv(
        all_long_path, index=False, compression="gzip", float_format="%.7g"
    )
    landmark.to_csv(lm_wide_path, index=False, float_format="%.7g")
    wide_to_long(landmark, selected).to_csv(
        lm_long_path, index=False, float_format="%.7g"
    )

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "geo_accession": GEO_ACCESSION,
        "perturbagen": "nintedanib",
        "pert_id": PERT_ID,
        "cell_id": EXPECTED_CELL,
        "exposure_time": EXPECTED_TIME,
        "dose_um": selected["dose_um"].tolist(),
        "signature_ids": selected["sig_id"].tolist(),
        "replicate_counts": selected["replicate_count"].astype(int).tolist(),
        "matrix_source_url": MATRIX_URL,
        "gctx_path": str(args.gctx),
        "gctx_sha256": sha256_file(args.gctx),
        "compressed_source_path": str(args.compressed_source) if args.compressed_source else None,
        "compressed_source_sha256": (
            sha256_file(args.compressed_source) if args.compressed_source else None
        ),
        "signature_info_path": str(args.sig_info),
        "signature_info_sha256": sha256_file(args.sig_info),
        "gene_info_path": str(args.gene_info),
        "gene_info_sha256": sha256_file(args.gene_info),
        "matrix_structure": matrix_structure,
        "extracted_all_gene_rows": int(len(wide)),
        "extracted_landmark_gene_rows": int(len(landmark)),
        "tidy_all_gene_rows": int(len(wide) * len(selected)),
        "tidy_landmark_gene_rows": int(len(landmark) * len(selected)),
        "quality_checks": {
            "all_selected_signatures_found": True,
            "expected_six_doses_found": True,
            "all_values_finite": True,
            "all_978_landmark_genes_retained": True,
        },
    }
    manifest_json_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    generated = [
        signature_metadata_path,
        all_wide_path,
        all_long_path,
        lm_wide_path,
        lm_long_path,
        manifest_json_path,
    ]
    write_checksum_manifest(generated, args.output_dir / "SHA256SUMS")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
