#!/usr/bin/env python3
"""Read-only verification of the extracted nintedanib-DILI release capsule.

This script intentionally uses only the Python standard library. It validates
package hashes, frozen structural requirements, notebook execution state, and a
small set of decisive scientific invariants. It does not download third-party
data or rerun computationally intensive upstream analyses.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(relative: str) -> list[dict[str, str]]:
    with (ROOT / relative).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "t", "yes"}


def close(observed: float, expected: float, tolerance: float = 5e-12) -> bool:
    return math.isclose(observed, expected, rel_tol=tolerance, abs_tol=tolerance)


def verify_hash_manifest() -> int:
    manifest = ROOT / "SHA256SUMS"
    assert manifest.is_file(), "SHA256SUMS is missing"
    manifest_entries: dict[str, str] = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        path = ROOT / relative
        assert path.is_file(), f"Manifest target missing: {relative}"
        observed = sha256(path)
        assert observed == expected, f"SHA256 mismatch: {relative}"
        assert relative not in manifest_entries, f"Duplicate SHA256SUMS entry: {relative}"
        manifest_entries[relative] = expected

    inventory_rows = read_csv("PACKAGE_INVENTORY.csv")
    inventory_entries = {row["relative_path"]: row for row in inventory_rows}
    assert len(inventory_entries) == len(inventory_rows), "Duplicate inventory path"
    assert set(inventory_entries) == set(manifest_entries), (
        "PACKAGE_INVENTORY.csv and SHA256SUMS path sets differ"
    )
    for relative, row in inventory_entries.items():
        path = ROOT / relative
        assert int(row["bytes"]) == path.stat().st_size, f"Inventory byte mismatch: {relative}"
        assert row["sha256"] == manifest_entries[relative], f"Inventory hash mismatch: {relative}"

    actual_files = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file()
    }
    symlinks = [
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_symlink()
    ]
    assert not symlinks, f"Release must not contain symlinks: {symlinks}"
    expected_files = set(manifest_entries) | {"PACKAGE_INVENTORY.csv", "SHA256SUMS"}
    assert actual_files == expected_files, (
        f"Actual file set differs from signed inventory; extra={sorted(actual_files - expected_files)}, "
        f"missing={sorted(expected_files - actual_files)}"
    )
    validation = json.loads((ROOT / "PACKAGE_VALIDATION.json").read_text())
    assert validation["final_file_count"] == len(actual_files)
    return len(manifest_entries)


def verify_structure() -> tuple[int, int]:
    validation = json.loads((ROOT / "PACKAGE_VALIDATION.json").read_text())
    assert validation["status"] == "PASS"
    assert all(validation["checks"].values())

    main_pdfs = sorted((ROOT / "Figures").glob("Figure_[1-6]_*.pdf"))
    supplementary_pdfs = sorted(
        (ROOT / "Supplementary_Figures").glob("Supplementary_Figure_S[1-7]_*.pdf")
    )
    assert len(main_pdfs) == 6, f"Expected six main-figure PDFs; found {len(main_pdfs)}"
    assert len(supplementary_pdfs) == 7, (
        f"Expected seven supplementary-figure PDFs; found {len(supplementary_pdfs)}"
    )
    for stem in [p.stem for p in main_pdfs]:
        for suffix in [".svg", ".png", ".tiff"]:
            assert (ROOT / "Figures" / f"{stem}{suffix}").is_file()
    for stem in [p.stem for p in supplementary_pdfs]:
        for suffix in [".svg", ".png", ".tiff"]:
            assert (ROOT / "Supplementary_Figures" / f"{stem}{suffix}").is_file()
    for suffix in [".pdf", ".svg", ".png", ".tiff"]:
        assert (ROOT / f"Graphical_Abstract_Public_Data_Target_Map{suffix}").is_file()
    return len(main_pdfs), len(supplementary_pdfs)


def verify_manuscript_and_workbook() -> tuple[int, int]:
    expected_title = (
        "Role-separated public-data triangulation prioritizes target–phenotype "
        "hypotheses for nintedanib-associated liver injury"
    )
    manuscript = (ROOT / "Main_Manuscript.md").read_text(encoding="utf-8")
    assert manuscript.splitlines()[0] == f"# {expected_title}"
    assert "FGFR1" in manuscript and "sole candidate meeting the prespecified exploratory T2" in manuscript
    assert "No T1 target was assigned" in manuscript
    assert "strict all-seed" in manuscript
    references = re.findall(r"(?m)^(\d+)\.\s", manuscript)
    assert references == [str(i) for i in range(1, 76)], "Expected 75 consecutive references"

    cited = read_csv("References/Cited_References_First_Appearance.csv")
    assert len(cited) == 75
    workbook = ROOT / "Supplementary_Tables.xlsx"
    assert workbook.is_file() and workbook.stat().st_size > 100_000
    with zipfile.ZipFile(workbook) as archive:
        table_parts = [
            name for name in archive.namelist()
            if re.fullmatch(r"xl/tables/table\d+\.xml", name)
        ]
    assert len(table_parts) == 107, f"Expected 107 structured workbook tables; found {len(table_parts)}"
    return len(references), len(table_parts)


def verify_notebook() -> tuple[int, int]:
    notebook = json.loads((ROOT / "Public_Data_Reproducibility_Notebook.ipynb").read_text())
    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
    cell_ids = [cell.get("id") for cell in notebook["cells"]]
    counts = [cell.get("execution_count") for cell in code_cells]
    assert counts == list(range(1, len(code_cells) + 1)), "Notebook is not sequentially executed"
    error_outputs = [
        output
        for cell in code_cells
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    assert not error_outputs, "Notebook contains error output"
    assert all(isinstance(cell_id, str) and cell_id for cell_id in cell_ids), "Notebook cell ID missing"
    assert len(cell_ids) == len(set(cell_ids)), "Notebook cell IDs are not unique"
    return len(notebook["cells"]), len(code_cells)


def verify_scientific_invariants() -> list[str]:
    faers = read_csv("Source_Data/FAERS/overall_active_comparator_signals.csv")
    by_endpoint = {row["endpoint"]: row for row in faers}
    assert close(float(by_endpoint["narrow"]["ROR"]), 3.9966057509880484)
    assert close(float(by_endpoint["broad"]["ROR"]), 3.3668703340517556)

    heterogeneity = read_csv("Source_Data/FAERS/year_effect_heterogeneity_tests.csv")
    assert len(heterogeneity) == 8
    heterogeneity_by_key = {
        (row["endpoint"], row["window"], row["test"]): row
        for row in heterogeneity
    }
    narrow_bd = heterogeneity_by_key[("narrow", "all_years", "Tarone-adjusted Breslow-Day")]
    broad_q = heterogeneity_by_key[("broad", "all_years", "Inverse-variance Cochran Q")]
    assert close(float(narrow_bd["statistic"]), 88.25434621529466)
    assert close(float(narrow_bd["p_value"]), 3.6637359812630166e-14)
    assert int(narrow_bd["degrees_of_freedom"]) == 11
    assert close(float(broad_q["statistic"]), 61.39384545688669)
    assert close(float(broad_q["p_value"]), 1.256004261089914e-8)
    assert close(float(broad_q["i2_percent"]), 80.45406683569463)

    temporal = read_csv("Source_Data/FAERS/temporal_influence_signals.csv")
    temporal_by_key = {(row["endpoint"], row["temporal_window"]): row for row in temporal}
    assert int(temporal_by_key[("narrow", "year_2022_only")]["pirfenidone_event"]) == 0
    assert close(float(temporal_by_key[("narrow", "exclude_2022_posthoc")]["ROR"]), 2.5503789028415342)

    specifications = read_csv("Source_Data/FAERS/report_level_specification_signals.csv")
    suspect_ids = {
        row["analysis"] for row in specifications
        if "suspect_drug_report_filter" in row["analysis"]
    }
    assert suspect_ids == {
        "suspect_drug_report_filter", "ild_plus_suspect_drug_report_filter"
    }
    assert all(
        'drugcharacterization:"1"' in row["filter_query"]
        for row in specifications
        if "suspect_drug_report_filter" in row["analysis"]
    )

    role_profiles = read_csv("Source_Data/FAERS/drug_role_and_indication_profiles.csv")
    role_lookup = {
        (row["endpoint"], row["drug"], row["variable"], row["level"]): int(row["reports"])
        for row in role_profiles
    }
    assert role_lookup[("narrow", "nintedanib", "matching_drug_role", "suspect")] == 330
    assert role_lookup[("narrow", "pirfenidone", "matching_drug_role", "suspect")] == 105
    assert role_lookup[("broad", "nintedanib", "matching_drug_ild_indication", "True")] == 1373
    assert role_lookup[("broad", "pirfenidone", "matching_drug_ild_indication", "True")] == 604

    proteomics = read_csv("Source_Data/Human_DILI_Proteomics/effect_estimates.csv")
    confirmatory = [
        row for row in proteomics
        if row["cohort"] == "confirmatory" and row["contrast"] == "DO_vs_NDO"
    ]
    fdr_hits = [row["protein"] for row in confirmatory if float(row["mann_whitney_q_bh"]) < 0.05]
    assert fdr_hits == ["FBP1"], f"Unexpected proteomic FDR set: {fdr_hits}"

    ranking = read_csv("Source_Data/Integration/candidate_consensus_ranking.csv")
    t2 = [
        row["gene_symbol"] for row in ranking
        if row["candidate_class"] == "exposure_proximal_target"
        and as_bool(row["primary_T2_gate"])
    ]
    strict = [row["gene_symbol"] for row in ranking if as_bool(row["strict_all_seed_T2_sensitivity_gate"])]
    assert t2 == ["FGFR1"], f"Unexpected exploratory T2-gate set: {t2}"
    assert strict == [], f"Strict all-seed T2 set must be empty: {strict}"

    blood = read_csv("Source_Data/Whole_Blood/target_longitudinal_summary.csv")
    assert len(blood) == 30 and all(float(row["q_bh"]) >= 0.05 for row in blood)

    return [
        "FAERS narrow and broad RORs",
        "FAERS year-effect heterogeneity, 2022 influence, suspect-role and ILD sensitivities",
        "FBP1-only confirmatory proteomic FDR set",
        "FGFR1-only exploratory T2-gate set, empty strict all-seed set, and no T1",
        "whole-blood negative-control FDR result",
    ]


def main() -> None:
    hash_count = verify_hash_manifest()
    main_count, supplementary_count = verify_structure()
    reference_count, workbook_table_count = verify_manuscript_and_workbook()
    notebook_cells, notebook_code_cells = verify_notebook()
    invariants = verify_scientific_invariants()
    print("PASS: portable release capsule verified")
    print(f"  package_root: {ROOT}")
    print(f"  sha256_entries: {hash_count}")
    print(f"  figures: {main_count} main, {supplementary_count} supplementary")
    print(f"  manuscript_references: {reference_count}")
    print(f"  workbook_structured_tables: {workbook_table_count}")
    print(f"  notebook: {notebook_cells} cells, {notebook_code_cells} executed code cells")
    print(f"  decisive_scientific_invariants: {len(invariants)}")
    for item in invariants:
        print(f"    - {item}")


if __name__ == "__main__":
    main()
