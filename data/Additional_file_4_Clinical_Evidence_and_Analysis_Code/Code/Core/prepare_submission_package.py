#!/usr/bin/env python3
"""Freeze and optionally archive the already assembled release capsule.

This script never recopies upstream development trees. It signs the current
portable package in place, validates its journal-facing structure and decisive
release invariants, and can create a deterministic ZIP beside the package.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INVENTORY = ROOT / "PACKAGE_INVENTORY.csv"
SUMS = ROOT / "SHA256SUMS"
VALIDATION = ROOT / "PACKAGE_VALIDATION.json"
ADMIN_FILES = {"PACKAGE_INVENTORY.csv", "SHA256SUMS", "PACKAGE_VALIDATION.json"}
TITLE = (
    "Role-separated public-data triangulation prioritizes target–phenotype "
    "hypotheses for nintedanib-associated liver injury"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def regular_files(*, exclude_admin: bool = False) -> list[Path]:
    files = [path for path in ROOT.rglob("*") if path.is_file() and not path.is_symlink()]
    if exclude_admin:
        files = [path for path in files if path.relative_to(ROOT).as_posix() not in ADMIN_FILES]
    return sorted(files, key=lambda path: path.relative_to(ROOT).as_posix())


def reference_count() -> int:
    text = (ROOT / "Main_Manuscript.md").read_text(encoding="utf-8")
    references = re.findall(r"(?m)^(\d+)\.\s", text)
    assert references == [str(i) for i in range(1, 76)], "References must be consecutively numbered 1–75"
    return len(references)


def workbook_table_count() -> int:
    with zipfile.ZipFile(ROOT / "Supplementary_Tables.xlsx") as archive:
        tables = [
            name for name in archive.namelist()
            if re.fullmatch(r"xl/tables/table\d+\.xml", name)
        ]
    return len(tables)


def notebook_counts() -> tuple[int, int]:
    notebook = json.loads((ROOT / "Public_Data_Reproducibility_Notebook.ipynb").read_text(encoding="utf-8"))
    code_cells = [cell for cell in notebook["cells"] if cell.get("cell_type") == "code"]
    cell_ids = [cell.get("id") for cell in notebook["cells"]]
    counts = [cell.get("execution_count") for cell in code_cells]
    errors = [
        output for cell in code_cells for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    assert counts == list(range(1, len(code_cells) + 1)), "Notebook execution counts are not sequential"
    assert not errors, "Notebook contains an error output"
    assert all(isinstance(cell_id, str) and cell_id for cell_id in cell_ids), "Notebook cell ID missing"
    assert len(cell_ids) == len(set(cell_ids)), "Notebook cell IDs are not unique"
    return len(notebook["cells"]), len(code_cells)


def figure_checks() -> dict[str, bool]:
    main_pdfs = sorted((ROOT / "Figures").glob("Figure_[1-6]_*.pdf"))
    supp_pdfs = sorted((ROOT / "Supplementary_Figures").glob("Supplementary_Figure_S[1-7]_*.pdf"))
    formats = (".pdf", ".svg", ".png", ".tiff")

    main_complete = len(main_pdfs) == 6 and all(
        all(path.with_suffix(suffix).is_file() for suffix in formats)
        for path in main_pdfs
    )
    supp_complete = len(supp_pdfs) == 7 and all(
        all(path.with_suffix(suffix).is_file() for suffix in formats)
        for path in supp_pdfs
    )
    panel_sets = []
    for pdf in main_pdfs:
        svg = pdf.with_suffix(".svg").read_text(encoding="utf-8")
        labels = {
            letter for letter in "ABCDEFGHIJ"
            if re.search(rf"(?:>\s*{letter}\s*<|<!--\s*{letter}\s*-->)", svg)
        }
        panel_sets.append(labels == set("ABCDEFGHIJ"))
    graphical_abstract = all(
        (ROOT / f"Graphical_Abstract_Public_Data_Target_Map{suffix}").is_file()
        for suffix in formats
    )
    return {
        "six_main_figure_sets_four_formats": main_complete,
        "each_main_figure_has_A_to_J": len(panel_sets) == 6 and all(panel_sets),
        "seven_supplementary_figure_sets_four_formats": supp_complete,
        "graphical_abstract_four_formats": graphical_abstract,
    }


def build_validation() -> dict:
    required = {
        "Main_Manuscript.md", "Main_Manuscript.docx", "Main_Manuscript.pdf",
        "Supplementary_Information.md", "Supplementary_Information.pdf",
        "Supplementary_Tables.xlsx", "Public_Data_Reproducibility_Notebook.ipynb",
        "Notebook_Validation_Report.md", "README.md", "Cover_Letter.md",
        "Data_and_Code_Availability.md", "Experimental_Validation_Protocol.md",
        "Checklists/PrePackage_Statistical_Editorial_Audit.json",
        "Checklists/PrePackage_Statistical_Editorial_Audit.md",
    }
    current = {path.relative_to(ROOT).as_posix() for path in regular_files(exclude_admin=True)}
    symlinks = [path.relative_to(ROOT).as_posix() for path in ROOT.rglob("*") if path.is_symlink()]
    forbidden_parts = {"__pycache__", "node_modules", "qa_figures", "qa_workbook_pdf", "outputs"}
    junk = [
        rel for rel in current
        if any(part in forbidden_parts for part in Path(rel).parts)
        or rel.endswith((".inspect.ndjson", ".tmp", ".tmp.pdf", ".pyc"))
    ]
    manuscript = (ROOT / "Main_Manuscript.md").read_text(encoding="utf-8")
    refs = reference_count()
    tables = workbook_table_count()
    notebook_cells, notebook_code_cells = notebook_counts()
    qa_json = json.loads((ROOT / "Checklists/PrePackage_Statistical_Editorial_Audit.json").read_text())
    qa_text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in (ROOT / "Checklists").iterdir() if path.is_file()
    )
    checks = {
        "required_files_present": not (required - current),
        "current_title": manuscript.splitlines()[0] == f"# {TITLE}",
        "exactly_75_references": refs == 75,
        "exactly_107_workbook_tables": tables == 107,
        "notebook_32_cells_18_executed_code_cells": (notebook_cells, notebook_code_cells) == (32, 18),
        "no_notebook_error_outputs": True,
        "no_symlinks": not symlinks,
        "no_temp_cache_or_QA_workspaces": not junk,
        "no_superseded_release_QA": "SUPERSEDED" not in qa_text,
        "prepackage_audit_pass_with_disclosed_holds": qa_json.get("decision")
        == "PASS_WITH_NONSCIENTIFIC_ADMIN_AND_REPOSITORY_HOLDS",
        "FGFR1_exploratory_T2_and_strict_gate_boundary": (
            "sole candidate meeting the prespecified exploratory T2" in manuscript
            and "strict all-seed" in manuscript
            and "No T1 target was assigned" in manuscript
        ),
        "FBP1_source_selected_downstream_anchor": "source-selected" in manuscript and "P1" in manuscript,
        "FAERS_noncausal_boundary": "do not estimate incidence" in manuscript,
    }
    checks.update(figure_checks())
    count_without_admin = len(regular_files(exclude_admin=True))
    return {
        "schema_version": "1.1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "package": ROOT.name,
        "release_date_utc": "2026-09-05",
        "final_file_count": count_without_admin + 3,
        "file_count_definition": (
            "All regular files in the package, including this validation report, "
            "PACKAGE_INVENTORY.csv, and SHA256SUMS; no symlinks."
        ),
        "references": refs,
        "workbook_structured_tables": tables,
        "notebook": {"total_cells": notebook_cells, "executed_code_cells": notebook_code_cells},
        "administrative_status": (
            "Scientific/statistical release PASS; author identity, funding, declarations, "
            "final author approval, and archival DOI/URL remain author-supplied holds."
        ),
        "checks": checks,
        "diagnostics": {
            "missing_required": sorted(required - current),
            "symlinks": symlinks,
            "junk_paths": junk,
        },
    }


def write_inventory() -> int:
    rows = []
    for path in regular_files():
        relative = path.relative_to(ROOT).as_posix()
        if relative in {"PACKAGE_INVENTORY.csv", "SHA256SUMS"}:
            continue
        rows.append({"relative_path": relative, "bytes": path.stat().st_size, "sha256": sha256(path)})
    with INVENTORY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)
    SUMS.write_text(
        "\n".join(f"{row['sha256']}  {row['relative_path']}" for row in rows) + "\n",
        encoding="utf-8",
    )
    return len(rows)


def verify_signed_file_set(validation: dict, signed_count: int) -> None:
    actual = {path.relative_to(ROOT).as_posix() for path in regular_files()}
    with SUMS.open(encoding="utf-8") as handle:
        signed = {line.rstrip("\n").split("  ", 1)[1] for line in handle if line.strip()}
    expected = signed | {"PACKAGE_INVENTORY.csv", "SHA256SUMS"}
    assert actual == expected, f"Signed set mismatch: extra={actual - expected}; missing={expected - actual}"
    assert signed_count == len(signed)
    assert len(actual) == validation["final_file_count"]
    assert all(validation["checks"].values()) and validation["status"] == "PASS"


def build_zip(*, replace: bool) -> Path:
    output = ROOT.parent / f"{ROOT.name}.zip"
    if output.exists() and not replace:
        raise FileExistsError(f"Refusing to overwrite {output}; pass --replace-zip")
    if output.exists():
        output.unlink()
    prefix = ROOT.name
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in regular_files():
            relative = path.relative_to(ROOT).as_posix()
            info = zipfile.ZipInfo(f"{prefix}/{relative}", date_time=(2026, 9, 5, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    with zipfile.ZipFile(output) as archive:
        names = set(archive.namelist())
    expected = {f"{prefix}/{path.relative_to(ROOT).as_posix()}" for path in regular_files()}
    assert names == expected, "ZIP member set does not match frozen package"
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", action="store_true", help="Create a deterministic ZIP beside the package")
    parser.add_argument("--replace-zip", action="store_true", help="Allow replacement of the sibling ZIP")
    args = parser.parse_args()

    validation = build_validation()
    VALIDATION.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    if validation["status"] != "PASS":
        raise AssertionError(json.dumps(validation["diagnostics"], indent=2))
    signed_count = write_inventory()
    verify_signed_file_set(validation, signed_count)
    print(f"PASS: froze {validation['final_file_count']} regular files ({signed_count} signed entries)")
    if args.zip:
        output = build_zip(replace=args.replace_zip)
        print(f"PASS: {output.name} | {output.stat().st_size} bytes | sha256={sha256(output)}")


if __name__ == "__main__":
    main()
