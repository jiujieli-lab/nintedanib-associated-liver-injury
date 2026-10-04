#!/usr/bin/env python3
"""Refresh release-facing hashes without rerunning upstream scientific models.

The analytical manifests preserve third-party/upstream provenance while adding
or rebasing hashes for the compact files actually distributed in this release.
No statistical table or figure content is modified by this script.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def require(relative: str) -> Path:
    path = ROOT / relative
    if not path.is_file():
        raise FileNotFoundError(relative)
    return path


def hash_map(paths: list[str]) -> dict[str, str]:
    return {relative: sha256(require(relative)) for relative in paths}


def refresh_integration() -> None:
    path = ROOT / "Source_Data/Integration/integration_analysis_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))

    available_inputs = {
        "fresh_analysis/proteomics/candidate_evidence_table.csv": "Source_Data/Human_DILI_Proteomics/candidate_evidence_table.csv",
        "fresh_analysis/liver_singlecell/target_expression_by_compartment.csv": "Source_Data/Liver_Single_Cell/target_expression_by_compartment.csv",
        "fresh_analysis/pharmacology/chembl_target_summary.csv": "Source_Data/Pharmacology/chembl_target_summary.csv",
        "fresh_analysis/virtual_ko/virtual_knockout_stability_summary.csv": "Source_Data/Virtual_Knockout/virtual_knockout_stability_summary.csv",
        "fresh_analysis/virtual_ko/virtual_knockout_seed_summary.csv": "Source_Data/Virtual_Knockout/virtual_knockout_seed_summary.csv",
        "fresh_analysis/virtual_ko/virtual_knockout_gene_results.csv.gz": "Source_Data/Virtual_Knockout/virtual_knockout_gene_results.csv.gz",
    }
    original_inputs = dict(manifest["input_hashes"])
    upstream_inputs = dict(manifest.get("upstream_input_hashes_not_redistributed", {}))
    manifest["input_hashes"] = hash_map(list(available_inputs.values()))
    upstream_inputs.update({
        old: digest
        for old, digest in original_inputs.items()
        if old not in available_inputs and old not in available_inputs.values()
    })
    manifest["upstream_input_hashes_not_redistributed"] = upstream_inputs

    output_paths = [
        f"Source_Data/Integration/{Path(old).name}"
        for old in manifest["output_table_hashes"]
    ]
    manifest["output_table_hashes"] = hash_map(output_paths)
    original_text_hashes = dict(manifest["text_artifact_hashes"])
    upstream_text_hashes = dict(manifest.get("upstream_text_artifact_hashes_not_redistributed", {}))
    text_paths = [
        "Figures/Figure_4_Liver_Localization_and_Virtual_Knockout_caption.md",
        "Figures/Figure_5_Multi_Evidence_Target_Prioritization_caption.md",
        "Source_Data/Integration/integration_methods_results.md",
    ]
    manifest["text_artifact_hashes"] = hash_map(text_paths)
    upstream_text_hashes.update({
        old: digest
        for old, digest in original_text_hashes.items()
        if Path(old).name == "input_audit.md"
    })
    manifest["upstream_text_artifact_hashes_not_redistributed"] = upstream_text_hashes

    figure_stems = {
        "figure4": "Figure_4_Liver_Localization_and_Virtual_Knockout",
        "figure5": "Figure_5_Multi_Evidence_Target_Prioritization",
    }
    manifest["figure_hashes"] = {
        key: {
            "svg": sha256(require(f"Figures/{stem}.svg")),
            "pdf": sha256(require(f"Figures/{stem}.pdf")),
            "png_600dpi": sha256(require(f"Figures/{stem}.png")),
            "tiff_600dpi": sha256(require(f"Figures/{stem}.tiff")),
        }
        for key, stem in figure_stems.items()
    }
    manifest["analysis_script"] = "Code/Core/build_virtual_integration_figures.py"
    manifest["analysis_script_sha256"] = sha256(require(manifest["analysis_script"]))
    manifest["validation_script"] = "Code/Core/validate_virtual_integration_outputs.py"
    manifest["validation_script_sha256"] = sha256(require(manifest["validation_script"]))
    manifest["release_manifest_refreshed_utc"] = "2026-09-05"
    manifest["release_scope"] = (
        "Hashes under input_hashes, output_table_hashes, text_artifact_hashes and figure_hashes "
        "refer to files distributed in this compact release. Upstream-only inputs retain their "
        "original recorded hashes under upstream_input_hashes_not_redistributed."
    )
    write_json(path, manifest)


def refresh_lincs() -> None:
    path = ROOT / "Source_Data/LINCS/analysis_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    source_map = {
        "sciplex_nintedanib_signature_long.csv": "Source_Data/Pharmacology/sciplex_nintedanib_signature_long.csv",
        "chembl_target_summary.csv": "Source_Data/Pharmacology/chembl_target_summary.csv",
        "chembl_human_single_protein_activities.csv": "Source_Data/Pharmacology/chembl_human_single_protein_activities.csv",
    }
    for entry in manifest["sources"]:
        basename = Path(entry["path"]).name
        if basename in source_map:
            relative = source_map[basename]
            entry.update(path=relative, sha256=sha256(require(relative)), available_in_release=True)
        else:
            entry["available_in_release"] = False

    for entry in manifest["outputs"]:
        old = entry["path"]
        basename = Path(old).name
        if old.startswith(("Figures/", "Supplementary_Figures/", "Source_Data/LINCS/")):
            relative = old
        elif "Figure_3_Nintedanib_HepG2_Dose_Perturbation" in basename:
            relative = f"Figures/{basename}"
        elif "Supplementary_Figure_S3_Nintedanib_Perturbation_Robustness" in basename:
            extension = Path(basename).suffix
            relative = (
                "Supplementary_Figures/"
                f"Supplementary_Figure_S4_Nintedanib_Perturbation_Robustness{extension}"
            )
        else:
            relative = f"Source_Data/LINCS/{basename}"
        entry.update(path=relative, sha256=sha256(require(relative)), available_in_release=True)

    manifest["release_analysis_script"] = "Code/Core/analyze_lincs_dose_and_figure3.py"
    manifest["release_analysis_script_sha256"] = sha256(require(manifest["release_analysis_script"]))
    manifest["release_manifest_refreshed_utc"] = "2026-09-05"
    write_json(path, manifest)


def preclinical_release_path(original: str) -> str | None:
    basename = Path(original).name
    if original.startswith("fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints"):
        return f"Figures/{basename}"
    if original == "fresh_analysis/preclinical/scripts/analyze_public_preclinical.py":
        return "Code/Public_Preclinical_Core/analyze_public_preclinical.py"
    if original.startswith("fresh_analysis/preclinical/results/"):
        return f"Source_Data/Public_Preclinical_Core/{basename}"
    if original.startswith("fresh_analysis/preclinical/") and "/raw/" not in original and "/metadata/" not in original:
        candidate = f"Source_Data/Public_Preclinical_Core/{basename}"
        return candidate if (ROOT / candidate).is_file() else None
    return None


def refresh_preclinical() -> None:
    path = ROOT / "Source_Data/Public_Preclinical_Core/analysis_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    for entry in manifest["files"]:
        release_path = preclinical_release_path(entry["path"])
        if release_path is None:
            entry["available_in_release"] = False
            entry.pop("release_path", None)
            continue
        file_path = require(release_path)
        entry.update(
            release_path=release_path,
            available_in_release=True,
            bytes=file_path.stat().st_size,
            sha256=sha256(file_path),
        )
    manifest["release_manifest_refreshed_utc"] = "2026-09-05"
    manifest["release_scope"] = (
        "release_path identifies the compact distributed counterpart when available; raw and "
        "metadata mirrors that are not redistributed retain their original provenance hashes."
    )
    write_json(path, manifest)


def main() -> None:
    refresh_integration()
    refresh_lincs()
    refresh_preclinical()
    print("PASS: release-facing Integration, LINCS and public-preclinical manifests refreshed")


if __name__ == "__main__":
    main()
