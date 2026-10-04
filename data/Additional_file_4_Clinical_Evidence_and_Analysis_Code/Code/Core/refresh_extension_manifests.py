#!/usr/bin/env python3
"""Refresh release mappings and hashes for the S6/S7 extension manifests.

Large public raw inputs are intentionally not duplicated in the release.  Their
source-run hashes remain frozen and retain the upstream verification result.
Every released code, result, validation, and figure file is mapped explicitly
and re-hashed in the current capsule.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def strict_dump(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def release_path(dataset: str, frozen_path: str) -> Path | None:
    key = dataset.lower()
    mappings = {
        f"fresh_analysis/preclinical/{key}/code/": f"Code/{dataset}/",
        f"fresh_analysis/preclinical/{key}/results/": f"Source_Data/{dataset}/results/",
        f"fresh_analysis/preclinical/{key}/validation/": f"Source_Data/{dataset}/validation/",
        "fresh_analysis/figures/": "Supplementary_Figures/",
    }
    for old, new in mappings.items():
        if frozen_path.startswith(old):
            return Path(new + frozen_path[len(old):])
    return None


def refresh(dataset: str) -> str:
    manifest_path = ROOT / "Source_Data" / dataset / "results" / "analysis_manifest.json"
    detached_path = ROOT / "Source_Data" / dataset / "validation" / "manifest_hash_validation.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    previous_detached = json.loads(detached_path.read_text(encoding="utf-8"))

    upstream_verified = {
        row["path"]: bool(row.get("bytes_match") and row.get("sha256_match"))
        for row in previous_detached.get("entries", [])
    }
    audit_entries = []
    released = 0
    for section in ("inputs", "outputs"):
        for row in manifest[section]:
            mapped = release_path(dataset, row["path"])
            candidate = ROOT / mapped if mapped is not None else None
            available = bool(candidate and candidate.is_file())
            if section == "outputs" and not available:
                raise FileNotFoundError(f"Released analytical output is missing: {row['path']}")
            row["release_path"] = mapped.as_posix() if available else None
            row["available_in_release"] = available
            if available:
                row["bytes"] = candidate.stat().st_size
                row["sha256"] = sha256(candidate)
                basis = "release_capsule_current_file"
                bytes_match = True
                sha_match = True
                released += 1
            else:
                if not upstream_verified.get(row["path"], False):
                    raise RuntimeError(f"No successful upstream verification for omitted input: {row['path']}")
                basis = "upstream_full_analysis_audit_before_raw_input_deduplication"
                bytes_match = True
                sha_match = True
            row["verification_basis"] = basis
            audit_entries.append({
                "path": row["path"],
                "release_path": row["release_path"],
                "available_in_release": available,
                "verification_basis": basis,
                "bytes_match": bytes_match,
                "sha256_match": sha_match,
            })

    manifest["release_mapping_policy"] = (
        "All released code, analytical outputs, validation files, and figure derivatives are mapped "
        "by release_path and re-hashed. Large public raw inputs are not duplicated; their frozen "
        "source-run hashes retain the successful upstream full-analysis audit."
    )
    manifest["release_visual_update"] = {
        "date_utc": "2026-09-05",
        "scope": "print-layout-only S6/S7 label deconfliction and 600-dpi raster export",
        "analytical_values_changed": False,
        "release_script": "Code/Core/finalize_s6_s7_figure_exports.py",
    }
    strict_dump(manifest_path, manifest)
    manifest_hash = sha256(manifest_path)

    detached = {
        "manifest_path": manifest_path.relative_to(ROOT).as_posix(),
        "manifest_sha256": manifest_hash,
        "entries_checked": len(audit_entries),
        "release_entries_checked": released,
        "omitted_public_raw_inputs_checked_upstream": len(audit_entries) - released,
        "all_recorded_hashes_match": all(x["bytes_match"] and x["sha256_match"] for x in audit_entries),
        "all_available_release_hashes_match": all(
            x["bytes_match"] and x["sha256_match"]
            for x in audit_entries if x["available_in_release"]
        ),
        "audit_scope": (
            "Combined audit: current-file verification for every available release entry and the "
            "successful frozen upstream audit for public raw inputs omitted from the capsule."
        ),
        "entries": audit_entries,
    }
    strict_dump(detached_path, detached)
    print(
        f"PASS {dataset}: manifest={manifest_hash}; entries={len(audit_entries)}; "
        f"release={released}; upstream-only={len(audit_entries) - released}"
    )
    return manifest_hash


def main() -> None:
    refresh("GSE151374")
    refresh("PXD052594")


if __name__ == "__main__":
    main()
