#!/usr/bin/env python3
"""Refresh release-facing FAERS hashes after an audited, frozen rerender.

The scientific payload is not recomputed here. This utility only reconciles
the existing FAERS manifest with the packaged code, tables, figures, captions,
panel source data, and final figure QA report.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "Source_Data/FAERS/faers_frozen_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def figure_path(name: str) -> Path:
    if name == "FAERS_figure_validation_report.json":
        return ROOT / "Checklists" / name
    if name.endswith("_panel_source_data.csv"):
        return ROOT / "Source_Data/Figure_Panels" / name
    if name.startswith("Supplementary_"):
        return ROOT / "Supplementary_Figures" / name
    return ROOT / "Figures" / name


def strict_json(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: strict_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [strict_json(item) for item in value]
    return value


def main() -> None:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    groups = [
        (payload["script_sha256"], lambda key: ROOT / key),
        (payload["result_sha256"], lambda key: ROOT / "Source_Data/FAERS" / key),
        (payload["figure_sha256"], figure_path),
    ]
    refreshed = 0
    for mapping, resolver in groups:
        for key in list(mapping):
            path = resolver(key)
            if not path.is_file():
                raise FileNotFoundError(f"Manifest target is absent: {path}")
            mapping[key] = sha256(path)
            refreshed += 1
    payload["release_manifest_refreshed_utc"] = "2026-09-05"
    MANIFEST.write_text(
        json.dumps(strict_json(payload), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"PASS: refreshed {refreshed} FAERS manifest hashes")


if __name__ == "__main__":
    main()
