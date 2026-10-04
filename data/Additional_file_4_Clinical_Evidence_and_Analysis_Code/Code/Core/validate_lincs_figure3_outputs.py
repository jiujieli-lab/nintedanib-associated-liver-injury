#!/usr/bin/env python3
"""Independent structural and numerical QA for the Figure 3 analysis layer."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Source_Data" / "LINCS"
FIG = ROOT / "Figures"
SUPP_FIG = ROOT / "Supplementary_Figures"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check(condition: bool, message: str, checks: list[dict]) -> None:
    checks.append({"check": message, "passed": bool(condition)})
    if not condition:
        raise AssertionError(message)


def main() -> None:
    checks: list[dict] = []
    manifest = json.loads((OUT / "analysis_manifest.json").read_text())
    for source in manifest["sources"]:
        if not source.get("available_in_release", True):
            continue
        path = ROOT / source["path"]
        check(path.exists(), f"source exists: {source['path']}", checks)
        check(sha256(path) == source["sha256"], f"source hash matches: {source['path']}", checks)

    metrics = pd.read_csv(OUT / "gene_dose_response_metrics.csv.gz")
    check(len(metrics) == 12_328, "12,328 unique L1000 genes", checks)
    check(metrics["gene_symbol"].nunique() == 12_328, "gene symbols are unique", checks)
    check(int(metrics["is_landmark"].sum()) == 978, "978 measured landmark genes", checks)
    check(metrics["bh_q_all_genes"].between(0, 1).all(), "all genome-wide q values valid", checks)
    check(int(metrics["bh_q_all_genes"].lt(0.05).sum()) == 0, "no genome-wide exact-trend FDR hit", checks)
    check(np.isclose(metrics["spearman_exact_p"].min(), 2 / 720), "minimum exact P equals 2/720", checks)
    check(metrics["direction_consistency"].between(0, 1).all(), "direction consistency bounded", checks)
    check(metrics["loo_sign_stability"].between(0, 1).all(), "LOO stability bounded", checks)

    anchors = pd.read_csv(OUT / "dili_anchor_metrics.csv")
    check(len(anchors) == 13, "13 prespecified DILI anchors audited", checks)
    check(int(anchors["present_in_l1000"].sum()) == 12, "12 DILI anchors represented", checks)
    missing = anchors.loc[~anchors["present_in_l1000"], "gene_symbol"].tolist()
    check(missing == ["DMGDH"], "DMGDH is the only absent DILI anchor", checks)

    modules = pd.read_csv(OUT / "mechanism_module_statistics.csv")
    check(len(modules) == 7, "seven prespecified mechanism modules", checks)
    check(modules["bh_q_magnitude_family"].between(0, 1).all(), "module q values valid", checks)

    sci = pd.read_csv(OUT / "sciplex_hepg2_concordance.csv")
    check(len(sci) == 8, "four Sci-Plex dose pairs in two analysis subsets", checks)
    check(set(sci["subset"]) == {"all_overlap", "measured_landmark_only"}, "Sci-Plex sensitivity subsets complete", checks)
    check(sci["sign_concordance"].between(0, 1).all(), "Sci-Plex concordance bounded", checks)

    cross = pd.read_csv(OUT / "chembl_transcript_crosswalk.csv")
    high = cross.loc[cross["multi_document_potent"].eq(True) & cross["present_in_l1000"]]
    check(len(high) == 16, "16 multi-document potent targets represented in L1000", checks)
    check(int(high["is_canonical_nintedanib_target"].sum()) == 8, "all eight canonical target families represented", checks)

    main_stem = FIG / "Figure_3_Nintedanib_HepG2_Dose_Perturbation"
    supp_stem = SUPP_FIG / "Supplementary_Figure_S4_Nintedanib_Perturbation_Robustness"
    for stem in [main_stem, supp_stem]:
        for suffix in [".svg", ".pdf", ".png", ".tiff"]:
            check(stem.with_suffix(suffix).exists(), f"figure exists: {stem.name}{suffix}", checks)
        ET.parse(stem.with_suffix(".svg"))
        check(True, f"SVG is well-formed: {stem.name}", checks)
        for suffix in [".png", ".tiff"]:
            with Image.open(stem.with_suffix(suffix)) as image:
                width, height = image.size
                dpi = image.info.get("dpi", (0, 0))
                image.verify()
            check(width >= 3500 and height >= 5000, f"high-resolution raster dimensions: {stem.name}{suffix}", checks)
            check(min(dpi) >= 590, f"approximately 600 dpi: {stem.name}{suffix}", checks)

    caption = (OUT / "Figure_3_caption.md").read_text()
    svg_text = main_stem.with_suffix(".svg").read_text(encoding="utf-8")
    for panel in "ABCDEFGHIJ":
        present = re.search(rf">\s*{panel}\s*<", svg_text) is not None
        check(present, f"main figure contains panel {panel}", checks)
    check("no panel demonstrates nintedanib-specific clinical DILI" in caption,
          "caption retains clinical-causality boundary", checks)
    methods = (OUT / "manuscript_ready_lincs_methods_results.md").read_text()
    check("## Methods" in methods and "## Results" in methods, "manuscript-ready Methods and Results present", checks)
    check("do not demonstrate nintedanib-specific clinical DILI" in methods,
          "manuscript text retains interpretation boundary", checks)

    report = {
        "status": "PASS",
        "n_checks": len(checks),
        "n_passed": sum(x["passed"] for x in checks),
        "checks": checks,
        "main_figure_panels": 10,
        "main_figure_formats": ["SVG", "PDF", "PNG 600 dpi", "TIFF 600 dpi"],
    }
    (OUT / "validation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": report["n_checks"]}, indent=2))


if __name__ == "__main__":
    main()
