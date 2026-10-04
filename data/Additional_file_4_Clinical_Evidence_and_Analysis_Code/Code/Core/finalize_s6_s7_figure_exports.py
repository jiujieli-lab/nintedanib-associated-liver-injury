#!/usr/bin/env python3
"""Apply print-layout-only fixes to frozen S6/S7 SVGs and export 600-dpi files.

The scientific geometry, axes, points, heatmaps, thresholds, and numerical
labels are unchanged.  This release step removes only nominal top-P gene names
that do not represent a multiplicity-controlled discovery, repositions the S6D
title below its two-row legend, and refreshes PDF/PNG/TIFF derivatives.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path

from lxml import etree
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "Supplementary_Figures"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def text_value(group: etree._Element) -> str:
    return "".join(group.itertext()).strip()


def remove_annotation_groups(tree: etree._ElementTree, labels: set[str], *, arrows: bool) -> int:
    removed = 0
    for group in list(tree.xpath("//*[local-name()='g' and starts-with(@id, 'text_')]")):
        if text_value(group) not in labels:
            continue
        parent = group.getparent()
        if parent is None:
            raise RuntimeError("Annotation group has no parent")
        if arrows:
            previous = group.getprevious()
            if previous is None or not previous.get("id", "").startswith("patch_"):
                raise RuntimeError(f"Expected arrow immediately before {text_value(group)!r}")
            arrow_markup = etree.tostring(previous, encoding="unicode")
            if "#777777" not in arrow_markup:
                raise RuntimeError(f"Unexpected annotation-arrow style for {text_value(group)!r}")
            parent.remove(previous)
        parent.remove(group)
        removed += 1
    return removed


def write_svg(tree: etree._ElementTree, path: Path) -> None:
    tree.write(
        str(path),
        encoding="utf-8",
        xml_declaration=True,
        doctype=tree.docinfo.doctype or None,
        pretty_print=False,
    )


def edit_s6(path: Path) -> None:
    tree = etree.parse(str(path))
    old_title = tree.xpath("//*[local-name()='text' and text()='Broad-lineage composition at pool level']")
    if len(old_title) != 1:
        raise RuntimeError(f"Expected one S6D title, found {len(old_title)}")
    old_title[0].text = "Broad-lineage composition"
    old_title[0].set("y", "246.000000")
    old_title[0].set("transform", "rotate(-0 267.195512 246.000000)")

    labels = {
        "Mknk1", "Clint1", "RP23-116A10.3", "Map3k10", "2310009B15Rik",
        "Bik", "Galm", "Plcb4", "Gm28875", "Vps35",
    }
    removed = remove_annotation_groups(tree, labels, arrows=True)
    if removed != len(labels):
        raise RuntimeError(f"Expected to remove {len(labels)} S6 labels, removed {removed}")
    write_svg(tree, path)


def edit_s7(path: Path) -> None:
    tree = etree.parse(str(path))
    labels = {"Dennd11", "Kctd9", "Nle1", "Kif3b", "Elmo1", "Dnmbp"}
    # Dennd11 also occurs as a heatmap row label; restrict removal to the
    # 5-point annotation type used in panel C.
    removed = 0
    for group in list(tree.xpath("//*[local-name()='g' and starts-with(@id, 'text_')]")):
        label = text_value(group)
        if label not in labels:
            continue
        text_nodes = group.xpath(".//*[local-name()='text']")
        if len(text_nodes) != 1 or "font-size: 5px" not in text_nodes[0].get("style", ""):
            continue
        parent = group.getparent()
        if parent is None:
            raise RuntimeError("Annotation group has no parent")
        parent.remove(group)
        removed += 1
    if removed != len(labels):
        raise RuntimeError(f"Expected to remove {len(labels)} S7 labels, removed {removed}")
    write_svg(tree, path)


def export_derivatives(svg: Path) -> None:
    inkscape = shutil.which("inkscape")
    if not inkscape:
        raise RuntimeError("Inkscape is required for deterministic SVG export")
    pdf = svg.with_suffix(".pdf")
    png = svg.with_suffix(".png")
    tiff = svg.with_suffix(".tiff")
    subprocess.run(
        [inkscape, str(svg), "--export-type=pdf", f"--export-filename={pdf}"],
        check=True,
    )
    subprocess.run(
        [inkscape, str(svg), "--export-type=png", "--export-dpi=600", f"--export-filename={png}"],
        check=True,
    )
    # Pillow's large-strip LZW writer produced a header-readable but
    # pixel-undecodable S7 file in one environment.  ImageMagick writes tiled,
    # standards-compliant TIFF data and is independently decoded below.
    convert = shutil.which("convert")
    if not convert:
        raise RuntimeError("ImageMagick convert is required for robust TIFF export")
    subprocess.run(
        [convert, str(png), "-units", "PixelsPerInch", "-density", "600", "-compress", "LZW", str(tiff)],
        check=True,
    )

    with Image.open(png) as image:
        image.load()
        dpi = image.info.get("dpi", (0, 0))
        if min(dpi) < 599 or image.width < 8000:
            raise RuntimeError(f"Unexpected PNG export metadata for {png.name}: {image.size}, {dpi}")
    with Image.open(tiff) as image:
        image.load()
        dpi = image.info.get("dpi", (0, 0))
        if min(dpi) < 599 or image.width < 8000 or image.info.get("compression") != "tiff_lzw":
            raise RuntimeError(f"Unexpected TIFF export metadata for {tiff.name}: {image.size}, {dpi}, {image.info}")

    print(
        f"PASS {svg.name}: "
        + ", ".join(f"{p.suffix[1:]}={p.stat().st_size}B/{sha256(p)[:12]}" for p in (svg, pdf, png, tiff))
    )


def main() -> None:
    s6 = FIGURES / "Supplementary_Figure_S6_GSE151374.svg"
    s7 = FIGURES / "Supplementary_Figure_S7_PXD052594.svg"
    edit_s6(s6)
    edit_s7(s7)
    export_derivatives(s6)
    export_derivatives(s7)


if __name__ == "__main__":
    main()
