#!/usr/bin/env python3
"""Structural and preset audit for the submission DOCX."""

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path

from lxml import etree


ROOT = Path(__file__).resolve().parents[2]
DOCX = ROOT / "Main_Manuscript.docx"
REPORT = ROOT / "Checklists" / "Main_Manuscript_DOCX_Audit.json"
REPORT.parent.mkdir(parents=True, exist_ok=True)
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def attr(node, name):
    return node.get(f"{{{NS['w']}}}{name}") if node is not None else None


def main() -> None:
    checks = []

    def check(name, passed, detail):
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    with zipfile.ZipFile(DOCX) as zf:
        document = etree.fromstring(zf.read("word/document.xml"))
        styles = etree.fromstring(zf.read("word/styles.xml"))
        numbering = etree.fromstring(zf.read("word/numbering.xml"))
        footer_names = [x for x in zf.namelist() if x.startswith("word/footer") and x.endswith(".xml")]
        footers = [etree.fromstring(zf.read(x)) for x in footer_names]

    text = "\n".join(document.xpath(".//w:t/text()", namespaces=NS))
    check("docx_exists_and_nontrivial", DOCX.exists() and DOCX.stat().st_size > 100_000, DOCX.stat().st_size)
    check("no_unresolved_citation_tokens", not re.search(r"\[(?:PMID|DOI|DailyMed):", text), "none")
    check("no_raw_markdown_emphasis", "**" not in text and "`" not in text, "none")
    raw_headings = re.findall(r"(?:^|\n)(#{2,6}\s+[^\n]+)", text)
    check("no_raw_markdown_headings", not raw_headings, raw_headings)
    ref_numbers = [int(x) for x in re.findall(r"(?m)^\s*(\d+)\.\s", text)]
    check("exactly_75_references", sorted(set(ref_numbers)) == list(range(1, 76)), len(set(ref_numbers)))
    residual_numeric_citations = re.findall(r"\[\d+(?:[,–—-]\d+)*\]", text)
    check("no_residual_bracketed_numeric_citations", not residual_numeric_citations, residual_numeric_citations)
    for phrase in [
        "sole candidate meeting the prespecified exploratory T2",
        "failed the strict all-seed gate",
        "No T1 target was assigned",
        "source-selected",
        "unsigned network displacement",
        "do not estimate incidence",
        "Figure 6",
        "Supplementary Figure S5",
    ]:
        check(f"required_claim::{phrase}", phrase.lower() in text.lower(), phrase)

    sect = document.find(".//w:sectPr", NS)
    pg_sz = sect.find("w:pgSz", NS)
    pg_mar = sect.find("w:pgMar", NS)
    geometry = {
        "width": attr(pg_sz, "w"), "height": attr(pg_sz, "h"),
        "top": attr(pg_mar, "top"), "right": attr(pg_mar, "right"),
        "bottom": attr(pg_mar, "bottom"), "left": attr(pg_mar, "left"),
        "header": attr(pg_mar, "header"), "footer": attr(pg_mar, "footer"),
    }
    check("letter_page_size", geometry["width"] == "12240" and geometry["height"] == "15840", geometry)
    check("one_inch_margins", all(geometry[x] == "1440" for x in ["top", "right", "bottom", "left"]), geometry)
    check("header_footer_distance", geometry["header"] in {"708", "709"} and geometry["footer"] in {"708", "709"}, geometry)

    expected_styles = {
        "Normal": {"size": "22", "after": "160", "line": "320"},
        "Heading1": {"size": "32", "before": "360", "after": "200"},
        "Heading2": {"size": "26", "before": "240", "after": "120"},
        "Heading3": {"size": "24", "before": "160", "after": "80"},
    }
    for sid, expected in expected_styles.items():
        node = styles.find(f".//w:style[@w:styleId='{sid}']", NS)
        rpr = node.find("w:rPr", NS) if node is not None else None
        ppr = node.find("w:pPr", NS) if node is not None else None
        size = attr(rpr.find("w:sz", NS), "val") if rpr is not None else None
        spacing = ppr.find("w:spacing", NS) if ppr is not None else None
        observed = {"size": size, "before": attr(spacing, "before"), "after": attr(spacing, "after"), "line": attr(spacing, "line")}
        passed = observed["size"] == expected["size"] and all(observed.get(k) == v for k, v in expected.items() if k != "size")
        check(f"preset_style::{sid}", passed, observed)

    tables = document.findall(".//w:tbl", NS)
    table_details = []
    table_ok = True
    for idx, table in enumerate(tables, 1):
        tbl_pr = table.find("w:tblPr", NS)
        width = attr(tbl_pr.find("w:tblW", NS), "w")
        indent = attr(tbl_pr.find("w:tblInd", NS), "w")
        layout = attr(tbl_pr.find("w:tblLayout", NS), "type")
        grid = [int(attr(x, "w")) for x in table.findall("w:tblGrid/w:gridCol", NS)]
        cell_widths = [int(attr(x, "w")) for x in table.findall("w:tr[1]/w:tc/w:tcPr/w:tcW", NS)]
        row_heights = table.findall(".//w:trHeight", NS)
        ok = width == "9360" and indent == "120" and layout == "fixed" and sum(grid) == 9360 and grid == cell_widths and not row_heights
        table_ok &= ok
        table_details.append({"table": idx, "width": width, "indent": indent, "layout": layout, "grid": grid, "first_row_cell_widths": cell_widths, "fixed_row_heights": len(row_heights), "passed": ok})
    check("all_tables_fixed_dxa_geometry", table_ok and bool(tables), table_details)

    num_paras = document.xpath(".//w:p[w:pPr/w:numPr]", namespaces=NS)
    abs_nums = numbering.findall(".//w:abstractNum", NS)
    check("real_word_numbering_used", len(num_paras) > 0 and len(abs_nums) > 0, {"numbered_paragraphs": len(num_paras), "abstract_numbering_definitions": len(abs_nums)})
    page_fields = sum(len(f.xpath(".//*[@w:instr='PAGE']", namespaces=NS)) for f in footers)
    check("page_number_field_present", page_fields >= 1, page_fields)
    line_number_nodes = document.findall(".//w:sectPr/w:lnNumType", NS)
    check("continuous_line_numbering_present", bool(line_number_nodes) and attr(line_number_nodes[0], "restart") == "continuous", len(line_number_nodes))
    superscript_runs = document.xpath(
        ".//w:r[w:rPr/w:vertAlign[@w:val='superscript']]",
        namespaces=NS,
    )
    check("numbered_citations_render_as_superscript", len(superscript_runs) >= 30, len(superscript_runs))

    headings = document.xpath(".//w:p[w:pPr/w:pStyle[starts-with(@w:val,'Heading')]]", namespaces=NS)
    heading_text = ["".join(x.xpath(".//w:t/text()", namespaces=NS)) for x in headings]
    required_sections = ["Summary", "Introduction", "Results", "Discussion", "Limitations of the study", "STAR Methods", "References", "Figure legends", "Supplemental information"]
    check("all_required_sections_present", all(x in heading_text for x in required_sections), heading_text)

    allowed_prefixes = ("[AUTHOR", "[AFFILIATIONS", "[FULL NAME", "[EMAIL", "[ARCHIVAL", "[AUTHOR-VERIFIED")
    placeholders = sorted(set(re.findall(r"\[[A-Z][A-Z0-9 –—/;(),.'-]+\]", text)))
    scientific_group_labels = {"[DO]", "[NDO]", "[HV]", "[DF]", "[BDL]"}
    unexpected_placeholders = [
        x for x in placeholders
        if not x.startswith(allowed_prefixes) and x not in scientific_group_labels
    ]
    check("only_administrative_placeholders_remain", not unexpected_placeholders, {
        "allowed": [x for x in placeholders if x.startswith(allowed_prefixes)],
        "scientific_group_labels": sorted(scientific_group_labels & set(placeholders)),
        "unexpected": unexpected_placeholders,
    })

    report = {
        "status": "PASS" if all(x["passed"] for x in checks) else "FAIL",
        "document": DOCX.name,
        "n_checks": len(checks),
        "n_failed": sum(not x["passed"] for x in checks),
        "checks": checks,
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
