#!/usr/bin/env python3
"""Audit PMID, DOI, URL and DailyMed tokens in all currently visible manuscript text.

Layer 1 (reference-record integrity) is performed by audit_references.py. This
script performs Layer 2: every citation/source token in manuscript-facing text
must resolve either to a numbered reference or to a named official data-source
registry entry. Unknown tokens are hard failures.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFDIR = Path(__file__).resolve().parent
REFCSV = REFDIR / "Nintedanib_DILI_Verified_References.csv"

REQUIRED_PMIDS = {
    "38384288", "40387033", "40388329", "42092598", "42377397",
    "41866841", "41826148", "35714115", "36227799", "39892726",
    "39012714", "30406699",
}
REQUIRED_DAILYMED_SETID = "da1c9f37-779e-4682-816f-93d0faa4cfc9"

# Exact official endpoints intentionally present in Methods, figure legends, or
# the manuscript consistency audit. These are source locators, not fabricated
# journal references.
RESOURCE_URLS = {
    "https://api.fda.gov/drug/event.json": (
        "RESOURCE:OPENFDA_DRUG_EVENT_API",
        "US FDA openFDA drug adverse-event API endpoint",
    ),
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/": (
        "RESOURCE:GSE70138_FTP",
        "NCBI GEO official GSE70138 supplementary-file directory",
    ),
    "https://maayanlab.cloud/Harmonizome/": (
        "RESOURCE:HARMONIZOME",
        "Harmonizome official data portal",
    ),
    "https://string-db.org/api/": (
        "RESOURCE:STRING_API",
        "STRING official API documentation endpoint",
    ),
    "https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=CHEMBL502835&limit=1000": (
        "RESOURCE:CHEMBL_NINTEDANIB_ACTIVITY_API",
        "EMBL-EBI ChEMBL official activity API query for CHEMBL502835",
    ),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469": (
        "RESOURCE:GSE115469_GEO",
        "NCBI GEO accession GSE115469",
    ),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299128": (
        "RESOURCE:GSE299128_GEO",
        "NCBI GEO accession GSE299128",
    ),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138": (
        "RESOURCE:GSE70138_GEO",
        "NCBI GEO accession GSE70138",
    ),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120679": (
        "RESOURCE:GSE120679_GEO", "NCBI GEO accession GSE120679",
    ),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120804": (
        "RESOURCE:GSE120804_GEO", "NCBI GEO accession GSE120804",
    ),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151374": (
        "RESOURCE:GSE151374_GEO", "NCBI GEO accession GSE151374",
    ),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278200": (
        "RESOURCE:GSE278200_GEO", "NCBI GEO accession GSE278200",
    ),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308578": (
        "RESOURCE:GSE308578_GEO", "NCBI GEO accession GSE308578",
    ),
    "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA635636": (
        "RESOURCE:PRJNA635636", "NCBI BioProject PRJNA635636 linked to GSE151374",
    ),
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120679/suppl/": (
        "RESOURCE:GSE120679_FTP", "NCBI GEO official GSE120679 supplementary-file directory",
    ),
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278200/suppl/GSE278200_NINT_raw_counts.txt.gz": (
        "RESOURCE:GSE278200_COUNTS", "NCBI GEO official GSE278200 raw-count matrix",
    ),
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308578/suppl/GSE308578_Gene_counts.csv.gz": (
        "RESOURCE:GSE308578_COUNTS", "NCBI GEO official GSE308578 gene-count matrix",
    ),
    "https://www.ebi.ac.uk/pride/archive/projects/PXD024058": (
        "RESOURCE:PXD024058_PRIDE", "PRIDE project PXD024058",
    ),
    "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024058": (
        "RESOURCE:PXD024058_PROTEOMEXCHANGE", "ProteomeXchange record PXD024058",
    ),
    "https://www.ebi.ac.uk/pride/archive/projects/PXD052594": (
        "RESOURCE:PXD052594_PRIDE", "PRIDE project PXD052594",
    ),
    "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD052594": (
        "RESOURCE:PXD052594_PROTEOMEXCHANGE", "ProteomeXchange record PXD052594",
    ),
    "https://ddbj.nig.ac.jp/resource/bioproject/PRJDB12477": (
        "RESOURCE:PRJDB12477_DDBJ", "DDBJ BioProject PRJDB12477",
    ),
    "https://ddbj.nig.ac.jp/resource/sra-submission/DRA012991": (
        "RESOURCE:DRA012991_DDBJ", "DDBJ Sequence Read Archive submission DRA012991",
    ),
    "https://www.ebi.ac.uk/ena/browser/view/PRJDB12477": (
        "RESOURCE:PRJDB12477_ENA", "ENA mirror of BioProject PRJDB12477",
    ),
    "https://open.fda.gov/apis/drug/event/": (
        "RESOURCE:OPENFDA_DRUG_EVENT_DOCUMENTATION", "US FDA openFDA drug adverse-event API documentation",
    ),
    "https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL502835": (
        "RESOURCE:CHEMBL502835", "ChEMBL nintedanib compound record CHEMBL502835",
    ),
    "https://string-db.org/": (
        "RESOURCE:STRING_PORTAL", "STRING protein-association database",
    ),
    "https://www.python.org/": (
        "RESOURCE:PYTHON", "Python language official website",
    ),
    "https://github.com/cailab-tamu/scTenifoldNet": (
        "RESOURCE:SCTENIFOLDNET_CODE", "Official scTenifoldNet source-code repository",
    ),
    "https://numpy.org/": ("RESOURCE:NUMPY", "NumPy official project website"),
    "https://pandas.pydata.org/": ("RESOURCE:PANDAS", "pandas official project website"),
    "https://scipy.org/": ("RESOURCE:SCIPY", "SciPy official project website"),
    "https://scikit-learn.org/": ("RESOURCE:SCIKIT_LEARN", "scikit-learn official project website"),
    "https://www.statsmodels.org/": ("RESOURCE:STATSMODELS", "statsmodels official project website"),
    "https://matplotlib.org/": ("RESOURCE:MATPLOTLIB", "Matplotlib official project website"),
    "https://seaborn.pydata.org/": ("RESOURCE:SEABORN", "seaborn official project website"),
    "https://networkx.org/": ("RESOURCE:NETWORKX", "NetworkX official project website"),
    "https://openpyxl.readthedocs.io/": ("RESOURCE:OPENPYXL", "openpyxl official documentation"),
    "https://requests.readthedocs.io/": ("RESOURCE:REQUESTS", "Requests official documentation"),
    "https://www.h5py.org/": ("RESOURCE:H5PY", "h5py official project website"),
    "https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=9925b305-c1fb-43b6-a5fd-e03f14f28e6c": (
        "RESOURCE:DAILYMED_GENERIC_NINTEDANIB_COMPARATOR",
        "Official DailyMed generic nintedanib label used only as a documented comparison in the consistency audit",
    ),
}

# Persistent data-object DOIs are official source locators, not numbered journal
# references. This registry prevents DOI-shaped repository links in data-
# availability prose from being misclassified as literature citations.
RESOURCE_DOIS = {
    "10.5281/zenodo.11395642": (
        "RESOURCE:ZENODO_11395642",
        "Zenodo record 11395642 containing processed PXD052594 source files",
    ),
}

PMID_RE = re.compile(r"(?i)\bPMID\s*[:#]?\s*(\d{6,9})\b")
DOI_RE = re.compile(r"(?i)(?<![A-Za-z0-9])10\.\d{4,9}/[-._;()/:A-Z0-9]+")
URL_RE = re.compile(r"https?://[^\s<>\]}`\"'）。，；：！？】]+")
DAILYMED_RE = re.compile(r"(?i)DailyMed:setid-([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})")


def norm_doi(value):
    value = value.strip().lower()
    value = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value)
    value = value.rstrip(".,;:")
    while value.endswith(")") and value.count(")") > value.count("("):
        value = value[:-1]
    return value


def norm_url(value):
    value = value.strip()
    while value and value[-1] in ".,;:!?)]}>'\"。，；：！？）】":
        value = value[:-1]
    return value


def manuscript_files():
    candidates = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".md", ".txt"}:
            continue
        if REFDIR in path.parents:
            continue
        name = path.name.lower()
        manuscript_facing = (
            "draft" in name
            or "manuscript_ready" in name
            or "methods_results" in name
            or "caption" in name
            or "data_and_code_availability" in name
            or name == "experimental_validation_protocol.md"
            or name == "manuscript_consistency_audit.md"
        )
        if manuscript_facing:
            candidates.append(path)
    return sorted(set(candidates))


def main():
    with REFCSV.open(encoding="utf-8-sig", newline="") as fh:
        refs = list(csv.DictReader(fh))
    by_pmid = {r["pmid"]: r for r in refs if r["pmid"]}
    by_doi = {r["doi"].lower(): r for r in refs if r["doi"]}
    by_url = {r["official_url"]: r for r in refs if r["official_url"]}

    required_library_failures = []
    for pmid in sorted(REQUIRED_PMIDS):
        if pmid not in by_pmid:
            required_library_failures.append(f"Required PMID absent from reference library: {pmid}")
    branded_url = (
        "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid="
        + REQUIRED_DAILYMED_SETID
    )
    branded_ref = by_url.get(branded_url)
    if not branded_ref or branded_ref.get("reference_type") != "WEB":
        required_library_failures.append("Required DailyMed branded label absent or not typed WEB")

    rows = []
    files = manuscript_files()

    def add(path, lineno, token_type, raw, normalized, mapping_type, mapping_id, title, status):
        rows.append({
            "file": str(path.relative_to(ROOT)),
            "line": lineno,
            "token_type": token_type,
            "raw_token": raw,
            "normalized_token": normalized,
            "mapping_type": mapping_type,
            "mapping_id": mapping_id,
            "mapped_title_or_source": title,
            "status": status,
        })

    for path in files:
        for lineno, line in enumerate(path.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
            for match in PMID_RE.finditer(line):
                pmid = match.group(1)
                ref = by_pmid.get(pmid)
                add(
                    path, lineno, "PMID", match.group(0), pmid,
                    "numbered_reference" if ref else "unmapped",
                    f"REF:{ref['number']}" if ref else "",
                    ref["title"] if ref else "",
                    "MAPPED" if ref else "MISSING",
                )
            for match in DOI_RE.finditer(line):
                doi = norm_doi(match.group(0))
                ref = by_doi.get(doi)
                resource = RESOURCE_DOIS.get(doi)
                if ref:
                    add(path, lineno, "DOI", match.group(0), doi, "numbered_reference", f"REF:{ref['number']}", ref["title"], "MAPPED")
                elif resource:
                    rid, title = resource
                    add(path, lineno, "DOI", match.group(0), doi, "official_resource", rid, title, "MAPPED")
                else:
                    add(path, lineno, "DOI", match.group(0), doi, "unmapped", "", "", "MISSING")
            for match in DAILYMED_RE.finditer(line):
                setid = match.group(1).lower()
                url = "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=" + setid
                ref = by_url.get(url)
                add(
                    path, lineno, "DAILYMED_SETID", match.group(0), setid,
                    "numbered_web_reference" if ref else "unmapped",
                    f"REF:{ref['number']}" if ref else "",
                    ref["title"] if ref else "",
                    "MAPPED" if ref else "MISSING",
                )
            for match in URL_RE.finditer(line):
                raw = match.group(0)
                url = norm_url(raw)
                ref = by_url.get(url)
                if not ref and url.lower().startswith("https://doi.org/"):
                    ref = by_doi.get(norm_doi(url))
                resource = None
                if not ref and url.lower().startswith("https://doi.org/"):
                    resource = RESOURCE_DOIS.get(norm_doi(url))
                if not ref:
                    pm = re.match(r"https?://pubmed\.ncbi\.nlm\.nih\.gov/(\d{6,9})/?$", url, flags=re.I)
                    if pm:
                        ref = by_pmid.get(pm.group(1))
                if ref:
                    add(path, lineno, "URL", raw, url, "numbered_reference", f"REF:{ref['number']}", ref["title"], "MAPPED")
                elif resource:
                    rid, title = resource
                    add(path, lineno, "URL", raw, url, "official_resource", rid, title, "MAPPED")
                elif url in RESOURCE_URLS:
                    rid, title = RESOURCE_URLS[url]
                    add(path, lineno, "URL", raw, url, "official_resource", rid, title, "MAPPED")
                else:
                    add(path, lineno, "URL", raw, url, "unmapped", "", "", "MISSING")

    # Preserve occurrences for auditability while removing exact duplicates that
    # can arise only from duplicate regex passes on the same line/token/type.
    unique_rows = []
    seen = set()
    for row in rows:
        key = (row["file"], row["line"], row["token_type"], row["raw_token"])
        if key not in seen:
            seen.add(key)
            unique_rows.append(row)
    rows = unique_rows
    missing = [r for r in rows if r["status"] == "MISSING"]

    with (REFDIR / "Nintedanib_DILI_Manuscript_Token_Map.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        fields = list(rows[0]) if rows else [
            "file", "line", "token_type", "raw_token", "normalized_token",
            "mapping_type", "mapping_id", "mapped_title_or_source", "status",
        ]
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(r["token_type"] for r in rows)
    mapped_counts = Counter(r["token_type"] for r in rows if r["status"] == "MAPPED")
    md = [
        "# Manuscript citation and source-token coverage audit",
        "",
        "## Layer-2 result",
        "",
        f"- Currently visible manuscript-facing text files scanned: **{len(files)}**",
        f"- Token occurrences audited: **{len(rows)}**",
        f"- Mapped token occurrences: **{len(rows) - len(missing)}**",
        f"- Missing token occurrences: **{len(missing)}**",
        f"- Required added PMID records present in library: **{len(REQUIRED_PMIDS) - sum(x.startswith('Required PMID') for x in required_library_failures)}/{len(REQUIRED_PMIDS)}**",
        f"- Required DailyMed branded label present and typed as WEB: **{'YES' if not any('DailyMed' in x for x in required_library_failures) else 'NO'}**",
        f"- Required-library failures: **{len(required_library_failures)}**",
        "",
        "A token is mapped only when it resolves to a numbered JOUR/WEB reference or to an explicitly named official resource endpoint. Unknown domains or unregistered endpoints fail the audit; they are not silently accepted.",
        "",
        "## Token counts",
        "",
        "| Type | Audited | Mapped | Missing |",
        "|---|---:|---:|---:|",
    ]
    for token_type in sorted(counts):
        md.append(f"| {token_type} | {counts[token_type]} | {mapped_counts[token_type]} | {counts[token_type] - mapped_counts[token_type]} |")
    md += ["", "## Files scanned", ""]
    md.extend(f"- `{path.relative_to(ROOT)}`" for path in files)
    md += ["", "## Required newly added records", "", "| Token | Mapping | Title/source |", "|---|---|---|"]
    for pmid in sorted(REQUIRED_PMIDS):
        ref = by_pmid.get(pmid)
        md.append(f"| PMID:{pmid} | {'REF:' + ref['number'] if ref else 'MISSING'} | {ref['title'] if ref else ''} |")
    md.append(
        f"| DailyMed:{REQUIRED_DAILYMED_SETID} | "
        f"{'REF:' + branded_ref['number'] if branded_ref else 'MISSING'} | "
        f"{branded_ref['title'] if branded_ref else ''} |"
    )
    md += ["", "## Official resource URL registry", "", "| URL | Mapping | Source |", "|---|---|---|"]
    for url, (rid, title) in sorted(RESOURCE_URLS.items()):
        md.append(f"| {url} | {rid} | {title} |")
    md += ["", "## Official data-object DOI registry", "", "| DOI | Mapping | Source |", "|---|---|---|"]
    for doi, (rid, title) in sorted(RESOURCE_DOIS.items()):
        md.append(f"| {doi} | {rid} | {title} |")
    md += ["", "## Missing tokens", ""]
    if missing:
        md += ["| File | Line | Type | Token |", "|---|---:|---|---|"]
        for row in missing:
            md.append(f"| {row['file']} | {row['line']} | {row['token_type']} | {row['raw_token']} |")
    else:
        md.append("None. All currently visible manuscript PMID/DOI/URL/DailyMed tokens are mapped.")
    if required_library_failures:
        md += ["", "## Required-library failures", ""] + [f"- {x}" for x in required_library_failures]
    md += [
        "",
        "## Scope rule for later drafts",
        "",
        "The script discovers manuscript-facing Markdown/text files dynamically by filename. It must be rerun after manuscript assembly or citation replacement; any new unregistered token will produce a non-zero exit code and a visible MISSING row.",
        "",
    ]
    (REFDIR / "Nintedanib_DILI_Manuscript_Token_Coverage_Audit.md").write_text("\n".join(md), encoding="utf-8")
    failures = {"required_library_failures": required_library_failures, "missing_tokens": missing}
    (REFDIR / "manuscript_token_coverage_failures.json").write_text(json.dumps(failures, indent=2), encoding="utf-8")
    print(json.dumps({
        "files": len(files), "token_occurrences": len(rows),
        "mapped": len(rows) - len(missing), "missing": len(missing),
        "required_library_failures": len(required_library_failures),
    }, indent=2))
    if missing or required_library_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
