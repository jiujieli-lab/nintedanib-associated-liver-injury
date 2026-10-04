#!/usr/bin/env python3
"""Build a clean-slate, verified reference library for the nintedanib-DILI study."""

from __future__ import annotations

import csv
import html
import json
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


OUTDIR = Path(__file__).resolve().parent


def s(identifier, category, section, why, reference_type="JOUR", metadata=None):
    record = {
        "identifier": identifier,
        "category": category,
        "section": section,
        "why": why,
        "reference_type": reference_type,
    }
    if metadata:
        record["metadata"] = metadata
    return record


SEEDS = [
    # ILD/IPF/PPF definitions, burden, pivotal trials and current unmet need
    s("PMID:35486072", "ILD guidelines and trials", "Introduction; Methods—cohort definition", "Current ATS/ERS/JRS/ALAT definitions of IPF and progressive pulmonary fibrosis."),
    s("PMID:30168753", "ILD guidelines and trials", "Introduction; Methods—cohort definition", "Diagnostic framework for IPF."),
    s("PMID:21471066", "ILD guidelines and trials", "Introduction", "Foundational international IPF guideline."),
    s("PMID:24836310", "ILD guidelines and trials", "Introduction; Discussion", "Pivotal INPULSIS evidence for nintedanib efficacy and adverse-event context."),
    s("PMID:31566307", "ILD guidelines and trials", "Introduction; Discussion", "Pivotal INBUILD evidence across progressive fibrosing ILDs."),
    s("PMID:31112379", "ILD guidelines and trials", "Introduction; Discussion", "Pivotal SENSCIS evidence in systemic-sclerosis ILD."),
    s("PMID:21992121", "ILD guidelines and trials", "Introduction", "Early TOMORROW nintedanib trial."),
    s("PMID:24836312", "ILD guidelines and trials", "Introduction; Discussion", "Active-comparator pirfenidone efficacy evidence."),
    s("PMID:21571362", "ILD guidelines and trials", "Introduction", "CAPACITY pirfenidone trials and antifibrotic treatment context."),
    s("PMID:30224318", "ILD guidelines and trials", "Introduction; Discussion", "Long-term INPULSIS-ON safety and tolerability."),
    s("PMID:34475231", "ILD guidelines and trials", "Introduction; Discussion", "Whole-trial INBUILD follow-up."),
    s("PMID:35199968", "ILD guidelines and trials", "Introduction; Discussion", "INBUILD autoimmune-disease subgroup."),
    s("PMID:34514739", "ILD guidelines and trials", "Introduction; Discussion", "SENSCIS subgroup heterogeneity."),
    s("DOI:10.1183/13993003.00085-2020", "ILD guidelines and trials", "Introduction", "Natural history and residual progression risk in fibrosing ILD."),
    s("DOI:10.1016/S2213-2600(21)00503-8", "ILD proteomics", "Introduction; Discussion", "Proteomic biomarkers of progressive fibrosing ILD."),
    s("DOI:10.1038/s41598-025-12952-1", "ILD proteomics", "Introduction; Discussion", "Longitudinal circulating proteins during nintedanib treatment."),
    s("PMID:42253602", "Recent unmet need", "Introduction; Discussion", "Recent synthesis of the evolving pulmonary-fibrosis therapeutic pipeline."),
    s("DOI:10.1056/NEJMoa2512911", "Recent unmet need", "Introduction; Discussion", "Recent phase 3 evidence illustrating ongoing need for better-tolerated/effective IPF therapies."),
    s("DOI:10.1038/s41467-026-75291-3", "Recent unmet need", "Introduction; Discussion", "Recent early-phase inhaled therapy illustrating active search for safer mechanisms."),

    # Nintedanib pharmacology, PK, safety and hepatotoxicity
    s("PMID:31016670", "Nintedanib pharmacology and safety", "Introduction; Methods—drug definition; Discussion", "Clinical PK/PD synthesis for nintedanib."),
    s("PMID:29106740", "Nintedanib pharmacology and safety", "Introduction; Discussion", "Impact of hepatic impairment on nintedanib exposure."),
    s("PMID:36927840", "Nintedanib pharmacology and safety", "Introduction; Discussion—mechanism", "Human enzymes involved in nintedanib metabolism."),
    s("PMID:26400368", "Nintedanib pharmacology and safety", "Introduction; Discussion", "Management and clinical interpretation of nintedanib adverse events."),
    s("PMID:29794129", "Nintedanib pharmacology and safety", "Introduction; Discussion", "Real-world/expanded-access safety and tolerability."),
    s("PMID:32767182", "Nintedanib pharmacology and safety", "Introduction; Discussion", "Global pharmacovigilance experience."),
    s("PMID:25474320", "Nintedanib pharmacology and safety", "Introduction; Discussion—mechanism", "Drug discovery, kinase target spectrum and translation."),
    s("DOI:10.1158/0008-5472.CAN-07-6307", "Nintedanib pharmacology and safety", "Introduction; Discussion—mechanism", "Foundational triple angiokinase pharmacology."),
    s("PMID:35195327", "Nintedanib DILI", "Introduction; Discussion", "Published clinical nintedanib-induced liver injury."),
    s("DOI:10.3390/ph15050645", "Nintedanib DILI", "Introduction; Discussion", "Integrated pharmacovigilance–pharmacokinetic appraisal of nintedanib liver injury."),
    s("PMID:40356972", "Nintedanib DILI", "Introduction; Discussion", "Recent FAERS/JADER comparison of pirfenidone and nintedanib safety."),
    s("PMID:38860173", "Nintedanib DILI", "Introduction; Discussion", "FAERS disproportionality comparison of antifibrotics."),
    s("PMID:38464722", "Nintedanib DILI", "Introduction; Discussion", "FAERS and VigiAccess antifibrotic safety study."),
    s("PMID:41675318", "Nintedanib DILI", "Introduction; Discussion", "Recent hospital cohort plus FAERS analysis defining the unresolved molecular gap."),
    s("PMID:39068856", "Nintedanib DILI", "Discussion—mechanism", "Reactive iminium metabolite formation by CYP3A4."),
    s("PMID:39466587", "Nintedanib pharmacology and safety", "Introduction; Discussion", "Post-marketing safety in fibrosing ILDs."),
    s("PMID:39714546", "Nintedanib pharmacology and safety", "Introduction; Discussion", "Japanese post-marketing surveillance."),

    # DILI definitions, causality, natural history and biomarkers/proteomics
    s("PMID:30926241", "DILI guidance and biomarkers", "Methods—outcome definition; Discussion", "EASL DILI definitions, diagnostic work-up and causality principles."),
    s("PMID:35899384", "DILI guidance and biomarkers", "Methods—outcome definition; Discussion", "AASLD practice guidance for drug-induced liver injury."),
    s("DOI:10.1038/s41467-023-36858-6", "DILI proteomics", "Introduction; Methods—proteomics; Results; Discussion", "Human longitudinal serum-proteomic DILI discovery and validation resource analyzed in this study."),
    s("PMID:18955056", "DILI guidance and biomarkers", "Introduction; Discussion", "Prospective DILIN clinical phenotypes and outcomes."),
    s("PMID:25754159", "DILI guidance and biomarkers", "Introduction; Discussion", "Features and outcomes of 899 patients in the US DILIN prospective cohort."),
    s("PMID:24681128", "DILI guidance and biomarkers", "Introduction; Discussion", "Substantial short-term morbidity and mortality after idiosyncratic DILI."),
    s("PMID:23419359", "DILI guidance and biomarkers", "Introduction; Discussion", "Population-based DILI incidence, presentation and outcomes."),
    s("PMID:20512999", "DILI guidance and biomarkers", "Methods—causality; Discussion", "Structured expert causality assessment versus RUCAM in DILIN."),
    s("PMID:23390034", "DILI guidance and biomarkers", "Introduction; Discussion—mechanism", "Mechanistic circulating miR-122, HMGB1 and keratin-18 biomarkers."),
    s("PMID:28341748", "DILI guidance and biomarkers", "Introduction; Discussion", "Advances in DILI diagnosis and risk assessment."),
    s("DOI:10.1038/s41572-019-0105-0", "DILI guidance and biomarkers", "Introduction; Discussion", "Comprehensive DILI disease-primer framework."),
    s("DOI:10.1038/s41573-019-0048-x", "DILI guidance and biomarkers", "Introduction; Discussion—mechanism", "Mechanisms and experimental test systems for DILI."),
    s("PMID:41233320", "DILI toxicogenomics", "Introduction; Discussion", "Recent large-scale human toxicogenomics resource for DILI prediction."),

    # Pharmacovigilance and READUS-PV
    s("PMID:11828828", "Pharmacovigilance methods", "Methods—disproportionality", "Proportional reporting ratio method."),
    s("PMID:9696956", "Pharmacovigilance methods", "Methods—disproportionality", "Bayesian confidence propagation neural-network method."),
    s("PMID:11998548", "Pharmacovigilance methods", "Methods—disproportionality", "Comparison of disproportionality measures."),
    s("PMID:23571771", "Pharmacovigilance methods", "Methods—disproportionality; Discussion", "Benchmarking signal-detection algorithms in FAERS."),
    s("PMID:27193236", "Pharmacovigilance methods", "Methods—FAERS preprocessing", "Curated standardized FAERS resource and normalization principles."),
    s("PMID:38713346", "Pharmacovigilance reporting", "Methods; Reporting checklist", "READUS-PV development and statement."),
    s("PMID:38713347", "Pharmacovigilance reporting", "Methods; Reporting checklist", "READUS-PV explanation and elaboration."),
    s("DOI:10.1002/pds.1001", "Pharmacovigilance methods", "Methods—disproportionality", "Reporting odds ratio and signal-detection framework."),
    s("PMID:21705438", "Pharmacovigilance methods", "Methods—Bayesian shrinkage", "Shrinkage observed-to-expected ratios for sparse safety data."),

    # Human liver atlases and open data resources
    s("PMID:30348985", "Human liver atlas", "Methods—single-cell atlas; Results; Discussion", "GSE115469 human-liver single-cell atlas used for cell-type localization."),
    s("PMID:31292543", "Human liver atlas", "Methods—single-cell atlas; Discussion", "Independent human liver cell atlas."),
    s("PMID:31597160", "Human liver atlas", "Methods—single-cell atlas; Discussion—mechanism", "Single-cell dissection of human liver fibrosis and cellular niches."),
    s("PMID:23193258", "Public data resources", "Methods—data availability", "NCBI GEO database update."),
    s("PMID:29195078", "Perturbation resources", "Methods—observed perturbation", "LINCS L1000 Connectivity Map resource."),
    s("PMID:17008526", "Perturbation resources", "Introduction; Methods—observed perturbation", "Foundational Connectivity Map concept."),
    s("PMID:31806696", "Perturbation resources", "Methods—observed perturbation", "Sci-Plex single-cell chemical perturbation platform."),
    s("PMID:27374120", "Public data resources", "Methods—observed perturbation", "Harmonizome integrative portal used to access standardized signatures."),
    s("DOI:10.1093/nar/gkad1004", "Drug–target resources", "Methods—target pharmacology", "Current ChEMBL database architecture and bioactivity content."),
    s("PMID:36370105", "Drug–target resources", "Methods—network analysis", "STRING protein-interaction database and network evidence."),
    s("PMID:29140462", "Perturbation resources", "Methods—observed perturbation", "LINCS data portal for integrated access to perturbation-response data."),
    s("DOI:10.1093/nar/gky1106", "Proteomics resources", "Methods—data availability", "PRIDE proteomics repository infrastructure."),

    # Virtual perturbation and single-cell machine learning
    s("DOI:10.1016/j.patter.2022.100434", "Virtual perturbation", "Methods—virtual knockout; Results; Discussion", "scTenifoldKnk virtual-knockout algorithm."),
    s("DOI:10.1016/j.patter.2020.100139", "Virtual perturbation", "Methods—gene-regulatory network", "scTenifoldNet single-cell network reconstruction and comparison."),
    s("DOI:10.1038/s41592-025-02909-7", "Virtual perturbation", "Methods—perturbation analysis", "Pertpy unified perturbation-analysis framework."),
    s("DOI:10.1038/s41592-025-02980-0", "Virtual perturbation", "Methods—benchmarking; Discussion", "Current systematic benchmark of single-cell perturbation prediction."),
    s("DOI:10.1038/s41592-019-0494-8", "Virtual perturbation", "Methods—sensitivity analysis; Discussion", "scGen generative prediction of single-cell perturbation responses."),
    s("DOI:10.15252/msb.202211517", "Virtual perturbation", "Methods—sensitivity analysis; Discussion", "Compositional perturbation autoencoder for drug-response prediction."),
    s("DOI:10.1038/s41587-023-01905-6", "Virtual perturbation", "Methods—sensitivity analysis; Discussion", "GEARS multigene perturbation prediction."),
    s("DOI:10.1038/s41592-018-0229-2", "Single-cell machine learning", "Methods—latent representation", "scVI probabilistic single-cell representation framework."),
    s("DOI:10.1038/s41592-023-01969-x", "Virtual perturbation", "Methods—sensitivity analysis", "Neural optimal-transport modeling of single-cell perturbation responses."),

    # Statistical robustness, interpretable ML and reporting
    s("DOI:10.1111/j.2517-6161.1995.tb02031.x", "Statistical methods", "Methods—multiplicity", "Benjamini–Hochberg false-discovery-rate control."),
    s("DOI:10.1111/j.1467-9868.2010.00740.x", "Statistical methods", "Methods—stability selection", "Stability selection for reproducible feature selection."),
    s("DOI:10.1038/s41562-020-0912-z", "Statistical methods", "Methods—specification curve; Results—robustness", "Specification-curve analysis across defensible analytic choices."),
    s("DOI:10.1214/aos/1176344552", "Statistical methods", "Methods—bootstrap uncertainty", "Foundational nonparametric bootstrap resampling framework."),
    s("DOI:10.1080/01621459.1997.10474007", "Statistical methods", "Methods—bootstrap validation", ".632+ bootstrap correction for prediction-error estimation."),
    s("PMID:16504092", "Statistical methods", "Methods—machine learning validation", "Nested cross-validation to limit model-selection bias."),
    s("PMID:25569120", "Prediction reporting", "Methods; Reporting checklist", "TRIPOD reporting guideline."),
    s("PMID:30596875", "Prediction reporting", "Methods; Reporting checklist", "PROBAST risk-of-bias and applicability framework."),
    s("PMID:32908275", "AI reporting", "Methods; Reporting checklist", "MI-CLAIM minimum information for clinical AI modeling."),
    s("DOI:10.1136/bmj-2023-078378", "AI reporting", "Methods; Reporting checklist", "TRIPOD+AI updated reporting guidance."),
    s("DOI:10.1038/s42256-019-0138-9", "Interpretable machine learning", "Methods—model interpretation", "Consistent individualized feature attribution for tree ensembles."),
    s("PMID:3203132", "Statistical methods", "Methods—ROC comparison", "DeLong nonparametric comparison of correlated ROC curves."),
    s("PMID:12883005", "Statistical methods", "Methods—multiplicity", "q-values and genome-wide false-discovery-rate interpretation."),
    s("PMID:17938396", "Observational reporting", "Methods; Reporting checklist", "STROBE reporting guidance."),
    s("PMID:26440803", "Observational reporting", "Methods; Reporting checklist", "RECORD guidance for routinely collected health data."),

    # Public preclinical datasets used for orthogonal antifibrotic validation
    s("PMID:30489156", "Preclinical public data", "Methods—GSE120679; Results—ex vivo validation; Discussion", "Precision-cut rat lung-slice transcriptomic model for evaluating antifibrotic drugs, including nintedanib."),
    s("DOI:10.1152/ajpgi.00281.2018", "Preclinical public data", "Methods—GSE120804; Results—hepatic-context sensitivity; Discussion", "Primary precision-cut rat liver-slice transcriptomic study linked to GSE120804."),
    s("PMID:33995084", "Preclinical public data", "Methods—PXD024058; Results—proteomic validation; Discussion", "Pulmonary-fibrosis proteomics dataset containing an orthogonal nintedanib comparator arm."),

    # Additional clinical/translational papers explicitly cited in the manuscript drafts
    s("PMID:38384288", "Nintedanib DILI", "Introduction; Discussion", "Active-comparator health-claims analysis of potential hepatotoxicity with pirfenidone or nintedanib."),
    s("PMID:40387033", "Recent unmet need", "Introduction; Discussion", "Phase 3 FIBRONEER-IPF trial of nerandomilast."),
    s("PMID:40388329", "Recent unmet need", "Introduction; Discussion", "Phase 3 FIBRONEER-ILD trial of nerandomilast in progressive pulmonary fibrosis."),
    s("PMID:42092598", "Translational delivery", "Discussion—dual-organ validation", "Recent inhalable liposomal nintedanib study designed to decouple pulmonary efficacy from systemic toxicity."),
    s("PMID:42377397", "Translational delivery", "Discussion—future directions", "Safety and pharmacokinetic evidence for inhaled nintedanib dry-powder administration."),
    s("PMID:41866841", "Preclinical public data", "Methods—GSE278200; Results—animal RNA-seq validation; Discussion", "Independent pulmonary-fibrosis animal RNA-seq dataset used for orthogonal validation."),
    s("PMID:41826148", "Preclinical public data", "Methods—GSE308578; Results—animal RNA-seq validation; Discussion", "Independent nintedanib-treated pulmonary-fibrosis animal RNA-seq dataset used for orthogonal validation."),
    s("PMID:35714115", "Preclinical public data", "Methods—public-data eligibility audit; Discussion—limitations", "DRA012991/PRJDB12477 nintedanib-treated iRA-ILD mouse-lung RNA-seq; retained as audit-only because multiple sequencing runs map to each animal and the source analysis used post hoc animal exclusions."),
    s("PMID:36227799", "Preclinical public data", "Methods—GSE151374; Results—pooled-library single-cell sensitivity; Discussion", "Primary report linked by GEO to GSE151374, used for the nintedanib-treated bleomycin mouse-lung pooled-library single-cell sensitivity analysis."),
    s("PMID:39892726", "Preclinical public data", "Discussion—macrophage mechanism context", "Secondary mechanistic report used to contextualize, but not define the design or sample mapping of, the accession-linked GSE151374 primary study."),
    s("PMID:39012714", "Preclinical public data", "Methods—PXD052594; Results—animal-level lung-proteome sensitivity; Discussion", "Primary report for PXD052594; its animal-level lung proteome and separate phosphoproteome subset define an orthogonal pulmonary-efficacy evidence boundary."),
    s(
        "URL:https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=da1c9f37-779e-4682-816f-93d0faa4cfc9",
        "Regulatory label",
        "Introduction; Methods—clinical safety context; Discussion",
        "Official branded OFEV US Prescribing Information supporting indications, DILI warning, hepatic monitoring and dose modification.",
        reference_type="WEB",
        metadata={
            "authors_list": ["Boehringer Ingelheim Pharmaceuticals, Inc."],
            "authors": "Boehringer Ingelheim Pharmaceuticals, Inc.",
            "title": "OFEV (nintedanib) capsules, for oral use: full prescribing information",
            "journal": "DailyMed",
            "year": "2025",
            "volume": "",
            "issue": "",
            "pages": "",
            "doi": "",
            "pmid": "",
            "publisher": "National Library of Medicine (US)",
            "revised_date": "2025-05",
            "accessed_date": "2026-09-04",
            "verification_status": "Verified against the official DailyMed Structured Product Label; set ID da1c9f37-779e-4682-816f-93d0faa4cfc9; revised May 2025; accessed 2026-09-04",
        },
    ),
]


# Crossref occasionally omits legacy pagination or supplies a non-Index-Medicus
# journal title. These narrow, manually verified overrides make the exported
# Vancouver strings submission-ready without changing the source identifiers.
METADATA_OVERRIDES = {
    "10.1214/aos/1176344552": {"journal": "Ann Stat", "pages": "1-26"},
    "10.1080/01621459.1997.10474007": {"journal": "J Am Stat Assoc"},
    "10.1111/j.2517-6161.1995.tb02031.x": {"journal": "J R Stat Soc Series B Stat Methodol"},
    "10.1111/j.1467-9868.2010.00740.x": {"journal": "J R Stat Soc Series B Stat Methodol"},
}


def get_json(url, attempts=4):
    headers = {"User-Agent": "nintedanib-dili-reference-builder/1.0 (academic use)"}
    err = None
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=50) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:  # network retry
            err = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Unable to retrieve {url}: {err}")


def norm_doi(value):
    if not value:
        return ""
    value = html.unescape(str(value)).strip().lower()
    value = re.sub(r"^https?://(dx\.)?doi\.org/", "", value)
    return value.rstrip(".")


def initials(given):
    parts = re.findall(r"[A-Za-zÀ-ÖØ-öø-ÿ]+", given or "")
    return "".join(part[0].upper() for part in parts if part)


def crossref_record(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    msg = get_json(url)["message"]
    authors = []
    for a in msg.get("author", []):
        name = (a.get("family", "") + " " + initials(a.get("given", ""))).strip()
        if name:
            authors.append(name)
    date = msg.get("published-print") or msg.get("published-online") or msg.get("issued") or {}
    year = str((date.get("date-parts") or [[""]])[0][0])
    return {
        "authors_list": authors,
        "authors": ", ".join(authors),
        "title": html.unescape((msg.get("title") or [""])[0]),
        "journal": html.unescape((msg.get("container-title") or [""])[0]),
        "year": year,
        "volume": str(msg.get("volume", "")),
        "issue": str(msg.get("issue", "")),
        "pages": str(msg.get("page") or msg.get("article-number") or ""),
        "doi": norm_doi(msg.get("DOI") or doi),
        "pmid": "",
        "verification_status": "Verified against Crossref DOI metadata",
    }


def europepmc_record(seed):
    ident = seed["identifier"]
    if ident.startswith("URL:"):
        rec = dict(seed.get("metadata") or {})
        if not rec:
            raise RuntimeError(f"Missing web-reference metadata for {ident}")
        return rec
    if ident.startswith("PMID:"):
        value = ident.split(":", 1)[1]
        query = f"EXT_ID:{value} AND SRC:MED"
    else:
        value = norm_doi(ident.split(":", 1)[1])
        query = f'DOI:"{value}"'
    params = urllib.parse.urlencode({"query": query, "format": "json", "resultType": "core", "pageSize": 5})
    data = get_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + params)
    results = data.get("resultList", {}).get("result", [])
    chosen = None
    if ident.startswith("PMID:"):
        chosen = next((x for x in results if str(x.get("pmid", "")) == value), None)
    else:
        chosen = next((x for x in results if norm_doi(x.get("doi")) == value), None)
    if not chosen:
        if ident.startswith("DOI:"):
            return crossref_record(value)
        raise RuntimeError(f"No Europe PMC match for {ident}")
    author_list = []
    for a in chosen.get("authorList", {}).get("author", []):
        name = a.get("collectiveName") or a.get("fullName") or ""
        if name:
            author_list.append(html.unescape(name))
    if not author_list and chosen.get("authorString"):
        author_list = [x.strip() for x in chosen["authorString"].split(",") if x.strip()]
    doi = norm_doi(chosen.get("doi"))
    journal_info = chosen.get("journalInfo") or {}
    journal_obj = journal_info.get("journal") or {}
    journal_name = (
        chosen.get("journalTitle")
        or journal_obj.get("medlineAbbreviation")
        or journal_obj.get("isoabbreviation")
        or journal_obj.get("title")
        or ""
    )
    status = "Verified against Europe PMC/PubMed"
    if doi:
        status += " (PMID and DOI metadata)"
    else:
        status += " (PMID metadata; DOI not indexed)"
    return {
        "authors_list": author_list,
        "authors": ", ".join(author_list),
        "title": html.unescape(re.sub(r"<[^>]+>", "", chosen.get("title", ""))),
        "journal": html.unescape(journal_name),
        "year": str(chosen.get("pubYear") or journal_info.get("yearOfPublication") or ""),
        "volume": str(chosen.get("journalVolume") or journal_info.get("volume") or ""),
        "issue": str(chosen.get("issue") or journal_info.get("issue") or ""),
        "pages": str(chosen.get("pageInfo", "")),
        "doi": doi,
        "pmid": str(chosen.get("pmid", "")),
        "verification_status": status,
    }


def vancouver(rec):
    if rec.get("reference_type") == "WEB":
        author_text = (rec["authors"] or rec["journal"]).rstrip(".")
        revised = rec.get("revised_date", "")
        accessed = rec.get("accessed_date", "")
        revised_text = ""
        if revised:
            year, month = (revised.split("-", 1) + [""])[:2]
            month_name = {
                "01": "Jan", "02": "Feb", "03": "Mar", "04": "Apr",
                "05": "May", "06": "Jun", "07": "Jul", "08": "Aug",
                "09": "Sep", "10": "Oct", "11": "Nov", "12": "Dec",
            }.get(month, month)
            revised_text = f"revised {year} {month_name}"
        accessed_text = ""
        if accessed:
            year, month, day = (accessed.split("-") + ["", ""])[:3]
            month_name = {
                "01": "Jan", "02": "Feb", "03": "Mar", "04": "Apr",
                "05": "May", "06": "Jun", "07": "Jul", "08": "Aug",
                "09": "Sep", "10": "Oct", "11": "Nov", "12": "Dec",
            }.get(month, month)
            accessed_text = f"cited {year} {month_name} {int(day) if day.isdigit() else day}"
        dates = "; ".join(x for x in (revised_text, accessed_text) if x)
        return (
            f"{author_text}. {rec['title']} [Internet]. Bethesda (MD): "
            f"{rec.get('publisher') or rec['journal']}; {rec['year']}"
            f" [{dates}]. Available from: {rec['official_url']}"
        )
    authors = rec["authors_list"]
    author_text = ", ".join(authors[:6])
    if len(authors) > 6:
        author_text += ", et al"
    parts = [f"{author_text}." if author_text else "", f"{rec['title'].rstrip('.') }.", rec["journal"] + "."]
    bibliographic = rec["year"]
    if rec["volume"]:
        bibliographic += ";" + rec["volume"]
        if rec["issue"]:
            bibliographic += f"({rec['issue']})"
        if rec["pages"]:
            bibliographic += ":" + rec["pages"]
    elif rec["pages"]:
        bibliographic += ":" + rec["pages"]
    bibliographic += "."
    parts.append(bibliographic)
    if rec["doi"]:
        parts.append("doi: " + rec["doi"] + ".")
    return " ".join(p for p in parts if p).replace("..", ".")


def ris_record(number, rec):
    lines = ["TY  - " + ("ELEC" if rec.get("reference_type") == "WEB" else "JOUR")]
    for author in rec["authors_list"]:
        lines.append("AU  - " + author)
    lines.append("TI  - " + rec["title"])
    if rec.get("reference_type") == "WEB":
        lines.append("T2  - " + rec["journal"])
        if rec.get("publisher"):
            lines.append("PB  - " + rec["publisher"])
    else:
        lines.append("JO  - " + rec["journal"])
    lines.append("PY  - " + rec["year"])
    if rec["volume"]:
        lines.append("VL  - " + rec["volume"])
    if rec["issue"]:
        lines.append("IS  - " + rec["issue"])
    if rec["pages"]:
        if "-" in rec["pages"]:
            start, end = rec["pages"].split("-", 1)
            lines += ["SP  - " + start, "EP  - " + end]
        else:
            lines.append("SP  - " + rec["pages"])
    if rec["doi"]:
        lines += ["DO  - " + rec["doi"], "UR  - https://doi.org/" + rec["doi"]]
    elif rec["pmid"]:
        lines.append("UR  - https://pubmed.ncbi.nlm.nih.gov/" + rec["pmid"] + "/")
    elif rec.get("official_url"):
        lines.append("UR  - " + rec["official_url"])
    if rec.get("revised_date"):
        lines.append("DA  - " + rec["revised_date"].replace("-", "/"))
    if rec.get("accessed_date"):
        lines.append("Y2  - " + rec["accessed_date"].replace("-", "/"))
    if rec["pmid"]:
        lines.append("AN  - PMID:" + rec["pmid"])
    lines += ["N1  - Reference number: " + str(number), "ER  - ", ""]
    return "\n".join(lines)


def main():
    records = [None] * len(SEEDS)
    failures = []
    with ThreadPoolExecutor(max_workers=12) as pool:
        future_map = {pool.submit(europepmc_record, seed): i for i, seed in enumerate(SEEDS)}
        for future in as_completed(future_map):
            i = future_map[future]
            try:
                rec = future.result()
                rec.update(SEEDS[i])
                records[i] = rec
            except Exception as exc:
                failures.append({"index": i + 1, "seed": SEEDS[i], "error": str(exc)})
    good = []
    seen = set()
    for rec in records:
        if not rec:
            continue
        key = rec["doi"] or ("pmid:" + rec["pmid"] if rec["pmid"] else rec["identifier"])
        if key in seen:
            failures.append({
                "stage": "unique-seed coverage",
                "seed": rec["identifier"],
                "error": f"Resolved identifier collided with an earlier canonical record: {key}",
            })
            continue
        seen.add(key)
        rec.update(METADATA_OVERRIDES.get(rec["doi"], {}))
        rec["number"] = len(good) + 1
        if rec["pmid"]:
            rec["official_url"] = "https://pubmed.ncbi.nlm.nih.gov/" + rec["pmid"] + "/"
        elif rec["doi"]:
            rec["official_url"] = "https://doi.org/" + rec["doi"]
        else:
            rec["official_url"] = rec["identifier"].split(":", 1)[1]
        rec["vancouver"] = vancouver(rec)
        good.append(rec)

    if len(good) != len(SEEDS) and not any(f.get("stage") == "unique-seed coverage" for f in failures):
        failures.append({
            "stage": "seed coverage",
            "error": f"Expected {len(SEEDS)} unique verified records but generated {len(good)}",
        })

    columns = [
        "number", "reference_type", "authors", "title", "journal", "year", "volume", "issue", "pages",
        "doi", "pmid", "official_url", "category", "section", "why",
        "publisher", "revised_date", "accessed_date", "verification_status", "vancouver", "identifier",
    ]
    with (OUTDIR / "Nintedanib_DILI_Verified_References.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(good)

    with (OUTDIR / "Nintedanib_DILI_Verified_References.ris").open("w", encoding="utf-8") as fh:
        fh.write("\n".join(ris_record(rec["number"], rec) for rec in good))

    md = [
        "# Verified references: nintedanib-associated drug-induced liver injury",
        "",
        f"Clean-slate library generated from {len(SEEDS)} prespecified records; {len(good)} unique records were resolved and verified.",
        "",
        "Verification hierarchy: Europe PMC/PubMed metadata first; Crossref DOI metadata when a PubMed record was unavailable. URLs link to PubMed when indexed, otherwise to the DOI landing page.",
        "",
        "## Numbered Vancouver reference list",
        "",
    ]
    for rec in good:
        md += [
            f"{rec['number']}. {rec['vancouver']}",
            f"   - Category: {rec['category']}",
            f"   - Suggested citation location: {rec['section']}",
            f"   - Use: {rec['why']}",
            f"   - Reference type: {rec['reference_type']}",
            f"   - Identifiers: DOI {rec['doi'] or 'not applicable'}; PMID {rec['pmid'] or 'not applicable'}",
            f"   - Verified: {rec['verification_status']}",
            f"   - Official record: {rec['official_url']}",
            "",
        ]
    md += ["## Unresolved prespecified records", ""]
    if failures:
        for f in failures:
            seed = f.get("seed", "")
            if isinstance(seed, dict):
                seed = seed.get("identifier", "")
            label = f"Seed {f.get('index', '?')} ({seed})" if seed else f.get("stage", "resolution")
            md.append(f"- {label}: {f['error']}")
    else:
        md.append("None.")
    md += [
        "",
        "## Citation-use guardrails",
        "",
        "- Cite FAERS disproportionality results as reporting associations or safety signals, not incidence, risk, prevalence, or proof of causation.",
        "- Cite virtual-knockout outputs as model-based network perturbation evidence, not an experimentally observed knockout phenotype.",
        "- Keep observed human proteomic effects, computational inference, and proposed wet-lab validation explicitly separated.",
        "- Recheck online-ahead-of-print pagination for 2026 articles immediately before submission.",
        "",
    ]
    (OUTDIR / "Nintedanib_DILI_Verified_References.md").write_text("\n".join(md), encoding="utf-8")
    (OUTDIR / "reference_resolution_failures.json").write_text(json.dumps(failures, indent=2), encoding="utf-8")
    print(json.dumps({"seed_count": len(SEEDS), "verified_unique": len(good), "failures": len(failures)}, indent=2))
    if failures or len(good) != len(SEEDS):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
