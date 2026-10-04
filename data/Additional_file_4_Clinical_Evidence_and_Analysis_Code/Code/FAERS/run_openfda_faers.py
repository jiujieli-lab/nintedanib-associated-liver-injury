#!/usr/bin/env python3
"""Reproducible openFDA/FAERS analysis of nintedanib-associated liver injury.

This is a report-level pharmacovigilance analysis, not an incidence study.
All endpoints, aliases, strata, exclusions, and estimands are frozen below.

Data source: official openFDA drug/event endpoint
Retrieval date: 2026-09-04
API data release discovered at runtime (expected last_updated: 2026-07-30)

Run from the project root with:
    .venv/bin/python fresh_analysis/faers/run_openfda_faers.py

Optional:
    OPENFDA_API_KEY=... .venv/bin/python ... --refresh
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable

from audit_year_heterogeneity import compute_year_effect_heterogeneity


BASE_URL = "https://api.fda.gov/drug/event.json"
RETRIEVAL_DATE = "2026-09-04"
START_DATE = "20141015"  # same-day US approvals of OFEV and ESBRIET
EXPECTED_DATA_END = "20260630"
RANDOM_SEED = 20260904

# Official openFDA/ICH E2B mapping for patient.drug.drugcharacterization:
#   1 = Suspect, 2 = Concomitant, 3 = Interacting.
# Source schema: https://open.fda.gov/fields/drugevent.yaml
DRUG_CHARACTERIZATION_LABELS = {
    "1": "suspect",
    "2": "concomitant",
    "3": "interacting",
}

HERE = Path(__file__).resolve().parent
RAW_DIR = HERE / "raw_cache"
OUT_DIR = HERE / "results"
RAW_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)


# Prespecified custom endpoints. They are not asserted to be licensed MedDRA SMQs.
NARROW_PTS = (
    "DRUG-INDUCED LIVER INJURY",
    "HEPATOTOXICITY",
    "HEPATITIS TOXIC",
    "HEPATOCELLULAR INJURY",
    "CHOLESTATIC LIVER INJURY",
    "MIXED LIVER INJURY",
    "HEPATIC FAILURE",
    "ACUTE HEPATIC FAILURE",
    "HEPATIC NECROSIS",
    "FULMINANT HEPATITIS",
)

BROAD_ONLY_PTS = (
    "ALANINE AMINOTRANSFERASE INCREASED",
    "ASPARTATE AMINOTRANSFERASE INCREASED",
    "TRANSAMINASES INCREASED",
    "HEPATIC ENZYME INCREASED",
    "LIVER FUNCTION TEST ABNORMAL",
    "LIVER FUNCTION TEST INCREASED",
    "BLOOD BILIRUBIN INCREASED",
    "HYPERBILIRUBINAEMIA",
    "JAUNDICE",
    "CHOLESTASIS",
    "HEPATITIS",
    "HEPATITIS ACUTE",
    "HEPATIC FUNCTION ABNORMAL",
    "LIVER DISORDER",
    "GAMMA-GLUTAMYLTRANSFERASE INCREASED",
    "BLOOD ALKALINE PHOSPHATASE INCREASED",
)
BROAD_PTS = NARROW_PTS + BROAD_ONLY_PTS

AMINOTRANSFERASE_PTS = (
    "ALANINE AMINOTRANSFERASE INCREASED",
    "ASPARTATE AMINOTRANSFERASE INCREASED",
    "TRANSAMINASES INCREASED",
)
BILIRUBIN_JAUNDICE_PTS = (
    "BLOOD BILIRUBIN INCREASED",
    "HYPERBILIRUBINAEMIA",
    "JAUNDICE",
)

# Transparent, frozen indication vocabulary used only for a report-level
# specification analysis. openFDA search operates on report documents; it does
# not expose a relational guarantee that an indication term and role code occur
# on the same array element as the queried drug. This limitation is carried into
# the output label and manuscript caption rather than hidden.
ILD_INDICATION_TERMS = (
    "IDIOPATHIC PULMONARY FIBROSIS",
    "PULMONARY FIBROSIS",
    "INTERSTITIAL LUNG DISEASE",
    "IDIOPATHIC INTERSTITIAL PNEUMONIA",
    "PROGRESSIVE FIBROSING INTERSTITIAL LUNG DISEASE",
    "FIBROSING INTERSTITIAL LUNG DISEASE",
    "SYSTEMIC SCLEROSIS ASSOCIATED INTERSTITIAL LUNG DISEASE",
    "SYSTEMIC SCLEROSIS PULMONARY",
)

GENERIC_QUERY = {
    "nintedanib": 'patient.drug.openfda.generic_name.exact:"NINTEDANIB"',
    "pirfenidone": 'patient.drug.openfda.generic_name.exact:"PIRFENIDONE"',
}

EXPANDED_QUERY = {
    "nintedanib": "(" + " OR ".join((
        'patient.drug.openfda.generic_name.exact:"NINTEDANIB"',
        'patient.drug.medicinalproduct.exact:"OFEV"',
        'patient.drug.medicinalproduct.exact:"NINTEDANIB"',
        'patient.drug.medicinalproduct.exact:"NINTEDANIB ESYLATE"',
    )) + ")",
    "pirfenidone": "(" + " OR ".join((
        'patient.drug.openfda.generic_name.exact:"PIRFENIDONE"',
        'patient.drug.medicinalproduct.exact:"ESBRIET"',
        'patient.drug.medicinalproduct.exact:"PIRFENIDONE"',
    )) + ")",
}

DRUG_ALIASES = {
    "nintedanib": {"NINTEDANIB", "NINTEDANIB ESYLATE", "OFEV"},
    "pirfenidone": {"PIRFENIDONE", "ESBRIET"},
}

SEX_STRATA = {
    "male": 'patient.patientsex:"1"',
    "female": 'patient.patientsex:"2"',
}

# Age analyses deliberately use only reports whose age unit is years (ICH code 801).
AGE_STRATA = {
    "18-64_years": '(patient.patientonsetageunit:"801" AND patient.patientonsetage:[18 TO 64])',
    "65plus_years": '(patient.patientonsetageunit:"801" AND patient.patientonsetage:[65 TO 150])',
}

SERIOUS_FIELDS = {
    "serious": "serious",
    "death": "seriousnessdeath",
    "hospitalization": "seriousnesshospitalization",
    "life_threatening": "seriousnesslifethreatening",
    "disabling": "seriousnessdisabling",
    "congenital_anomaly": "seriousnesscongenitalanomali",
    "other_serious": "seriousnessother",
}


@dataclass(frozen=True)
class QuerySpec:
    label: str
    params: dict[str, Any]


class OpenFDAClient:
    """Cached, rate-conscious client with a complete audit log."""

    def __init__(self, refresh: bool = False, max_workers: int = 3):
        self.refresh = refresh
        self.max_workers = max_workers
        self.api_key = os.environ.get("OPENFDA_API_KEY")
        self.log: list[dict[str, Any]] = []

    @staticmethod
    def _cache_path(params: dict[str, Any]) -> Path:
        payload = json.dumps(params, sort_keys=True, separators=(",", ":"))
        return RAW_DIR / f"{hashlib.sha256(payload.encode()).hexdigest()}.json"

    def get(self, spec: QuerySpec) -> dict[str, Any]:
        params = dict(spec.params)
        if self.api_key:
            params["api_key"] = self.api_key
        cache_params = {k: v for k, v in params.items() if k != "api_key"}
        cache_path = self._cache_path(cache_params)
        started = datetime.now(timezone.utc).isoformat()
        from_cache = cache_path.exists() and not self.refresh
        status = 200
        elapsed = 0.0
        url = BASE_URL + "?" + urlencode(cache_params)

        if from_cache:
            raw = cache_path.read_bytes()
            data = json.loads(raw)
        else:
            last_error: Exception | None = None
            for attempt in range(7):
                t0 = time.monotonic()
                try:
                    response = requests.get(BASE_URL, params=params, timeout=240)
                    elapsed += time.monotonic() - t0
                    status = response.status_code
                    url = response.url.replace(f"api_key={self.api_key}", "api_key=REDACTED") if self.api_key else response.url
                    if response.status_code == 200:
                        raw = response.content
                        data = response.json()
                        cache_path.write_bytes(raw)
                        break
                    if response.status_code == 404:
                        # openFDA uses 404 for a valid query with zero matching records.
                        raw = response.content
                        data = {"meta": {"results": {"total": 0}}, "results": []}
                        cache_path.write_bytes(json.dumps(data, indent=2).encode())
                        break
                    if response.status_code in {429, 500, 502, 503, 504}:
                        time.sleep(min(60, 3 * (2**attempt)) + random.random())
                        continue
                    response.raise_for_status()
                except Exception as exc:  # network retries are recorded, not hidden
                    last_error = exc
                    if attempt == 6:
                        raise
                    time.sleep(min(60, 3 * (2**attempt)) + random.random())
            else:
                raise RuntimeError(f"openFDA query failed: {spec.label}: {last_error}")

        meta = data.get("meta", {})
        total = meta.get("results", {}).get("total")
        sha = hashlib.sha256(cache_path.read_bytes()).hexdigest()
        self.log.append({
            "label": spec.label,
            "retrieval_date": RETRIEVAL_DATE,
            "request_started_utc": started,
            "url": url,
            "status_code": status,
            "from_cache": from_cache,
            "elapsed_seconds": round(elapsed, 3),
            "meta_last_updated": meta.get("last_updated"),
            "meta_total": total,
            "cache_file": cache_path.name,
            "response_sha256": sha,
        })
        return data

    def many(self, specs: Iterable[QuerySpec]) -> dict[str, dict[str, Any]]:
        specs = list(specs)
        out: dict[str, dict[str, Any]] = {}
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(self.get, spec): spec for spec in specs}
            for future in as_completed(futures):
                spec = futures[future]
                out[spec.label] = future.result()
        return out

    def save_log(self) -> None:
        pd.DataFrame(self.log).sort_values(["label", "url"]).to_csv(
            OUT_DIR / "api_query_log.csv", index=False
        )


def q_or_terms(pts: Iterable[str]) -> str:
    return "(" + " OR ".join(
        f'patient.reaction.reactionmeddrapt.exact:"{pt}"' for pt in pts
    ) + ")"


def q_or_indications(indications: Iterable[str]) -> str:
    # drugindication is free text and its `.exact` keyword field is case-sensitive;
    # the analyzed phrase field is required to recover capitalization variants
    # that are visible in the downloaded report JSON.
    return "(" + " OR ".join(
        f'patient.drug.drugindication:"{term}"' for term in indications
    ) + ")"


def q_or_indications_exact(indications: Iterable[str]) -> str:
    """Case-sensitive keyword variant retained only as a field-sensitivity audit."""
    return "(" + " OR ".join(
        f'patient.drug.drugindication.exact:"{term}"' for term in indications
    ) + ")"


def exclusive_drug_query(drug: str, aliases: dict[str, str] = GENERIC_QUERY) -> str:
    other = "pirfenidone" if drug == "nintedanib" else "nintedanib"
    return f"({aliases[drug]}) AND NOT ({aliases[other]})"


def in_window(query: str, start: str, end: str) -> str:
    return f"({query}) AND receivedate:[{start} TO {end}]"


def total_from(data: dict[str, Any]) -> int:
    return int(data.get("meta", {}).get("results", {}).get("total", 0))


def fetch_all(client: OpenFDAClient, label: str, search: str, expected: int) -> list[dict[str, Any]]:
    """Fetch every report in deterministic skip pages and verify completeness."""
    page_size = 100  # no-key-safe maximum
    specs = [
        QuerySpec(f"{label}__page_{skip:05d}", {"search": search, "limit": page_size, "skip": skip})
        for skip in range(0, expected, page_size)
    ]
    pages = client.many(specs)
    records: list[dict[str, Any]] = []
    for spec in specs:
        records.extend(pages[spec.label].get("results", []))

    # Retain the highest safetyreportversion for a duplicated report ID, if encountered.
    by_id: dict[str, dict[str, Any]] = {}
    for record in records:
        rid = str(record.get("safetyreportid", ""))
        version = int(record.get("safetyreportversion", 0) or 0)
        current = by_id.get(rid)
        current_version = int(current.get("safetyreportversion", 0) or 0) if current else -1
        if current is None or version > current_version:
            by_id[rid] = record

    if len(records) != expected:
        raise RuntimeError(f"Pagination mismatch for {label}: fetched={len(records)}, expected={expected}")
    if len(by_id) != expected:
        raise RuntimeError(
            f"Report-ID duplication detected for {label}: unique={len(by_id)}, expected={expected}. "
            "Stop rather than silently changing the denominator."
        )
    return [by_id[k] for k in sorted(by_id, key=lambda x: int(x) if x.isdigit() else x)]


def _normalized_names(drug_obj: dict[str, Any]) -> set[str]:
    names: set[str] = set()
    for key in ("medicinalproduct",):
        value = drug_obj.get(key)
        if value:
            names.add(str(value).strip().upper())
    active = drug_obj.get("activesubstance", {}).get("activesubstancename")
    if active:
        names.add(str(active).strip().upper())
    ofda = drug_obj.get("openfda") or {}
    for key in ("generic_name", "brand_name", "substance_name"):
        for value in ofda.get(key, []) or []:
            names.add(str(value).strip().upper())
    return names


def _is_matching_drug(drug_obj: dict[str, Any], drug: str) -> bool:
    aliases = DRUG_ALIASES[drug]
    names = _normalized_names(drug_obj)
    return any(any(alias == name or alias in name for alias in aliases) for name in names)


def _parse_date(value: Any) -> pd.Timestamp | None:
    if value is None:
        return None
    text = str(value).strip()
    if not re.fullmatch(r"\d{8}", text):
        return None
    parsed = pd.to_datetime(text, format="%Y%m%d", errors="coerce")
    return None if pd.isna(parsed) else parsed


def _event_date(record: dict[str, Any]) -> pd.Timestamp | None:
    text = ((record.get("patient") or {}).get("summary") or {}).get("narrativeincludeclinical", "")
    match = re.search(r"CASE EVENT DATE:\s*(\d{8})(?!\d)", str(text))
    return _parse_date(match.group(1)) if match else None


def _age_years(patient: dict[str, Any]) -> float | None:
    try:
        age = float(patient.get("patientonsetage"))
    except (TypeError, ValueError):
        return None
    factor = {
        "800": 10.0,       # decade
        "801": 1.0,        # year
        "802": 1 / 12,     # month
        "803": 7 / 365.25, # week
        "804": 1 / 365.25, # day
        "805": 1 / (365.25 * 24),
    }.get(str(patient.get("patientonsetageunit")))
    if factor is None:
        return None
    value = age * factor
    return value if 0 <= value <= 120 else None


def tidy_report(record: dict[str, Any], drug: str) -> dict[str, Any]:
    patient = record.get("patient") or {}
    reactions = patient.get("reaction") or []
    pts = sorted({str(r.get("reactionmeddrapt", "")).upper() for r in reactions if r.get("reactionmeddrapt")})
    hepatic_pts = sorted(set(pts).intersection(BROAD_PTS))
    matched_drugs = [d for d in patient.get("drug", []) or [] if _is_matching_drug(d, drug)]
    roles = [str(d.get("drugcharacterization", "")) for d in matched_drugs]
    role_rank = {"1": 1, "2": 2, "3": 3}
    best_role_code = min(roles, key=lambda x: role_rank.get(x, 9)) if roles else ""
    role_label = DRUG_CHARACTERIZATION_LABELS.get(best_role_code, "unknown")

    starts: list[pd.Timestamp] = []
    for d in matched_drugs:
        if str(d.get("drugstartdateformat", "")) == "102":
            date = _parse_date(d.get("drugstartdate"))
            if date is not None:
                starts.append(date)
    start = min(starts) if starts else None
    event = _event_date(record)
    tto = (event - start).days if start is not None and event is not None else None

    indications = sorted({str(d.get("drugindication", "")).strip().upper() for d in matched_drugs if d.get("drugindication")})
    ild_keywords = ("PULMONARY FIBROSIS", "INTERSTITIAL LUNG", "FIBROSING INTERSTITIAL", "SYSTEMIC SCLEROSIS")
    ild_indication = any(any(keyword in indication for keyword in ild_keywords) for indication in indications)
    sex = {"1": "male", "2": "female"}.get(str(patient.get("patientsex", "")), "unknown")
    age_unit_code = str(patient.get("patientonsetageunit", ""))
    try:
        age_value_reported = float(patient.get("patientonsetage"))
    except (TypeError, ValueError):
        age_value_reported = np.nan

    result: dict[str, Any] = {
        "drug": drug,
        "safetyreportid": str(record.get("safetyreportid", "")),
        "safetyreportversion": int(record.get("safetyreportversion", 0) or 0),
        "receivedate": record.get("receivedate"),
        "receiptdate": record.get("receiptdate"),
        "report_year": int(str(record.get("receivedate"))[:4]),
        "primarysourcecountry": record.get("primarysourcecountry"),
        "occurcountry": record.get("occurcountry"),
        "reporter_qualification": (record.get("primarysource") or {}).get("qualification"),
        "sex": sex,
        "age_years": _age_years(patient),
        # Retain the unconverted source fields so the years-only age analysis can
        # enforce exactly the same ICH-unit and range restrictions in its numerator
        # as in the API-derived denominator.
        "age_value_reported": age_value_reported,
        "age_unit_code": age_unit_code,
        "all_pts": "|".join(pts),
        "hepatic_pts": "|".join(hepatic_pts),
        "narrow_case": bool(set(pts).intersection(NARROW_PTS)),
        "broad_case": bool(set(pts).intersection(BROAD_PTS)),
        "aminotransferase_plus_bilirubin_jaundice_proxy": bool(
            set(pts).intersection(AMINOTRANSFERASE_PTS)
            and set(pts).intersection(BILIRUBIN_JAUNDICE_PTS)
        ),
        "matching_drug_role": role_label,
        "matching_drug_suspect_or_concomitant": best_role_code in {"1", "2"},
        "matching_drug_indications": "|".join(indications),
        "matching_drug_ild_indication": ild_indication,
        "drug_start_date": start.strftime("%Y-%m-%d") if start is not None else None,
        "case_event_date": event.strftime("%Y-%m-%d") if event is not None else None,
        "time_to_onset_days": tto,
        "time_to_onset_0_365_valid": tto is not None and 0 <= tto <= 365,
        "fatal_hepatic_reaction": any(
            str(r.get("reactionoutcome")) == "5"
            and str(r.get("reactionmeddrapt", "")).upper() in BROAD_PTS
            for r in reactions
        ),
    }
    for label, field in SERIOUS_FIELDS.items():
        result[label] = str(record.get(field, "")) == "1"
    return result


def contingency_metrics(a: int, drug_total: int, c: int, comparator_total: int) -> dict[str, Any]:
    """Active-comparator 2x2 metrics at report level."""
    b = drug_total - a
    d = comparator_total - c
    if min(a, b, c, d) < 0:
        raise ValueError(f"Invalid 2x2 table: {(a, b, c, d)}")

    # When neither arm contains an endpoint report, the exposure-event contrast
    # is non-identifiable. A continuity correction would manufacture an arbitrary
    # point estimate, so retain the counts and mark all comparative estimands NA.
    if a + c == 0:
        return {
            "nintedanib_event": a,
            "nintedanib_nonevent": b,
            "pirfenidone_event": c,
            "pirfenidone_nonevent": d,
            "nintedanib_total": drug_total,
            "pirfenidone_total": comparator_total,
            "continuity_corrected": False,
            "ROR": np.nan,
            "ROR_low95": np.nan,
            "ROR_high95": np.nan,
            "PRR": np.nan,
            "PRR_low95": np.nan,
            "PRR_high95": np.nan,
            "IC_median": np.nan,
            "IC_low95": np.nan,
            "IC_high95": np.nan,
            "fisher_exact_p": np.nan,
            "pearson_yates_chi2": np.nan,
            "pearson_yates_p": np.nan,
        }
    corrected = any(x == 0 for x in (a, b, c, d))
    aa, bb, cc, dd = (x + 0.5 for x in (a, b, c, d)) if corrected else (a, b, c, d)
    log_ror = math.log((aa * dd) / (bb * cc))
    se_log_ror = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    ror = math.exp(log_ror)
    ror_low, ror_high = math.exp(log_ror - 1.96 * se_log_ror), math.exp(log_ror + 1.96 * se_log_ror)

    risk1, risk0 = aa / (aa + bb), cc / (cc + dd)
    prr = risk1 / risk0
    se_log_prr = math.sqrt(max(0.0, 1 / aa - 1 / (aa + bb) + 1 / cc - 1 / (cc + dd)))
    prr_low, prr_high = math.exp(math.log(prr) - 1.96 * se_log_prr), math.exp(math.log(prr) + 1.96 * se_log_prr)

    odds, fisher_p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    try:
        chi2, chi2_p, _, _ = stats.chi2_contingency([[a, b], [c, d]], correction=True)
    except ValueError:
        # A calendar stratum with no endpoint reports in either arm has a zero
        # expected event column; retain the row but mark chi-square undefined.
        chi2, chi2_p = np.nan, np.nan

    # Bayesian information component: Jeffreys-Dirichlet posterior over the 2x2 table.
    # This explicitly stated active-comparator implementation avoids presenting it as
    # the proprietary WHO database's production BCPNN estimator.
    rng = np.random.default_rng(RANDOM_SEED + a + c)
    draws = rng.dirichlet(np.asarray([a, b, c, d], dtype=float) + 0.5, size=100_000)
    p11, p10, p01, _ = draws.T
    ic_draws = np.log2(p11 / ((p11 + p10) * (p11 + p01)))
    ic_med, ic_low, ic_high = np.quantile(ic_draws, [0.5, 0.025, 0.975])

    return {
        "nintedanib_event": a,
        "nintedanib_nonevent": b,
        "pirfenidone_event": c,
        "pirfenidone_nonevent": d,
        "nintedanib_total": drug_total,
        "pirfenidone_total": comparator_total,
        "continuity_corrected": corrected,
        "ROR": ror,
        "ROR_low95": ror_low,
        "ROR_high95": ror_high,
        "PRR": prr,
        "PRR_low95": prr_low,
        "PRR_high95": prr_high,
        "IC_median": ic_med,
        "IC_low95": ic_low,
        "IC_high95": ic_high,
        "fisher_exact_p": fisher_p,
        "pearson_yates_chi2": chi2,
        "pearson_yates_p": chi2_p,
    }


def endpoint_count(cases: pd.DataFrame, drug: str, endpoint: str) -> int:
    subset = cases[cases["drug"] == drug]
    if endpoint == "narrow":
        return int(subset["narrow_case"].sum())
    if endpoint == "broad":
        return int(subset["broad_case"].sum())
    if endpoint == "hy_proxy":
        return int(subset["aminotransferase_plus_bilirubin_jaundice_proxy"].sum())
    raise KeyError(endpoint)


def denominator_specs(end_date: str, aliases: dict[str, str] = GENERIC_QUERY, prefix: str = "primary") -> list[QuerySpec]:
    specs: list[QuerySpec] = []
    for drug in ("nintedanib", "pirfenidone"):
        dq = exclusive_drug_query(drug, aliases)
        specs.append(QuerySpec(f"{prefix}__total__{drug}", {"search": in_window(dq, START_DATE, end_date), "limit": 1}))
        for year in range(2014, int(end_date[:4]) + 1):
            start = max(START_DATE, f"{year}0101")
            end = min(end_date, f"{year}1231")
            specs.append(QuerySpec(f"{prefix}__year_{year}__{drug}", {"search": in_window(dq, start, end), "limit": 1}))
        for label, sq in SEX_STRATA.items():
            specs.append(QuerySpec(f"{prefix}__sex_{label}__{drug}", {"search": in_window(f"({dq}) AND {sq}", START_DATE, end_date), "limit": 1}))
        for label, aq in AGE_STRATA.items():
            specs.append(QuerySpec(f"{prefix}__age_{label}__{drug}", {"search": in_window(f"({dq}) AND {aq}", START_DATE, end_date), "limit": 1}))
    return specs


def mh_pool(year_df: pd.DataFrame, endpoint: str) -> dict[str, Any]:
    values = []
    for _, row in year_df[year_df["endpoint"] == endpoint].iterrows():
        a, b = row["nintedanib_event"], row["nintedanib_nonevent"]
        c, d = row["pirfenidone_event"], row["pirfenidone_nonevent"]
        n = a + b + c + d
        if n and a + c > 0:
            values.append((a, b, c, d, n))
    num = sum(a * d / n for a, b, c, d, n in values)
    den = sum(b * c / n for a, b, c, d, n in values)
    mh = num / den if den else np.nan
    # Canonical fixed-stratum Mantel-Haenszel confidence interval using the
    # Robins-Breslow-Greenland large-sample variance implemented by statsmodels.
    tables = np.stack(
        [np.asarray([[a, b], [c, d]], dtype=float) for a, b, c, d, _ in values],
        axis=2,
    )
    stratified = StratifiedTable(tables)
    mh_rbg = float(stratified.oddsratio_pooled)
    mh_rbg_low, mh_rbg_high = map(float, stratified.oddsratio_pooled_confint())
    # Robust uncertainty by parametric multinomial bootstrap within year.
    rng = np.random.default_rng(RANDOM_SEED + len(values))
    boots = []
    for _ in range(30_000):
        bn = bd = 0.0
        for a, b, c, d, n in values:
            draw = rng.multinomial(int(n), np.asarray([a, b, c, d]) / n)
            aa, bb, cc, dd = draw
            bn += (aa + 0.5) * (dd + 0.5) / (n + 2)
            bd += (bb + 0.5) * (cc + 0.5) / (n + 2)
        if bd > 0:
            boots.append(bn / bd)
    low, high = np.quantile(boots, [0.025, 0.975]) if boots else (np.nan, np.nan)
    return {
        "endpoint": endpoint,
        "MH_ROR_year_adjusted": mh,
        "MH_ROR_RBG_check": mh_rbg,
        "MH_low95_RBG": mh_rbg_low,
        "MH_high95_RBG": mh_rbg_high,
        "MH_low95_bootstrap_sensitivity": low,
        "MH_high95_bootstrap_sensitivity": high,
        "n_year_strata": len(values),
    }


def main(refresh: bool = False, workers: int = 3) -> None:
    client = OpenFDAClient(refresh=refresh, max_workers=workers)

    latest = client.get(QuerySpec("dataset_latest_receivedate", {"sort": "receivedate:desc", "limit": 1}))
    end_date = str(latest["results"][0]["receivedate"])
    last_updated = latest.get("meta", {}).get("last_updated")
    if end_date != EXPECTED_DATA_END:
        raise RuntimeError(
            f"Data release changed (end={end_date}, expected={EXPECTED_DATA_END}). "
            "Review and update the frozen analysis window before rerunning."
        )

    endpoint_rows = []
    for name, pts in (("narrow", NARROW_PTS), ("broad", BROAD_PTS)):
        for pt in pts:
            endpoint_rows.append({"endpoint": name, "preferred_term": pt, "prespecified": True})
    pd.DataFrame(endpoint_rows).to_csv(OUT_DIR / "prespecified_meddra_pt_sets.csv", index=False)

    # Primary denominator queries and aliases sensitivity denominator queries.
    denom_data = client.many(denominator_specs(end_date, GENERIC_QUERY, "primary"))
    expanded_data = client.many(denominator_specs(end_date, EXPANDED_QUERY, "expanded"))

    # Fetch all broad-endpoint cases for the primary query; narrow is a strict subset.
    broad_rx = q_or_terms(BROAD_PTS)
    records: list[dict[str, Any]] = []
    expected_by_drug: dict[str, int] = {}
    for drug in ("nintedanib", "pirfenidone"):
        search = in_window(f"({exclusive_drug_query(drug)}) AND {broad_rx}", START_DATE, end_date)
        count_data = client.get(QuerySpec(f"primary__broad_count__{drug}", {"search": search, "limit": 1}))
        expected = total_from(count_data)
        expected_by_drug[drug] = expected
        fetched = fetch_all(client, f"primary__broad_cases__{drug}", search, expected)
        records.extend(tidy_report(record, drug) for record in fetched)

    cases = pd.DataFrame(records)
    cases.to_csv(OUT_DIR / "broad_liver_event_reports_deduplicated.csv", index=False)

    # Completeness invariants.
    for drug, expected in expected_by_drug.items():
        observed = int((cases["drug"] == drug).sum())
        if observed != expected:
            raise RuntimeError(f"Case completeness failed for {drug}: {observed} != {expected}")
    if cases["safetyreportid"].duplicated().any():
        # Two-drug overlap was excluded in the query, so duplicates are not expected even across groups.
        raise RuntimeError("Duplicate safetyreportid remained after active-comparator exclusion")

    totals = {
        drug: total_from(denom_data[f"primary__total__{drug}"])
        for drug in ("nintedanib", "pirfenidone")
    }

    # Overall endpoint signal estimates.
    overall = []
    for endpoint in ("narrow", "broad", "hy_proxy"):
        row = contingency_metrics(
            endpoint_count(cases, "nintedanib", endpoint), totals["nintedanib"],
            endpoint_count(cases, "pirfenidone", endpoint), totals["pirfenidone"],
        )
        row.update({"analysis": "primary_generic_exclusive", "endpoint": endpoint, "start_date": START_DATE, "end_date": end_date})
        overall.append(row)
    overall_df = pd.DataFrame(overall)
    overall_df.to_csv(OUT_DIR / "overall_active_comparator_signals.csv", index=False)

    # Brand/product-alias sensitivity uses API aggregate counts rather than reusing primary records.
    sensitivity_specs: list[QuerySpec] = []
    for drug in ("nintedanib", "pirfenidone"):
        dq = exclusive_drug_query(drug, EXPANDED_QUERY)
        for endpoint, pts in (("narrow", NARROW_PTS), ("broad", BROAD_PTS)):
            sensitivity_specs.append(QuerySpec(
                f"expanded__{endpoint}__{drug}",
                {"search": in_window(f"({dq}) AND {q_or_terms(pts)}", START_DATE, end_date), "limit": 1},
            ))
    sensitivity_data = client.many(sensitivity_specs)
    expanded_totals = {
        drug: total_from(expanded_data[f"expanded__total__{drug}"])
        for drug in ("nintedanib", "pirfenidone")
    }
    sensitivity_rows = []
    for endpoint in ("narrow", "broad"):
        row = contingency_metrics(
            total_from(sensitivity_data[f"expanded__{endpoint}__nintedanib"]), expanded_totals["nintedanib"],
            total_from(sensitivity_data[f"expanded__{endpoint}__pirfenidone"]), expanded_totals["pirfenidone"],
        )
        row.update({"analysis": "expanded_exact_aliases_exclusive", "endpoint": endpoint, "start_date": START_DATE, "end_date": end_date})
        sensitivity_rows.append(row)
    pd.DataFrame(sensitivity_rows).to_csv(OUT_DIR / "alias_sensitivity_signals.csv", index=False)

    # Report-level clinical-context specifications. These are useful checks for
    # non-ILD use of nintedanib, but the flattened report query cannot prove that
    # role/indication and name belong to the same drug-array element. Accordingly,
    # these remain sensitivity specifications rather than replacing the primary.
    specification_filters = {
        "ild_indication_report_filter": q_or_indications(ILD_INDICATION_TERMS),
        "ild_indication_exact_keyword_sensitivity": q_or_indications_exact(ILD_INDICATION_TERMS),
        "suspect_drug_report_filter": 'patient.drug.drugcharacterization:"1"',
        "ild_plus_suspect_drug_report_filter": (
            f"({q_or_indications(ILD_INDICATION_TERMS)}) AND "
            'patient.drug.drugcharacterization:"1"'
        ),
    }
    specification_specs: list[QuerySpec] = []
    for specification, context_query in specification_filters.items():
        for drug in ("nintedanib", "pirfenidone"):
            dq = exclusive_drug_query(drug)
            denominator_search = in_window(
                f"({dq}) AND ({context_query})", START_DATE, end_date
            )
            specification_specs.append(QuerySpec(
                f"specification__{specification}__total__{drug}",
                {"search": denominator_search, "limit": 1},
            ))
            for endpoint, pts in (("narrow", NARROW_PTS), ("broad", BROAD_PTS)):
                specification_specs.append(QuerySpec(
                    f"specification__{specification}__{endpoint}__{drug}",
                    {
                        "search": in_window(
                            f"({dq}) AND ({context_query}) AND {q_or_terms(pts)}",
                            START_DATE,
                            end_date,
                        ),
                        "limit": 1,
                    },
                ))
    specification_data = client.many(specification_specs)
    specification_rows = []
    for specification in specification_filters:
        specification_totals = {
            drug: total_from(specification_data[
                f"specification__{specification}__total__{drug}"
            ])
            for drug in ("nintedanib", "pirfenidone")
        }
        for endpoint in ("narrow", "broad"):
            row = contingency_metrics(
                total_from(specification_data[
                    f"specification__{specification}__{endpoint}__nintedanib"
                ]),
                specification_totals["nintedanib"],
                total_from(specification_data[
                    f"specification__{specification}__{endpoint}__pirfenidone"
                ]),
                specification_totals["pirfenidone"],
            )
            row.update({
                "analysis": specification,
                "endpoint": endpoint,
                "filter_query": specification_filters[specification],
                "interpretation_boundary": (
                    "report-level openFDA filter; role/indication cannot be proven "
                    "to belong to the same drug-array element as the indexed product"
                ),
                "start_date": START_DATE,
                "end_date": end_date,
            })
            specification_rows.append(row)
    pd.DataFrame(specification_rows).to_csv(
        OUT_DIR / "report_level_specification_signals.csv", index=False
    )

    # Year-specific active-comparator estimates.
    annual_rows = []
    for year in range(2014, int(end_date[:4]) + 1):
        year_cases = cases[cases["report_year"] == year]
        year_totals = {
            drug: total_from(denom_data[f"primary__year_{year}__{drug}"])
            for drug in ("nintedanib", "pirfenidone")
        }
        for endpoint in ("narrow", "broad"):
            a = endpoint_count(year_cases, "nintedanib", endpoint)
            c = endpoint_count(year_cases, "pirfenidone", endpoint)
            if year_totals["nintedanib"] == 0 or year_totals["pirfenidone"] == 0:
                continue
            row = contingency_metrics(a, year_totals["nintedanib"], c, year_totals["pirfenidone"])
            row.update({"year": year, "endpoint": endpoint})
            annual_rows.append(row)
    annual_df = pd.DataFrame(annual_rows)
    annual_df.to_csv(OUT_DIR / "annual_active_comparator_signals.csv", index=False)
    pd.DataFrame([mh_pool(annual_df, endpoint) for endpoint in ("narrow", "broad")]).to_csv(
        OUT_DIR / "year_adjusted_mantel_haenszel_signals.csv", index=False
    )
    compute_year_effect_heterogeneity(annual_df).to_csv(
        OUT_DIR / "year_effect_heterogeneity_tests.csv", index=False
    )

    # Influence analysis for the conspicuous 2022 narrow-endpoint reporting
    # cluster. This is explicitly labelled post hoc and complements, rather than
    # replaces, the full-window primary analysis.
    temporal_windows = {
        "all_years_primary": lambda d: d,
        "exclude_2022_posthoc": lambda d: d[d["year"] != 2022],
        "pre_2022": lambda d: d[d["year"] <= 2021],
        "year_2022_only": lambda d: d[d["year"] == 2022],
        "post_2022": lambda d: d[d["year"] >= 2023],
    }
    temporal_rows = []
    for endpoint in ("narrow", "broad"):
        endpoint_annual = annual_df[annual_df["endpoint"] == endpoint]
        for window_label, selector in temporal_windows.items():
            selected = selector(endpoint_annual)
            if selected.empty:
                continue
            row = contingency_metrics(
                int(selected["nintedanib_event"].sum()),
                int(selected["nintedanib_total"].sum()),
                int(selected["pirfenidone_event"].sum()),
                int(selected["pirfenidone_total"].sum()),
            )
            row.update({
                "endpoint": endpoint,
                "temporal_window": window_label,
                "first_year": int(selected["year"].min()),
                "last_year": int(selected["year"].max()),
                "prespecification": (
                    "primary" if window_label == "all_years_primary" else "post_hoc_influence_analysis"
                ),
            })
            temporal_rows.append(row)
    pd.DataFrame(temporal_rows).to_csv(
        OUT_DIR / "temporal_influence_signals.csv", index=False
    )

    # Descriptive audit of the 2022 narrow-endpoint cluster. Country and month
    # counts have no exposure denominators and must not be interpreted as rates.
    cluster = cases[(cases["report_year"] == 2022) & cases["narrow_case"]].copy()
    cluster["country_for_audit"] = cluster["primarysourcecountry"].fillna("MISSING")
    cluster["receive_month"] = pd.to_datetime(
        cluster["receivedate"].astype(str), format="%Y%m%d", errors="coerce"
    ).dt.strftime("%Y-%m")
    country_rows = (
        cluster.groupby(["drug", "country_for_audit"], dropna=False)
        .size().rename("reports").reset_index()
    )
    country_rows["drug_2022_narrow_reports"] = country_rows.groupby("drug")["reports"].transform("sum")
    country_rows["percent"] = 100 * country_rows["reports"] / country_rows["drug_2022_narrow_reports"]
    country_rows.to_csv(OUT_DIR / "narrow_2022_country_audit.csv", index=False)
    month_rows = (
        cluster.groupby(["drug", "receive_month"], dropna=False)
        .size().rename("reports").reset_index()
    )
    month_rows.to_csv(OUT_DIR / "narrow_2022_month_audit.csv", index=False)

    # Sex and age-unit-years-only stratum estimates.
    stratum_rows = []
    for dimension, strata in (("sex", SEX_STRATA), ("age", AGE_STRATA)):
        for stratum in strata:
            if dimension == "sex":
                subset = cases[cases["sex"] == stratum]
            elif stratum == "18-64_years":
                subset = cases[
                    (cases["age_unit_code"].astype(str) == "801")
                    & cases["age_value_reported"].between(18, 64, inclusive="both")
                ]
            else:
                subset = cases[
                    (cases["age_unit_code"].astype(str) == "801")
                    & cases["age_value_reported"].between(65, 150, inclusive="both")
                ]
            stotals = {
                drug: total_from(denom_data[f"primary__{dimension}_{stratum}__{drug}"])
                for drug in ("nintedanib", "pirfenidone")
            }
            for endpoint in ("narrow", "broad"):
                row = contingency_metrics(
                    endpoint_count(subset, "nintedanib", endpoint), stotals["nintedanib"],
                    endpoint_count(subset, "pirfenidone", endpoint), stotals["pirfenidone"],
                )
                row.update({"dimension": dimension, "stratum": stratum, "endpoint": endpoint})
                stratum_rows.append(row)
    pd.DataFrame(stratum_rows).to_csv(OUT_DIR / "sex_age_stratified_signals.csv", index=False)

    # PT-specific report counts and signals, preserving overlap across PTs.
    pt_rows = []
    for pt in BROAD_PTS:
        counts = {}
        for drug in ("nintedanib", "pirfenidone"):
            subset = cases[cases["drug"] == drug]
            counts[drug] = int(subset["all_pts"].fillna("").str.split("|").apply(lambda xs: pt in xs).sum())
        row = contingency_metrics(counts["nintedanib"], totals["nintedanib"], counts["pirfenidone"], totals["pirfenidone"])
        row.update({"preferred_term": pt, "endpoint_membership": "narrow" if pt in NARROW_PTS else "broad_only"})
        pt_rows.append(row)
    pd.DataFrame(pt_rows).sort_values("ROR", ascending=False).to_csv(OUT_DIR / "preferred_term_signals.csv", index=False)

    # Serious outcomes are descriptive within endpoint cases; categories may overlap.
    outcome_rows = []
    for endpoint in ("narrow", "broad"):
        endpoint_col = f"{endpoint}_case"
        for drug in ("nintedanib", "pirfenidone"):
            subset = cases[(cases["drug"] == drug) & cases[endpoint_col]]
            for outcome in list(SERIOUS_FIELDS) + ["fatal_hepatic_reaction"]:
                n = int(subset[outcome].sum())
                outcome_rows.append({
                    "endpoint": endpoint,
                    "drug": drug,
                    "outcome": outcome,
                    "events": n,
                    "endpoint_reports": len(subset),
                    "percent": 100 * n / len(subset) if len(subset) else np.nan,
                })
    pd.DataFrame(outcome_rows).to_csv(OUT_DIR / "serious_outcome_profiles.csv", index=False)

    # Matching-drug role and ILD-indication sensitivity descriptors.
    role_rows = []
    for endpoint in ("narrow", "broad"):
        endpoint_col = f"{endpoint}_case"
        for drug in ("nintedanib", "pirfenidone"):
            subset = cases[(cases["drug"] == drug) & cases[endpoint_col]]
            for variable in ("matching_drug_role", "matching_drug_ild_indication"):
                counts = subset[variable].value_counts(dropna=False)
                for level, count in counts.items():
                    role_rows.append({
                        "endpoint": endpoint, "drug": drug, "variable": variable,
                        "level": level, "reports": int(count), "endpoint_reports": len(subset),
                        "percent": 100 * count / len(subset) if len(subset) else np.nan,
                    })
    pd.DataFrame(role_rows).to_csv(OUT_DIR / "drug_role_and_indication_profiles.csv", index=False)

    # Time-to-onset: only exact YYYYMMDD dates and 0-365 days in the primary summary.
    tto_rows = []
    tto_summary = []
    tto_exclusion_audit = []
    for endpoint in ("narrow", "broad"):
        endpoint_col = f"{endpoint}_case"
        for drug in ("nintedanib", "pirfenidone"):
            subset = cases[(cases["drug"] == drug) & cases[endpoint_col]].copy()
            exact_pairs = subset[subset["time_to_onset_days"].notna()].copy()
            valid = subset[subset["time_to_onset_0_365_valid"]].copy()
            vals = valid["time_to_onset_days"].astype(float).to_numpy()
            tto_exclusion_audit.append({
                "endpoint": endpoint,
                "drug": drug,
                "endpoint_reports": len(subset),
                "exact_date_pairs_before_window": len(exact_pairs),
                "negative_intervals_excluded": int((exact_pairs["time_to_onset_days"] < 0).sum()),
                "intervals_0_365_included": len(valid),
                "intervals_over_365_excluded": int((exact_pairs["time_to_onset_days"] > 365).sum()),
                "exact_pair_inclusion_percent": (
                    100 * len(valid) / len(exact_pairs) if len(exact_pairs) else np.nan
                ),
            })
            beta = scale = np.nan
            if len(vals) >= 3:
                # +0.5 permits same-day reports while fixing location at zero.
                beta, _, scale = stats.weibull_min.fit(vals + 0.5, floc=0)
            tto_summary.append({
                "endpoint": endpoint, "drug": drug,
                "endpoint_reports": len(subset), "exact_date_pairs_0_365": len(vals),
                "date_pair_completeness_percent": 100 * len(vals) / len(subset) if len(subset) else np.nan,
                "median_days": np.median(vals) if len(vals) else np.nan,
                "q1_days": np.quantile(vals, 0.25) if len(vals) else np.nan,
                "q3_days": np.quantile(vals, 0.75) if len(vals) else np.nan,
                "weibull_shape_beta": beta, "weibull_scale_days": scale,
            })
            for _, row in valid.iterrows():
                tto_rows.append({
                    "endpoint": endpoint, "drug": drug,
                    "safetyreportid": row["safetyreportid"],
                    "time_to_onset_days": int(row["time_to_onset_days"]),
                    "matching_drug_role": row["matching_drug_role"],
                })
    pd.DataFrame(tto_rows).to_csv(OUT_DIR / "time_to_onset_case_values.csv", index=False)
    pd.DataFrame(tto_summary).to_csv(OUT_DIR / "time_to_onset_summary.csv", index=False)
    pd.DataFrame(tto_exclusion_audit).to_csv(
        OUT_DIR / "time_to_onset_exclusion_audit.csv", index=False
    )

    # Frozen configuration and integrity summary.
    config = {
        "retrieval_date": RETRIEVAL_DATE,
        "api_endpoint": BASE_URL,
        "api_last_updated": last_updated,
        "analysis_start_receivedate": START_DATE,
        "analysis_end_receivedate": end_date,
        "primary_drug_queries": GENERIC_QUERY,
        "expanded_alias_queries": EXPANDED_QUERY,
        "coexposure_rule": "reports containing both nintedanib and pirfenidone are excluded from both active-comparator arms",
        "narrow_pts": NARROW_PTS,
        "broad_only_pts": BROAD_ONLY_PTS,
        "ild_indication_terms_for_report_level_specification": ILD_INDICATION_TERMS,
        "drugcharacterization_schema": {
            "field": "patient.drug.drugcharacterization",
            "official_mapping": DRUG_CHARACTERIZATION_LABELS,
            "source": "https://open.fda.gov/fields/drugevent.yaml",
        },
        "age_rule": "age-stratified denominators require ICH age unit 801 (years); bands 18-64 and 65-150",
        "time_to_onset_rule": "exact YYYYMMDD drug start and CASE EVENT DATE only; primary window 0-365 days",
        "random_seed": RANDOM_SEED,
        "raw_record_counts": expected_by_drug,
        "primary_exclusive_drug_report_totals": totals,
        "response_cache_directory": "fresh_analysis/faers/raw_cache",
    }
    (OUT_DIR / "analysis_configuration.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    client.save_log()
    print(json.dumps({
        "status": "complete",
        "last_updated": last_updated,
        "receivedate_end": end_date,
        "drug_report_totals": totals,
        "broad_case_counts": expected_by_drug,
        "outputs": str(OUT_DIR),
    }, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true", help="ignore cached API responses")
    parser.add_argument("--workers", type=int, default=3, choices=range(1, 5))
    args = parser.parse_args()
    main(refresh=args.refresh, workers=args.workers)
