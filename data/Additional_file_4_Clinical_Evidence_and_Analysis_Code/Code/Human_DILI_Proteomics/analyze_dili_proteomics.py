#!/usr/bin/env python3
"""Reproducible re-analysis of the Ravindra et al. DILI source-data workbook.

This script deliberately treats every protein block as a separate de-identified
series.  The source workbook does not expose a stable participant identifier
linking proteins, and the row orders differ between figure-source blocks.
Consequently, no multivariable participant-level model is fitted.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests


SEED = 20260904
N_BOOT = 2000
RNG = np.random.default_rng(SEED)

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "fresh_data_dili_source.xlsx"
OUT = Path(__file__).resolve().parent

CANDIDATES = [
    "ACO1", "ALDOB", "ASS1", "CES1", "CPS1", "DMGDH", "FAH",
    "FBP1", "GSTA1", "HPD", "LECT2", "OTC", "PCK2",
]
COMPARATORS = ["GLDH", "CK18"]

DISCOVERY_N = {"HV": 10, "DO": 10, "DF": 10, "NDO": 5, "NDF": 5, "NAFLD": 10}
CONFIRMATORY_N = {"HV": 60, "DO": 82, "DF": 77, "NDO": 34, "NDF": 22}

CONTRASTS_DISCOVERY = [
    ("DO", "HV"), ("DO", "NDO"), ("DO", "DF"),
    ("NDO", "HV"), ("NDO", "NDF"), ("DO", "NAFLD"),
]
CONTRASTS_CONFIRMATORY = [
    ("DO", "HV"), ("DO", "NDO"), ("DO", "DF"),
    ("NDO", "HV"), ("NDO", "NDF"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def parse_group(value: object) -> str:
    text = clean_text(value)
    match = re.match(r"([A-Z]+)", text)
    return match.group(1) if match else text


def parse_value(value: object) -> tuple[float, bool, bool, str]:
    """Return numeric value, ND flag, starred flag, and original text."""
    raw = clean_text(value)
    if not raw:
        return np.nan, False, False, raw
    if raw.upper() in {"ND", "NOT DETECTED", "N/A", "NA"}:
        return np.nan, True, False, raw
    starred = "*" in raw
    normalized = raw.replace("*", "").replace(",", "")
    try:
        return float(normalized), False, starred, raw
    except ValueError:
        return np.nan, False, starred, raw


def invalid_text_flag(value: float, nd: bool, raw: str) -> bool:
    return bool(raw and not nd and not np.isfinite(value))


def marker_rows(df: pd.DataFrame, marker: str) -> list[int]:
    target = marker.strip().lower()
    return [i for i, v in df.iloc[:, 0].items() if clean_text(v).lower() == target]


def find_header_row(df: pd.DataFrame, start: int, end: int) -> int:
    for idx in range(start, min(end, start + 10)):
        first = clean_text(df.iat[idx, 0]).lower()
        if first.startswith("subgroup") or first == "biomarker":
            return idx
    raise ValueError(f"No header row found from {start} to {end}")


def find_protein(df: pd.DataFrame, start: int, header: int) -> str:
    text = " | ".join(clean_text(x).upper() for x in df.iloc[start:header].to_numpy().ravel())
    hits = [p for p in CANDIDATES + COMPARATORS if re.search(rf"\b{re.escape(p)}\b", text)]
    if not hits:
        raise ValueError(f"No protein label found in rows {start}:{header}")
    return hits[0]


def parse_column_block(
    df: pd.DataFrame,
    header_row: int,
    protein: str,
    cohort: str,
    expected_n: dict[str, int],
    source_sheet: str,
    assay: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for col in range(1, df.shape[1]):
        group = parse_group(df.iat[header_row, col])
        if group not in expected_n:
            continue
        for j in range(expected_n[group]):
            source_row = header_row + 1 + j
            value, nd, starred, raw = parse_value(df.iat[source_row, col])
            rows.append(
                {
                    "cohort": cohort,
                    "protein": protein,
                    "group": group,
                    "source_row_index_within_group": j + 1,
                    "value": value,
                    "is_not_detected": nd,
                    "is_invalid_nonnumeric_text": invalid_text_flag(value, nd, raw),
                    "is_starred_source_value": starred,
                    "source_value_text": raw,
                    "assay": assay,
                    "source_sheet": source_sheet,
                    "source_excel_row": source_row + 1,
                    "source_excel_column": col + 1,
                }
            )
    return pd.DataFrame(rows)


def parse_observed_column_block(
    df: pd.DataFrame,
    header_row: int,
    block_end: int,
    protein: str,
    cohort: str,
    allowed_groups: set[str],
    source_sheet: str,
    assay: str,
) -> pd.DataFrame:
    """Parse every supplied cell through the last nonblank cell in each group column.

    This is used for supplementary follow-up blocks whose observed column lengths
    disagree with the publication-level sample counts. Internal blank cells are
    retained as missing; trailing structural padding is not emitted.
    """
    rows: list[dict[str, object]] = []
    for col in range(1, df.shape[1]):
        group = parse_group(df.iat[header_row, col])
        if group not in allowed_groups:
            continue
        raw_cells = [clean_text(df.iat[r, col]) for r in range(header_row + 1, block_end)]
        populated = [j for j, text in enumerate(raw_cells) if text]
        if not populated:
            continue
        last = max(populated)
        for j in range(last + 1):
            source_row = header_row + 1 + j
            value, nd, starred, raw = parse_value(df.iat[source_row, col])
            rows.append(
                {
                    "cohort": cohort,
                    "protein": protein,
                    "group": group,
                    "source_row_index_within_group": j + 1,
                    "value": value,
                    "is_not_detected": nd,
                    "is_invalid_nonnumeric_text": invalid_text_flag(value, nd, raw),
                    "is_starred_source_value": starred,
                    "source_value_text": raw,
                    "assay": assay,
                    "source_sheet": source_sheet,
                    "source_excel_row": source_row + 1,
                    "source_excel_column": col + 1,
                }
            )
    return pd.DataFrame(rows)


def parse_discovery() -> pd.DataFrame:
    df = pd.read_excel(SOURCE, sheet_name="Figure 2C", header=None, dtype=object)
    starts = marker_rows(df, "Fig 2C")
    frames = []
    for k, start in enumerate(starts):
        end = starts[k + 1] if k + 1 < len(starts) else len(df)
        header = find_header_row(df, start, end)
        protein = find_protein(df, start, header)
        frames.append(
            parse_column_block(
                df, header, protein, "discovery", DISCOVERY_N,
                "Figure 2C", "TMT VSN-normalized relative intensity",
            )
        )

    df_lect2 = pd.read_excel(SOURCE, sheet_name="Supp Fig 1-LECT2", header=None, dtype=object)
    header = find_header_row(df_lect2, 0, len(df_lect2))
    frames.append(
        parse_column_block(
            df_lect2, header, "LECT2", "discovery", DISCOVERY_N,
            "Supp Fig 1-LECT2", "TMT VSN-normalized relative intensity",
        )
    )
    out = pd.concat(frames, ignore_index=True)
    out["analysis_sample_key"] = (
        "DISC_" + out["group"] + "_" +
        out["source_row_index_within_group"].astype(int).astype(str).str.zfill(3)
    )
    return out


def parse_supp_figure2() -> pd.DataFrame:
    df = pd.read_excel(SOURCE, sheet_name="Supp Fig 2", header=None, dtype=object)
    starts = marker_rows(df, "Supp Fig 2")
    frames = []
    for k, start in enumerate(starts):
        end = starts[k + 1] if k + 1 < len(starts) else len(df)
        header = find_header_row(df, start, end)
        protein = find_protein(df, start, header)
        frames.append(parse_observed_column_block(
            df, header, end, protein, "confirmatory", {"DO", "DF", "NDO", "NDF"},
            "Supp Fig 2", "SureQuant light/heavy relative intensity",
        ))
    return pd.concat(frames, ignore_index=True)


def parse_figure3() -> pd.DataFrame:
    df = pd.read_excel(SOURCE, sheet_name="Figure 3", header=None, dtype=object)
    starts = marker_rows(df, "Fig 3")
    frames = []
    for k, start in enumerate(starts):
        end = starts[k + 1] if k + 1 < len(starts) else len(df)
        header = find_header_row(df, start, end)
        protein = find_protein(df, start, header)
        frames.append(
            parse_column_block(
                df, header, protein, "confirmatory", {"HV": 60, "DO": 82, "NDO": 34},
                "Figure 3", "SureQuant light/heavy relative intensity",
            )
        )
    return pd.concat(frames, ignore_index=True)


def parse_pck2() -> pd.DataFrame:
    df = pd.read_excel(SOURCE, sheet_name="Supp Fig 3-PCK2", header=None, dtype=object)
    header = find_header_row(df, 0, len(df))
    out = parse_observed_column_block(
        df, header, len(df), "PCK2", "confirmatory", set(CONFIRMATORY_N),
        "Supp Fig 3-PCK2", "PCK2 ELISA (source units not stated in workbook)",
    )
    return out


def compare_repeated_series(fig3: pd.DataFrame, supp2: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein in sorted(set(fig3.protein) & set(supp2.protein)):
        for group in ["DO", "NDO"]:
            a = fig3.query("protein == @protein and group == @group").sort_values(
                "source_row_index_within_group"
            ).value.to_numpy()
            b = supp2.query("protein == @protein and group == @group").sort_values(
                "source_row_index_within_group"
            ).value.to_numpy()
            a_finite = a[np.isfinite(a)]
            b_finite = b[np.isfinite(b)]
            min_n = min(len(a), len(b))
            pair_complete = np.isfinite(a[:min_n]) & np.isfinite(b[:min_n])
            same_numeric_multiset = (
                len(a_finite) == len(b_finite)
                and np.allclose(np.sort(a_finite), np.sort(b_finite), atol=1e-12, rtol=1e-12)
            )
            rows.append(
                {
                    "protein": protein,
                    "group": group,
                    "n_figure3_rows": len(a),
                    "n_supp_figure2_rows": len(b),
                    "n_pair_complete": int(pair_complete.sum()),
                    "n_missingness_discordant_in_overlapping_positions": int(
                        np.sum(np.isfinite(a[:min_n]) != np.isfinite(b[:min_n]))
                    ),
                    "n_value_discordant_at_tolerance_1e-12": int(
                        np.sum(~np.isclose(
                            a[:min_n][pair_complete], b[:min_n][pair_complete], atol=1e-12, rtol=1e-12
                        ))
                    ),
                    "max_absolute_difference": (
                        float(np.max(np.abs(a[:min_n][pair_complete] - b[:min_n][pair_complete])))
                        if pair_complete.any() else np.nan
                    ),
                    "same_numeric_value_multiset": bool(same_numeric_multiset),
                }
            )
    return pd.DataFrame(rows)


def assemble_confirmatory() -> tuple[pd.DataFrame, pd.DataFrame]:
    supp2 = parse_supp_figure2()
    fig3 = parse_figure3()
    pck2 = parse_pck2()
    repeated_qc = compare_repeated_series(fig3, supp2)

    # Figure 3 is the canonical onset source because it declares n for HV/DO/NDO.
    # Supp Fig 2 contributes only DF/NDF; its onset copies are retained for QC.
    onset = fig3.query("group in ['HV', 'DO', 'NDO']").copy()
    followup = supp2.query("group in ['DF', 'NDF']").copy()
    pck2_use = pck2.query("group in ['HV', 'DO', 'NDO', 'DF', 'NDF']").copy()
    out = pd.concat([onset, followup, pck2_use], ignore_index=True)
    out["analysis_sample_key"] = (
        "CONF_" + out["group"] + "_" +
        out["source_row_index_within_group"].astype(int).astype(str).str.zfill(3)
    )
    return out, repeated_qc


def parse_zonal_scores() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_excel(SOURCE, sheet_name="Figure 6")
    out = df.iloc[:, :3].copy()
    out.columns = ["zone", "score", "group"]
    out = out[out.zone.isin(["Zone 1", "Zone 2", "Zone 3"])].copy()
    out["score"] = pd.to_numeric(out["score"], errors="coerce")
    out["group"] = out["group"].astype(str).str.strip()
    out["source_sequence"] = out.groupby(["zone", "group"]).cumcount() + 1

    duplicates = []
    cleaned = []
    for zone, zdf in out.groupby("zone", sort=False):
        hv = zdf[zdf.group == "HV"].copy()
        first = hv.iloc[:10].score.to_numpy()
        second = hv.iloc[10:20].score.to_numpy()
        exact_repeat = len(hv) == 20 and np.array_equal(first, second)
        duplicates.append(
            {
                "zone": zone,
                "hv_rows_in_source": len(hv),
                "second_ten_exact_duplicate_of_first_ten": exact_repeat,
                "hv_rows_retained": 10 if exact_repeat else len(hv),
            }
        )
        if exact_repeat:
            drop_index = hv.iloc[10:].index
            zdf = zdf.drop(index=drop_index)
        cleaned.append(zdf)
    out = pd.concat(cleaned, ignore_index=True)
    out["source_row_index_within_group"] = out.groupby(["zone", "group"]).cumcount() + 1
    out["analysis_sample_key"] = (
        "ZONE_" + out.zone.str.replace(" ", "", regex=False) + "_" + out.group + "_" +
        out.source_row_index_within_group.astype(str).str.zfill(3)
    )
    out["source_sheet"] = "Figure 6"
    return out, pd.DataFrame(duplicates)


def parse_correlation_sheet(sheet: str, marker: str) -> pd.DataFrame:
    df = pd.read_excel(SOURCE, sheet_name=sheet, header=None, dtype=object)
    starts = marker_rows(df, marker)
    frames = []
    for k, start in enumerate(starts):
        end = starts[k + 1] if k + 1 < len(starts) else len(df)
        header = find_header_row(df, start, end)
        protein = find_protein(df, start, header)
        labels = [clean_text(x).upper() for x in df.iloc[header].tolist()]
        alt_cols = [i for i, x in enumerate(labels) if "ALT" in x]
        if not alt_cols:
            raise ValueError(f"No ALT column in {sheet}, protein {protein}")
        alt_col = alt_cols[0]
        group_cols = [i for i, x in enumerate(labels) if x in {"SUBGROUP", "TYPE"}]
        if not group_cols:
            raise ValueError(f"No group column in {sheet}, protein {protein}")
        group_col = group_cols[-1] if labels[group_cols[-1]] == "TYPE" else group_cols[0]
        biomarker_cols = [
            i for i, x in enumerate(labels)
            if i != alt_col and (protein in x or "BIOMARKER VALUE" in x)
        ]
        if not biomarker_cols:
            raise ValueError(f"No biomarker column in {sheet}, protein {protein}")
        biomarker_col = biomarker_cols[0]

        records = []
        for source_row in range(header + 1, end):
            group = clean_text(df.iat[source_row, group_col]).upper()
            if group not in {"HV", "DO"}:
                continue
            alt, alt_nd, alt_star, alt_raw = parse_value(df.iat[source_row, alt_col])
            val, nd, starred, raw = parse_value(df.iat[source_row, biomarker_col])
            records.append(
                {
                    "protein": protein,
                    "group": group,
                    "log2_biomarker": val,
                    "log2_alt": alt,
                    "biomarker_is_not_detected": nd,
                    "biomarker_is_invalid_nonnumeric_text": invalid_text_flag(val, nd, raw),
                    "biomarker_is_starred_source_value": starred,
                    "source_biomarker_text": raw,
                    "source_alt_text": alt_raw,
                    "source_sheet": sheet,
                    "source_excel_row": source_row + 1,
                }
            )
        block = pd.DataFrame(records)
        block["source_row_index_within_group_and_protein"] = block.groupby("group").cumcount() + 1
        frames.append(block)
    return pd.concat(frames, ignore_index=True)


def parse_alt_correlations() -> pd.DataFrame:
    fig5 = parse_correlation_sheet("Figure 5", "Fig 5")
    supp7 = parse_correlation_sheet("Supp Fig 7", "Supp Fig 7")
    out = pd.concat([fig5, supp7], ignore_index=True)

    # Three Supp Fig 7 blocks numerically substitute values for observations
    # labeled ND in the canonical boxplot sheets, without documenting the rule.
    hp_otc_placeholder = np.isclose(out.log2_biomarker, -14.2877123795494, atol=1e-10)
    pck2_placeholder = np.isclose(out.log2_biomarker, -2.867752, atol=1e-5)
    out["is_likely_numeric_ND_substitution"] = (
        (out.protein.isin(["HPD", "OTC"]) & hp_otc_placeholder)
        | ((out.protein == "PCK2") & pck2_placeholder)
    )
    out["is_ck18_reported_lower_limit"] = (
        (out.protein == "CK18")
        & np.isclose(out.log2_biomarker, math.log2(100), atol=1e-10)
    )
    out["log2_biomarker_primary"] = out.log2_biomarker.mask(
        out.is_likely_numeric_ND_substitution
    )
    return out


def bootstrap_two_sample(
    x: np.ndarray,
    y: np.ndarray,
    statistic,
    n_boot: int = N_BOOT,
) -> tuple[float, float]:
    vals = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        xb = RNG.choice(x, size=len(x), replace=True)
        yb = RNG.choice(y, size=len(y), replace=True)
        try:
            vals[i] = statistic(xb, yb)
        except Exception:
            vals[i] = np.nan
    finite = vals[np.isfinite(vals)]
    if len(finite) < max(100, n_boot // 2):
        return np.nan, np.nan
    return tuple(np.percentile(finite, [2.5, 97.5]))


def bootstrap_effect_bundle(
    x: np.ndarray,
    y: np.ndarray,
    include_ratio: bool,
    n_boot: int = N_BOOT,
) -> dict[str, tuple[float, float]]:
    """Vectorized stratified bootstrap for all reported two-sample effects."""
    ix = RNG.integers(0, len(x), size=(n_boot, len(x)))
    iy = RNG.integers(0, len(y), size=(n_boot, len(y)))
    xb, yb = x[ix], y[iy]
    med_x, med_y = np.median(xb, axis=1), np.median(yb, axis=1)
    med_diff = med_x - med_y

    nx, ny = len(x), len(y)
    vx, vy = np.var(xb, axis=1, ddof=1), np.var(yb, axis=1, ddof=1)
    pooled = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    correction = 1 - 3 / (4 * (nx + ny) - 9)
    g = correction * (np.mean(xb, axis=1) - np.mean(yb, axis=1)) / pooled

    delta = np.empty(n_boot, dtype=float)
    batch = 100
    for start in range(0, n_boot, batch):
        stop = min(start + batch, n_boot)
        a, b = xb[start:stop], yb[start:stop]
        gt = (a[:, :, None] > b[:, None, :]).mean(axis=(1, 2))
        lt = (a[:, :, None] < b[:, None, :]).mean(axis=(1, 2))
        delta[start:stop] = gt - lt

    def ci(values: np.ndarray) -> tuple[float, float]:
        finite = values[np.isfinite(values)]
        return tuple(np.percentile(finite, [2.5, 97.5])) if len(finite) else (np.nan, np.nan)

    out = {
        "median_difference": ci(med_diff),
        "hedges_g": ci(g),
        "cliffs_delta": ci(delta),
        "auc": ci((delta + 1) / 2),
    }
    if include_ratio:
        ratio = np.divide(med_x, med_y, out=np.full_like(med_x, np.nan), where=med_y > 0)
        out["median_ratio"] = ci(ratio)
    return out


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    u = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").statistic
    return float(2 * u / (len(x) * len(y)) - 1)


def hedges_g(x: np.ndarray, y: np.ndarray) -> float:
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    pooled = math.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    if pooled == 0:
        return np.nan
    d = (np.mean(x) - np.mean(y)) / pooled
    correction = 1 - 3 / (4 * (nx + ny) - 9)
    return float(correction * d)


def half_min_positive(values: pd.Series) -> float:
    positive = pd.to_numeric(values, errors="coerce")
    positive = positive[positive > 0]
    return float(positive.min() / 2) if len(positive) else 0.5


def analyze_contrasts(
    data: pd.DataFrame,
    cohort: str,
    contrasts: list[tuple[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    results = []
    missingness = []
    sensitivity = []
    is_discovery = cohort == "discovery"

    for protein, pdf in data.groupby("protein"):
        pseudo = None if is_discovery else half_min_positive(pdf.value)
        for group1, group2 in contrasts:
            d1 = pdf[pdf.group == group1]
            d2 = pdf[pdf.group == group2]
            x = d1.value.dropna().to_numpy(float)
            y = d2.value.dropna().to_numpy(float)
            if len(x) < 2 or len(y) < 2:
                continue

            mw = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
            delta = cliffs_delta(x, y)
            raw_boot = bootstrap_effect_bundle(x, y, include_ratio=not is_discovery)
            delta_lo, delta_hi = raw_boot["cliffs_delta"]
            med_diff = float(np.median(x) - np.median(y))
            md_lo, md_hi = raw_boot["median_difference"]

            if is_discovery:
                xt, yt = x, y
                transform = "source VSN scale"
                median_ratio = median_ratio_lo = median_ratio_hi = np.nan
            else:
                xt = np.log2(x + pseudo)
                yt = np.log2(y + pseudo)
                transform = f"log2(value + {pseudo:.8g})"
                denom = float(np.median(y))
                median_ratio = float(np.median(x) / denom) if denom > 0 else np.nan
                median_ratio_lo, median_ratio_hi = raw_boot["median_ratio"]

            welch = stats.ttest_ind(xt, yt, equal_var=False, nan_policy="omit")
            g = hedges_g(xt, yt)
            transformed_boot = (
                raw_boot if is_discovery
                else bootstrap_effect_bundle(xt, yt, include_ratio=False)
            )
            g_lo, g_hi = transformed_boot["hedges_g"]

            labels = np.r_[np.ones(len(x)), np.zeros(len(y))]
            scores = np.r_[x, y]
            auc = float(roc_auc_score(labels, scores))
            auc_lo, auc_hi = raw_boot["auc"]

            results.append(
                {
                    "cohort": cohort,
                    "protein": protein,
                    "contrast": f"{group1}_vs_{group2}",
                    "group1": group1,
                    "group2": group2,
                    "n1_source_rows": len(d1),
                    "n1_observed": len(x),
                    "n2_source_rows": len(d2),
                    "n2_observed": len(y),
                    "group1_median": float(np.median(x)),
                    "group1_q1": float(np.percentile(x, 25)),
                    "group1_q3": float(np.percentile(x, 75)),
                    "group2_median": float(np.median(y)),
                    "group2_q1": float(np.percentile(y, 25)),
                    "group2_q3": float(np.percentile(y, 75)),
                    "median_difference_group1_minus_group2": med_diff,
                    "median_difference_ci_low": md_lo,
                    "median_difference_ci_high": md_hi,
                    "median_ratio_group1_over_group2": median_ratio,
                    "median_ratio_ci_low": median_ratio_lo,
                    "median_ratio_ci_high": median_ratio_hi,
                    "cliffs_delta_group1_minus_group2": delta,
                    "cliffs_delta_ci_low": delta_lo,
                    "cliffs_delta_ci_high": delta_hi,
                    "mann_whitney_p": float(mw.pvalue),
                    "secondary_transform": transform,
                    "hedges_g_group1_minus_group2": g,
                    "hedges_g_ci_low": g_lo,
                    "hedges_g_ci_high": g_hi,
                    "welch_t_p": float(welch.pvalue),
                    "roc_auc_group1_positive": auc,
                    "roc_auc_ci_low": auc_lo,
                    "roc_auc_ci_high": auc_hi,
                    "roc_auc_separation": max(auc, 1 - auc),
                    "higher_group": group1 if delta > 0 else group2 if delta < 0 else "tie",
                }
            )

            # Differential missingness/ND is itself a possible assay artifact.
            table = np.array(
                [
                    [int(d1.value.notna().sum()), int(d1.value.isna().sum())],
                    [int(d2.value.notna().sum()), int(d2.value.isna().sum())],
                ]
            )
            fisher = stats.fisher_exact(table, alternative="two-sided")
            missingness.append(
                {
                    "cohort": cohort,
                    "protein": protein,
                    "contrast": f"{group1}_vs_{group2}",
                    "group1_observed": int(table[0, 0]),
                    "group1_missing_or_ND": int(table[0, 1]),
                    "group2_observed": int(table[1, 0]),
                    "group2_missing_or_ND": int(table[1, 1]),
                    "fisher_odds_ratio_observed_vs_missing": float(fisher.statistic),
                    "fisher_p": float(fisher.pvalue),
                }
            )

            # For confirmatory data, evaluate a standard left-censor substitution.
            if not is_discovery:
                x_imp = d1.value.fillna(pseudo).to_numpy(float)
                y_imp = d2.value.fillna(pseudo).to_numpy(float)
                mw_imp = stats.mannwhitneyu(x_imp, y_imp, alternative="two-sided", method="asymptotic")
                sensitivity.append(
                    {
                        "cohort": cohort,
                        "protein": protein,
                        "contrast": f"{group1}_vs_{group2}",
                        "substitution_value_half_min_positive": pseudo,
                        "complete_case_cliffs_delta": delta,
                        "half_min_substitution_cliffs_delta": cliffs_delta(x_imp, y_imp),
                        "complete_case_mann_whitney_p": float(mw.pvalue),
                        "half_min_substitution_mann_whitney_p": float(mw_imp.pvalue),
                    }
                )

    result_df = pd.DataFrame(results)
    missing_df = pd.DataFrame(missingness)
    sens_df = pd.DataFrame(sensitivity)
    for _, idx in result_df.groupby(["cohort", "contrast"]).groups.items():
        result_df.loc[idx, "mann_whitney_q_bh"] = multipletests(
            result_df.loc[idx, "mann_whitney_p"], method="fdr_bh"
        )[1]
        result_df.loc[idx, "welch_t_q_bh"] = multipletests(
            result_df.loc[idx, "welch_t_p"], method="fdr_bh"
        )[1]
    for _, idx in missing_df.groupby(["cohort", "contrast"]).groups.items():
        missing_df.loc[idx, "fisher_q_bh"] = multipletests(
            missing_df.loc[idx, "fisher_p"], method="fdr_bh"
        )[1]
    if not sens_df.empty:
        for _, idx in sens_df.groupby(["cohort", "contrast"]).groups.items():
            sens_df.loc[idx, "half_min_substitution_q_bh"] = multipletests(
                sens_df.loc[idx, "half_min_substitution_mann_whitney_p"], method="fdr_bh"
            )[1]
    return result_df, missing_df, sens_df


def bootstrap_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    n = len(x)
    idx = RNG.integers(0, n, size=(N_BOOT, n))
    xr = stats.rankdata(x[idx], axis=1)
    yr = stats.rankdata(y[idx], axis=1)
    xr = xr - xr.mean(axis=1, keepdims=True)
    yr = yr - yr.mean(axis=1, keepdims=True)
    denom = np.sqrt(np.sum(xr * xr, axis=1) * np.sum(yr * yr, axis=1))
    values = np.divide(
        np.sum(xr * yr, axis=1), denom,
        out=np.full(N_BOOT, np.nan), where=denom > 0,
    )
    values = values[np.isfinite(values)]
    return tuple(np.percentile(values, [2.5, 97.5])) if len(values) else (np.nan, np.nan)


def analyze_alt_correlations(data: pd.DataFrame) -> pd.DataFrame:
    results = []
    for protein, pdf in data.groupby("protein"):
        for input_rule, biomarker_col in [
            ("primary_exclude_undocumented_numeric_ND_substitutions", "log2_biomarker_primary"),
            ("sensitivity_include_values_as_published", "log2_biomarker"),
        ]:
            for stratum in ["HV_plus_DO", "DO", "HV"]:
                sdf = pdf if stratum == "HV_plus_DO" else pdf[pdf.group == stratum]
                sdf = sdf.dropna(subset=[biomarker_col, "log2_alt"])
                if len(sdf) < 5:
                    continue
                corr = stats.spearmanr(sdf[biomarker_col], sdf.log2_alt)
                lo, hi = bootstrap_spearman(
                    sdf[biomarker_col].to_numpy(float), sdf.log2_alt.to_numpy(float)
                )
                results.append(
                    {
                        "protein": protein,
                        "correlation_input_rule": input_rule,
                        "stratum": stratum,
                        "n_complete_pairs": len(sdf),
                        "n_excluded_undocumented_numeric_ND_substitutions": int(
                            (pdf if stratum == "HV_plus_DO" else pdf[pdf.group == stratum])
                            .is_likely_numeric_ND_substitution.sum()
                        ),
                        "spearman_rho": float(corr.statistic),
                        "spearman_ci_low": lo,
                        "spearman_ci_high": hi,
                        "spearman_p": float(corr.pvalue),
                    }
                )
    out = pd.DataFrame(results)
    for _, idx in out.groupby(["correlation_input_rule", "stratum"]).groups.items():
        out.loc[idx, "spearman_q_bh"] = multipletests(out.loc[idx, "spearman_p"], method="fdr_bh")[1]
    return out


def analyze_zones(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries = (
        data.groupby(["zone", "group"])
        .score.agg(n="count", mean="mean", sd="std", median="median", q1=lambda x: x.quantile(.25), q3=lambda x: x.quantile(.75))
        .reset_index()
    )
    summaries["se"] = summaries.sd / np.sqrt(summaries.n)
    rows = []
    contrasts = [("DO", "HV"), ("NDO", "HV"), ("DO", "NDO"), ("DO", "DF"), ("NDO", "NDF")]
    for zone, zdf in data.groupby("zone"):
        for g1, g2 in contrasts:
            x = zdf.loc[zdf.group == g1, "score"].dropna().to_numpy(float)
            y = zdf.loc[zdf.group == g2, "score"].dropna().to_numpy(float)
            if len(x) < 2 or len(y) < 2:
                continue
            p = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").pvalue
            d = cliffs_delta(x, y)
            lo, hi = bootstrap_effect_bundle(x, y, include_ratio=False)["cliffs_delta"]
            rows.append(
                {
                    "zone": zone,
                    "contrast": f"{g1}_vs_{g2}",
                    "n1": len(x), "n2": len(y),
                    "cliffs_delta_group1_minus_group2": d,
                    "cliffs_delta_ci_low": lo,
                    "cliffs_delta_ci_high": hi,
                    "mann_whitney_p": float(p),
                }
            )
    effects = pd.DataFrame(rows)
    for _, idx in effects.groupby("contrast").groups.items():
        effects.loc[idx, "mann_whitney_q_bh"] = multipletests(
            effects.loc[idx, "mann_whitney_p"], method="fdr_bh"
        )[1]
    return summaries, effects


def cross_cohort_concordance(all_effects: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    details = []
    summaries = []
    for contrast in ["DO_vs_HV", "DO_vs_NDO", "DO_vs_DF"]:
        disc = all_effects.query("cohort == 'discovery' and contrast == @contrast")[
            ["protein", "cliffs_delta_group1_minus_group2"]
        ].rename(columns={"cliffs_delta_group1_minus_group2": "discovery_delta"})
        conf = all_effects.query("cohort == 'confirmatory' and contrast == @contrast")[
            ["protein", "cliffs_delta_group1_minus_group2"]
        ].rename(columns={"cliffs_delta_group1_minus_group2": "confirmatory_delta"})
        merged = disc.merge(conf, on="protein", how="inner")
        # HPD is excluded because its 50-cell discovery series is an exact
        # duplicate of ALDOB in the source workbook (unresolved source anomaly).
        merged = merged[merged.protein != "HPD"].copy()
        merged["contrast"] = contrast
        merged["direction_concordant"] = np.sign(merged.discovery_delta) == np.sign(merged.confirmatory_delta)
        details.append(merged)
        corr = stats.spearmanr(merged.discovery_delta, merged.confirmatory_delta)
        matches = int(merged.direction_concordant.sum())
        summaries.append(
            {
                "contrast": contrast,
                "n_proteins": len(merged),
                "n_direction_concordant": matches,
                "proportion_direction_concordant": matches / len(merged),
                "exact_binomial_p_vs_0.5": float(stats.binomtest(matches, len(merged), 0.5).pvalue),
                "spearman_rho_across_proteins": float(corr.statistic),
                "spearman_p": float(corr.pvalue),
            }
        )
    return pd.concat(details, ignore_index=True), pd.DataFrame(summaries)


def build_candidate_evidence(
    effects: pd.DataFrame,
    correlations: pd.DataFrame,
    concordance_detail: pd.DataFrame,
) -> pd.DataFrame:
    def select_effect(cohort: str, contrast: str) -> pd.DataFrame:
        cols = [
            "protein", "cliffs_delta_group1_minus_group2", "cliffs_delta_ci_low",
            "cliffs_delta_ci_high", "mann_whitney_q_bh", "roc_auc_group1_positive",
            "roc_auc_ci_low", "roc_auc_ci_high", "roc_auc_separation",
        ]
        subset = effects.query("cohort == @cohort and contrast == @contrast")[cols].copy()
        return subset.rename(columns={c: f"{cohort}_{contrast}_{c}" for c in cols if c != "protein"})

    out = pd.DataFrame({"protein": CANDIDATES})
    for cohort, contrast in [
        ("confirmatory", "DO_vs_HV"),
        ("confirmatory", "DO_vs_NDO"),
        ("confirmatory", "DO_vs_DF"),
        ("discovery", "DO_vs_HV"),
        ("discovery", "DO_vs_NDO"),
        ("discovery", "DO_vs_DF"),
    ]:
        out = out.merge(select_effect(cohort, contrast), on="protein", how="left")

    primary_correlations = correlations.query(
        "correlation_input_rule == 'primary_exclude_undocumented_numeric_ND_substitutions'"
    )
    pooled_corr = primary_correlations.query("stratum == 'HV_plus_DO'")[
        ["protein", "spearman_rho", "spearman_ci_low", "spearman_ci_high", "spearman_q_bh"]
    ].rename(columns={c: f"pooled_alt_{c}" for c in ["spearman_rho", "spearman_ci_low", "spearman_ci_high", "spearman_q_bh"]})
    do_corr = primary_correlations.query("stratum == 'DO'")[
        ["protein", "spearman_rho", "spearman_ci_low", "spearman_ci_high", "spearman_q_bh"]
    ].rename(columns={c: f"within_DO_alt_{c}" for c in ["spearman_rho", "spearman_ci_low", "spearman_ci_high", "spearman_q_bh"]})
    out = out.merge(pooled_corr, on="protein", how="left").merge(do_corr, on="protein", how="left")

    for contrast in ["DO_vs_HV", "DO_vs_NDO", "DO_vs_DF"]:
        cc = concordance_detail.query("contrast == @contrast")[["protein", "direction_concordant"]]
        cc = cc.rename(columns={"direction_concordant": f"{contrast}_discovery_confirmatory_direction_concordant"})
        out = out.merge(cc, on="protein", how="left")

    q_spec = out["confirmatory_DO_vs_NDO_mann_whitney_q_bh"]
    d_spec = out["confirmatory_DO_vs_NDO_cliffs_delta_group1_minus_group2"].abs()
    q_detect = out["confirmatory_DO_vs_HV_mann_whitney_q_bh"]
    d_detect = out["confirmatory_DO_vs_HV_cliffs_delta_group1_minus_group2"].abs()
    q_recovery = out["confirmatory_DO_vs_DF_mann_whitney_q_bh"]
    d_recovery = out["confirmatory_DO_vs_DF_cliffs_delta_group1_minus_group2"].abs()

    out["dili_vs_non_dili_specificity_signal"] = (q_spec < 0.05) & (d_spec >= 0.20)
    out["acute_liver_injury_detection_signal"] = (q_detect < 0.05) & (d_detect >= 0.33)
    out["onset_vs_followup_recovery_signal"] = (q_recovery < 0.05) & (d_recovery >= 0.33)
    out["specificity_replicates_discovery_direction"] = out[
        "DO_vs_NDO_discovery_confirmatory_direction_concordant"
    ].fillna(False)
    out["within_DO_severity_correlation_signal"] = (
        (out["within_DO_alt_spearman_q_bh"] < 0.05) &
        (out["within_DO_alt_spearman_rho"].abs() >= 0.30)
    )

    def tier(row: pd.Series) -> str:
        if (
            row.dili_vs_non_dili_specificity_signal
            and row.specificity_replicates_discovery_direction
            and row.onset_vs_followup_recovery_signal
        ):
            return "A: etiology-discriminating, discovery-concordant, recovery-associated"
        if row.dili_vs_non_dili_specificity_signal and row.specificity_replicates_discovery_direction:
            return "B: etiology-discriminating and discovery-concordant"
        if row.acute_liver_injury_detection_signal and row.onset_vs_followup_recovery_signal:
            return "C: general injury/recovery marker, not etiology-specific"
        return "D: limited or inconsistent evidence in available source data"

    out["evidence_tier"] = out.apply(tier, axis=1)
    out["interpretation_boundary"] = (
        "Circulating protein evidence anchors a downstream injury phenotype; it does not establish "
        "an upstream causal susceptibility target or nintedanib specificity."
    )
    out["discovery_source_integrity_flag"] = np.where(
        out.protein == "HPD",
        "HPD discovery series is an exact duplicate of ALDOB; discovery concordance not credited",
        "none identified for this protein",
    )
    return out


def discovery_duplicate_series_qc(discovery: pd.DataFrame) -> pd.DataFrame:
    wide = discovery.pivot_table(
        index=["group", "source_row_index_within_group"], columns="protein", values="value",
        aggfunc="first", dropna=False,
    )
    rows = []
    proteins = sorted(wide.columns)
    for i, p1 in enumerate(proteins):
        for p2 in proteins[i + 1:]:
            a, b = wide[p1].to_numpy(float), wide[p2].to_numpy(float)
            same_missing = np.array_equal(np.isnan(a), np.isnan(b))
            finite = np.isfinite(a) & np.isfinite(b)
            exact = same_missing and finite.any() and np.array_equal(a[finite], b[finite])
            if exact:
                rows.append(
                    {
                        "protein_1": p1,
                        "protein_2": p2,
                        "n_positions": len(a),
                        "n_joint_numeric": int(finite.sum()),
                        "exact_duplicate_including_missingness": True,
                        "interpretation": "Probable source workbook copy/label anomaly; do not count as independent discovery evidence",
                    }
                )
    return pd.DataFrame(rows)


def pck2_outlier_sensitivity(confirmatory: pd.DataFrame) -> pd.DataFrame:
    pdf = confirmatory[confirmatory.protein == "PCK2"].copy()
    rows = []
    for rule, filtered in [
        ("retain_starred_3077.0266_value", pdf),
        ("exclude_starred_3077.0266_value", pdf[~pdf.is_starred_source_value]),
    ]:
        for g1, g2 in [("DO", "HV"), ("DO", "NDO"), ("DO", "DF")]:
            x = filtered.loc[filtered.group == g1, "value"].dropna().to_numpy(float)
            y = filtered.loc[filtered.group == g2, "value"].dropna().to_numpy(float)
            if len(x) < 2 or len(y) < 2:
                continue
            pseudo = half_min_positive(pd.concat([
                filtered.loc[filtered.group == g1, "value"],
                filtered.loc[filtered.group == g2, "value"],
            ]))
            xt, yt = np.log2(x + pseudo), np.log2(y + pseudo)
            mw = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
            welch = stats.ttest_ind(xt, yt, equal_var=False)
            rows.append(
                {
                    "analysis_rule": rule,
                    "contrast": f"{g1}_vs_{g2}",
                    "n1": len(x), "n2": len(y),
                    "group1_median": float(np.median(x)),
                    "group2_median": float(np.median(y)),
                    "cliffs_delta_group1_minus_group2": cliffs_delta(x, y),
                    "mann_whitney_p": float(mw.pvalue),
                    "hedges_g_log2": hedges_g(xt, yt),
                    "welch_t_p_log2": float(welch.pvalue),
                }
            )
    return pd.DataFrame(rows)


def row_order_linkage_qc(confirmatory: pd.DataFrame, correlations: pd.DataFrame) -> pd.DataFrame:
    """Demonstrate why protein blocks cannot be positionally joined for ML."""
    rows = []
    for protein in sorted(set(CANDIDATES) & set(correlations.protein)):
        for group in ["HV", "DO"]:
            box = confirmatory.query("protein == @protein and group == @group").sort_values(
                "source_row_index_within_group"
            ).value.to_numpy(float)
            corr = correlations.query("protein == @protein and group == @group").sort_values(
                "source_row_index_within_group_and_protein"
            ).log2_biomarker.to_numpy(float)
            box_log = np.full_like(box, np.nan, dtype=float)
            positive = box > 0
            box_log[positive] = np.log2(box[positive])
            finite = np.isfinite(box_log) & np.isfinite(corr)
            order_equal = bool(np.allclose(box_log[finite], corr[finite], atol=1e-9, rtol=1e-9)) if finite.any() else False
            sorted_box = np.sort(box_log[np.isfinite(box_log)])
            sorted_corr = np.sort(corr[np.isfinite(corr)])
            same_multiset = (
                len(sorted_box) == len(sorted_corr)
                and np.allclose(sorted_box, sorted_corr, atol=1e-9, rtol=1e-9)
            )
            rows.append(
                {
                    "protein": protein,
                    "group": group,
                    "n_boxplot_values": len(box),
                    "n_correlation_values": len(corr),
                    "same_value_multiset_after_log2": bool(same_multiset),
                    "same_row_order_after_log2": order_equal,
                    "participant_linkage_interpretation": (
                        "Same values but different order; row number is not a stable participant identifier"
                        if same_multiset and not order_equal
                        else "See source-level QC"
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_data_dictionary() -> pd.DataFrame:
    rows = [
        ("discovery_tidy.csv", "cohort", "Discovery cohort; source values are TMT VSN-normalized relative intensities."),
        ("discovery_tidy.csv", "protein", "HGNC-style protein/gene symbol supplied by the source workbook."),
        ("discovery_tidy.csv", "group", "HV healthy volunteer; DO DILI onset; DF DILI follow-up; NDO non-DILI acute liver injury onset; NDF non-DILI follow-up; NAFLD chronic control."),
        ("discovery_tidy.csv", "source_row_index_within_group", "Sequential position inside one protein/group block; not a participant identifier across proteins."),
        ("discovery_tidy.csv", "value", "Numeric source value; ND and blank entries are represented as missing."),
        ("discovery_tidy.csv", "is_not_detected", "True only when the workbook explicitly labels the cell ND/not detected."),
        ("*_tidy.csv", "source_value_text", "Exact trimmed cell text retained for auditability."),
        ("*_tidy.csv", "is_starred_source_value", "True for a value carrying an unexplained source asterisk, currently the PCK2 3077.0266 value."),
        ("*_tidy.csv", "source_sheet/source_excel_row/source_excel_column", "One-based workbook provenance for the parsed cell."),
        ("confirmatory_tidy.csv", "assay", "SureQuant light/heavy ratio for most candidates; PCK2 was measured by ELISA and workbook units are unstated."),
        ("confirmatory_tidy.csv", "analysis_sample_key", "Convenience key for a protein-specific row series. It must not be used to join proteins."),
        ("confirmatory_tidy.csv", "is_invalid_nonnumeric_text", "Nonblank source text that could not be parsed as numeric and was not explicitly ND."),
        ("alt_correlation_pairs.csv", "log2_biomarker", "Workbook-provided log2 biomarker value paired with ALT within a single protein scatterplot block."),
        ("alt_correlation_pairs.csv", "log2_alt", "Workbook-provided log2 ALT value paired with the biomarker in the same row."),
        ("effect_estimates.csv", "cliffs_delta_group1_minus_group2", "P(group1 > group2) minus P(group1 < group2); range -1 to 1."),
        ("effect_estimates.csv", "n1_source_rows/n2_source_rows", "Rows present in the selected source series, including missing/ND rows."),
        ("effect_estimates.csv", "n1_observed/n2_observed", "Finite numeric observations included in the complete-case primary test."),
        ("effect_estimates.csv", "group*_median/q1/q3", "Group-specific median and interquartile range on the original workbook scale."),
        ("effect_estimates.csv", "median_difference_*", "Group1 median minus group2 median, with percentile bootstrap CI."),
        ("effect_estimates.csv", "median_ratio_*", "Group1 median divided by group2 median for confirmatory positive-scale values; not defined for discovery VSN values."),
        ("effect_estimates.csv", "mann_whitney_q_bh", "Benjamini-Hochberg FDR within cohort and contrast across candidate proteins."),
        ("effect_estimates.csv", "hedges_g_*", "Bias-corrected standardized mean difference on the documented secondary transform, with bootstrap CI."),
        ("effect_estimates.csv", "welch_t_p/welch_t_q_bh", "Secondary Welch-test p value and within-contrast BH FDR."),
        ("effect_estimates.csv", "roc_auc_group1_positive", "Univariate apparent AUC with group1 as the positive class; values below 0.5 indicate lower values in group1."),
        ("effect_estimates.csv", "roc_auc_separation", "max(AUC, 1-AUC), reported only as direction-free separation magnitude."),
        ("left_censor_sensitivity.csv", "substitution_value_half_min_positive", "One-half of the smallest positive value for that protein, used only in sensitivity analysis."),
        ("missingness_stats.csv", "fisher_*", "Fisher exact test of observed versus missing/ND proportions and BH FDR."),
        ("alt_correlation_pairs.csv", "is_likely_numeric_ND_substitution", "QC flag for undocumented numeric replacements of canonical ND values in Supp Fig 7."),
        ("alt_correlation_pairs.csv", "log2_biomarker_primary", "Biomarker value after setting undocumented numeric ND substitutions to missing."),
        ("correlation_stats.csv", "stratum", "HV_plus_DO is pooled and susceptible to group separation; DO and HV are within-group sensitivity analyses."),
        ("correlation_stats.csv", "correlation_input_rule", "Whether undocumented numeric ND substitutions were excluded (primary) or retained (sensitivity)."),
        ("correlation_stats.csv", "spearman_*", "Spearman rho, bootstrap CI, two-sided p value, and BH FDR within input rule and stratum."),
        ("candidate_evidence_table.csv", "evidence_tier", "Rule-based descriptive tier from FDR-controlled univariate evidence, not a causal target rank."),
        ("zonal_scores_tidy.csv", "score", "Workbook-provided ssGSEA pathway/zone signature score."),
        ("source_block_counts.csv", "declared_n", "Sample count printed in a block header, where available."),
        ("source_block_counts.csv", "rows_through_last_nonblank", "Rows from block header through the last supplied cell; retains internal blanks."),
        ("row_order_linkage_qc.csv", "same_value_multiset_after_log2", "Whether boxplot and ALT-correlation blocks contain the same numeric values regardless of order."),
        ("row_order_linkage_qc.csv", "same_row_order_after_log2", "Whether those same values occur in the same order; failures demonstrate unsafe positional joining."),
    ]
    return pd.DataFrame(rows, columns=["file", "field", "definition"])


def audit_column_blocks(sheet: str, marker: str) -> pd.DataFrame:
    df = pd.read_excel(SOURCE, sheet_name=sheet, header=None, dtype=object)
    candidate_starts = marker_rows(df, marker)
    starts = []
    for start in candidate_starts:
        try:
            find_header_row(df, start, min(len(df), start + 10))
            starts.append(start)
        except ValueError:
            continue
    rows = []
    for k, start in enumerate(starts):
        end = starts[k + 1] if k + 1 < len(starts) else len(df)
        header = find_header_row(df, start, end)
        protein = find_protein(df, start, header)
        for col in range(1, df.shape[1]):
            header_text = clean_text(df.iat[header, col])
            group = parse_group(header_text)
            if group not in set(DISCOVERY_N) | set(CONFIRMATORY_N):
                continue
            declared_match = re.search(r"n\s*=\s*(\d+)", header_text, flags=re.I)
            declared_n = int(declared_match.group(1)) if declared_match else np.nan
            raw = [clean_text(df.iat[r, col]) for r in range(header + 1, end)]
            populated = [j for j, text in enumerate(raw) if text]
            last = max(populated) if populated else -1
            within = raw[:last + 1]
            parsed = [parse_value(x) for x in within]
            rows.append(
                {
                    "source_sheet": sheet,
                    "protein": protein,
                    "group": group,
                    "declared_n": declared_n,
                    "rows_through_last_nonblank": len(within),
                    "n_nonblank": sum(bool(x) for x in within),
                    "n_numeric": sum(np.isfinite(x[0]) for x in parsed),
                    "n_explicit_ND": sum(x[1] for x in parsed),
                    "n_invalid_text": sum(invalid_text_flag(x[0], x[1], x[3]) for x in parsed),
                    "n_internal_blank": sum(not bool(x) for x in within),
                    "declared_vs_rows_difference": (
                        len(within) - declared_n if np.isfinite(declared_n) else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def qc_summary(
    discovery: pd.DataFrame,
    confirmatory: pd.DataFrame,
    correlations: pd.DataFrame,
    row_order_qc: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for name, data in [("discovery", discovery), ("confirmatory", confirmatory)]:
        for protein, pdf in data.groupby("protein"):
            for group, gdf in pdf.groupby("group"):
                rows.append(
                    {
                        "dataset": name,
                        "protein": protein,
                        "group": group,
                        "n_rows": len(gdf),
                        "n_numeric": int(gdf.value.notna().sum()),
                        "n_missing_total": int(gdf.value.isna().sum()),
                        "n_explicit_ND": int(gdf.is_not_detected.sum()),
                        "n_starred": int(gdf.is_starred_source_value.sum()),
                        "minimum": float(gdf.value.min()) if gdf.value.notna().any() else np.nan,
                        "maximum": float(gdf.value.max()) if gdf.value.notna().any() else np.nan,
                    }
                )
    out = pd.DataFrame(rows)
    linkage = {
        "correlation_blocks_checked": int(len(row_order_qc)),
        "same_multiset_but_different_order": int(
            (row_order_qc.same_value_multiset_after_log2 & ~row_order_qc.same_row_order_after_log2).sum()
        ),
    }
    (OUT / "linkage_qc_summary.json").write_text(json.dumps(linkage, indent=2), encoding="utf-8")
    return out


def save_csv(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, index=False, na_rep="NA")


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)
    OUT.mkdir(parents=True, exist_ok=True)

    discovery = parse_discovery()
    confirmatory, repeated_qc = assemble_confirmatory()
    zones, zone_duplicate_qc = parse_zonal_scores()
    correlations = parse_alt_correlations()
    row_order_qc = row_order_linkage_qc(confirmatory, correlations)
    source_block_counts = pd.concat(
        [
            audit_column_blocks("Figure 2C", "Fig 2C"),
            audit_column_blocks("Figure 3", "Fig 3"),
            audit_column_blocks("Supp Fig 2", "Supp Fig 2"),
            audit_column_blocks("Supp Fig 3-PCK2", "Supp Fig 3-PCK2"),
        ],
        ignore_index=True,
    )
    discovery_duplicates = discovery_duplicate_series_qc(discovery)

    disc_effects, disc_missing, _ = analyze_contrasts(
        discovery, "discovery", CONTRASTS_DISCOVERY
    )
    conf_effects, conf_missing, conf_sensitivity = analyze_contrasts(
        confirmatory, "confirmatory", CONTRASTS_CONFIRMATORY
    )
    effects = pd.concat([disc_effects, conf_effects], ignore_index=True)
    missingness = pd.concat([disc_missing, conf_missing], ignore_index=True)
    correlation_stats = analyze_alt_correlations(correlations)
    zone_summary, zone_effects = analyze_zones(zones)
    concordance_detail, concordance_summary = cross_cohort_concordance(effects)
    evidence = build_candidate_evidence(effects, correlation_stats, concordance_detail)
    pck2_sensitivity = pck2_outlier_sensitivity(confirmatory)
    qc = qc_summary(discovery, confirmatory, correlations, row_order_qc)

    # Protein-specific source values are long-format only.  A wide participant
    # matrix is intentionally not emitted because cross-protein linkage is absent.
    model_feasibility = pd.DataFrame(
        [
            {
                "proposed_analysis": "Multivariable or nested-CV classifier across proteins",
                "status": "Not performed",
                "reason": (
                    "The workbook provides individual case values within each protein block but no common "
                    "participant identifier. Cross-figure QC demonstrates that protein-specific source rows "
                    "are reordered. Positional joining would fabricate participant-level linkage."
                ),
                "minimum_additional_data_needed": (
                    "De-identified participant ID linked across targeted protein measurements and outcomes, "
                    "or an official participant-by-protein matrix from PXD034882/Panorama."
                ),
            }
        ]
    )

    save_csv(discovery, "discovery_tidy.csv")
    save_csv(confirmatory, "confirmatory_tidy.csv")
    save_csv(zones, "zonal_scores_tidy.csv")
    save_csv(correlations, "alt_correlation_pairs.csv")
    save_csv(effects, "effect_estimates.csv")
    save_csv(missingness, "missingness_stats.csv")
    save_csv(conf_sensitivity, "left_censor_sensitivity.csv")
    save_csv(correlation_stats, "correlation_stats.csv")
    save_csv(zone_summary, "zonal_score_summary.csv")
    save_csv(zone_effects, "zonal_effect_estimates.csv")
    save_csv(concordance_detail, "cross_cohort_concordance_detail.csv")
    save_csv(concordance_summary, "cross_cohort_concordance_summary.csv")
    save_csv(evidence, "candidate_evidence_table.csv")
    save_csv(qc, "qc_counts.csv")
    save_csv(repeated_qc, "repeated_series_qc.csv")
    save_csv(zone_duplicate_qc, "zonal_duplicate_qc.csv")
    save_csv(row_order_qc, "row_order_linkage_qc.csv")
    save_csv(source_block_counts, "source_block_counts.csv")
    save_csv(discovery_duplicates, "discovery_duplicate_series_qc.csv")
    save_csv(pck2_sensitivity, "pck2_outlier_sensitivity.csv")
    save_csv(model_feasibility, "multivariable_model_feasibility.csv")
    save_csv(build_data_dictionary(), "data_dictionary.csv")

    manifest = {
        "analysis": "Fresh re-analysis of Nature Communications 2023 DILI source-data workbook",
        "source_file": SOURCE.name,
        "source_sha256": sha256(SOURCE),
        "source_article_doi": "10.1038/s41467-023-36858-6",
        "random_seed": SEED,
        "bootstrap_replicates": N_BOOT,
        "python": sys.version,
        "platform": platform.platform(),
        "package_versions": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit-learn": importlib.metadata.version("scikit-learn"),
            "statsmodels": importlib.metadata.version("statsmodels"),
            "openpyxl": importlib.metadata.version("openpyxl"),
        },
        "key_guardrail": (
            "No cross-protein positional join and no multivariable model, because the source workbook "
            "does not provide stable participant linkage."
        ),
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Wrote {len(list(OUT.glob('*.csv')))} CSV files to {OUT}")
    print(evidence[[
        "protein", "dili_vs_non_dili_specificity_signal",
        "acute_liver_injury_detection_signal", "onset_vs_followup_recovery_signal",
        "evidence_tier",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
