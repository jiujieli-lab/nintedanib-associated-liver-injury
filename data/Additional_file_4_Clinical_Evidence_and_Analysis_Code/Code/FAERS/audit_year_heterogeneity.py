#!/usr/bin/env python3
"""Audit calendar-year heterogeneity in the frozen FAERS 2x2 tables.

The Tarone-adjusted Breslow-Day test evaluates equality of odds ratios across
estimable annual strata. Cochran's Q is calculated from inverse-variance
weighted annual log reporting odds ratios (RORs); a 0.5 continuity correction
is applied to all four cells only for an otherwise estimable zero-cell stratum,
matching the frozen annual ROR calculation. Double-zero endpoint strata are
non-estimable and are excluded from both tests.

These tests determine whether the Mantel-Haenszel value can be interpreted as a
stable common effect. They do not convert spontaneous-report disproportionality
into incidence or a causal effect.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.contingency_tables import StratifiedTable


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
DEFAULT_INPUT = PACKAGE_ROOT / "Source_Data" / "FAERS" / "annual_active_comparator_signals.csv"
DEFAULT_OUTPUT = PACKAGE_ROOT / "Source_Data" / "FAERS" / "year_effect_heterogeneity_tests.csv"

CELL_COLUMNS = (
    "nintedanib_event",
    "nintedanib_nonevent",
    "pirfenidone_event",
    "pirfenidone_nonevent",
)


def _estimable_strata(frame: pd.DataFrame) -> pd.DataFrame:
    """Return annual strata with an identifiable exposure-event contrast."""
    cells = frame.loc[:, CELL_COLUMNS].apply(pd.to_numeric, errors="raise")
    if (cells < 0).any().any():
        raise ValueError("Annual 2x2 cells must be non-negative")
    # A year with no endpoint report in either comparator arm has no odds-ratio
    # information. Retain single-zero strata; Breslow-Day supports them and the
    # log-ROR calculation below uses the prespecified 0.5 correction.
    keep = cells["nintedanib_event"].add(cells["pirfenidone_event"]).gt(0)
    return frame.loc[keep].copy()


def _tarone_breslow_day(frame: pd.DataFrame) -> tuple[float, int, float]:
    tables = np.stack(
        [
            np.asarray(
                [
                    [row.nintedanib_event, row.nintedanib_nonevent],
                    [row.pirfenidone_event, row.pirfenidone_nonevent],
                ],
                dtype=float,
            )
            for row in frame.itertuples(index=False)
        ],
        axis=2,
    )
    result = StratifiedTable(tables, shift_zeros=False).test_equal_odds(adjust=True)
    return float(result.statistic), int(len(frame) - 1), float(result.pvalue)


def _cochran_q(frame: pd.DataFrame) -> tuple[float, int, float, float, int]:
    log_rors: list[float] = []
    variances: list[float] = []
    corrected = 0
    for row in frame.itertuples(index=False):
        cells = np.asarray(
            [
                row.nintedanib_event,
                row.nintedanib_nonevent,
                row.pirfenidone_event,
                row.pirfenidone_nonevent,
            ],
            dtype=float,
        )
        if np.any(cells == 0):
            cells = cells + 0.5
            corrected += 1
        a, b, c, d = cells
        log_rors.append(float(np.log((a * d) / (b * c))))
        variances.append(float(np.sum(1.0 / cells)))
    weights = 1.0 / np.asarray(variances)
    effects = np.asarray(log_rors)
    common = float(np.sum(weights * effects) / np.sum(weights))
    statistic = float(np.sum(weights * np.square(effects - common)))
    df = int(len(frame) - 1)
    p_value = float(chi2.sf(statistic, df))
    i2 = float(max(0.0, (statistic - df) / statistic) * 100.0) if statistic > 0 else 0.0
    return statistic, df, p_value, i2, corrected


def compute_year_effect_heterogeneity(annual: pd.DataFrame) -> pd.DataFrame:
    required = {"year", "endpoint", *CELL_COLUMNS}
    missing = sorted(required.difference(annual.columns))
    if missing:
        raise ValueError(f"Annual table is missing columns: {missing}")

    rows: list[dict[str, object]] = []
    windows = {
        "all_years": lambda frame: frame,
        "exclude_2022": lambda frame: frame.loc[frame["year"].ne(2022)],
    }
    for endpoint in ("narrow", "broad"):
        endpoint_frame = annual.loc[annual["endpoint"].eq(endpoint)].copy()
        for window, select in windows.items():
            selected = select(endpoint_frame)
            estimable = _estimable_strata(selected)
            bd_stat, bd_df, bd_p = _tarone_breslow_day(estimable)
            rows.append(
                {
                    "endpoint": endpoint,
                    "window": window,
                    "test": "Tarone-adjusted Breslow-Day",
                    "statistic": bd_stat,
                    "degrees_of_freedom": bd_df,
                    "p_value": bd_p,
                    "i2_percent": np.nan,
                    "calendar_strata_total": int(len(selected)),
                    "estimable_strata": int(len(estimable)),
                    "continuity_corrected_strata": 0,
                    "interpretation": "tests equality of annual odds ratios",
                }
            )
            q_stat, q_df, q_p, i2, corrected = _cochran_q(estimable)
            rows.append(
                {
                    "endpoint": endpoint,
                    "window": window,
                    "test": "Inverse-variance Cochran Q",
                    "statistic": q_stat,
                    "degrees_of_freedom": q_df,
                    "p_value": q_p,
                    "i2_percent": i2,
                    "calendar_strata_total": int(len(selected)),
                    "estimable_strata": int(len(estimable)),
                    "continuity_corrected_strata": corrected,
                    "interpretation": "heterogeneity of estimable annual log-RORs",
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    annual = pd.read_csv(args.input)
    result = compute_year_effect_heterogeneity(annual)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
