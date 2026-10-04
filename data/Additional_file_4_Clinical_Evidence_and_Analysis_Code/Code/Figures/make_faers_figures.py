#!/usr/bin/env python3
"""Validate the frozen openFDA/FAERS results and render publication figures.

The main figure contains exactly ten data panels (A-J) and no workflow panel.
All proportions use spontaneous-report denominators and are not incidence rates.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAERS = ROOT / "Code" / "FAERS"
RESULTS = ROOT / "Source_Data" / "FAERS"
MAIN_FIGURES = ROOT / "Figures"
SUPPLEMENTARY_FIGURES = ROOT / "Supplementary_Figures"
PANEL_DATA = ROOT / "Source_Data" / "Figure_Panels"
CHECKLISTS = ROOT / "Checklists"

NINT = "#D55E00"       # Okabe-Ito vermillion
PIRF = "#0072B2"       # Okabe-Ito blue
NARROW = "#CC79A7"     # reddish purple
BROAD = "#009E73"      # bluish green
PROXY = "#6F6F6F"
GOLD = "#E69F00"
SKY = "#56B4E9"
BLACK = "#222222"
MID_GREY = "#777777"
LIGHT_GREY = "#D9D9D9"


mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.0,
    "axes.titlesize": 7.8,
    "axes.titleweight": "bold",
    "axes.labelsize": 7.0,
    "axes.linewidth": 0.6,
    "xtick.labelsize": 6.4,
    "ytick.labelsize": 6.4,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.fontsize": 6.3,
    "legend.frameon": False,
    "lines.linewidth": 1.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": "white",
})


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def ror_ci(a: int, n1: int, c: int, n0: int) -> tuple[float, float, float]:
    b, d = n1 - a, n0 - c
    if a + c == 0:
        return np.nan, np.nan, np.nan
    vals = [a, b, c, d]
    if any(v == 0 for v in vals):
        aa, bb, cc, dd = [v + 0.5 for v in vals]
    else:
        aa, bb, cc, dd = vals
    log_est = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(sum(1 / v for v in (aa, bb, cc, dd)))
    return math.exp(log_est), math.exp(log_est - 1.96 * se), math.exp(log_est + 1.96 * se)


def style_axis(ax: mpl.axes.Axes, grid: str | None = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis=grid, color="#ECECEC", linewidth=0.5, zorder=0)
    ax.tick_params(direction="out", pad=1.5)


def panel_label(ax: mpl.axes.Axes, letter: str, x: float = -0.17, y: float = 1.10) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10.5, fontweight="bold",
            va="top", ha="left", color=BLACK)


def save_all(fig: mpl.figure.Figure, stem: str) -> list[Path]:
    destination = SUPPLEMENTARY_FIGURES if stem.startswith("Supplementary_") else MAIN_FIGURES
    paths = []
    for suffix in ("svg", "pdf"):
        path = destination / f"{stem}.{suffix}"
        fig.savefig(path, bbox_inches="tight")
        paths.append(path)
    png = destination / f"{stem}.png"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    paths.append(png)
    tiff = destination / f"{stem}.tiff"
    fig.savefig(tiff, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    paths.append(tiff)
    return paths


def validate_inputs(tables: dict[str, pd.DataFrame]) -> dict[str, object]:
    overall = tables["overall"]
    cases = tables["cases"]
    annual = tables["annual"]
    strata = tables["strata"]
    temporal = tables["temporal"]
    mh = tables["mh"]
    heterogeneity = tables["heterogeneity"]

    assert cases["safetyreportid"].astype(str).is_unique
    assert dict(cases.groupby("drug").size()) == {"nintedanib": 1773, "pirfenidone": 717}
    expected = {
        "narrow": (340, 30233, 112, 39467),
        "broad": (1773, 30233, 717, 39467),
        "hy_proxy": (19, 30233, 12, 39467),
    }
    checks = []
    for endpoint, counts in expected.items():
        row = overall.loc[overall["endpoint"].eq(endpoint)].iloc[0]
        observed = (int(row.nintedanib_event), int(row.nintedanib_total),
                    int(row.pirfenidone_event), int(row.pirfenidone_total))
        assert observed == counts
        est, low, high = ror_ci(*counts)
        assert np.allclose([row.ROR, row.ROR_low95, row.ROR_high95], [est, low, high], rtol=1e-12)
        checks.append({"endpoint": endpoint, "counts": observed, "ROR": est, "low95": low, "high95": high})

    # The years-only numerator must exactly use source ICH unit 801.
    unit = pd.to_numeric(cases["age_unit_code"], errors="coerce")
    age = pd.to_numeric(cases["age_value_reported"], errors="coerce")
    age_expected = {
        "18-64_years": (99, 359, 14, 71),
        "65plus_years": (193, 1088, 57, 302),
    }
    for band, (nn, nb, pn, pb) in age_expected.items():
        if band == "18-64_years":
            mask = unit.eq(801) & age.between(18, 64)
        else:
            mask = unit.eq(801) & age.between(65, 150)
        subset = cases[mask]
        observed = (
            int(((subset.drug == "nintedanib") & subset.narrow_case).sum()),
            int(((subset.drug == "nintedanib") & subset.broad_case).sum()),
            int(((subset.drug == "pirfenidone") & subset.narrow_case).sum()),
            int(((subset.drug == "pirfenidone") & subset.broad_case).sum()),
        )
        assert observed == (nn, nb, pn, pb)
        for endpoint, ne, pe in (("narrow", nn, pn), ("broad", nb, pb)):
            row = strata.query("dimension == 'age' and stratum == @band and endpoint == @endpoint").iloc[0]
            assert int(row.nintedanib_event) == ne and int(row.pirfenidone_event) == pe

    # Annual 2x2 cells must reconstitute the full-window table.
    for endpoint in ("narrow", "broad"):
        d = annual[annual.endpoint.eq(endpoint)]
        row = overall[overall.endpoint.eq(endpoint)].iloc[0]
        assert int(d.nintedanib_event.sum()) == int(row.nintedanib_event)
        assert int(d.pirfenidone_event.sum()) == int(row.pirfenidone_event)
        assert int(d.nintedanib_total.sum()) == int(row.nintedanib_total)
        assert int(d.pirfenidone_total.sum()) == int(row.pirfenidone_total)

    # The double-zero 2014 narrow contrast must remain explicitly non-estimable.
    zero = annual.query("year == 2014 and endpoint == 'narrow'").iloc[0]
    assert zero.nintedanib_event == 0 and zero.pirfenidone_event == 0 and pd.isna(zero.ROR)

    # Independent checks of the 2022 influence estimates.
    excluded = temporal.query("endpoint == 'narrow' and temporal_window == 'exclude_2022_posthoc'").iloc[0]
    assert np.allclose(
        [excluded.ROR, excluded.ROR_low95, excluded.ROR_high95],
        ror_ci(int(excluded.nintedanib_event), int(excluded.nintedanib_total),
               int(excluded.pirfenidone_event), int(excluded.pirfenidone_total)),
        rtol=1e-12,
    )
    post = temporal.query("endpoint == 'narrow' and temporal_window == 'post_2022'").iloc[0]
    assert np.allclose(
        [post.ROR, post.ROR_low95, post.ROR_high95],
        ror_ci(int(post.nintedanib_event), int(post.nintedanib_total),
               int(post.pirfenidone_event), int(post.pirfenidone_total)),
        rtol=1e-12,
    )

    # Canonical year-adjusted MH values and RBG confidence intervals.
    mh_expect = {
        "narrow": (4.511587193291437, 3.6056881245985086, 5.645085847500401),
        "broad": (4.048094766549167, 3.6866077228046237, 4.445027101092307),
    }
    for endpoint, vals in mh_expect.items():
        row = mh[mh.endpoint.eq(endpoint)].iloc[0]
        assert np.allclose([row.MH_ROR_year_adjusted, row.MH_low95_RBG, row.MH_high95_RBG], vals, rtol=1e-12)

    # A common annual effect is rejected for both endpoints, even after 2022 is
    # removed. Therefore the MH values are descriptive year-stratified summaries.
    heterogeneity_expect = {
        ("narrow", "all_years", "Tarone-adjusted Breslow-Day"): (88.25434621529466, 11, 3.6637359812630166e-14, np.nan),
        ("narrow", "exclude_2022", "Tarone-adjusted Breslow-Day"): (33.92521640584324, 10, 1.9016298099105988e-4, np.nan),
        ("broad", "all_years", "Tarone-adjusted Breslow-Day"): (63.325369250923984, 12, 5.552200810221564e-9, np.nan),
        ("broad", "exclude_2022", "Tarone-adjusted Breslow-Day"): (40.28344742835485, 11, 3.199346023374794e-5, np.nan),
        ("narrow", "all_years", "Inverse-variance Cochran Q"): (43.11266305705881, 11, 1.0381310309410039e-5, 74.48545457411967),
        ("broad", "all_years", "Inverse-variance Cochran Q"): (61.39384545688669, 12, 1.2560042610899139e-8, 80.45406683569463),
    }
    for (endpoint, window, test), expected_values in heterogeneity_expect.items():
        row = heterogeneity.query(
            "endpoint == @endpoint and window == @window and test == @test"
        ).iloc[0]
        observed = (row.statistic, int(row.degrees_of_freedom), row.p_value, row.i2_percent)
        assert np.allclose(observed[:3], expected_values[:3], rtol=1e-12)
        if np.isfinite(expected_values[3]):
            assert np.isclose(observed[3], expected_values[3], rtol=1e-12)

    italy = tables["country"].query("drug == 'nintedanib' and country_for_audit == 'IT'").iloc[0]
    july = tables["month"].query("drug == 'nintedanib' and receive_month == '2022-07'").iloc[0]
    assert int(italy.reports) == 104 and int(july.reports) == 96
    exclusion_expected = {
        "nintedanib": (783, 125, 600, 58),
        "pirfenidone": (199, 4, 180, 15),
    }
    for drug, expected_counts in exclusion_expected.items():
        row = tables["tto_exclusion"].query(
            "endpoint == 'broad' and drug == @drug"
        ).iloc[0]
        observed_counts = (
            int(row.exact_date_pairs_before_window),
            int(row.negative_intervals_excluded),
            int(row.intervals_0_365_included),
            int(row.intervals_over_365_excluded),
        )
        assert observed_counts == expected_counts

    return {
        "status": "pass",
        "unit": "latest-version safetyreportid (report)",
        "overall_independent_recalculation": checks,
        "age_unit_801_event_count_check": age_expected,
        "double_zero_contrast": "2014 narrow marked non-estimable",
        "mh_rbg_check": mh_expect,
        "year_effect_heterogeneity_check": {
            "status": "common annual effect rejected; MH retained as descriptive summary",
            "expected": {" | ".join(key): value for key, value in heterogeneity_expect.items()},
        },
        "2022_cluster": {"nintedanib_narrow_reports": 123, "Italy": 104, "July_2022": 96},
        "broad_tto_exclusion_audit": exclusion_expected,
    }


def forest_points(ax: mpl.axes.Axes, frame: pd.DataFrame, y: np.ndarray,
                  color: str | list[str], marker: str | list[str] = "o",
                  size: float = 4.0, zorder: int = 3) -> None:
    colors = [color] * len(frame) if isinstance(color, str) else color
    markers = [marker] * len(frame) if isinstance(marker, str) else marker
    for (_, row), yy, cc, mm in zip(frame.iterrows(), y, colors, markers):
        if pd.isna(row.ROR):
            continue
        ax.plot([row.ROR_low95, row.ROR_high95], [yy, yy], color=cc, lw=0.85, zorder=zorder)
        ax.plot(row.ROR, yy, marker=mm, ms=size, color=cc, markeredgecolor=cc,
                markerfacecolor="white" if mm == "o" and row.get("open_marker", False) else cc,
                markeredgewidth=0.75, zorder=zorder + 1)


def figure_one(t: dict[str, pd.DataFrame]) -> tuple[mpl.figure.Figure, pd.DataFrame]:
    fig, axes = plt.subplots(5, 2, figsize=(7.20, 10.20), layout="constrained")
    fig.suptitle(
        "Figure 1 | Active-comparator pharmacovigilance identifies disproportionate hepatic-event reporting",
        x=0.02, ha="left", fontsize=9.2, fontweight="bold",
    )
    axs = axes.ravel()
    for ax, letter in zip(axs, "ABCDEFGHIJ"):
        panel_label(ax, letter)

    overall = t["overall"].copy()

    # A: fractions of indexed reports meeting each endpoint, with exact counts.
    ax = axs[0]
    order = ["narrow", "broad", "hy_proxy"]
    labels = ["Narrow", "Broad", "ALT/AST +\nbilirubin/jaundice"]
    x = np.arange(len(order))
    width = 0.34
    for j, (drug, color, offset) in enumerate((("nintedanib", NINT, -width / 2), ("pirfenidone", PIRF, width / 2))):
        events = []
        totals = []
        for endpoint in order:
            row = overall[overall.endpoint.eq(endpoint)].iloc[0]
            events.append(int(row[f"{drug}_event"]))
            totals.append(int(row[f"{drug}_total"]))
        pct = 100 * np.asarray(events) / np.asarray(totals)
        bars = ax.bar(x + offset, pct, width=width, color=color, edgecolor="white", linewidth=0.5,
                      label=drug.capitalize(), zorder=2)
        for bar, n, N, p in zip(bars, events, totals, pct):
            ax.text(bar.get_x() + bar.get_width() / 2, p + max(0.09, 0.018 * p),
                    f"{n:,}", ha="center", va="bottom", fontsize=6.2, color=BLACK)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Endpoint reports / indexed reports (%)")
    ax.set_ylim(0, 7.25)
    ax.set_title("Report-set composition")
    ax.text(0.01, 0.985, "Indexed reports: N=30,233; P=39,467",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.2, color=MID_GREY)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.80), ncol=1,
              handlelength=1.2, labelspacing=0.3)
    style_axis(ax, "y")

    # B: overall active-comparator RORs.
    ax = axs[1]
    d = overall.set_index("endpoint").loc[["broad", "narrow", "hy_proxy"]].reset_index()
    yy = np.arange(len(d))[::-1]
    colors = [BROAD, NARROW, PROXY]
    forest_points(ax, d, yy, colors, "o", 4.3)
    for (_, row), yv in zip(d.iterrows(), yy):
        ax.text(row.ROR_high95 * 1.035, yv, f"{row.ROR:.2f}", va="center", fontsize=6.2)
    ax.axvline(1, color=BLACK, lw=0.65, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(0.75, 6.3)
    ax.set_xticks([1, 2, 4, 6], ["1", "2", "4", "6"])
    ax.set_yticks(yy, ["Broad", "Narrow", "Co-reporting proxy"])
    ax.set_xlabel("Reporting odds ratio (95% CI)")
    ax.set_title("Primary active-comparator signals")
    style_axis(ax, "x")

    # C: specification curve, including phrase/keyword and 2022 influence checks.
    ax = axs[2]
    primary = overall[overall.endpoint.isin(["narrow", "broad"])].copy()
    primary["spec"] = "Primary generic"
    alias = t["alias"].copy(); alias["spec"] = "Expanded aliases"
    specs = t["specs"].copy()
    spec_names = {
        "ild_indication_report_filter": "ILD phrase filter",
        "ild_indication_exact_keyword_sensitivity": "ILD exact-keyword filter",
        "suspect_drug_report_filter": "Any suspect drug in report",
        "ild_plus_suspect_drug_report_filter": "ILD phrase + suspect drug",
    }
    specs["spec"] = specs.analysis.map(spec_names)
    exc = t["temporal"].query("temporal_window == 'exclude_2022_posthoc'").copy()
    exc["spec"] = "Exclude 2022 (post hoc)"
    s = pd.concat([
        primary[["spec", "endpoint", "ROR", "ROR_low95", "ROR_high95"]],
        alias[["spec", "endpoint", "ROR", "ROR_low95", "ROR_high95"]],
        specs[["spec", "endpoint", "ROR", "ROR_low95", "ROR_high95"]],
        exc[["spec", "endpoint", "ROR", "ROR_low95", "ROR_high95"]],
    ], ignore_index=True)
    spec_order = ["Primary generic", "Expanded aliases", "ILD phrase filter",
                  "ILD exact-keyword filter", "Any suspect drug in report",
                  "ILD phrase + suspect drug", "Exclude 2022 (post hoc)"]
    ybase = np.arange(len(spec_order))[::-1]
    for endpoint, color, delta, marker in (("broad", BROAD, -0.11, "s"), ("narrow", NARROW, 0.11, "o")):
        dd = s[s.endpoint.eq(endpoint)].set_index("spec").loc[spec_order].reset_index()
        forest_points(ax, dd, ybase + delta, color, marker, 3.7)
    ax.axvline(1, color=BLACK, lw=0.65, ls="--")
    ax.set_xscale("log"); ax.set_xlim(0.9, 6.2)
    ax.set_xticks([1, 2, 4, 6], ["1", "2", "4", "6"])
    ax.set_yticks(ybase, spec_order)
    ax.set_xlabel("ROR (95% CI)")
    ax.set_title("Specification and influence analyses")
    ax.plot([], [], "s", color=BROAD, ms=4, label="Broad")
    ax.plot([], [], "o", color=NARROW, ms=4, label="Narrow")
    ax.legend(loc="upper center", bbox_to_anchor=(0.54, -0.15), ncol=2,
              columnspacing=0.8, handletextpad=0.3)
    style_axis(ax, "x")

    # D/E: annual signals; partial years are open circles and 2022 is highlighted.
    for ax, endpoint, title, color, ylim, ticks in (
        (axs[3], "broad", "Annual broad-endpoint ROR", BROAD, (0.3, 30), [0.5, 1, 2, 4, 8, 20]),
        (axs[4], "narrow", "Annual narrow-endpoint ROR", NARROW, (0.08, 8000), [0.1, 1, 10, 100, 1000]),
    ):
        d = t["annual"].query("endpoint == @endpoint").copy()
        d = d[d.ROR.notna()]
        ax.axvspan(2021.72, 2022.28, color="#F4D6C7", alpha=0.8, lw=0, zorder=0)
        ax.axhline(1, color=BLACK, lw=0.65, ls="--")
        for _, row in d.iterrows():
            open_marker = int(row.year) in (2014, 2026)
            ax.plot([row.year, row.year], [row.ROR_low95, row.ROR_high95], color=color, lw=0.7, zorder=2)
            ax.plot(row.year, row.ROR, "o", ms=3.3, markeredgecolor=color, markeredgewidth=0.75,
                    markerfacecolor="white" if open_marker else color, zorder=3)
        ax.plot(d.year, d.ROR, color=color, lw=0.7, alpha=0.65, zorder=1)
        ax.set_yscale("log"); ax.set_ylim(*ylim)
        ax.set_yticks(ticks, [str(v) for v in ticks])
        ax.set_xlim(2013.6, 2026.4); ax.set_xticks(range(2014, 2027, 2))
        ax.set_xlabel("Report receipt year")
        ax.set_ylabel("ROR (95% CI)")
        ax.set_title(title)
        style_axis(ax, "y")
    axs[3].text(2022, 9.4, "2022", ha="center", va="bottom", fontsize=6.2, color=NINT)
    axs[4].text(2022, 700, "123 vs 0 reports", ha="center", va="bottom", fontsize=6.2, color=NINT)
    axs[4].text(2014.0, 0.115, "2014: 0 vs 0; NE", fontsize=6.2, color=MID_GREY, ha="left")

    # F: crude and descriptive year-stratified Mantel-Haenszel summaries.
    ax = axs[5]
    ybase = np.array([1, 0])
    for endpoint, yv, color in (("narrow", 1, NARROW), ("broad", 0, BROAD)):
        raw = overall[overall.endpoint.eq(endpoint)].iloc[0]
        mh = t["mh"][t["mh"].endpoint.eq(endpoint)].iloc[0]
        ax.plot([raw.ROR_low95, raw.ROR_high95], [yv + 0.10] * 2, color=color, lw=0.85)
        ax.plot(raw.ROR, yv + 0.10, "o", color=color, ms=4, label="Crude" if endpoint == "narrow" else None)
        ax.plot([mh.MH_low95_RBG, mh.MH_high95_RBG], [yv - 0.10] * 2, color=color, lw=0.85)
        ax.plot(mh.MH_ROR_year_adjusted, yv - 0.10, "D", color=color, ms=3.6,
                label="Descriptive year-stratified MH" if endpoint == "narrow" else None)
    ax.axvline(1, color=BLACK, lw=0.65, ls="--")
    ax.set_xscale("log"); ax.set_xlim(0.9, 6.5)
    ax.set_xticks([1, 2, 4, 6], ["1", "2", "4", "6"])
    ax.set_yticks(ybase, ["Narrow", "Broad"])
    ax.set_xlabel("ROR (95% CI)")
    ax.set_title("Crude and descriptive year-stratified")
    ax.text(
        0.98, 0.48,
        "Tarone–BD heterogeneity\n"
        "Narrow: P=3.7×10⁻¹⁴\nBroad: P=5.6×10⁻⁹",
        transform=ax.transAxes, ha="right", va="center", fontsize=5.5, color=MID_GREY,
    )
    ax.legend(loc="upper left", ncol=1, handletextpad=0.3)
    style_axis(ax, "x")

    # G/H: sex and years-only age strata.
    for ax, dimension, labels_map, title in (
        (axs[6], "sex", {"male": "Male", "female": "Female"}, "Sex-stratified signals"),
        (axs[7], "age", {"18-64_years": "18–64 years", "65plus_years": "≥65 years"}, "Age-stratified signals (years only)"),
    ):
        d = t["strata"].query("dimension == @dimension").copy()
        order_strata = list(labels_map)
        ybase = np.arange(len(order_strata))[::-1]
        for endpoint, color, delta, marker in (("broad", BROAD, -0.10, "s"), ("narrow", NARROW, 0.10, "o")):
            dd = d[d.endpoint.eq(endpoint)].set_index("stratum").loc[order_strata].reset_index()
            forest_points(ax, dd, ybase + delta, color, marker, 3.8)
        ax.axvline(1, color=BLACK, lw=0.65, ls="--")
        ax.set_xscale("log"); ax.set_xlim(0.9, 10)
        ax.set_xticks([1, 2, 4, 8], ["1", "2", "4", "8"])
        ax.set_yticks(ybase, [labels_map[s] for s in order_strata])
        ax.set_xlabel("ROR (95% CI)")
        ax.set_title(title)
        ax.plot([], [], "s", color=BROAD, ms=4, label="Broad")
        ax.plot([], [], "o", color=NARROW, ms=4, label="Narrow")
        ax.legend(loc="upper left", ncol=1, handletextpad=0.3)
        style_axis(ax, "x")

    # I: serious outcomes among broad-endpoint reports (overlapping categories).
    ax = axs[8]
    out = t["outcomes"].query("endpoint == 'broad'").copy()
    outcome_order = ["serious", "hospitalization", "death", "life_threatening", "fatal_hepatic_reaction"]
    outcome_labels = ["Any serious", "Hospitalization", "Death", "Life-threatening", "Fatal hepatic PT"]
    y = np.arange(len(outcome_order))[::-1]
    h = 0.34
    for drug, color, delta in (("nintedanib", NINT, h / 2), ("pirfenidone", PIRF, -h / 2)):
        dd = out[out.drug.eq(drug)].set_index("outcome").loc[outcome_order]
        bars = ax.barh(y + delta, dd.percent, height=h, color=color, edgecolor="white", linewidth=0.4,
                       label=drug.capitalize(), zorder=2)
        for bar, (_, row) in zip(bars, dd.iterrows()):
            ax.text(min(row.percent + 1.1, 96), bar.get_y() + bar.get_height() / 2,
                    f"{int(row.events)}/{int(row.endpoint_reports)}", va="center", fontsize=6.2)
    ax.set_yticks(y, outcome_labels)
    ax.set_xlim(0, 105); ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Broad-endpoint reports (%)")
    ax.set_title("Serious outcomes within endpoint reports")
    ax.legend(loc="lower right", ncol=2, columnspacing=0.8, handlelength=1.0)
    style_axis(ax, "x")

    # J: time-to-onset empirical cumulative distributions for broad cases.
    ax = axs[9]
    vals = t["tto_values"].query("endpoint == 'broad'").copy()
    summ = t["tto_summary"].query("endpoint == 'broad'").set_index("drug")
    for drug, color in (("nintedanib", NINT), ("pirfenidone", PIRF)):
        arr = np.sort(vals.loc[vals.drug.eq(drug), "time_to_onset_days"].astype(float).to_numpy())
        ecdf = np.arange(1, len(arr) + 1) / len(arr)
        med = float(summ.loc[drug, "median_days"])
        ax.step(arr, ecdf * 100, where="post", color=color,
                label=f"{drug.capitalize()} (n={len(arr)}, median {med:g} d)")
        ax.axvline(med, color=color, lw=0.7, ls=":")
    ax.set_xlim(0, 365); ax.set_ylim(0, 101)
    ax.set_xticks([0, 30, 90, 180, 365])
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Days from recorded drug start to case event")
    ax.set_ylabel("Cumulative valid date pairs (%)")
    ax.set_title("Time-to-onset among exact date pairs")
    ax.legend(loc="lower right", handlelength=2.2)
    style_axis(ax, "both")

    # Ensure no accidental panel count drift.
    assert len(axs) == 10
    fig.align_ylabels()

    panel_sources = []
    source_map = {
        "A-B-F": ("overall_active_comparator_signals.csv", overall),
        "C_alias": ("alias_sensitivity_signals.csv", t["alias"]),
        "C_context": ("report_level_specification_signals.csv", t["specs"]),
        "C_temporal": ("temporal_influence_signals.csv", t["temporal"]),
        "D-E": ("annual_active_comparator_signals.csv", t["annual"]),
        "F_MH": ("year_adjusted_mantel_haenszel_signals.csv", t["mh"]),
        "F_heterogeneity": ("year_effect_heterogeneity_tests.csv", t["heterogeneity"]),
        "G-H": ("sex_age_stratified_signals.csv", t["strata"]),
        "I": ("serious_outcome_profiles.csv", t["outcomes"]),
        "J_values": ("time_to_onset_case_values.csv", t["tto_values"]),
        "J_summary": ("time_to_onset_summary.csv", t["tto_summary"]),
    }
    for panel, (source, frame) in source_map.items():
        copy = frame.copy(); copy.insert(0, "source_file", source); copy.insert(0, "panel", panel)
        panel_sources.append(copy)
    return fig, pd.concat(panel_sources, ignore_index=True, sort=False)


def supplementary_figure(t: dict[str, pd.DataFrame]) -> tuple[mpl.figure.Figure, pd.DataFrame]:
    fig = plt.figure(figsize=(7.35, 10.8), layout="constrained")
    gs = fig.add_gridspec(4, 2, height_ratios=[1.35, 1.35, 1, 1])
    ax_a = fig.add_subplot(gs[0:2, :])
    ax_b = fig.add_subplot(gs[2, 0])
    ax_c = fig.add_subplot(gs[2, 1])
    ax_d = fig.add_subplot(gs[3, 0])
    nested = gs[3, 1].subgridspec(1, 2, wspace=0.35)
    ax_e = fig.add_subplot(nested[0, 0])
    ax_f = fig.add_subplot(nested[0, 1])
    axes = [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f]
    for ax, letter in zip(axes[:4], "ABCD"):
        panel_label(ax, letter, -0.10 if ax is ax_a else -0.20, 1.04 if ax is ax_a else 1.11)
    panel_label(ax_e, "E", -0.08, 1.25)
    panel_label(ax_f, "F", -0.08, 1.25)

    # A: PT-level report signals (overlapping terms).
    pt = t["pt"].copy()
    pt["events_total"] = pt.nintedanib_event + pt.pirfenidone_event
    pt = pt[pt.ROR.notna() & pt.events_total.ge(4)].sort_values("ROR")
    pretty = {
        "DRUG-INDUCED LIVER INJURY": "Drug-induced liver injury",
        "ALANINE AMINOTRANSFERASE INCREASED": "ALT increased",
        "ASPARTATE AMINOTRANSFERASE INCREASED": "AST increased",
        "GAMMA-GLUTAMYLTRANSFERASE INCREASED": "GGT increased",
        "BLOOD ALKALINE PHOSPHATASE INCREASED": "Alkaline phosphatase increased",
        "BLOOD BILIRUBIN INCREASED": "Blood bilirubin increased",
        "LIVER FUNCTION TEST INCREASED": "Liver function test increased",
        "LIVER FUNCTION TEST ABNORMAL": "Liver function test abnormal",
        "HEPATIC FUNCTION ABNORMAL": "Hepatic function abnormal",
        "HEPATIC ENZYME INCREASED": "Hepatic enzyme increased",
    }
    labels = [pretty.get(x, x.title()) for x in pt.preferred_term]
    y = np.arange(len(pt))
    colors = [NARROW if x == "narrow" else BROAD for x in pt.endpoint_membership]
    forest_points(ax_a, pt, y, colors, "o", 3.4)
    ax_a.axvline(1, color=BLACK, lw=0.65, ls="--")
    ax_a.set_xscale("log"); ax_a.set_xlim(0.25, 420)
    ax_a.set_xticks([0.5, 1, 2, 5, 10, 50, 200], ["0.5", "1", "2", "5", "10", "50", "200"])
    ax_a.set_yticks(y, labels)
    ax_a.set_xlabel("Report-level ROR (95% CI)")
    ax_a.set_title("Preferred-term signals (terms overlap within reports)")
    ax_a.plot([], [], "o", color=NARROW, ms=4, label="Narrow endpoint term")
    ax_a.plot([], [], "o", color=BROAD, ms=4, label="Broad-only term")
    ax_a.legend(loc="lower right", ncol=2)
    style_axis(ax_a, "x")

    # B: matching-drug role profile among broad endpoint reports.
    role = t["roles"].query("endpoint == 'broad' and variable == 'matching_drug_role'").copy()
    roles = ["suspect", "concomitant", "interacting", "unknown"]
    role_colors = [NINT, GOLD, SKY, LIGHT_GREY]
    y = [1, 0]
    for iy, drug in zip(y, ["nintedanib", "pirfenidone"]):
        dd = role[role.drug.eq(drug)].set_index("level").percent
        left = 0
        for level, color in zip(roles, role_colors):
            value = float(dd.get(level, 0))
            ax_b.barh(iy, value, left=left, color=color, height=0.55, edgecolor="white", lw=0.4)
            left += value
    ax_b.set_yticks(y, ["Nintedanib", "Pirfenidone"])
    ax_b.set_xlim(0, 100); ax_b.set_xticks([0, 50, 100])
    ax_b.set_xlabel("Broad reports (%)")
    ax_b.set_title("Matching-drug role")
    ax_b.legend([Patch(color=c) for c in role_colors[:3]],
                ["Suspect", "Concomitant", "Interacting"], loc="lower center", ncol=3,
                bbox_to_anchor=(0.5, -0.43), columnspacing=0.6, handlelength=0.8)
    style_axis(ax_b, "x")

    # C: recorded ILD indication in the matched-drug object among endpoint reports.
    ind = t["roles"][(t["roles"].variable == "matching_drug_ild_indication") &
                     (t["roles"].level.astype(str).str.lower() == "true")].copy()
    x = np.arange(2); width = 0.34
    for drug, color, delta in (("nintedanib", NINT, -width / 2), ("pirfenidone", PIRF, width / 2)):
        dd = ind[ind.drug.eq(drug)].set_index("endpoint").loc[["narrow", "broad"]]
        bars = ax_c.bar(x + delta, dd.percent, width, color=color, label=drug.capitalize())
        for bar, value in zip(bars, dd.percent):
            ax_c.text(bar.get_x() + bar.get_width()/2, value + 1, f"{value:.0f}%", ha="center", fontsize=5.4)
    ax_c.set_xticks(x, ["Narrow", "Broad"]); ax_c.set_ylim(0, 105)
    ax_c.set_ylabel("Endpoint reports with ILD indication (%)")
    ax_c.set_title("Matched-drug indication field")
    ax_c.legend(loc="lower left", ncol=2, columnspacing=0.7, handlelength=0.9)
    style_axis(ax_c, "y")

    # D: exact-date completeness with median/IQR annotations.
    summ = t["tto_summary"].copy()
    key_order = [("narrow", "nintedanib"), ("narrow", "pirfenidone"),
                 ("broad", "nintedanib"), ("broad", "pirfenidone")]
    labels = ["Narrow / N", "Narrow / P", "Broad / N", "Broad / P"]
    y = np.arange(4)[::-1]
    colors = [NINT, PIRF, NINT, PIRF]
    for yy, key, color in zip(y, key_order, colors):
        row = summ[(summ.endpoint == key[0]) & (summ.drug == key[1])].iloc[0]
        ax_d.barh(yy, row.date_pair_completeness_percent, color=color, height=0.58)
        ax_d.text(row.date_pair_completeness_percent + 0.8, yy,
                  f"{int(row.exact_date_pairs_0_365)}/{int(row.endpoint_reports)}; "
                  f"median {row.median_days:g} d", va="center", fontsize=5.1)
    ax_d.set_yticks(y, labels); ax_d.set_xlim(0, 56)
    ax_d.set_xlabel("Exact 0–365-day date pairs (%)")
    ax_d.set_title("Time-to-onset completeness")
    style_axis(ax_d, "x")

    # E/F: denominator-free audit of the 2022 narrow-endpoint cluster.
    country = t["country"].query("drug == 'nintedanib'").copy().sort_values("reports", ascending=False)
    top = country.head(5).copy()
    other = int(country.iloc[5:].reports.sum())
    labels = list(top.country_for_audit) + (["Other"] if other else [])
    values = list(top.reports.astype(int)) + ([other] if other else [])
    y = np.arange(len(labels))[::-1]
    bars = ax_e.barh(y, values, color=[NINT] + ["#C6C6C6"] * (len(values) - 1), height=0.65)
    for bar, value in zip(bars, values):
        ax_e.text(value + 1, bar.get_y() + bar.get_height()/2, str(value), va="center", fontsize=5.0)
    ax_e.set_yticks(y, labels); ax_e.set_xlim(0, 122)
    ax_e.set_xlabel("Reports")
    ax_e.set_title("2022 narrow: country")
    style_axis(ax_e, "x")

    month = t["month"].query("drug == 'nintedanib'").copy()
    month["month_num"] = month.receive_month.str[-2:].astype(int)
    ax_f.plot(month.month_num, month.reports, color=NINT, marker="o", ms=2.8)
    ax_f.fill_between(month.month_num, month.reports, color=NINT, alpha=0.12)
    ax_f.text(7, 96, "96", ha="center", va="bottom", fontsize=5.3, color=NINT)
    ax_f.set_xticks([1, 4, 7, 10, 12]); ax_f.set_xlim(0.5, 12.5)
    ax_f.set_ylim(0, 105); ax_f.set_xlabel("Receipt month")
    ax_f.set_ylabel("Reports")
    ax_f.set_title("2022 narrow: month")
    style_axis(ax_f, "both")

    source_map = {
        "A": ("preferred_term_signals.csv", t["pt"]),
        "B-C": ("drug_role_and_indication_profiles.csv", t["roles"]),
        "D_summary": ("time_to_onset_summary.csv", t["tto_summary"]),
        "D_exclusions": ("time_to_onset_exclusion_audit.csv", t["tto_exclusion"]),
        "E": ("narrow_2022_country_audit.csv", t["country"]),
        "F": ("narrow_2022_month_audit.csv", t["month"]),
    }
    panel_sources = []
    for panel, (source, frame) in source_map.items():
        copy = frame.copy(); copy.insert(0, "source_file", source); copy.insert(0, "panel", panel)
        panel_sources.append(copy)
    return fig, pd.concat(panel_sources, ignore_index=True, sort=False)


def write_captions() -> None:
    main = """# Figure 1. Active-comparator pharmacovigilance identifies disproportionate hepatic-event reporting with a marked 2022 narrow-endpoint cluster.

**(A)** Proportion of mutually exclusive nintedanib- and pirfenidone-indexed reports meeting the prespecified narrow, broad, and exploratory aminotransferase-plus-bilirubin/jaundice co-reporting definitions; bar labels give endpoint-report counts and the annotation gives each indexed-report denominator. **(B)** Full-window active-comparator reporting odds ratios (RORs) with Wald 95% confidence intervals (CIs). The co-reporting constellation is exploratory and is not Hy's law because laboratory values, upper limits of normal, alkaline-phosphatase exclusion, and clinical adjudication were unavailable. **(C)** Specification curve comparing the primary generic-name analysis with expanded product aliases, ILD-indication phrase filtering, the case-sensitive exact-keyword field-behavior sensitivity check, a report-level filter requiring any suspect drug (official openFDA `drugcharacterization` code 1), their joint report-level filter, and a post hoc influence analysis excluding 2022. The openFDA indication and role filters operate at report level and cannot guarantee that the filtered field belongs to the same drug-array element as the indexed product. **(D,E)** Calendar-year broad- and narrow-endpoint RORs. Open symbols denote partial calendar years (15 October–31 December 2014 and 1 January–30 June 2026); 2014 narrow was non-estimable (0 versus 0 endpoint reports). Shading marks 2022, when the narrow estimate was driven by 123 nintedanib versus 0 pirfenidone reports, demonstrating marked temporal variation rather than stability. **(F)** Crude RORs and descriptive calendar-year-stratified Mantel–Haenszel (MH) summaries; MH CIs use the Robins–Breslow–Greenland variance. Equality of annual odds ratios was rejected by the Tarone-adjusted Breslow–Day test for the narrow endpoint (χ²=88.25, 11 df, P=3.66×10⁻¹⁴) and broad endpoint (χ²=63.33, 12 df, P=5.55×10⁻⁹), and remained rejected after excluding 2022 (narrow: χ²=33.93, 10 df, P=1.90×10⁻⁴; broad: χ²=40.28, 11 df, P=3.20×10⁻⁵). Inverse-variance Cochran Q likewise showed full-window heterogeneity (narrow: Q=43.11, 11 df, P=1.04×10⁻⁵, I²=74.5%; broad: Q=61.39, 12 df, P=1.26×10⁻⁸, I²=80.5%). Accordingly, MH values summarize year-stratified reporting descriptively and are not stable common effects. **(G,H)** Sex- and age-stratified RORs. Age strata include only reports with ICH age-unit code 801 (years). **(I)** Serious outcomes among broad-endpoint reports; seriousness fields describe the whole report, outcome categories can overlap, and fatal hepatic preferred term denotes reaction outcome code 5 on a hepatic term. **(J)** Empirical cumulative distributions of time to onset among broad-endpoint reports with an exact YYYYMMDD matching-drug start date, an exact report-level CASE EVENT DATE, and an interval of 0–365 days. The event date need not refer specifically to the hepatic reaction; these selected-report timings are not hazards or causal latency estimates. Data are latest-version spontaneous reports received from 15 October 2014 through 30 June 2026 in openFDA (release updated 30 July 2026; retrieved 4 September 2026). The unit is a safety report, not a patient or treated-person denominator. RORs quantify disproportionate reporting relative to pirfenidone and do not estimate incidence, absolute risk, or causality. N, nintedanib; P, pirfenidone; NE, non-estimable; MH, Mantel–Haenszel.
"""
    supp = """# Supplementary Figure S1. Preferred-term, reporting-context, date-completeness, and 2022-cluster audits.

**(A)** Preferred-term-specific RORs with 95% CIs for terms with at least four combined endpoint reports and an estimable contrast. Preferred terms overlap within reports and were not subjected to multiplicity-adjusted confirmatory inference. **(B)** Role of the matching drug among broad-endpoint reports, parsed directly from downloaded report objects using the official openFDA mapping: code 1, suspect; code 2, concomitant; code 3, interacting. openFDA does not subdivide suspect drugs into primary and secondary categories. **(C)** Fraction of endpoint reports in which the matching-drug indication text contained a prespecified ILD-related phrase. **(D)** Availability of exact 0–365-day time-to-onset date pairs; annotations give valid pairs/endpoint reports and the median among valid pairs. Before windowing, 783 nintedanib broad-endpoint reports had exact date pairs (125 negative intervals and 58 intervals >365 days excluded) and 199 pirfenidone reports had exact pairs (4 negative and 15 >365 days excluded), leaving 600 and 180 plotted values, respectively. **(E,F)** Denominator-free country and receipt-month counts for the 2022 nintedanib narrow-endpoint cluster. Of 123 reports, 104 (84.6%) listed Italy as the primary-source country and 96 (78.0%) were received in July. These distributions demonstrate reporting concentration and are not geographic or temporal rates. All panels use latest-version openFDA safety reports and must not be interpreted as incidence or causal effects.
"""
    (MAIN_FIGURES / "Figure_1_FAERS_caption.md").write_text(main, encoding="utf-8")
    (SUPPLEMENTARY_FIGURES / "Supplementary_Figure_S1_FAERS_caption.md").write_text(supp, encoding="utf-8")


def write_manifest(figure_paths: list[Path], validation: dict[str, object]) -> None:
    config = json.loads((RESULTS / "analysis_configuration.json").read_text(encoding="utf-8"))
    result_paths = sorted(
        p for p in RESULTS.iterdir()
        if p.is_file() and p.name != "faers_frozen_manifest.json"
    )
    scripts = [
        FAERS / "run_openfda_faers.py",
        FAERS / "audit_year_heterogeneity.py",
        FAERS / "migrate_frozen_faers_semantics.py",
        HERE / "make_faers_figures.py",
    ]
    manifest = {
        "status": "frozen_after_complete_rerun_and_independent_validation",
        "retrieval_date": config["retrieval_date"],
        "openfda_last_updated": config["api_last_updated"],
        "analysis_receivedate_window": [config["analysis_start_receivedate"], config["analysis_end_receivedate"]],
        "unit": "one latest-version openFDA safetyreportid (spontaneous report)",
        "primary_exposure_definition": config["primary_drug_queries"],
        "coexposure_rule": config["coexposure_rule"],
        "indication_phrase_sensitivity": {
            "primary_field": "patient.drug.drugindication (analyzed phrase; capitalization-insensitive)",
            "sensitivity_field": "patient.drug.drugindication.exact (case-sensitive keyword)",
            "reason": "the exact keyword field missed capitalization variants visible in downloaded JSON",
            "terms": config["ild_indication_terms_for_report_level_specification"],
            "relational_boundary": "report-level filtering cannot prove same-array-element linkage to the indexed drug",
        },
        "validation": validation,
        "script_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in scripts},
        "result_sha256": {p.name: sha256(p) for p in result_paths},
        "figure_sha256": {p.name: sha256(p) for p in sorted(figure_paths)},
    }
    path = RESULTS / "faers_frozen_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> None:
    tables = {
        "overall": read("overall_active_comparator_signals.csv"),
        "alias": read("alias_sensitivity_signals.csv"),
        "specs": read("report_level_specification_signals.csv"),
        "annual": read("annual_active_comparator_signals.csv"),
        "mh": read("year_adjusted_mantel_haenszel_signals.csv"),
        "heterogeneity": read("year_effect_heterogeneity_tests.csv"),
        "strata": read("sex_age_stratified_signals.csv"),
        "outcomes": read("serious_outcome_profiles.csv"),
        "tto_values": read("time_to_onset_case_values.csv"),
        "tto_summary": read("time_to_onset_summary.csv"),
        "tto_exclusion": read("time_to_onset_exclusion_audit.csv"),
        "pt": read("preferred_term_signals.csv"),
        "roles": read("drug_role_and_indication_profiles.csv"),
        "cases": read("broad_liver_event_reports_deduplicated.csv"),
        "temporal": read("temporal_influence_signals.csv"),
        "country": read("narrow_2022_country_audit.csv"),
        "month": read("narrow_2022_month_audit.csv"),
    }
    validation = validate_inputs(tables)
    validation_path = CHECKLISTS / "FAERS_figure_validation_report.json"
    validation_path.write_text(json.dumps(validation, indent=2), encoding="utf-8")

    fig1, source1 = figure_one(tables)
    fig1_paths = save_all(fig1, "Figure_1_FAERS")
    plt.close(fig1)
    source1_path = PANEL_DATA / "Figure_1_FAERS_panel_source_data.csv"
    source1.to_csv(source1_path, index=False)

    figs, sources = supplementary_figure(tables)
    supp_paths = save_all(figs, "Supplementary_Figure_S1_FAERS")
    plt.close(figs)
    source_s_path = PANEL_DATA / "Supplementary_Figure_S1_FAERS_panel_source_data.csv"
    sources.to_csv(source_s_path, index=False)

    write_captions()
    caption_paths = [MAIN_FIGURES / "Figure_1_FAERS_caption.md", SUPPLEMENTARY_FIGURES / "Supplementary_Figure_S1_FAERS_caption.md"]
    all_paths = fig1_paths + supp_paths + [source1_path, source_s_path, validation_path] + caption_paths
    write_manifest(all_paths, validation)

    print(json.dumps({
        "status": "complete",
        "main_figure_panels": list("ABCDEFGHIJ"),
        "main_figure_files": [str(p) for p in fig1_paths],
        "supplementary_figure_files": [str(p) for p in supp_paths],
        "validation": validation["status"],
        "manifest": str(RESULTS / "faers_frozen_manifest.json"),
    }, indent=2))


if __name__ == "__main__":
    main()
