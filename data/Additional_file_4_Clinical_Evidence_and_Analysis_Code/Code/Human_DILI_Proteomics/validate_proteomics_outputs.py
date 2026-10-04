#!/usr/bin/env python3
"""Deterministic integrity checks for the DILI proteomics re-analysis outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(condition), "detail": detail})
        if not condition:
            raise AssertionError(f"{name}: {detail}")

    manifest = json.loads((HERE / "analysis_manifest.json").read_text())
    source = ROOT / "fresh_data_dili_source.xlsx"
    check(
        "source_hash_matches_manifest",
        file_hash(source) == manifest["source_sha256"],
        manifest["source_sha256"],
    )

    discovery = pd.read_csv(HERE / "discovery_tidy.csv", keep_default_na=False, na_values=["NA"])
    confirmatory = pd.read_csv(HERE / "confirmatory_tidy.csv", keep_default_na=False, na_values=["NA"])
    effects = pd.read_csv(HERE / "effect_estimates.csv")
    correlations = pd.read_csv(HERE / "correlation_stats.csv")
    duplicate_qc = pd.read_csv(HERE / "discovery_duplicate_series_qc.csv")

    check("thirteen_discovery_candidates", discovery.protein.nunique() == 13, str(discovery.protein.nunique()))
    check("thirteen_confirmatory_candidates", confirmatory.protein.nunique() == 13, str(confirmatory.protein.nunique()))
    check(
        "no_duplicate_long_keys_discovery",
        not discovery.duplicated(["protein", "group", "source_row_index_within_group"]).any(),
        "protein-group-position keys are unique",
    )
    check(
        "no_duplicate_long_keys_confirmatory",
        not confirmatory.duplicated(["protein", "group", "source_row_index_within_group"]).any(),
        "protein-group-position keys are unique",
    )

    onset = confirmatory[confirmatory.group.isin(["HV", "DO", "NDO"])]
    counts = onset.groupby(["protein", "group"]).size().unstack()
    check("canonical_onset_rows", bool((counts.HV == 60).all() and (counts.DO == 82).all() and (counts.NDO == 34).all()), str(counts.to_dict()))

    for col in ["mann_whitney_q_bh", "welch_t_q_bh", "roc_auc_group1_positive"]:
        check(f"valid_range_{col}", bool(effects[col].dropna().between(0, 1).all()), col)
    check("valid_correlation_q", bool(correlations.spearman_q_bh.dropna().between(0, 1).all()), "spearman_q_bh")

    fbp1 = effects.query(
        "cohort == 'confirmatory' and protein == 'FBP1' and contrast == 'DO_vs_NDO'"
    ).iloc[0]
    check(
        "fbp1_specificity_effect_reproduced",
        np.isclose(fbp1.cliffs_delta_group1_minus_group2, -0.505345394736842, atol=1e-12),
        str(fbp1.cliffs_delta_group1_minus_group2),
    )
    check(
        "fbp1_specificity_fdr_reproduced",
        np.isclose(fbp1.mann_whitney_q_bh, 0.0004687724004871, atol=1e-12),
        str(fbp1.mann_whitney_q_bh),
    )
    check(
        "discovery_duplicate_detected",
        bool(((duplicate_qc.protein_1 == "ALDOB") & (duplicate_qc.protein_2 == "HPD")).any()),
        "ALDOB/HPD exact duplicate is explicitly flagged",
    )
    check(
        "no_wide_multivariable_matrix_emitted",
        not (HERE / "confirmatory_wide.csv").exists(),
        "Cross-protein positional joins are prohibited",
    )

    report = {
        "status": "PASS",
        "n_checks": len(checks),
        "checks": checks,
    }
    (HERE / "validation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
