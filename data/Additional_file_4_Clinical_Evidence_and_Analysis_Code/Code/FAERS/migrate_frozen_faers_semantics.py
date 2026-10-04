#!/usr/bin/env python3
"""Correct drugcharacterization labels in the frozen FAERS release tables.

The original API queries and counts were correct, but the human-readable role
labels were mapped incorrectly. This deterministic migration changes labels and
one derived column name only; it does not alter report membership or counts.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from typing import Callable


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
DEFAULT_RESULTS = PACKAGE_ROOT / "Source_Data" / "FAERS"

ROLE_MAP = {
    "primary_suspect": "suspect",
    "secondary_suspect": "concomitant",
    "concomitant": "interacting",
}
ANALYSIS_MAP = {
    "primary_suspect_report_filter": "suspect_drug_report_filter",
    "ild_plus_primary_suspect_report_filter": "ild_plus_suspect_drug_report_filter",
}


def _rewrite_csv(path: Path, transform: Callable[[dict[str, str]], dict[str, str]],
                 rename: dict[str, str] | None = None) -> int:
    rename = rename or {}
    temporary = path.with_name(f".{path.name}.semantic-correction.tmp")
    changed = 0
    with path.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {path}")
        fieldnames = [rename.get(field, field) for field in reader.fieldnames]
        with temporary.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(destination, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            for original in reader:
                row = {rename.get(key, key): value for key, value in original.items()}
                corrected = transform(row)
                changed += int(corrected != row)
                writer.writerow(corrected)
    os.replace(temporary, path)
    return changed


def migrate(results: Path) -> dict[str, int]:
    # The renamed Boolean column is an unambiguous migration-state marker. This
    # guard is essential because the old code-3 label ("concomitant") is also a
    # valid corrected code-2 label and therefore cannot be remapped twice.
    cases = results / "broad_liver_event_reports_deduplicated.csv"
    with cases.open("r", encoding="utf-8", newline="") as source:
        header = next(csv.reader(source))
    if "matching_drug_suspect_or_concomitant" in header:
        validate(results)
        return {
            name: 0
            for name in (
                cases.name,
                "time_to_onset_case_values.csv",
                "drug_role_and_indication_profiles.csv",
                "report_level_specification_signals.csv",
                "api_query_log.csv",
            )
        }
    if "matching_drug_primary_or_secondary_suspect" not in header:
        raise AssertionError("Cannot determine frozen-table semantic version")

    changes: dict[str, int] = {}

    def case_row(row: dict[str, str]) -> dict[str, str]:
        row = row.copy()
        row["matching_drug_role"] = ROLE_MAP.get(row["matching_drug_role"], row["matching_drug_role"])
        return row

    changes[cases.name] = _rewrite_csv(
        cases,
        case_row,
        rename={
            "matching_drug_primary_or_secondary_suspect":
                "matching_drug_suspect_or_concomitant"
        },
    )

    onset = results / "time_to_onset_case_values.csv"
    changes[onset.name] = _rewrite_csv(onset, case_row)

    def profile_row(row: dict[str, str]) -> dict[str, str]:
        row = row.copy()
        if row["variable"] == "matching_drug_role":
            row["level"] = ROLE_MAP.get(row["level"], row["level"])
        return row

    profile = results / "drug_role_and_indication_profiles.csv"
    changes[profile.name] = _rewrite_csv(profile, profile_row)

    def specification_row(row: dict[str, str]) -> dict[str, str]:
        row = row.copy()
        row["analysis"] = ANALYSIS_MAP.get(row["analysis"], row["analysis"])
        return row

    specifications = results / "report_level_specification_signals.csv"
    changes[specifications.name] = _rewrite_csv(specifications, specification_row)

    def query_log_row(row: dict[str, str]) -> dict[str, str]:
        row = row.copy()
        for old, new in ANALYSIS_MAP.items():
            row["label"] = row["label"].replace(old, new)
        return row

    query_log = results / "api_query_log.csv"
    changes[query_log.name] = _rewrite_csv(query_log, query_log_row)
    return changes


def validate(results: Path) -> None:
    prohibited = (
        "primary_suspect",
        "secondary_suspect",
        "matching_drug_primary_or_secondary_suspect",
    )
    for name in (
        "broad_liver_event_reports_deduplicated.csv",
        "time_to_onset_case_values.csv",
        "drug_role_and_indication_profiles.csv",
        "report_level_specification_signals.csv",
        "api_query_log.csv",
    ):
        text = (results / name).read_text(encoding="utf-8")
        present = [term for term in prohibited if term in text]
        if present:
            raise AssertionError(f"Uncorrected labels in {name}: {present}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if not args.validate_only:
        for name, count in migrate(args.results).items():
            print(f"{name}: {count} corrected rows")
    validate(args.results)
    print("FAERS drugcharacterization semantics: PASS")


if __name__ == "__main__":
    main()
