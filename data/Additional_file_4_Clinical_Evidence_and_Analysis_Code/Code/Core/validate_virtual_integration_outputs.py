#!/usr/bin/env python3
"""Programmed QA for the frozen Figure 4/5 integration deliverables.

The report intentionally contains exactly 44 high-level checks. Hash checks are
grouped by artifact family so methodological checks can be added without
changing the agreed QA count.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None


ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "Figures"
INT = ROOT / "Source_Data" / "Integration"
PHARMACOLOGY = ROOT / "Source_Data" / "Pharmacology"
LINCS = ROOT / "Source_Data" / "LINCS"
VIRTUAL_KO = ROOT / "Source_Data" / "Virtual_Knockout"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def add(checks: list[dict[str, object]], name: str, passed: bool, detail: object) -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": detail})


def svg_text(path: Path) -> list[str]:
    root = ET.parse(path).getroot()
    return ["".join(node.itertext()).strip() for node in root.iter() if node.tag.endswith("text")]


def hash_family(entries: dict[str, str]) -> tuple[bool, dict[str, object]]:
    detail: dict[str, object] = {}
    passed = True
    for rel, expected in entries.items():
        path = ROOT / rel
        observed = sha256(path) if path.exists() else None
        item_ok = observed == expected
        passed &= item_ok
        detail[rel] = {"passed": item_ok, "expected": expected, "observed": observed}
    return passed, detail


def figure_hash_family(fig_key: str, entries: dict[str, str]) -> tuple[bool, dict[str, object]]:
    pattern = "Figure_4*" if fig_key == "figure4" else "Figure_5*"
    suffixes = {
        "svg": ".svg",
        "pdf": ".pdf",
        "png_600dpi": ".png",
        "tiff_600dpi": ".tiff",
    }
    detail: dict[str, object] = {}
    passed = True
    for fmt, expected in entries.items():
        matches = [p for p in FIG.glob(pattern) if p.name.endswith(suffixes[fmt])]
        observed = sha256(matches[0]) if len(matches) == 1 else None
        item_ok = len(matches) == 1 and observed == expected
        passed &= item_ok
        detail[fmt] = {
            "passed": item_ok,
            "expected": expected,
            "observed": observed,
            "matches": [p.name for p in matches],
        }
    return passed, detail


def main() -> None:
    checks: list[dict[str, object]] = []
    manifest_path = INT / "integration_analysis_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    add(
        checks,
        "manifest_exists",
        manifest_path.exists(),
        manifest_path.relative_to(ROOT).as_posix(),
    )
    script = ROOT / manifest["analysis_script"]
    validation_script = ROOT / manifest["validation_script"]
    scripts_ok = (
        sha256(script) == manifest["analysis_script_sha256"]
        and sha256(validation_script) == manifest["validation_script_sha256"]
    )
    add(
        checks,
        "analysis_and_validation_script_hashes",
        scripts_ok,
        {"analysis": sha256(script), "validation": sha256(validation_script)},
    )

    ok, detail = hash_family(manifest["input_hashes"])
    add(checks, "all_input_hashes_match", ok, detail)
    output_entries = manifest["output_table_hashes"] | manifest["text_artifact_hashes"]
    ok, detail = hash_family(output_entries)
    add(checks, "all_output_and_text_artifact_hashes_match", ok, detail)
    for fig_key in ["figure4", "figure5"]:
        ok, detail = figure_hash_family(fig_key, manifest["figure_hashes"][fig_key])
        add(checks, f"{fig_key}_all_format_hashes_match", ok, detail)

    for number, stem in [
        (4, "Figure_4_Liver_Localization_and_Virtual_Knockout"),
        (5, "Figure_5_Multi_Evidence_Target_Prioritization"),
    ]:
        texts = svg_text(FIG / f"{stem}.svg")
        counts = {letter: texts.count(letter) for letter in "ABCDEFGHIJ"}
        add(checks, f"figure_{number}_exactly_ten_panel_letters", all(v == 1 for v in counts.values()), counts)
        for ext in ["png", "tiff"]:
            path = FIG / f"{stem}.{ext}"
            with Image.open(path) as im:
                dpi = tuple(float(v) for v in im.info.get("dpi", (0, 0)))
                add(
                    checks,
                    f"figure_{number}_{ext}_600dpi",
                    min(dpi) >= 599,
                    {"dpi": dpi, "pixels": tuple(int(v) for v in im.size)},
                )

    fig5_text = " ".join(svg_text(FIG / "Figure_5_Multi_Evidence_Target_Prioritization.svg"))
    add(checks, "figure5e_title_correct", "FBP1 is the sole cross-context anchor overlap" in fig5_text, "requested title")
    add(checks, "figure5e_global_rho_annotations", "rho=0.09" in fig5_text and "rho=0.14" in fig5_text, "rho=0.09/0.14")
    fig4_text = " ".join(svg_text(FIG / "Figure_4_Liver_Localization_and_Virtual_Knockout.svg"))
    add(checks, "figure4i_exploratory_label", "Exploratory matched-null" in fig4_text, "explicit exploratory label")

    sci = pd.read_csv(PHARMACOLOGY / "sciplex_dose_concordance.csv")
    both = sci.loc[sci["observed_in_both_doses"].eq(True) & sci["direction_concordant"].eq(True), "gene_symbol"]
    anchors = {"ACO1", "ASS1", "FAH", "CPS1", "ALDOB", "HPD", "OTC", "DMGDH", "GSTA1", "FBP1", "PCK2", "CES1", "LECT2"}
    overlap = sorted(set(both) & anchors)
    add(checks, "fbp1_sole_anchor_in_both_sciplex_dose_lists", overlap == ["FBP1"], overlap)

    concordance = pd.read_csv(LINCS / "sciplex_hepg2_concordance.csv")
    c1 = concordance.loc[(concordance["dose_pair"].eq("HepG2 1.11 µM vs MCF7 1uM")) & concordance["subset"].eq("all_overlap")].iloc[0]
    c10 = concordance.loc[(concordance["dose_pair"].eq("HepG2 10 µM vs MCF7 10uM")) & concordance["subset"].eq("all_overlap")].iloc[0]
    add(
        checks,
        "global_overlap_correlations_match",
        round(c1.spearman_rho, 2) == 0.09 and round(c10.spearman_rho, 2) == 0.14,
        {"rho_1uM": c1.spearman_rho, "n_1uM": int(c1.n_overlap), "rho_10uM": c10.spearman_rho, "n_10uM": int(c10.n_overlap)},
    )

    downstream = pd.read_csv(VIRTUAL_KO / "stable_top_downstream_genes.csv")
    drug_targets = pd.read_csv(VIRTUAL_KO / "virtual_knockout_stability_summary.csv").query("ko_class == 'drug_target'")["ko_label"]
    anchor_hits = downstream.loc[
        downstream["gene"].isin(anchors)
        & downstream["top20_seed_frequency"].ge(2)
        & downstream["ko_label"].isin(drug_targets)
    ]
    add(checks, "no_stable_anchor_top20_after_drug_target_ko", anchor_hits.empty, int(len(anchor_hits)))

    matched = pd.read_csv(INT / "virtual_ko_expression_matched_null_stability.csv")
    add(checks, "no_matched_null_bh_signal_repeats_two_seeds", not matched["n_seed_bh_q_lt_0_05"].ge(2).any(), int(matched["n_seed_bh_q_lt_0_05"].max()))

    ranking = pd.read_csv(INT / "candidate_consensus_ranking.csv")
    target = ranking.loc[ranking["candidate_class"].eq("exposure_proximal_target")].copy()
    phenotype = ranking.loc[ranking["candidate_class"].eq("downstream_DILI_anchor")].copy()
    add(checks, "candidate_roles_separate", set(ranking["candidate_class"]) == {"downstream_DILI_anchor", "exposure_proximal_target"}, sorted(ranking["candidate_class"].unique()))
    add(checks, "no_tier1_claim", not ranking["evidence_tier"].astype(str).str.startswith(("T1", "Tier 1")).any(), "none assigned")

    caption4 = (FIG / "Figure_4_Liver_Localization_and_Virtual_Knockout_caption.md").read_text()
    caption5 = (FIG / "Figure_5_Multi_Evidence_Target_Prioritization_caption.md").read_text()
    add(checks, "virtual_ko_unsigned_guardrail", "unsigned network displacement" in caption4, "caption statement")
    add(checks, "phenotype_and_target_tiers_explicitly_separated", "P labels describe phenotype evidence and are not drug-target tiers" in caption5, "caption statement")
    add(checks, "same_compartment_no_splicing_stated", "values from different compartments are never spliced" in caption5, "caption statement")

    tier_sets = {
        "P1": set(phenotype.loc[phenotype["evidence_tier"].str.startswith("P1"), "gene_symbol"]),
        "P2": set(phenotype.loc[phenotype["evidence_tier"].str.startswith("P2"), "gene_symbol"]),
        "P3": set(phenotype.loc[phenotype["evidence_tier"].str.startswith("P3"), "gene_symbol"]),
    }
    expected_tiers = {
        "P1": {"FBP1"},
        "P2": {"ACO1", "ASS1", "FAH", "CPS1", "ALDOB", "HPD", "OTC", "DMGDH", "GSTA1"},
        "P3": {"PCK2", "CES1", "LECT2"},
    }
    for tier in ["P1", "P2", "P3"]:
        add(checks, f"{tier.lower()}_phenotype_set", tier_sets[tier] == expected_tiers[tier], sorted(tier_sets[tier]))

    key_ranks = phenotype.set_index("gene_symbol").loc[["GSTA1", "FBP1", "OTC"], ["consensus_rank", "loso_rank_best", "loso_rank_worst"]]
    expected_key = {"GSTA1": (3, 1, 8), "FBP1": (7, 5, 9), "OTC": (8, 7, 9)}
    observed_key = {gene: tuple(int(v) for v in key_ranks.loc[gene]) for gene in key_ranks.index}
    add(checks, "p_track_key_ranks_and_loso_ranges", observed_key == expected_key, observed_key)

    target_key_ranks = target.set_index("gene_symbol").loc[["FGFR1", "FLT1", "FGFR3"], ["consensus_rank", "loso_rank_best", "loso_rank_worst"]]
    expected_target_key = {"FGFR1": (3, 2, 7), "FLT1": (7, 5, 8), "FGFR3": (9, 5, 13)}
    observed_target_key = {gene: tuple(int(v) for v in target_key_ranks.loc[gene]) for gene in target_key_ranks.index}
    add(checks, "t_track_key_ranks_after_same_compartment_rerun", observed_target_key == expected_target_key, observed_target_key)

    t2 = set(target.loc[target["evidence_tier"].str.startswith("T2"), "gene_symbol"])
    add(checks, "exploratory_t2_gate_set", t2 == {"FGFR1"}, sorted(t2))
    downgraded = target.set_index("gene_symbol").loc[["FLT1", "FGFR3"], "evidence_tier"].astype(str)
    add(checks, "flt1_and_fgfr3_downgraded_to_t3", downgraded.str.startswith("T3").all(), downgraded.to_dict())
    strict = set(target.loc[target["strict_all_seed_T2_sensitivity_gate"].eq(True), "gene_symbol"])
    add(checks, "strict_all_seed_t2_is_empty", not strict, sorted(strict))

    target_indexed = target.set_index("gene_symbol")
    add(checks, "t2_has_same_compartment_gate", bool(target_indexed.loc["FGFR1", "same_compartment_primary_gate"]), target_indexed.loc["FGFR1", "joint_context_compartment"])
    add(checks, "t2_has_strong_pharmacology", bool(target_indexed.loc["FGFR1", "strong_pharmacology_gate"]), target_indexed.loc["FGFR1", "pharmacology_tier"])

    primary_recomputed = target["strong_pharmacology_gate"].eq(True) & (
        target["same_compartment_primary_gate"].eq(True) | target["hepg2_alternative_gate"].eq(True)
    )
    add(checks, "exploratory_t2_gate_recomputes_exactly", primary_recomputed.equals(target["primary_T2_gate"].astype(bool)), int(primary_recomputed.sum()))
    strict_recomputed = target["strong_pharmacology_gate"].eq(True) & (
        target["same_compartment_all_seed_sensitivity_gate"].eq(True) | target["hepg2_alternative_gate"].eq(True)
    )
    add(checks, "strict_t2_gate_recomputes_exactly", strict_recomputed.equals(target["strict_all_seed_T2_sensitivity_gate"].astype(bool)), int(strict_recomputed.sum()))

    valid = target["hepg2_measurement_valid"].eq(True)
    valid_expected = target["is_landmark"].eq(1) & target["constant_profile"].eq(False) & target["spearman_rho"].notna()
    add(checks, "hepg2_measurement_validity_recomputes_exactly", valid.equals(valid_expected), sorted(target.loc[valid, "gene_symbol"]))
    alternative_expected = valid & target["spearman_rho"].abs().ge(0.70)
    add(checks, "hepg2_alternative_gate_recomputes_exactly", target["hepg2_alternative_gate"].eq(True).equals(alternative_expected), sorted(target.loc[alternative_expected, "gene_symbol"]))
    strong_alt = set(target.loc[target["strong_pharmacology_gate"].eq(True) & target["hepg2_alternative_gate"].eq(True), "gene_symbol"])
    add(checks, "hepg2_alternative_adds_no_strong_pharmacology_target", not strong_alt, sorted(strong_alt))

    pairs = pd.read_csv(INT / "target_same_compartment_context_all_pairs.csv")
    pair_gate_expected = pairs["mean_detection_fraction"].ge(0.10) & pairs["n_seeds"].eq(3) & pairs["vko_within_compartment_rank"].ge(0.75)
    add(checks, "all_pair_primary_gates_recompute_exactly", pair_gate_expected.equals(pairs["same_compartment_primary_gate"].astype(bool)), int(pair_gate_expected.sum()))

    selected = pairs.sort_values(
        ["ko_label", "same_compartment_primary_gate", "same_compartment_context_score", "mean_detection_fraction"],
        ascending=[True, False, False, False],
    ).drop_duplicates("ko_label").set_index("ko_label")
    joint_match = all(
        target_indexed.loc[g, "joint_context_compartment"] == selected.loc[g, "compartment"]
        and np.isclose(target_indexed.loc[g, "same_compartment_context_score"], selected.loc[g, "same_compartment_context_score"], equal_nan=True)
        for g in target_indexed.index
    )
    add(checks, "joint_context_matches_best_single_compartment_pair", joint_match, {g: selected.loc[g, "compartment"] for g in ["FGFR1", "FLT1", "FGFR3"]})
    add(checks, "liver_context_equals_same_compartment_maximin", np.allclose(target["liver_context"], target["same_compartment_context_score"], equal_nan=True), "all targets")

    shortlist = pd.read_csv(INT / "tier2_candidate_shortlist.csv")
    add(checks, "tier2_shortlist_contains_fgfr1_only", shortlist["gene_symbol"].tolist() == ["FGFR1"], shortlist["gene_symbol"].tolist())

    gate_table = pd.read_csv(INT / "target_tier_gate_sensitivity.csv")
    phenotype_table = pd.read_csv(INT / "phenotype_anchor_tiers.csv")
    domain_table = pd.read_csv(INT / "candidate_domain_scores_long.csv")
    loso_table = pd.read_csv(INT / "candidate_leave_one_source_out_ranks.csv")
    aligned = (
        set(gate_table["gene_symbol"]) == set(target["gene_symbol"])
        and set(phenotype_table["gene_symbol"]) == set(phenotype["gene_symbol"])
        and set(domain_table["gene_symbol"]) == set(ranking["gene_symbol"])
        and set(loso_table["gene_symbol"]) == set(ranking["gene_symbol"])
    )
    add(checks, "ranking_shortlist_domain_loso_candidate_sets_aligned", aligned, {"ranking": len(ranking), "domains": domain_table["gene_symbol"].nunique(), "loso": loso_table["gene_symbol"].nunique()})

    if len(checks) != 44:
        raise RuntimeError(f"QA contract requires exactly 44 checks; assembled {len(checks)}")
    report = {"status": "PASS" if all(c["passed"] for c in checks) else "FAIL", "n_checks": len(checks), "checks": checks}
    (INT / "integration_validation_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
