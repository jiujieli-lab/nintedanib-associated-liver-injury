# Fresh analysis of public human DILI serum-proteomics source data

This directory contains an independent, from-scratch re-analysis of the source-data workbook accompanying Ravindra et al., *Nature Communications* (2023), DOI: [10.1038/s41467-023-36858-6](https://doi.org/10.1038/s41467-023-36858-6). The input SHA-256 is recorded in `analysis_manifest.json`.

## Scientific role in the nintedanib–ILD project

These cohorts contain adjudicated DILI caused by multiple drugs, not nintedanib-specific cases. They therefore provide a **human distal liver-injury phenotype anchor**. They do not, by themselves, identify a nintedanib-specific susceptibility mechanism or causal intervention target. Drug specificity must come from the project’s pharmacovigilance, direct nintedanib perturbation, target-pharmacology and liver-cell network layers.

## Primary conclusions from the available workbook

- FBP1 was the only one of 13 candidates that remained different between DILI onset (DO) and acute non-DILI liver injury onset (NDO) after Benjamini–Hochberg correction: Cliff’s delta for DO minus NDO = −0.505 (95% bootstrap CI −0.724 to −0.288), q = 4.69×10⁻⁴. With NDO as the positive direction, the apparent univariate AUC was 0.753 (95% bootstrap CI 0.644–0.862). The half-minimum substitution sensitivity analysis remained significant (q = 0.00156).
- FBP1 also separated DO from healthy volunteers (HV): Cliff’s delta = 0.922 (95% CI 0.847–0.978), q = 4.02×10⁻¹⁵; apparent univariate AUC = 0.961 (95% CI 0.923–0.989). This comparison is partly influenced by differential non-detection, which is explicitly reported in `missingness_stats.csv`.
- The discovery and confirmatory directions agreed descriptively for 11 of 12 evaluable proteins in DO versus NDO; HPD was excluded because its entire discovery series is an exact copy of ALDOB in the source workbook. Across-protein direction agreement and the effect-rank correlation (Spearman rho = 0.561) are descriptive because proteins are correlated features drawn from the same source cohorts, not independent Bernoulli or sampling units.
- ACO1, ALDOB, ASS1, CPS1, DMGDH, FAH, GSTA1, HPD and OTC showed strong DO–HV separation and lower values at follow-up in the available unpaired source series. These are best interpreted as general hepatocellular injury/recovery markers, not DILI-etiology-specific markers.
- Within DO cases, ALT correlated most strongly with ASS1 (rho = 0.873), ACO1 (0.800), ALDOB (0.798), HPD (0.750) and FAH (0.740), all FDR-significant. FBP1 retained a moderate within-DO association (rho = 0.425, q = 1.59×10⁻⁴). The pooled HV+DO correlations are larger and should not be interpreted as within-patient severity relationships because between-group separation contributes to them.
- No liver-zone score comparison survived FDR correction (all q ≥ 0.182). The zone layer is descriptive/exploratory in this workbook.

## Statistical specification

- Primary two-group inference: two-sided Mann–Whitney test; Benjamini–Hochberg FDR within each cohort/contrast across proteins.
- Primary effect size: Cliff’s delta with 2,000 stratified bootstrap resamples and percentile 95% CI.
- Secondary estimates: median difference; confirmatory median ratio; Welch test and Hedges’ g after a documented log2 transform; univariate apparent ROC AUC with bootstrap CI.
- Non-detects: excluded in complete-case primary analysis, tested by differential-missingness Fisher tests, and re-analysed using half of the protein-specific minimum positive value. Numeric ND substitutions in Supplementary Figure 7 are excluded from primary correlations and included only as sensitivity analysis.
- Correlations: Spearman rho in pooled HV+DO, DO-only and HV-only strata, with bootstrap CI and BH correction.
- Follow-up comparisons are unpaired because no subject identifiers were supplied.
- Random seed: 20260904. The analysis is rerunnable with `../../.venv/bin/python analyze_dili_proteomics.py` from this directory or `.venv/bin/python fresh_analysis/proteomics/analyze_dili_proteomics.py` from the project root.

## Critical data-integrity findings

1. The discovery ALDOB and HPD series are exact duplicates across all 50 positions, including missingness; HPD discovery evidence is not counted as independent.
2. Figure 3 labels DMGDH HV as n=60 but provides 61 numeric values. The declared first 60 were retained for the primary table and the discrepancy is recorded.
3. Supplementary Figure 2 follow-up column lengths disagree with the cohort counts stated in the paper and contain an invalid `NDF` text value in the HPD NDF numeric column. All supplied follow-up cells are retained as source observations, but recovery evidence is explicitly exploratory.
4. PCK2 contains an unexplained `3077.0266*` value. Rank-based conclusions were stable when it was excluded; parametric estimates were more sensitive (`pck2_outlier_sensitivity.csv`).
5. Several protein arrays contain explicit ND values. Supplementary Figure 7 silently replaces HPD/OTC/PCK2 NDs with numeric constants in log2 space; these are excluded from primary correlations.
6. Cross-protein row order is not a participant identifier. Direct comparison of the same protein’s boxplot and correlation blocks shows identical value multisets but different order for multiple proteins. A cross-protein participant matrix, multivariable classifier, nested CV and paired onset–follow-up tests would therefore fabricate linkage and were not performed.

## Output map

| File | Purpose |
|---|---|
| `discovery_tidy.csv` | Long-form discovery values for 13 candidates |
| `confirmatory_tidy.csv` | Long-form confirmatory values; canonical onset plus available follow-up series |
| `effect_estimates.csv` | All univariate contrasts, effect sizes, CIs, p values, FDR and AUCs |
| `left_censor_sensitivity.csv` | Half-minimum substitution sensitivity analyses |
| `missingness_stats.csv` | Detection/missingness contingency analyses |
| `alt_correlation_pairs.csv` | Within-block biomarker–ALT pairs, with ND-substitution flags |
| `correlation_stats.csv` | Primary and sensitivity Spearman results |
| `candidate_evidence_table.csv` | Descriptive, rule-based evidence table for downstream integration |
| `cross_cohort_concordance_*.csv` | Discovery-to-confirmatory directional/effect concordance |
| `zonal_scores_tidy.csv`, `zonal_*csv` | Deduplicated zone-score values and exploratory summaries |
| `source_block_counts.csv` | Declared/observed sample-count and invalid-cell audit |
| `row_order_linkage_qc.csv` | Evidence that positional cross-protein joins are unsafe |
| `discovery_duplicate_series_qc.csv` | Exact duplicated discovery series |
| `pck2_outlier_sensitivity.csv` | Retain/exclude-starred-value comparison |
| `multivariable_model_feasibility.csv` | Formal record of why participant-level ML was not run |
| `data_dictionary.csv` | Field definitions |
| `xlsx_audit.md` | Full sheet-by-sheet structural audit |

## Interpretation boundary

The strongest defensible proteomic conclusion is that FBP1 recapitulates a circulating phenotype reported by the source study and distinguishes adjudicated DILI from alternative acute liver injury among its 13 source-selected candidates. It is a candidate **phenotype anchor**, not a nintedanib-specific causal target. Intervention priority requires independent evidence that a candidate is engaged by nintedanib or lies downstream of its hepatic target network, localization to a relevant human liver cell state, and graded bidirectional functional testing with rescue; unsigned virtual perturbation alone cannot assign direction.
