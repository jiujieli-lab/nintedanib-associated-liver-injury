# Public preclinical nintedanib perturbation module

This directory contains the frozen derived tables, audit records, and manuscript-ready text supporting main **Figure 6**. The inferential hierarchy is fixed:

1. **Independent animals:** GSE278200 and GSE308578 support animal-level gene and module inference.
2. **Pooled-tissue slices:** GSE120804 and GSE120679 support exploratory within-pool slice reproducibility only.
3. **Ambiguous replicate provenance:** PXD024058 supports descriptive direction checks only.

None of the lung datasets is interpreted as liver-toxicity evidence. GSE120804 is a cholestatic liver-context perturbation, not a nintedanib DILI model.

## Verify and regenerate

Verify the complete extracted release from the package root without modifying it:

```bash
python Code/Core/verify_release_capsule.py
```

The original analysis used deterministic seed 20260904. The frozen module script is supplied at `Code/Public_Preclinical_Core/analyze_public_preclinical.py`, but full raw-data regeneration requires downloading the public inputs and recreating the development-relative `raw/` and `metadata/` layout recorded in the manifest. The compact package does not redistribute those large third-party mirrors and does not represent this module as a zero-configuration raw rebuild. No sample was removed after inspecting the reanalysis outcome.

## Primary deliverables

- `Methods_Results_Preclinical.md`: manuscript-ready methods, results, source links, and evidence boundaries.
- `Figure_6_caption.md`: submission-ready A–J legend.
- `dataset_screening_ledger.csv`: systematic search and inclusion/exclusion audit.
- `Code/Public_Preclinical_Core/analyze_public_preclinical.py`: frozen analysis and plotting source code; external inputs and path reconstruction are required for raw regeneration.
- `results/validation_report.json`: automated computational QA; `validation_passed` must be `true`.
- `results/analysis_manifest.json`: versions, seed, file sizes, and SHA-256 hashes; the manifest correctly excludes its own mutable self-hash.
- `Figures/Figure_6_Public_Preclinical_Constraints.{pdf,svg,png,tiff}`: vector and high-resolution submission formats at package root.

## Core result tables

- `results/GSE278200_differential_expression.csv.gz`, `GSE278200_candidate_results.csv`, `GSE278200_module_contrasts.csv`, and `GSE278200_reversal_geometry.csv.gz`.
- `results/GSE308578_differential_expression.csv.gz`, `GSE308578_candidate_results.csv`, `GSE308578_module_contrasts.csv`, `GSE308578_reversal_geometry.csv.gz`, and `GSE308578_leave_one_animal_out_summary.csv`.
- `results/GSE120804_candidate_results.csv` and `GSE120804_module_contrasts.csv`; all p/q fields are slice-level exploratory quantities.
- `results/GSE120679_candidate_results.csv` and `GSE120679_module_contrasts.csv`; all p/q fields are slice-level exploratory quantities.
- `results/PXD024058_protein_effects_descriptive.csv.gz`; no inferential p/q values are supplied.
- `results/cross_model_candidate_module_effects.csv`: exact source table for Figure 6J.

## Source-level adjudications

- GSE308578: GEO SOFT and the source ARRIVE Supplementary Figure S1 agree on CTRL 10 / vehicle 16 / nintedanib 17. The paper's main Figure 6 caption swaps the latter two counts. No sample was relabelled.
- GSE120804 and GSE120679: original Methods explicitly state that tissues or PCLSs from three rats were pooled before randomized slice treatment. Slice labels are not animal identifiers.
- PXD024058: the Methods describe three selected lung tissues per group, whereas the Results call the corresponding columns technical replicates. The public files cannot resolve this conflict.

## Interpretation lock

The larger-magnitude GSE278200 pulmonary response and lower-magnitude GSE308578 response are intentionally retained together. Their difference is a model-response boundary, not a reason to select the favorable dataset. The liver PCLS shifts demonstrate hepatic context dependence but cannot confirm or refute clinical DILI causality.
