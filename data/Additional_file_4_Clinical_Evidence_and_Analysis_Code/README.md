# Supplementary data and code

This archive accompanies “Clinical proteomics guides calibrated perturbation analysis of nintedanib-associated liver injury”. It contains the frozen source tables, revised analysis inputs and outputs, and analysis code. No new wet-laboratory experiments are claimed.

## Result locations

- `Source_Data/`: source-specific frozen data and original analysis tables. `Source_Data/Integration/` is HISTORICAL comparator material, not current target inference. Old tier labels in these original machine-readable inputs must not be used as current conclusions.
- `Analysis/network_diffusion/`: donor-balanced phenotype-seeded random walk, matched nulls, parameter grid, donor omission and covariance/degree sensitivities.
- `Analysis/integrated_analysis/primary/`: final six-method comparison and candidate ranking. PDGFRA and FLT4 are experimental priorities; no DILI-causal target is confirmed.
- `Analysis/integrated_analysis/vko_*`: numerical run audit. In `target_protein_projection_by_seed.csv`, filter `numerically_nonzero == True` before interpretation; the uncorrected workflow's numerical-zero rank results are invalid for biological claims.
- `Analysis/independent_liver/`: final deposited-matrix GSE125975 analysis. Raw CEL normalization was not executed. Only the 13,392 deposited features support reported expression results; unmeasured features are unavailable.
- `Analysis/heldout_rank_association/deposited_matrix/`: exact 92,378-animal-allocation evaluation of frozen algorithm ranks versus hepatic response magnitude.
- `Analysis/pulmonary_sensitivity/`: raw-feature shared-reference and disjoint-reference pulmonary analyses.
- `Source_Data/Revised_Figure_Panels/`: new panel-source maps and plotted values. Map aliases Network, Integration, Hepatic and Pulmonary refer to the corresponding Analysis directories.
- `Code/`: archived source-specific analysis and figure code. Historical integration/figure-generation scripts reconstruct the original comparator and are not entry points for the revised selection.

## Reproduce the revised analyses

Use a separate working copy of this archive, Python with the scientific packages recorded in `Analysis/revision_runtime_versions.json`, and run from the archive root. Module-specific original manifests remain authoritative for historical calculations. Some complete raw public mirrors are not redistributed; accession and download manifests identify their sources. Starting from the included frozen processed inputs reproduces the revised analyses below.

```bash
python Analysis/integrated_analysis/audit_virtual_knockout_numerics.py --input-root .
python Analysis/network_diffusion/run_network_diffusion.py --package . --output Analysis/network_diffusion
python Analysis/integrated_analysis/run_primary_integration.py --input-root . --network-dir Analysis/network_diffusion
python Analysis/integrated_analysis/build_candidate_bridge_map.py --input-root . --network-dir Analysis/network_diffusion
python Analysis/integrated_analysis/validate_primary_integration.py --input-root . --network-dir Analysis/network_diffusion
python Analysis/pulmonary_sensitivity/analyze_shared_reference.py
python Analysis/independent_liver/analyze_independent_liver.py
python Analysis/independent_liver/batch_sensitivity.py
python Analysis/heldout_rank_association/analyze_frozen_rank_association.py --external-dir Analysis/independent_liver --expression-file Analysis/independent_liver/expression_RMA_symbols.csv.gz --metadata-file Analysis/independent_liver/sample_metadata.csv --ranking-dir Analysis/integrated_analysis/primary --output-dir Analysis/heldout_rank_association/deposited_matrix
```

Run the held-out rank script with `--help` for its explicit paths to the frozen ranking, hepatic matrix, and output directory. Do not refit integration using these external outcomes. The full numerical projection comparator can be reconstructed with `run_integrated_benchmark.py --input-root .`; this diagnostic script flags numerically zero runs and its uncorrected fusion output is not the primary result. Reproducing a script may overwrite the corresponding derived files, so preserve the frozen copy and compare results separately.

No omitted dataset was removed because it disagreed with a preferred target. No gene or cell count is used as an independent animal, donor, or patient denominator. Prospective experimental confirmation is described separately in Supplementary Information.

## Data rights and provenance

Original data remain attributable to their repository depositors and published source studies. Accession links, download URLs, and original source manifests accompany the tables. This archive does not grant rights beyond the original terms. `SHA256SUMS.txt` lists the packaged byte-level checksums; it is an integrity record, not a claim of biological validation.
