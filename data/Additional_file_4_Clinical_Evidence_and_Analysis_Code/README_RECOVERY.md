# Computational archive provenance

This release preserves every scientific path in the original 578-entry integrity record. Of these files, 570 retain their original SHA-256 hash. Eight auxiliary paths were reconstructed or updated; their sources, historical hashes, current hashes and verification details are listed in RECOVERY_PROVENANCE.json and CURRENT_FILE_PROVENANCE.json. The frozen scientific result tables have not been replaced by replay outputs.

## Independent hepatic analysis

The GEO series matrix and SOFT records, sample metadata, expression matrices, feature annotation and original result tables retain their archived values. The publicly downloaded mouse gene-annotation file differs from the earlier compressed file, but all 13392 feature-to-symbol mappings used by the analysis are identical. Re-execution of the original hepatic script produces all 12 checked result CSV files byte-identically, including genome-wide contrasts, fixed candidates, module analyses, PCA, sample QC and disjoint-reference sensitivity. Reconstructed source-reference, numerical-QA and provenance records document this comparison.

## Pulmonary sensitivity analysis

Both deposited count files and all three archived sensitivity-result CSV files retain their original hashes. The two log-CPM matrices, reconstruction script and analysis summary were rebuilt from the deposited counts and the recorded analysis specification. The script writes replay outputs to a separate user-specified directory and leaves archived results unchanged.

A complete replay verified 506 disjoint-reference splits and 1070 whole-animal allocations. All observed statistics and the six reported P and q values match the archived outputs. All split assignments and values match within the documented numerical tolerance. Two rat permutation Spearman coefficients differ by at most 5.49 × 10⁻⁸ because of floating-point tie handling; restoration differences are at most 3.56 × 10⁻¹⁵. These differences do not alter any reported inferential result. The original permutation results remain the publication source.

The recovery records are an integrity and reproducibility description. They do not add biological replicates, improve external validation or establish a causal target.
