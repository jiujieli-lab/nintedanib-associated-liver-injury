# Journal of Translational Medicine supplementary data reconstruction

Extract Additional files 4–8 into the same empty directory, preserving paths. Each is a standalone ZIP archive; do not concatenate the archives. Additional file 3 contains the unchanged supplementary workbook. Additional file 9 contains the vector figure masters and composition code.

## Integrity and provenance

This release retains all 578 original scientific data and code paths. Most files are byte-identical to the archived analysis. The eight reconstructed auxiliary paths are identified in CURRENT_FILE_PROVENANCE.json and RECOVERY_PROVENANCE.json. They comprise an updated public gene-annotation download, associated provenance and numerical-verification metadata, the two restored pulmonary log-CPM matrices, and their reconstruction script and summary. Archived hepatic and pulmonary result tables remain unchanged.

SHA256SUMS.txt contains the current hashes for all 578 scientific paths. Run `sha256sum -c SHA256SUMS.txt` after extraction. The earlier integrity record is preserved in Original_Integrity/SHA256SUMS_20260921.txt; it is historical and will flag the explicitly documented reconstructed paths. JTM_PACKAGE_PART_4.json through JTM_PACKAGE_PART_8.json record each archive's contents and hashes. README_RECOVERY.md records numerical recovery checks and distinguishes byte identity from numerical equivalence.

The current mouse annotation retains every original feature-to-symbol mapping used by the liver analysis. Regeneration from the unchanged GEO matrix reproduces the archived hepatic result tables. The pulmonary reconstruction starts from the unchanged deposited count files and preserves the archived result tables; use a separate output directory for replay. No new biological observations, independent replications or functional-rescue experiments are introduced by this release.

## Analysis entry points

Read the original README.md for the study design and historical analysis organization. For the reconstructed pulmonary replay, use:

```
python Analysis/pulmonary_sensitivity/analyze_shared_reference.py --output-dir replay/pulmonary_sensitivity
```

Keep replay outputs separate from the frozen tables. The clinical, donor-network, integration, independent-liver and external rank-association scripts retain their source-specific environment requirements and statistical units. Original source manifests may refer to historical checksums; CURRENT_FILE_PROVENANCE.json and RECOVERY_PROVENANCE.json resolve the eight reconstructed paths.

## Archive contents

- Additional file 4: Clinical reporting, serum proteins, pharmacology, drug-response data, analysis code, calibrated network/integration results and integrity records.
- Additional file 5: Healthy human liver data from five donors and pooled-library mouse pulmonary immune-cell data at days 7 and 14. The time-point groups contain different animals.
- Additional file 6: Public hepatic/pulmonary response, proteome/phosphoproteome and whole-blood context data.
- Additional file 7: Independent GSE125975 hepatic pharmacodynamic evaluation and frozen-rank association testing. This experiment does not establish drug-induced liver injury.
- Additional file 8: Deposited pulmonary count inputs and shared-reference/disjoint-reference sensitivity analyses.

Repository terms remain applicable. The results prioritize exposure-controlled perturbation and rescue experiments; they do not establish a causal DILI mediator.
