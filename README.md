# Additional file 6. Analysis code and source data

This single attachment consolidates the supplied targeted-analysis and screening code with the recovered original analysis modules and historical figure-construction sources. It contains 113 unique original source-code files: 109 Python files and four JavaScript modules. Five byte-identical source-code copies are mapped to their retained canonical files. All retained original scripts are byte-identical to their source archives.

## Start here

Run the offline, read-only integrity check from the extracted package root:

```bash
python verify_package.py
```

The verifier uses the Python standard library. It checks packaged file hashes, the original-code inventory, Python syntax and archived output counts. It does not import analysis modules, modify source data, run statistical models, launch docking or call external services. Integrity checks do not establish biological validity or demonstrate that the full analysis has been rerun.

## Directory map

| Location | Contents |
| --- | --- |
| `data/Additional_file_4_Clinical_Evidence_and_Analysis_Code/` | Original clinical proteomics, pharmacovigilance, perturbation and target-integration source code, data, network analyses and validation records. |
| `data/Additional_file_5_Liver_and_Lung_Single_Cell_Data/` | Frozen liver/lung single-cell inputs and derived tables. |
| `data/Additional_file_6_Preclinical_Response_and_Proteomics_Data/` | Frozen preclinical, lung-proteomics and whole-blood inputs/results. |
| `data/Additional_file_7_External_Hepatic_Evaluation/` | Independent hepatic evaluation and held-out rank-association code and data. |
| `data/Additional_file_8_Pulmonary_Sensitivity_Data/` | Pulmonary shared-reference sensitivity code, raw counts and results. |
| `new_virtual/` | Target-pair graph perturbation, calibrations, tables and plotting source. |
| `binding_analysis/` | Nintedanib binding-evidence reanalysis, frozen activity records and raw property response. |
| `jcm_revision/candidate_library/` | Compound identities, assay audit, source records and clinical-context tables. |
| `jcm_revision/screening/` and `jcm_revision/boltz-experiments/` | Provider responses, all screening dispositions, exported cropped pocket models and summarization scripts. |
| `jcm_revision/docking/` | Receptor preparation, search code, configurations, scores, poses, logs and validation. |
| `jcm_revision/figures/` and `jcm_revision/tables/` | Screening figure, molecular-panel and workbook construction scripts. |
| `source_indices/` | Archived mapping between analytical source panels and final manuscript panels. |
| `provenance/author_archive/` | Historical authoring, reference-audit, figure-assembly and packaging sources; these are distinct from scientific analysis entry points. |
| `documentation/` | Code inventory, source-to-package hash map, module coverage, runtime notes, path audit and consolidation records. |

In the supplementary PDF, the path `Code/GSE151374/analyze_gse151374.py` is relative to the module root `data/Additional_file_4_Clinical_Evidence_and_Analysis_Code/`.

The `Additional_file_4` through `Additional_file_8` labels inside `data/` preserve original directory names required by source scripts. They are internal provenance labels, not separate current submission attachments. Historical README files, output filenames, figure labels and archived statuses retain their original scope. Use this README for the consolidated package and the current manuscript for publication numbering. In particular, an early binding-analysis manifest stating that docking or provider inference was not performed describes that earlier module only; the subsequently completed screening and docking outputs are supplied under `jcm_revision/`.

## Scope and execution limits

All supplied and recovered source code is retained, including distinct historical versions. The package includes all contents of the five original scientific source-data archives after extraction. Some full public raw-data mirrors and a referenced third-party R script were not present in those archives. Historical scripts can retain development-layout assumptions; the presence of source code does not mean every script is a zero-configuration entry point. See `documentation/Runtime_Notes.md` and `documentation/Module_Coverage.csv` before running anything.

Use a separate working copy for any reproduction. Numerous original scripts write outputs beside their inputs, and retrieval scripts can replace frozen records. The consolidated frozen inputs and outputs have not been recalculated. Provider re-inference requires external accounts and can incur charges. Docking launchers can consume substantial CPU time. No production molecular-dynamics simulation is represented by this archive; `jcm_revision/md/` contains coordinate/contact quality checks only.

## What was consolidated or omitted

- Seven identical shared-path files were stored once, including overlapping frozen inputs and the nintedanib property response.
- Five identical code copies were mapped to their canonical retained files; no distinct source-code version was removed.
- Eleven zero-byte incomplete-download placeholders were removed.
- Five nested source ZIPs were expanded into `data/` rather than nested again.
- The nested original-panel archive was inspected. Its three scripts, layout specifications and audit records are retained. Its 198 raster/PDF assets were excluded from this code attachment to avoid duplicating figure assets. Their names, sizes and hashes are listed in `documentation/Omitted_Panel_Assets.csv`; figure-assembly code requiring them needs the original panel assets. The supplied final submission figures remain separate publication files.

`documentation/Source_File_Map.csv` records every source member and its disposition. `SHA256SUMS.txt` is the integrity list for this consolidated attachment. Older checksum files are preserved as provenance and must not be used as the top-level check for the merged archive.
