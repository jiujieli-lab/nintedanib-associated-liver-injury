# Runtime and input notes

## 1. Current upstream analyses

The five original scientific modules retain their original package boundaries under `data/`. Their current `Analysis/` scripts expect a combined working root containing `Analysis/` and `Source_Data/`. To reproduce them, copy those directories from the five modules into a new working directory, preserving relative paths. Static inspection found no conflicting paths across those five source modules. Do not merge their separately named package metadata files or overwrite the frozen attachment.

The combined data include the clinical effect tables, candidate pharmacology, historical numerical audit inputs, three liver network matrices and cell metadata, independent hepatic matrix/annotation/metadata, and two pulmonary count files. The network and integration scripts expose explicit input/output arguments. The independent hepatic scripts read and write beside their own files and therefore must run from a disposable working copy. Commands retained in old README files may refer to a formerly combined root.

For the pulmonary sensitivity script, the current interface requires a distinct output directory:

```bash
python Analysis/pulmonary_sensitivity/analyze_shared_reference.py --source-root . --output-dir replay/pulmonary_sensitivity
```

Its frozen results are not overwritten. Network diffusion, primary integration, donor sensitivity, held-out association and historical numerical checks have separate entry points documented in their original module README files; running them can be computationally substantial.

## 2. Historical source-specific workflows

The original `Code/` files are archived source code. Several were written for a development tree containing `fresh_data`, `fresh_analysis`, `preclinical/raw` and `preclinical/metadata`; they cannot all run directly from the consolidated storage locations. `Path_Audit.csv` lists static root assumptions and possible network calls.

Examples:

- `Code/Human_DILI_Proteomics/analyze_dili_proteomics.py` expects the source workbook `fresh_data_dili_source.xlsx` at its original project root and writes derived outputs beside the script. The frozen processed clinical tables are supplied; the workbook itself was not present in the recovered source modules.
- `Code/GSE151374/analyze_gse151374.py` and `Code/PXD052594/analyze_pxd052594.py` assume a historical dataset-specific `fresh_analysis/preclinical/<dataset>/code/` layout. Reconstruct that layout in a working copy before using their raw-data path logic.
- `Code/Public_Preclinical_Core/analyze_public_preclinical.py` requires the source `preclinical/raw/` and `preclinical/metadata/` arrangement documented in its manifests.
- Several `Code/Core/` scripts build historical manuscripts, supplementary files or release packages. They are retained for provenance and should not be used as scientific entry points or to overwrite the current submission documents.

Original script bytes were preserved; no analysis algorithm, numerical parameter, seed or result table was altered to remove these path assumptions.

## 3. Inputs not redistributed in the recovered modules

Full raw regeneration of historical source-specific modules requires some external inputs listed in the original source manifests, including the Ravindra source workbook; full GSE115469 expression/cell-type files; original GSE299128 counts; certain original ChEMBL activity JSON and Sci-Plex pages; the complete GSE70138 GCTX and annotations; the 18 original GSE151374 H5 files; the PXD052594 source workbook; several original GTF and GSE120804/GSE120679 inputs; and original STRING TSV files. Frozen derived data and completed outputs are retained where supplied. These are input-availability limits, not evidence that the reported analysis was rerun during consolidation.

The referenced third-party `R_Script_Lauer_et_al_2024.R` is not included. `reprocess_CEL_standard_RMA.R` is listed in the historical exclusion record; the accompanying analysis status states that raw CEL RMA reprocessing was not executed. No such script or analysis has been fabricated. `Public_Data_Reproducibility_Notebook.ipynb` is not included; its original generator is supplied.

## 4. Targeted perturbation and binding modules

The merged package keeps the exact paths used by `new_virtual/run_target_pair_perturbation.py` and `binding_analysis/analyze_binding_evidence.py`. Their principal frozen inputs are included. `new_virtual/validate_extension.py` requires the generated Figure S11 TIFF and nine PNG panels, so `make_figure_s11.py` must be run first in a working copy if those assets are needed. The historical Figure S11 source name is preserved; publication-panel numbering is recorded separately.

## 5. Screening, docking and figure modules

- `candidate_library/retrieve_*.py` performs online retrieval; two scripts import `retrieve_library.py`, whose top-level retrieval code also executes. Review these scripts before importing them. The frozen source records allow inspection without retrieval.
- `screening/code/summarize_boltz.py` and `summarize_inductive.py` read completed local provider responses. They do not submit new prediction jobs.
- `docking/run_docking.py` returns an existing `result.json` instead of recomputing that run. A fresh calculation requires a deliberately prepared working copy and review of the original code.
- Primary docking has 96 completed searches; the separate redocking control and three apo-PDGFRA sensitivity searches remain separate. Search seeds are computational repeats.
- `render_molecular_panels.py` and `render_boltzmol_panels.py` create intermediate PNGs required by the screening figure assemblers. Those intermediate figure exports were not in the original screening source ZIP. `build_main_screening.py` also requires the path to `Boltz_All32_Dispositions.csv` as its command-line argument.
- `build_screening_tables.mjs` depends on `@oai/artifact-tool`, a runtime-specific package that is not vendored here. The source CSV files and submitted workbook are available independently.

## 6. Environments

Do not merge the historical requirements files into a single claimed tested environment. `Code/requirements-analysis.txt` and `Code/requirements-virtual-knockout.txt` describe different recorded environments. The revised network analysis runtime is documented in `Analysis/revision_runtime_versions.json`; its source imports `threadpoolctl`, which is absent from the older analysis requirements list but recorded in that runtime manifest.

The targeted extension's recorded numerical/plotting versions are retained under `provenance/original_package_metadata/Additional_file_6_Targeted_Analysis_Code/ENVIRONMENT.json`. The docking protocol records AutoDock Vina 1.2.7, Meeko 0.8.0 and RDKit 2026.03.6. Receptor repair additionally uses PDBFixer/OpenMM; that preparation step is not production molecular dynamics. The code inventory and imported-module table provide a static dependency map. No fresh dependency installation or full execution test was performed for consolidation.
