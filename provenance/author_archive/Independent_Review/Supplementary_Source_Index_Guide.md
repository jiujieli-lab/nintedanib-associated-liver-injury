# Final figure routes for retained supplementary source indices

Additional file 3 is the source workbook dated 19 September 2026. Its historical figure identifiers, artifact filenames, source paths and hashes preserve the provenance of the analyses. They are not all final manuscript figure numbers. This guide maps those labels to the completed revision without changing the scientific values or source identifiers.

`Additional3_Final_Figure_Index_Map.csv` supplies all 89 populated data rows from the workbook's `Figure_Source_Map` sheet, including exact worksheet cells, original labels, original artifact/source paths and current reader routes. Header rows 1 and 56 and blank rows 53–55 are excluded. Original paths in columns C, D and H remain explicitly labelled as historical source locators.

| Workbook cells | Retained source label | Final reader route |
|---|---|---|
| Figure_Source_Map!A3:B3 | Figure 1, A–J | Main Figure 1, B–K. Final panel A is the study workflow. |
| Figure_Source_Map!A4:B4 | Figure 2, A–J | Main Figure 2, A–J. |
| Figure_Source_Map!A5:B5 | Figure 3, A–J | Main Figure 3, A–J. |
| Figure_Source_Map!A6:B6 | Figure 4, A–J | **Supplementary Figure S12, B–K**. Final S12A is the localization/calibration workflow. |
| Figure_Source_Map!A7:B16 | Figure 5, A–J | Main Figure 5, B–K, respectively. |
| Figure_Source_Map!A57:B58 | Figure 5, K–L | Main Figure 5, L–M, respectively. |
| Figure_Source_Map!A83:B94 | Figure 6, A–L | Main Figure 6, B–M, respectively. Final panel A is the hepatic-evaluation workflow. |
| Figure_Source_Map rows 17–26, 27–51 and 59–82 | Supplementary Figures S1–S10 and their indexed panels | Retain the corresponding S1–S10 identifiers and panel scopes in Additional file 2. |
| Figure_Source_Map row 52 | PXD052594 phosphoproteome, table-only sensitivity | Remains table-only in Additional file 3; no figure panel is assigned. |
| Figure_Source_Map row 2 | Historical graphical abstract | Source architecture schematic; no final numbered quantitative figure is assigned by this index. |

For the former Figure 4 specifically, `Figure_Source_Map!C6` contains `Figures/Figure_4_Liver_Localization_and_Perturbation_QC.pdf`, and `H6` contains `Main_Manuscript.docx, Figure 4 legend`. Both are historical source locators. The current image is Supplementary Figure S12 in Additional file 2, with the current legend in Additional file 4. `Data_Provenance!B36` contains `fresh_analysis/figures/Figure_4_caption.md`; retain this source filename and its associated provenance rather than renaming it as if it had originally been S12.

The **current main Figure 4** is the new endpoint-resolved biochemical figure. Its complete numerical records are in Additional file 5, with extended methods/results in Additional file 4. It is not the liver-localization figure indexed as source Figure 4 in Additional file 3.

## New compound-screening tables and figures

The new screening outputs have their own index in Additional file 8 and raw inputs/outputs in Additional file 9. They do not overwrite the historical source workbook.

| Supplementary table/sheet | Contents and source file | Final figure use |
|---|---|---|
| N20_Library | Compound identity and clinical context; Candidate_Library_Evidence.csv | Figure 7 and screening Methods/context |
| N21_Boltz | All 32 AI-screening dispositions; Boltz_All32_Dispositions.csv | Figure 7B,C; pocket structures for 7E,F are in Additional file 9 |
| N22_Docking_Runs | All 96 primary Vina searches; Docking_All_Completed_Runs.csv | Figure S13A,B |
| N23_Docking_Summary | All 32 complete compound–target summaries; Docking_Summary.csv | Figure S13C,D and docking Results |
| N24_Properties | All 48 Inductive predictions; Inductive_Properties_All16.csv | Figure S14A–C |
| N25_Assay_Disposition | All 191 experimental activity records; All_Target_Activity_Disposition.csv | Eligibility context for Figure S13E and experimental summaries |
| N26_Activity_Summary | Endpoint-specific eligible biochemical summaries; Endpoint_Separated_Experimental_Summary.csv | Figure S13E |
| N27_Matched_Kd | Davis source-matched binding records with censoring and duplicate flags; Davis2011_Matched_Kd_Benchmark.csv | Figure 7A |

Figure S13F's 30 retained primary PDGFRA nintedanib poses are in `Nintedanib_Pose_QC_All_Searches.csv` within Additional file 9; that complete file also includes 30 FLT4 poses. Figure S14D,E uses the seed-104729 docking structures in the same archive. Figure S14F's six top poses are in `PDGFRA_Nintedanib_Receptor_Sensitivity.csv`; all 53 retained poses across the primary and apo PDGFRA analyses are in `PDGFRA_Nintedanib_Sensitivity_All_Poses.csv`. The three additional apo searches are separate from the 96 primary rows in N22.

A reader-facing workbook update can append this map or add a final-figure/panel column to the index. Preserve the original labels alongside the new route, and retain original filenames, numerical data, source IDs and hashes unchanged. This review made no edits to the source workbook.
