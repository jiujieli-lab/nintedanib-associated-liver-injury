# Independent screening/docking manuscript audit

Reviewed `manuscript/JCM_Manuscript_Text.txt`, `manuscript/new_sections.json`, the original 42-page render and the final 40-page PDF/PNG render. No main-manuscript or source-workbook edits were made by this review.

## Numerical findings

The reported 16 compounds, 32 compound–target combinations, 96 primary docking searches, three fixed seeds, 191 experimental records, 12 scored Boltz compounds per target and four filtered compounds per target agree with the terminal outputs. All named Boltz ranking scores agree at the reported precision. The filtered compounds are nintedanib, sunitinib, ponatinib and MAZ51 on both targets.

The matched experimental comparison agrees with the retained Davis records: PDGFRA/FLT4 Kd 16/95 nM for nintedanib, 31/>10,000 nM for imatinib, and 0.51/170 nM for axitinib. The manuscript appropriately retains censoring and source-duplicate qualifiers.

The 48 Inductive predictions, model versions, logD range 1.376–4.604 and nintedanib logD 2.523 agree with the actual responses. Eight pKa results carry out-of-domain flags (four acidic, four basic); three acidic-pKa results carry low-confidence flags. No logD result carries either flag. The nintedanib physicochemical predictions and their absence of applicability flags also agree with the source.

The separate preflight imatinib redocking RMSD is 0.57427456 Å, correctly reported as 0.574 Å. The three primary-library control RMSDs are 0.55649154, 0.77616255 and 0.71178899 Å. The WT receptor, AlphaFold database version, 222-pair core alignment, 1.074 Å alignment RMSD and pocket mean pLDDT 85.88 agree with the preparation records.

The newly added nintedanib paragraph is numerically correct: primary medians −8.520 and −10.081 kcal/mol; maximum top-pose RMSDs 9.351 and 1.788 Å; apo-PDGFRA median −6.495 kcal/mol and maximum RMSD 5.094 Å. The receptor-dependent interpretation and lack of calibrated cross-target affinity are appropriate.

The displayed Boltz pocket residue counts 199/362 and 195/329 match the exported model provenance. Imatinib–PDGFRA Cys677 N proximity is 2.930026810 Å; tivozanib–FLT4 Cys930 N proximity is 2.974576227 Å. The initial text's 2.98 Å for the latter was a double-rounding error; the final PDF now correctly reports **2.930 and 2.975 Å**.

The initial Methods overview's Supplementary Figures S1–S12 reference was outdated. The final PDF uses **S1–S14** and retains the new screening/sensitivity sources.

## Render review

The original pages 22–42 and corresponding final pages 22–40 show no clipped panels, cut tables or detached figure captions. In the final render, Figures 5, 6 and 7 remain with complete captions on pages 23, 25 and 28. The corrected screening, docking and apo-sensitivity results fit cleanly on page 29. References remain within the page margins. There is substantial whitespace on final page 22 before the following figure; it does not hide or omit content and is not a blocking layout defect.

## Retained source-index routing

The Additional file 3 workbook contains historical figure/panel identifiers. The exact routing guide and 89-row cell-level crosswalk are `Supplementary_Source_Index_Guide.md` and `Additional3_Final_Figure_Index_Map.csv`. Source Figure 4 maps to final Figure S12, and workflow insertion shifts the quantitative panels by one. Original filenames, source identifiers and hashes should remain historical; the current reader route belongs in a companion map or separate final-label columns.

No further quantitative correction was identified in the audited screening/docking claims after the two corrections above were applied. This audit does not claim to revalidate every earlier pharmacovigilance, network or transcriptomic calculation.
