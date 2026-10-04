# Four workflow replacements

Completed files are in `final/`, named `Current_Figure_1`, `Current_Figure_5`, `Current_Figure_6`, and `Current_Figure_7`, each as TIFF, PNG, and PDF. Numbering remains provisional for manuscript assembly.

The replacements use a 380-pixel-high workflow band. Quantitative content below the band retains its native pixels, dimensions, and horizontal position. Its only geometric change is a uniform vertical translation of 182, 164, 164, and 166 pixels, respectively. All outputs retain 300 dpi; the quantitative panels have not been upsampled. PDFs embed the same native raster content at its 300-dpi physical size.

Forty-four lower panels were validated against their original source crops, with the source extraction's documented neighboring-text exclusions preserved. The complete lower canvas of each figure is also pixel-identical to the corresponding original reassembled canvas after translation. No data chart or measurement was generated or redrawn.

TIFFs were encoded in memory, atomically written, and reopened after the writer process exited. All four have nonzero TIFF IFD offsets, decode successfully, and match their PNG exports pixel for pixel. PDF pages were rendered with Poppler and visually inspected. The workflow bands were inspected at native size; no text clipping, icon clipping, or overlaps remain.

Generated biomedical illustrations are schematic, not molecular structures, microscopy, or experimental results. The built-in image-generation tool created one coordinated icon sheet. The exact prompt is saved in `assets/imagegen_prompt.txt`. All labels, arrows, panel letters, layouts, and composite geometry were produced deterministically by `build_workflow_figures.py`.

QA records: `composite_validation.json`, `export_reopen_validation.json`. Review previews: `final/Workflow_Panels_Overview.png` and `pdf_qa/PDF_Contact_Sheet.png`.

Scientific safeguards: the 29 nonzero target runs are a subset, not the remainder of all 210 runs after numerical-zero exclusion; prioritization is not causal DILI validation; the independent liver study is a pharmacodynamic comparison, not an observed rescue experiment. Figure 2A was left untouched because it is a quantitative plot.
