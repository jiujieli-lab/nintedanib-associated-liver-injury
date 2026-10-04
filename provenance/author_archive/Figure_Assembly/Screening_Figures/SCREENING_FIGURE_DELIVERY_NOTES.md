# Final screening figure inventory

All three figure sets are complete and use final source data. Each set contains an RGB TIFF, a PNG and a PDF under `final/`. New quantitative plots are 600 dpi; coordinate-based molecular panels are embedded at 600 ppi. TIFFs were encoded in memory, written atomically, reopened after the writer process exited, and verified pixel-identical to the corresponding PNG. Rendered PDFs were visually checked.

| Figure | File stem | Pixels | Evidence |
| --- | --- | --- | --- |
| Main 7 | Figure_7_Compound_Comparison | 4320 × 4050 | 12 matched-assay Kd records; all 32 BoltzMol dispositions (24 scored, 8 filtered); experimental redocking control; two selected verified predicted pocket models |
| S13 | Figure_S13_Docking_and_Evidence | 4320 × 6090 | All 96 primary searches and all 32 compound–target summaries; endpoint evidence coverage; all 30 retained PDGFRA nintedanib poses |
| S14 | Figure_S14_Properties_and_Receptor_Sensitivity | 4320 × 4410 | All 48 property predictions; selected exact-coordinate nintedanib poses; six receptor-sensitivity searches (three primary plus three additional apo searches), with all 53 retained poses in source data |

Final captions are `Figure_7_Legend.txt`, `Figure_S13_Legend.txt` and `Figure_S14_Legend.txt`. Reproducible scripts are `build_main_screening.py`, `build_supplement_docking.py`, `build_supplement_properties.py` and `plot_fixed_evidence.py`; molecular rendering scripts and provenance are in `molecular/`. Per-figure source manifests and export QA JSON files document inputs and final checks. PDF previews are in `pdf_qa/`.

Interpretation constraints are retained in labels/captions: provider binding confidence is dimensionless and uncalibrated; filtered compounds are not assigned zero scores; exported BoltzMol models are cropped pockets; docking scores are assessed within each receptor; seeds are computational searches, not biological replicates; model bounds are not confidence intervals; structural proximity is not a hydrogen-bond assignment. No molecular dynamics results or safety claims are introduced.

The earlier four workflow composites remain separately named `../final/Current_Figure_{1,5,6,7}` pending root-controlled narrative numbering. Their lower original quantitative panels remain pixel-identical at the original native resolution.
