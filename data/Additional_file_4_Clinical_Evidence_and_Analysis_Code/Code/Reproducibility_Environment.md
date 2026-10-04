# Reproducibility environments

The analyses were frozen on 4 September 2026. Two Python environments are listed because the virtual-knockout stack and the journal-artifact/standard statistical stack have different compatible dependency sets.

## Standard analysis and artifact environment

Use Python 3.12 and `requirements-analysis.txt` for FAERS, LINCS, pharmacology, animal/cell reanalysis, integration figures, manuscript assembly, DOCX generation, and PDF generation. Individual frozen manifests remain authoritative for the exact software versions used by each analysis module.

## Virtual-knockout and notebook environment

Use Python 3.12 and `requirements-virtual-knockout.txt` for scTenifold virtual knockout and for the human DILI serum-proteomics reanalysis, whose manifest records the same NumPy, pandas, SciPy, scikit-learn, and statsmodels versions plus `openpyxl 3.1.5`. The installed distribution was `sctenifoldpy 0.4.0`; the imported `scTenifold` module reports internal version `0.3.0`. This package/module version distinction is recorded here rather than silently normalized.

The portable reproducibility notebook is a frozen-output audit rather than an upstream scTenifold or proteomics rebuild. It may be executed in the standard analysis environment; its exact arithmetic is also exercised by the deterministic fallback runner. The two requirement files intentionally pin incompatible numerical stacks and must be installed into separate environments.

## External data access

Raw third-party archives are not redistributed in the compact submission package. Accession numbers, persistent URLs, retrieval dates, source filenames, and available SHA-256 hashes are provided in `Data_and_Code_Availability.md`, the supplementary workbook, and the layer-specific manifests. The package includes processed source data sufficient to verify all reported figures and tables.

## Notebook execution note

Jupyter and `nbformat` were not installed in the final release-build runtime. The checked notebook was therefore generated with a standard-library nbformat-v4 writer and executed with the supplied deterministic fallback runner, which executes code cells sequentially in a shared namespace, captures outputs, and fails on any exception or stderr. The notebook remains a valid nbformat-v4 artifact and can be opened and rerun in a conventional Jupyter environment. The final validation report records zero errors and continuous execution counts.
