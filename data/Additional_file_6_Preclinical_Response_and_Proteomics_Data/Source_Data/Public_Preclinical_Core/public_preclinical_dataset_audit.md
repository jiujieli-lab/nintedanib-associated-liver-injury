# Public preclinical nintedanib dataset audit

**Freeze date:** 4 September 2026  
**Scope:** publicly accessible cell, tissue-slice, organoid, animal, proteomic, phosphoproteomic, transcriptomic, metabolomic, high-throughput toxicology, and morphology experiments containing an exact nintedanib/BIBF 1120 exposure or a directly relevant repository inventory search.

## Search and eligibility rule

GEO/SRA/BioProject, DDBJ/ENA, ArrayExpress/BioStudies, PRIDE/ProteomeXchange, MetaboLights, Open TG-GATEs, DrugMatrix/CEBS, ToxCast/Tox21 invitroDB, LINCS, and Broad Cell Painting/JUMP were searched with Boolean combinations of `nintedanib`, `BIBF 1120`, `lung`, `pulmonary fibrosis`, `bleomycin`, `liver`, `hepatic`, `hepatocyte`, `toxicity`, `proteome`, `phosphoproteome`, `transcriptome`, `single cell`, and `metabolome`. Eligibility required an exact exposure, a concurrent control, recoverable treatment timing, downloadable sample-level measurements, and an experimental unit that could be reconstructed without treating cells, slices, technical channels, or sequencing runs as independent animals.

The machine-readable record is `public_preclinical_dataset_audit.csv`; the broader search ledger, including negative inventories, is `dataset_screening_ledger.csv`.

## Freeze decisions

| Resource | Experimental unit | Decision | Permitted interpretation |
|---|---|---|---|
| GSE151374 / PRJNA635636 / SRP265115 | Independently prepared pooled 10x library; 3 libraries/condition/time, 3 mice pooled/library | Reanalyzed in Figure S6 | Pulmonary macrophage/composition sensitivity only. Cells and constituent mice never contribute degrees of freedom; not hepatic-toxicity evidence. |
| PXD052594 / Zenodo 11395642 | Individual mouse for the processed proteome: 10 nintedanib, 13 vehicle; separate 5-versus-5 phosphoproteome subset | Reanalyzed in Figure S7 and source tables | Unconditioned pulmonary proteome sensitivity only. Radiomic response clusters were not used. Phosphosites are a table-only unmoderated Welch–BH screen and were not assessed by an exact small-sample test or adjusted for matched total-protein abundance/site occupancy. |
| DRA012991 / PRJDB12477 | Individual mouse after aggregating 2–4 runs/animal; 5 nintedanib, 3 vehicle | Audit-only | Raw reconstruction is possible, but the deposited 26 runs are not 26 independent samples and the source post hoc animal exclusions cannot be inherited. |
| GSE278200 and GSE308578 | Individual rat/mouse | Included in Figure 6 | Animal-level pulmonary-response/model-boundary evidence, never liver-safety evidence. |
| GSE120679 and GSE120804 | Slice observation from tissue pooled across three source rats | Exploratory Figure 6 context | Within-pool slice reproducibility only, not animal-population inference. GSE120804 is cholestatic hepatic context, not nintedanib-induced liver injury. |
| PXD024058 | Three quantitative columns/group with conflicting biological-versus-technical provenance | Descriptive Figure 6 context | Protein direction only; no inferential p/q values or independent-animal claim. |

No eligible nintedanib-monotherapy animal-liver or primary-hepatocyte DILI omics experiment was identified. The closest public hepatic resources were pooled BDL liver slices (GSE120804) and transformed HepG2 L1000 signatures (GSE70138); neither constitutes an animal or primary-hepatocyte model of nintedanib-associated DILI.

## GSE151374 evidence boundary

The 18 deposited libraries span NaCl, BLM, and BLM+nintedanib at days 7 and 14. Each condition/time cell contains three independently prepared libraries, and each library pools three mice (54 mice represented). Nintedanib was given at 50 mg/kg on days 0–6 for the day-7 endpoint or days 7–13 for the day-14 endpoint. The GEO record does not provide an unambiguous BLM dose/route, so those details are not imputed. Figure S6 uses library-level QC summaries, broad marker-rule lineages, pooled-library composition, and macrophage pseudobulk. The exhaustive 3-versus-3 two-sided permutation grid has a minimum attainable P value of 0.10; no genome-wide macrophage differential-expression result or prespecified composition/module/candidate family survived its stated multiplicity control.

The accession-linked primary report is Watson et al., *American Journal of Respiratory Cell and Molecular Biology* 2023;68:366–380 (PMID 36227799; DOI 10.1165/rcmb.2022-0021OC). Zhong et al., *Cellular Signalling* 2025;128:111635 (PMID 39892726; DOI 10.1016/j.cellsig.2025.111635) is retained as secondary mechanistic context and is not used to define the accession's sample mapping.

## PXD052594 evidence boundary

The source study used female C57BL/6J mice, BLM 2 U/kg intratracheally, and double-blind randomized nintedanib 60 mg/kg orally once daily on days 7–20, with tissue collected on day 21. Six postrandomization exclusions preceded the final imaging cohort (10 nintedanib and 14 vehicle), and one additional vehicle animal (M29) was excluded upstream from comparative-proteome sample preparation. Figure S7 therefore compares all animals present in the processed matrix (10 nintedanib versus 13 vehicle) but is not presented as an unbiased causal effect. No animal was selected, excluded, stratified, or adjusted on the basis of the post-treatment radiomic clusters.

All 7,006 proteins were null after protein-wide BH correction (minimum q=0.607966); the five prespecified modules were null after permutation and within-family BH control (minimum q=0.480045); and the detected candidate-track family was null. The separate deposited phosphoproteome contains 20,043 complete rows from a random 5-versus-5 animal subset, of which 994 met site-wise BH q<0.05 in an unmoderated Welch–BH screen. Because only processed intensities were available, this was not an exact small-sample test; no matched total-protein adjustment or site-occupancy model was applied, so those site shifts cannot be separated from underlying protein-abundance shifts and do not establish target engagement or a signaling-state mechanism. The primary report is Lauer et al., *JCI Insight* 2024;9:e181757 (PMID 39012714; DOI 10.1172/jci.insight.181757).

## DRA012991 animal/run map and exclusion audit

Mikami et al. studied D1CC×D1BC induced rheumatoid-arthritis-associated ILD mice treated with nintedanib 90 mg/kg orally once daily during weeks 35–43 (PLoS One 2022;17:e0270056; PMID 35714115; DOI 10.1371/journal.pone.0270056). The raw deposit contains 26 MiSeq runs:

| Animal/source label | Run accessions |
|---|---|
| Vehicle VA | DRR326942, DRR326943, DRR326944, DRR326945, DRR326946, DRR326947 |
| Vehicle VB | DRR326948, DRR326949 |
| Vehicle VC | DRR326950, DRR326951 |
| Control run | DRR326952 |
| Nintedanib NA | DRR326953, DRR326954 |
| Nintedanib NB | DRR326955, DRR326956 |
| Nintedanib NC | DRR326957, DRR326958, DRR326959 |
| Nintedanib ND | DRR326960, DRR326961, DRR326962 |
| Nintedanib NE | DRR326963, DRR326964, DRR326965, DRR326966, DRR326941 |

The publication excluded ND, NE, and VC after inspecting the clustering result. A non-circular reanalysis would instead aggregate runs within each named animal, retain all five treated and all three vehicle animals, adjudicate the separate `Control` run before analysis, and never use the 26 runs as the sample size. Because only raw FASTQ files were deposited and this reconstruction was outside the frozen extension, DRA012991 remains audit-only.

- DDBJ submission: https://ddbj.nig.ac.jp/resource/sra-submission/DRA012991
- ENA BioProject mirror: https://www.ebi.ac.uk/ena/browser/view/PRJDB12477

## Negative and context-mismatched inventories

- Open TG-GATEs: no exact nintedanib/BIBF 1120 exposure in the official 170-compound inventory.
- DrugMatrix/CEBS: no exact exposure record was identified.
- ToxCast/Tox21 invitroDB v4.3: a directly measured nintedanib record could not be confirmed under the accessible interface; this is explicitly not an assertion of absolute absence.
- Broad cpg0004-LINCS: BRD-K49075727 exists in A549 Cell Painting (0.041–10 μM, 48 h; five technical wells/dose), but the morphology wells are neither a fibrotic-lung tissue model nor a hepatocyte/DILI model.
- JUMP cpg0016: no exact free/salt InChIKey match in the audited compound/source metadata.
- MetaboLights: no exact-nintedanib cell/animal metabolomics accession in the prespecified lung/liver contexts.
- PXD044051 and GSE218575 contain exact nintedanib exposures but were excluded because gastric-cancer cells and bladder-cancer xenografts do not address the prespecified lung-fibrosis or hepatic-injury questions.

## Reproducibility pointers

- Figure S6 analysis: `gse151374/code/analyze_gse151374.py`; manifest `gse151374/results/analysis_manifest.json`; detached hash audit `gse151374/validation/manifest_hash_validation.json`.
- Figure S7 analysis: `pxd052594/code/analyze_pxd052594.py`; manifest `pxd052594/results/analysis_manifest.json`; detached hash audit `pxd052594/validation/manifest_hash_validation.json`.
- Captions: `gse151374/docs/Supplementary_Figure_S6_caption.md` and `pxd052594/docs/Supplementary_Figure_S7_caption.md`.
