# Actual local docking protocol

This folder contains an executed PDGFRA redocking validation and the completed, prespecified comparative docking campaign: 96 successful searches across 16 compounds, two targets and three fixed seeds. The 32 compound–target summaries each contain all three searches. All ligand MMFF minimizations converged; 955 retained poses parse successfully and agree with their raw result scores. `Final_Primary_Docking_QC.json` records these checks. The separate nintedanib-only apo-receptor sensitivity experiment is described below and does not replace any primary result.

## Receptors and binding sites

**PDGFRA:** PDB 6JOL, human wild-type kinase in complex with imatinib, X-ray resolution 1.90 Å. All 279 coordinate-bearing residues were checked against the canonical P16234 sequence: no mismatches; T674 and D842 retain their wild-type identities. The experimental construct deletes part of the kinase insert and has unresolved loops. Chain A protein atoms were retained, water and ligand removed, and Meeko 0.8.0 templates assigned receptor chemistry and polar hydrogen positions. Alternate-location preference is A. No bad residue deletion was allowed. The imatinib-bound conformation is an inactive-state receptor; its ranking can favor ligands compatible with this conformation. This receptor and box were fixed for the full primary library. A later nintedanib-only apo-receptor sensitivity analysis was motivated by receptor-conformation bias and the observed ambiguity of the primary nintedanib poses. Its protocol was recorded before its own docking scores were generated, and its results remain separate from the primary campaign.

The PDGFRA box is centered at (-38.4125, 157.049, 0.794) Å with dimensions (22.0, 28.722, 25.84) Å. It is the deposited imatinib heavy-atom bounding box plus 6 Å on each face, with a minimum dimension of 22 Å.

**FLT4:** no experimentally determined kinase coordinates were identified by an exact current RCSB search for P35916. Entries 4BSJ/4BSK and newer 21EI/21ES cover extracellular regions. The receptor is the canonical AlphaFold DB model AF-P35916-F1 version 6, cropped to UniProt residues 845–1173. This is a predicted rigid receptor. The ATP-site box was located by aligning the experimental VEGFR2–nintedanib complex 3C7Q to the FLT4 kinase, using BLOSUM62 sequence alignment and an iterative Cα core superposition. The alignment retains 222 core Cα pairs with 1.074 Å RMSD (269 matched pairs before core selection). The transferred nintedanib coordinates define the box only; they are not a FLT4 docking result. The pocket residues within 6 Å have mean pLDDT 85.88, with several glycine-loop residues below 70. No FLT4 self-redocking control is available.

The FLT4 box is centered at (-11.70297, 3.17197, -11.88550) Å with dimensions (22.0, 22.75049, 29.37968) Å. It uses the same 6 Å-per-face rule and minimum size. The 3C7Q template contains modified/phosphorylated protein residues, but none were copied into FLT4: only the geometrical ligand box was transferred.

## Ligands, scoring, and repeated searches

The exact 16-compound approved library is read from `../candidate_library/Candidate_Screening_Library.csv`. It includes the index drug, mechanistic kinase comparators, and two clinical nonkinase references. The latter are not validated biochemical nonbinders. Supplied parent SMILES, tautomer/protonation states and stereochemistry are used unchanged apart from RDKit canonicalization and hydrogen addition. This is not a pH-dependent microspecies enumeration. Protonation-state and receptor-state dependence limit interpretation.

Each search independently embeds a ligand using ETKDGv3 and minimizes it with MMFF (up to 2,000 iterations); its convergence status is recorded. Meeko creates ligand and receptor PDBQT. AutoDock Vina 1.2.7 runs the Vina scoring function, exhaustiveness 32, maximum 10 retained poses, minimum pose separation 1 Å, energy range 5 kcal/mol, default grid spacing 0.375 Å. Three fixed search seeds were prespecified: 104729, 130363, 155921. Four processes each use two CPUs, matching the runtime's eight-CPU quota. Actual run-level scores, seeds, input conformers, configurations, poses, logs and package versions are retained.

The three seeds assess stochastic search consistency. They are not independent biological replicates or molecular-dynamics trajectories. Summary medians and ranges are descriptive; no significance test treats these searches as experiments. Scores from different receptor structures are not calibrated measures of biochemical selectivity.

## Completed method validation

Imatinib was removed from 6JOL and independently embedded from the authoritative CCD STI SMILES, then redocked with the stated box and exhaustiveness. The top-ranked pose has a symmetry-aware heavy-atom RMSD of **0.574 Å** from the crystal pose in the unchanged receptor coordinate frame (no ligand superposition), and Vina score **−12.407 kcal/mol**. This establishes a successful PDGFRA control for this ligand and receptor configuration. It does not validate every compound, the FLT4 model, or clinical mechanism. Full control output is in `runs/PDGFRA_Imatinib_redocking_104729/`.

In addition, the three imatinib searches within the actual library campaign independently reproduce the crystal pose at **0.556, 0.776 and 0.712 Å**, respectively. These use the library's canonical atom ordering for initial embedding and the campaign's two-CPU setting; the separate preflight control used the CCD SMILES ordering and four CPUs. The three campaign controls provide the directly matched validation range **0.56–0.78 Å**. All coordinates and the RMSD audit are retained; no ligand alignment was performed before RMSD calculation.

## Output interpretation

`Docking_All_Completed_Runs.csv` contains actual per-search scores. `Docking_Summary.csv` contains only completed compound–target groups, with completed-search counts, score medians/ranges and pairwise top-pose RMSD. `Nintedanib_Docking_Contacts.csv` uses a simple protein–ligand heavy-atom distance below 4 Å; these are proximity contacts, not a hydrogen-bond classification. The seed-104729 best-pose complex PDBs, when created, are docking structures and require additional chemistry/structure preparation before any MD.

Rigid docking omits solvent, induced fit, entropy and much receptor dynamics. A more favorable score is not a measured free energy, Kd, IC50, target activation/inhibition state, lower DILI risk, or therapeutic recommendation. The mechanistic objective is comparative target engagement. Binding or inhibition of a target implicated in a DILI pathway does not itself establish a beneficial intervention.

The seed-104729 nintedanib pose ranking has a concrete binding-mode limitation: its best PDGFRA pose occupies the back-pocket region and is 12.29 Å from hinge Cys677 (minimum heavy-atom distance). Pose 3 has a score only 0.161 kcal/mol higher and approaches gatekeeper Thr674 to 3.38 Å, but remains 6.63 Å from Cys677. Proximity to the gatekeeper must not be mislabeled a canonical hinge contact. Neither score ranking alone nor this alternative pose resolves canonical ATP-site engagement.

Across all three primary nintedanib searches, the median Vina scores are −8.520 kcal/mol for PDGFRA (range −8.613 to −8.502) and −10.081 kcal/mol for FLT4 (range −10.117 to −9.781). Nintedanib ranks 14th and 3rd, respectively, within the corresponding 16-compound receptor-specific score lists. These separate lists are not a cross-target selectivity calculation. The maximum pairwise top-pose heavy-atom RMSD is 9.351 Å for PDGFRA versus 1.788 Å for FLT4, without ligand superposition. Consequently, the narrow PDGFRA score range must not be presented as evidence of a stable or unique binding mode. `Nintedanib_Pose_QC_All_Searches.csv` retains all 60 primary nintedanib poses with separate gatekeeper and hinge proximity measurements.

## Separately specified apo-receptor sensitivity

The sensitivity receptor is wild-type apo PDGFRA 8PQJ, chain A, 1.82 Å resolution. All 329 coordinate-bearing residues match canonical P16234, including Thr674 and Asp842. This is an apo structure; the structure paper's title must not be used to label it as an avapritinib-bound complex. Core alignment to 6JOL gives 0.785 Å RMSD over 271 retained Cα pairs (279 initially matched); therefore this single alternative structure does not sample a full receptor-state ensemble.

The 6JOL imatinib reference was transferred through that alignment only to define the search box. The same bounding-box-plus-6-Å-per-face rule, minimum dimension 22 Å, yields center (26.94732, 4.68286, 21.49816) Å and dimensions (25.22713, 27.43338, 22.0) Å. No docking score was used to choose the box. Nintedanib alone was searched with the same three seeds, ligand preparation, exhaustiveness 32 and two CPUs per search. The 16-compound primary rankings remain unchanged.

Eleven residues in deposited 8PQJ have incomplete side chains. PDBFixer 1.12.0/OpenMM 8.4.0 added 32 missing side-chain heavy atoms and a terminal OXT, with missing-residue reconstruction disabled. Every deposited atom retained its original coordinates (maximum displacement 0 Å). Newly built Asn780 and Asp785 side-chain atoms were subsequently rotated about their CA–CB axes in 15° increments to remove local clashes; only newly added atoms were moved. Those sites are more than 17 Å from the transferred reference ligand. No loop was rebuilt, and this geometry repair is not a molecular-dynamics simulation. Full atom-addition and clash-repair audits, the pre-scoring plan, and both original and prepared receptors are retained.

Sensitivity output is reported in `PDGFRA_Nintedanib_Receptor_Sensitivity.csv`, `PDGFRA_Nintedanib_Sensitivity_All_Poses.csv` and `PDGFRA_8PQJ_Sensitivity_Summary.json`. The plan was specified after observing the primary binding-mode limitation, so this is a transparent follow-up sensitivity analysis, not an independently prespecified second primary receptor.

All three apo-receptor searches completed, retaining 23 poses in total. Their median Vina score is −6.495 kcal/mol (range −6.884 to −6.481), and the maximum pairwise top-pose RMSD is 5.094 Å. The three top poses approach Cys677 to 4.729, 3.578 and 4.311 Å, respectively, while remaining 8.238–9.320 Å from gatekeeper Thr674. These distances show receptor-dependent proximity changes; they do not establish hydrogen bonds or a reproducible canonical binding mode. The six top poses and all 53 retained poses from the two PDGFRA nintedanib receptor analyses are reported without favorable-score selection. The combined execution audit comprises 96 primary searches, three apo sensitivity searches and one separate preflight redocking search (988 retained poses in total).

## Provenance

- RCSB 6JOL: https://www.rcsb.org/structure/6JOL
- RCSB 3C7Q: https://www.rcsb.org/structure/3C7Q
- AlphaFold FLT4: https://alphafold.ebi.ac.uk/entry/P35916
- AutoDock Vina: https://doi.org/10.1021/acs.jcim.1c00203
- Meeko: https://github.com/forlilab/Meeko
- All downloaded coordinate and metadata files, RCSB FLT4 query response, receptor hashes and software versions are retained in this folder and `../proto`.
