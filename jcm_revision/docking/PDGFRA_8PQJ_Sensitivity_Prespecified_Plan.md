# Separately reported receptor-conformation sensitivity analysis

This plan is recorded before any docking score is generated for 8PQJ.

Rationale: the primary imatinib-bound PDGFRA structure 6JOL is an inactive conformation. Independent examination of its nintedanib outputs found no canonical hinge-Cys677 proximity in the leading pose. An alternative wild-type conformation is therefore a targeted robustness check. It will not replace the primary receptor or be selected on the basis of a more favorable score.

Use the deposited wild-type apo PDGFRA coordinates in **8PQJ**, chain A, X-ray resolution **1.82 Å**, released 27 December 2023. RCSB reports no point mutation; compare every coordinate-bearing amino acid with P16234 and explicitly check Thr674 and Asp842. The construct spans UniProt 550–973 with kinase-insert deletion and unresolved loops; preserve the actual coordinates and disclose missing regions. This is an experimental apo kinase construct, not a full-length native receptor.

Transfer the crystallographic imatinib coordinates from the primary PDGFRA 6JOL receptor to 8PQJ by a sequence-matched Cα structural superposition. Use an iterative conserved-core fit, retain the transform and matched-residue map, and report core size and RMSD. The transferred ligand defines the same style of search box: 6 Å on each face, minimum dimension 22 Å. It is not an inferred binding pose of imatinib in 8PQJ.

Dock only the exact approved parent nintedanib SMILES, with the unchanged Vina 1.2.7/Meeko 0.8.0/RDKit 2026.03.6 protocol, exhaustiveness 32, seeds **104729, 130363, 155921**, 10 poses, 1 Å minimum pose separation and 5 kcal/mol pose window. Retain all three scores, all poses, ligand preparation logs and receptor files regardless of outcome. Report the 8PQJ results alongside, and separately from, the primary 96-search panel.

Prespecified outputs: each top score; median and range over the three searches; top-pose consistency; distance to gatekeeper Thr674 and hinge Cys677 for every retained pose; direct comparison with the primary nintedanib search results. These are descriptive docking-sensitivity metrics and do not establish measured affinity or binding mode. No score threshold is used to discard a completed run. Do not update the primary ranking based on this sensitivity experiment.

Associated structure paper: Teuber A. et al. *Avapritinib-based SAR studies unveil a binding pocket in KIT and PDGFRA*. Nature Communications **2024**, **15**, 63. DOI https://doi.org/10.1038/s41467-023-44376-8 . The paper includes multiple apo and ligand-bound, wild-type and mutant structures; **8PQJ itself is wild-type apo**, not an avapritinib complex.

## Preparation record completed before scoring

The WT identity check passed for all 329 coordinate-bearing residues. The 6JOL-to-8PQJ alignment uses 279 common Cα atoms, of which 271 form the retained core (RMSD 0.785 Å). The transferred reference defines center (26.94732, 4.68286, 21.49816) Å and size (25.22713, 27.43338, 22.0) Å. The small core RMSD means this is an apo-versus-ligand-bound structural sensitivity check, not evidence that the structures span every biologically relevant kinase activation state.

PDBFixer 1.12.0/OpenMM 8.4.0 restored 32 absent side-chain atoms in 11 residues and one terminal OXT. No absent loop was reconstructed and no deposited atom moved (verified maximum displacement 0 Å). Rebuilt Asn780 and Asp785 side chains initially clashed; only their newly added atoms were rotated around Cα–Cβ bonds on a 15° grid to eliminate short nonbonded contacts. These residues were more than 17 Å from the transferred reference ligand before adjustment. The simple geometric repair is not an affinity calculation, force-field validation or MD. Its choices and atom-level audit are retained in `8PQJ_heavy_atom_repair_QC.json` and `8PQJ_rebuilt_sidechain_clash_repair.json`. Meeko then parameterized the entire receptor successfully with bad-residue deletion disabled.
