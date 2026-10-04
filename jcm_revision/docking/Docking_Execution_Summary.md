# Completed docking execution: PDGFRA and FLT4

The approved 16-compound library was actually docked with AutoDock Vina 1.2.7 against the fixed primary PDGFRA and FLT4 receptors using three fixed stochastic seeds. **All 96 primary searches completed successfully**, with all 32 compound–target groups containing three searches. The separately approved nintedanib-only apo-PDGFRA sensitivity completed all three searches. A separate imatinib preflight control also completed: **100 actual searches and 988 retained poses overall** (955 primary, 23 apo sensitivity, 10 preflight). No synthetic, missing or substituted scores appear in the tables.

## Nintedanib findings

| Analysis | Receptor | Median score (kcal/mol) | Range (kcal/mol) | Maximum top-pose pairwise RMSD (Å) | Primary within-target median-score position |
|---|---|---:|---|---:|---|
| Primary PDGFRA | WT 6JOL, imatinib-bound | −8.520 | −8.613 to −8.502 | 9.351 | 14 of 16 |
| Primary FLT4 | Predicted AF-P35916-F1 v6, residues 845–1173 | −10.081 | −10.117 to −9.781 | 1.788 | 3 of 16 |
| Separate PDGFRA sensitivity | WT apo 8PQJ | −6.495 | −6.884 to −6.481 | 5.094 | Not ranked; nintedanib only |

The PDGFRA primary score range is narrow while its top poses vary substantially. For seed 104729, the top pose is 12.29 Å from hinge Cys677. A near-degenerate alternative (rank 3; 0.161 kcal/mol above the best score) approaches **gatekeeper Thr674**, not the canonical hinge: Thr674 is 3.38 Å away but Cys677 remains 6.63 Å away. These observations must not be depicted as verified hinge hydrogen bonding or a unique canonical ATP-site pose.

In the apo sensitivity searches, top-pose Cys677 distances become 4.729, 3.578 and 4.311 Å, but pose variability persists and scores are less favorable. This supplies a receptor-dependent contact hypothesis, not a validated pose. Both receptor datasets are retained in full; the follow-up does not replace the primary receptor or its ranking.

The FLT4 receptor is predicted and has no target-matched experimental self-redocking control. Its docking box was transferred geometrically from the experimental **VEGFR2/KDR**, not FLT4, nintedanib complex 3C7Q. Mean pocket pLDDT is 85.88, with several glycine-loop residues below 70. Scores across these different receptors are not calibrated measures of biochemical selectivity. None of these calculations establishes efficacy, target activation/inhibition in vivo, a favorable DILI intervention, or measured binding affinity.

## Validation and execution

The three imatinib searches within the primary 6JOL campaign reproduce the crystallographic ligand with symmetry-aware heavy-atom RMSDs of **0.556, 0.776 and 0.712 Å** in the unchanged receptor frame (no ligand superposition). The separate preflight control gives 0.574 Å and −12.407 kcal/mol. This validates imatinib redocking in that receptor configuration; it does not validate every compound or the FLT4 model.

All prepared ligand conformers converged under MMFF. All 988 exported SDF poses parse, have finite scores and match the raw Vina result arrays. Run-level configurations, input ligand conformers, receptor coordinates, poses, logs and engine outputs are retained. `Final_Docking_Execution_QC.json` and `Final_Docking_Run_Inventory.csv` are the final completeness audit. The primary library SHA-256 is `8fd844b0a6abf4033ee0eaf5ce190a5b4c293de890b6b78044d10f09fac4812e`.

The main protocol uses ETKDGv3, MMFF, Meeko 0.8.0 and Vina 1.2.7; RDKit is 2026.03.6. Seeds are 104729, 130363 and 155921; exhaustiveness is 32, maximum retained poses 10, minimum pose separation 1 Å, energy window 5 kcal/mol, two CPU threads per search. The exact supplied parent states were used without a pH-dependent microspecies enumeration. These are repeated stochastic docking searches, not biological replicates or MD trajectories.

The 8PQJ plan was fixed before its own scores were generated, after the primary nintedanib pose limitation was observed. Deposited WT residue identities were verified. Missing side-chain atoms were reconstructed with PDBFixer, and two newly rebuilt side chains underwent local clash correction; every deposited atom retained its original coordinates. No missing loop was built and no receptor MD was performed as part of docking preparation. Full details and repair audits are in the Methods and sensitivity-plan files.

## Files for manuscript and figures

- `Docking_Methods_and_Limits.md`: full actual protocol, receptor provenance, interpretation limits and completed sensitivity result.
- `Verified_Docking_References.md`: verified method/software/structure references, including the direct 6JOL dataset citation because no associated published paper was verified.
- `Docking_All_Completed_Runs.csv`: 96 primary search-level scores.
- `Docking_Summary.csv`: 32 complete compound–target summaries with medians, ranges and pose variability.
- `Nintedanib_Pose_QC_All_Searches.csv`: all 60 retained nintedanib primary poses with separate gatekeeper and hinge measurements.
- `Nintedanib_Docking_Contacts.csv`: <4 Å protein–ligand heavy-atom proximity contacts, not H-bond classifications.
- `PDGFRA_Nintedanib_Receptor_Sensitivity.csv`: all six PDGFRA nintedanib top poses (three primary plus three apo).
- `PDGFRA_Nintedanib_Sensitivity_All_Poses.csv`: all 53 retained PDGFRA nintedanib poses across both receptors.
- `PDGFRA_8PQJ_Sensitivity_Summary.json`: completed apo versus primary summaries.
- `runs/`: exact raw inputs, configurations and outputs for every executed search.

`PDGFRA_nintedanib_best_docked_complex.pdb`, `FLT4_nintedanib_best_docked_complex.pdb` and `PDGFRA_8PQJ_nintedanib_best_docked_complex.pdb` use the prespecified seed 104729 top pose. They are representative docking complexes, not experimentally determined complexes or MD-ready systems.

## Proto and MD boundaries

Proto successfully retrieved the actual AlphaFold structures. A fresh account check after docking still reported no linked Modal workspace and zero deployed tools. Its deployable Vina/Boltz tools were not run; no paid deployment or duplicate paid prediction was initiated. Exact tool schemas, examples, retrieval outputs and the latest workspace response are in `../proto/`.

This docking campaign did not produce a ligand–protein MD trajectory. A separate team's ligand-free FLT4 OpenMM feasibility benchmark is documented independently under `../md/`; it must not be represented as nintedanib-complex MD or evidence of stable binding.

To regenerate the completed tables, run `env/bin/python summarize_docking.py`, `env/bin/python summarize_8pqj_sensitivity.py` and `env/bin/python finalize_docking_audit.py` from this folder. Archive the scripts and scientifically relevant inputs/outputs; exclude the locally installed `env/` directory.
