FOCUSED COMPOUND SCREENING: SOURCE DATA AND CODE

This archive contains completed computational outputs for 16 compounds against PDGFRA and FLT4, with all scored, filtered and retained results. There are 96 primary Vina searches, three separately specified apo-PDGFRA sensitivity searches and one separate crystallographic redocking control. The three seeds are computational searches, not biological replicates. No target–ligand production molecular-dynamics trajectories are included or claimed.

Directory layout
source_indices: complete89-row mapping of analysis-source figure IDs to final manuscript panels.
jcm_revision/candidate_library: exact compound identities, label sources, all191 experimental records, exclusions, endpoint summaries and the source-matched Kd panel.
jcm_revision/boltz-experiments: both completed provider jobs, all24 exported pocket models and raw metrics. Filtered candidates have no returned model score. Exported coordinates are cropped pocket structures.
jcm_revision/screening: complete32-input disposition, 48 molecular-property predictions, exact provider responses and summarization scripts.
jcm_revision/docking: receptor structures, preparation, grids, molecular inputs, all search logs/poses and terminal validation. Python environment binaries are omitted; resolved software versions are preserved in methods and manifests.
jcm_revision/proto: actual AlphaFold retrieval and tool-execution provenance. No unexecuted Proto prediction is represented as a result.
jcm_revision/md: Boltz coordinate/contact quality checks only; the directory name does not indicate an MD simulation.
jcm_revision/figures/screening: plotting scripts and underlying numerical/coordinate records for Figure7 and FiguresS13–S14.
binding_analysis/Inductive_Bio_Raw_Predictions.json: the matched-version raw nintedanib predictions used with the15 further compounds.

Reproduction
Run `python verify_archive.py` from this directory for an offline integrity and count check (Python standard library only). All SHA-256 entries must match. Analysis scripts use paths relative to their own locations unless an archived runtime path is explicitly retained as provenance. Summarization requires RDKit2026.03.6 and the standard numerical stack. Docking used AutoDockVina1.2.7, Meeko0.8.0 and RDKit2026.03.6. See docking/Docking_Methods_and_Limits.md and Additionalfile7 for exact preparation, seeds, grid and exhaustiveness. Rerunning provider inference requires the respective accounts and may incur new charges; no such rerun is necessary to inspect these archived results. Do not execute the campaign launcher merely to review results.

Interpretation
Binding confidence is a target-specific model output, not a measured dissociation constant. The returned_smiles_identity_verified field checks returned molecular metadata against input stereospecific InChIKeys; it is not an independent check of coordinate chirality. Contact distances do not prove hydrogen bonds. Rigid docking scores do not measure affinity, clinical benefit or hepatic safety. Source record inequalities and duplicate flags are retained.
