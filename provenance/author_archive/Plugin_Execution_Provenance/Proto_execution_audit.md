# Proto execution audit — 2 October 2026

The selected Proto connector is callable. Live `workspace_info` was queried again after the user's reconnection message: it still reported `modal.connected=false`, `tools_deployed=0`, and no HuggingFace token. This is the remote service response, not an assumption based on earlier configuration. The connector itself works: two authorized, free `alphafold-db-fetch` runs completed and their returned AlphaFold structures were downloaded.

## Completed work

- Read the entire 156-tool catalogue, deployed availability, and relevant docking/dynamics/prediction schemas and known-good examples.
- Ran `alphafold-db-fetch` for human PDGFRA P16234 and FLT4 P35916. Both returned model version 6, created 1 August 2025. Raw run outputs, full PDBs, kinase-domain excerpts, and confidence summaries are retained here.
- Inspected all callable tool names/descriptions for general-purpose Python, custom code, GPU, OpenMM, GROMACS and Modal runners. No custom scientific compute runner is exposed. `proto_deploy_tool` accepts a registered `tool_key`, not arbitrary code.
- No deploy or paid model run was submitted. Schema/example discovery calls are not biological model executions.
- Because remote Vina requires deployment, a separate task-local CPU Vina environment was installed. It completed 96 primary library searches, three apo-receptor sensitivity searches and one separate preflight control, with all raw results retained in `../docking`. Consult that folder's final execution summary and audit for actual docking evidence, not this capability audit. A final fresh Proto workspace query still reports Modal disconnected and zero deployed tools; its exact response is retained as `workspace_status_final_recheck_raw.json`.

## Actual available capabilities

| Tool | Engine/output | Device flag | Live status | Required input |
|---|---|---|---|---|
| `vina-docking` | Rigid-receptor AutoDock Vina or Vinardo search; ranked SDF/PDBQT poses, score and RMSD bounds | CPU | `needs_deploy` | Receptor PDB/mmCIF; one or more ligand Fragments/SMILES; explicit box or coordinate-bearing reference ligand |
| `boltz2-prediction` | Joint biomolecular structure prediction; one selected structure per complex | GPU | `needs_deploy` | `complexes: [{chains: [...]}]`; optional per-complex MSAs |
| `boltz2-affinity` | Boltz-2 protein–ligand structure, predicted log10 IC50 in μM and binder probability as described by wrapper | GPU | `needs_deploy` | Same complex input plus optional `binder_chain` selection |
| `bioemu-sample` | Generative protein conformational ensemble sampling | GPU | `needs_deploy` | Protein complex/sequence and optional MSA; does not provide a conventional ligand–protein MD trajectory |
| `alphafold-db-fetch` | Existing public monomer prediction retrieval | No GPU needed | `ready_to_run` | `uniprot_id` and optional isoform |
| `pdb-fetch-entry` | PDB title, method, resolution, source URL | CPU/network | `ready_to_run` | `pdb_id` |
| `pdb-fetch-fasta` | PDB entity/chain sequences | CPU/network | `ready_to_run` | `pdb_id` |
| `sequence-fetch` | Batched name/ID resolution and sequence/structure metadata from NCBI/UniProt/PDB | CPU/network | `ready_to_run` | Requests with `target_name`, `organism`, `sequence_types`; optional `pdb_id` |

No OpenMM or GROMACS production MD tool exists in the current catalogue. OpenMM occurs only in `freebindcraft-design` protein relaxation. PyRosetta FastRelax, BioEmu ensembles, Boltz diffusion samples and docking searches must not be described as time-resolved molecular dynamics.

## Exact input contract and reproducibility

Full JSON Schemas, canonical examples, implementation URLs, method citations and licenses are retained in each `*_metadata.json` file. These are the authoritative detailed contracts. Before any new remote run, call `get_tool_schema(tool_key)`, then `get_tool_example(tool_key)`, then construct exactly one `run_tool` request.

For Vina, `inputs.receptor` is a Structure object with raw `structure` PDB/mmCIF content; `inputs.ligands` is a list of Fragment objects containing `smiles`; `inputs.search_box` can be `{mode: "coordinates", center: [x,y,z], size: [sx,sy,sz]}` or `{mode: "reference_ligand", reference_ligand: {structure: "..."}, padding: 4}`. A remote worker cannot assume it can read this task's local file path; pass actual content or a supported accessible URL.

Vina config fields include `scoring_function` (`vina`/`vinardo`), `exhaustiveness`, `num_poses`, `energy_range`, `min_rmsd`, `max_evaluations`, `cpu`, `grid_spacing`, `allow_bad_residues`, positive integer `seed`, and `timeout`. Default exhaustiveness is 8; default poses 9; default energy range 3 kcal/mol. Keep `allow_bad_residues=false` so missing receptor chemistry cannot silently delete residues. Output includes each ligand's actual seed, conformer optimization status, warnings and pose SDF/PDBQT.

Two fully populated, schema-derived nintedanib–kinase Proto Boltz2-affinity payloads are prepared but explicitly NOT SUBMITTED. They preserve the independently verified kinase sequences and neutral-parent Z-isomer SMILES from the existing official Boltz API inputs. They are alternatives to the official Boltz API route; do not submit duplicate campaigns through both services.

## Remote blocker and concrete next action

Proto's account settings endpoint is https://proto.evodesign.org/settings/account/modal . Linking the Proto plugin is distinct from linking a Modal workspace. Once Modal is truly linked, re-query availability. If Vina remains `needs_deploy`, deploying `vina-docking` builds a paid Modal app. The tool explicitly requires user confirmation before the build. No exact build/runtime price or estimator is exposed by the callable Proto metadata; do not invent one. This account setup is unnecessary for the task-local CPU Vina route now being used.

For genuine MD, a separate installed OpenMM/GROMACS environment and suitable compute is needed. The current general executor has no detected NVIDIA GPU and originally lacked the MD packages. A parallel task subsequently installed OpenMM/PDBFixer in `../md/openmm_env` and completed a brief ligand-free FLT4 feasibility benchmark (0.44 ps total, as reported by that task). That benchmark is not ligand–target MD and does not replace a prespecified replicated production simulation or convergence assessment. See the MD task's own logs for actual parameters and status. Proto still exposes no production MD runner. Preserve any future actual trajectory and logs; report run length and completed replicates exactly.

## Interpretation limits

The downloaded AlphaFold structures are fallback receptor scaffolds, not ligand-complex predictions or docking/MD results. PDGFRA kinase-domain mean pLDDT is 74.34; FLT4 77.66. Kinase insertions contribute low-confidence residues. Confidence must be assessed at the binding pocket rather than using a global average as a guarantee. Docking scores and learned binding predictions cannot establish measured binding, target selectivity or DILI causality.
