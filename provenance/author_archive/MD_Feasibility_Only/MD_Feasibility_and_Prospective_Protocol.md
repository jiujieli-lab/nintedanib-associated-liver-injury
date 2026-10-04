# MD feasibility and prospective protocol

Audit date: 2 October 2026. This is an author-facing resource and preparation assessment. **No PDGFRA–nintedanib or FLT4–nintedanib production MD has been performed.** The completed 0.44 ps ligand-free runtime test is excluded from manuscript binding results and figures.

Production MD is technically possible with OpenMM, but the current shared executor is not an appropriate place to complete the proposed replicated campaign. It has an eight-core CPU quota, 8 GiB RAM and no exposed GPU. The docking campaign uses the available CPU quota. A small isolated OpenMM installation and a one-thread, explicit-solvent protein benchmark succeeded; the full OpenFF/AmberTools environment did not finish solving within available memory. Standard ligand parameters have consequently **not** been generated or verified.

## What actually ran

| Item | Observed result |
|---|---|
| Isolated environment | `openmm_env/bin/python`; Python 3.12.14, OpenMM 8.4.0, PDBFixer 1.12.0, NumPy 2.2.6 |
| Available platforms | Reference and CPU; no CUDA/OpenCL GPU device exposed |
| Ligand parameter environment attempt | Micromamba 2.3.3, conda-forge only; requested Python 3.11 / OpenMM 8.4.0 / OpenFF Toolkit 0.17.1 / AmberTools 24.8 / OpenMMForceFields 0.15.0 / PDBFixer 1.11 |
| Full environment outcome | Dry-run solve killed with exit 137; cgroup recorded an OOM kill. No compatible full environment was established. It was not retried. |
| Benchmark input | Existing ligand-free AlphaFold FLT4 fragment, UniProt residues 845–1173; unmodified original coordinate source |
| Benchmark preparation | Added terminal OXT and hydrogens at nominal pH 7.4, Amber ff14SB/TIP3P, neutralized with 150 mM salt, rectangular box with 1.2 nm padding per side |
| Benchmark size | 5,193 protein atoms including H; 57,629 total atoms; box 9.991 × 8.356 × 7.350 nm |
| Dynamics settings | CPU, one thread; PME, 1.0 nm cutoff; H-bond constraints, rigid water; LangevinMiddle 310 K, 1/ps friction, 2 fs timestep; no barostat in this runtime-only test |
| Runtime | Setup 10.54 s; limited minimization 175.27 s; 20 warm-up and 200 timed steps; 0.44 ps integrated in total |
| Timed throughput | 0.5602 and 0.6183 ns/day in two 100-step blocks; median 0.5893 ns/day |
| Integrity | Finite potential energies in both timed blocks; script, raw log, input SHA-256 and machine-readable manifest retained |

The throughput measurements are short, share CPU with docking, and use a ligand-free scaffold with uncertain loops. They support only an order-of-magnitude resource decision. They are not a statistically precise performance estimate, an equilibration test, a validated protein conformation, or evidence for ligand binding. No trajectory or ligand parameter file was produced. The saved `prepared_protein_only_benchmark.pdb` contains the solvated **starting** coordinates, not a final equilibrated structure.

## Starting structures still need scientific preparation

The exact measurements and source hashes are in `Structure_Preparation_QA.json` and its accompanying report.

* PDGFRA/6JOL coordinates contain A583–693 followed by A790–957, with a 6.423 Å C693–N790 separation and no TER separator. Native residues 697–768 were deliberately deleted from the crystallographic construct; 694–696 and 769–789 are unresolved. The construct deletion and unresolved residues must not be conflated. Directly inferring a peptide bond across the coordinate gap is invalid. A deletion-construct model needs the 24 unresolved residues reconstructed; a native-domain model instead needs its complete sequence and uncertainty documented.
* The current PDGFRA rank-1 docking pose is 8.428 Å from gatekeeper Thr674 and 10.741 Å from hinge residues 675–677. Rank 3 is closer (3.380 Å to Thr674; 5.266 Å to the hinge) and only 0.161 kcal/mol worse by docking score. Neither distance establishes a canonical hinge hydrogen bond. Review both hypotheses by chemical and structural criteria; selecting rank 1 alone is not justified by that small score difference. Correct duplicated protein/ligand atom serials before topology construction.
* The FLT4 AlphaFold fragment has continuous peptide geometry but artificial cropped termini and 45 residues with pLDDT below 50. Several ligand-contact residues in the glycine-rich region have pLDDT below 70. Low-confidence loops cannot simply be deleted without defining new termini and changing the model.
* Complete Boltz kinase-domain models could avoid the 6JOL construct gap only if the predicted sequence, backbone, pocket geometry and confidence pass independent checks. All 12 returned PDGFRA and all 12 FLT4 CIF exports contain only 199/362 and 195/329 declared residues, respectively; the shared cropping patterns are documented in Boltz_All24_Structure_QA.json. These cropped pocket models do not solve the MD continuity problem. Their atom B fields are all 100.000 and must not be interpreted as local pLDDT. The current focused-library Boltz screen filtered nintedanib and does **not** supply a nintedanib-bound prediction; its reference-ligand input is not an output complex. No candidate model is accepted for MD solely because Boltz produced it.
* Preserve nintedanib's selected Z keto-enamine identity and its exact atom mapping. PDB-only ligand records do not encode sufficient stereochemistry or bond orders; use the validated pose SDF and explicit connectivity. Review protonation and tautomer assignment. Existing pKa predictions indicate that the neutral docking microstate cannot automatically be treated as the only relevant pH 7.4 species. Investigate site-specific protonation and a plausible +1 microstate before choosing a primary simulation, and record the choice and sensitivity analysis.
* Define terminal chemistry, histidine states, salt bridges, disulfides if present, missing side chains, retained crystallographic waters and any cofactors. Explicitly inspect DFG, HRD, αC and hinge geometry. Geometry repair must preserve the intended receptor construct and ligand coordinate frame.

An experimental nintedanib–KDR/VEGFR2 complex may help validate the preparation/analysis pipeline. It cannot stand in for PDGFRA or FLT4 target evidence. An isolated kinase-domain simulation also cannot determine full-length receptor activation, lymphatic function, or hepatoprotection.

## Force-field route and software limits

The preferred conventional route is Amber ff14SB protein, GAFF2 2.2.20 ligand with validated AM1-BCC charges, compatible TIP3P/ion parameters, and OpenMM. An explicitly pinned OpenFF Sage release is an alternative ligand force field, with its charge method and compatibility separately checked. Either route requires all atoms to receive documented parameters, finite charges summing to the intended formal charge, preserved stereochemistry, and geometry/energy inspection. Partial charges generated for docking must not be reused as an undocumented MD charge model.

OpenFF recommends an isolated conda-forge environment; its free AM1-BCC workflow normally uses AmberTools. The unsuccessful local solve does not show these packages are intrinsically incompatible. It shows that installation and parameterization were not established on this memory-limited executor. On the intended compute host, resolve and test a compatible pinned environment, then save its exact lock/export and the ligand parameter cache. Do not label the failed requested pin set as a validated environment. Current OpenMMForceFields releases have different OpenMM requirements from older releases, so package versions must be checked together rather than mixed silently. [1–3]

Espaloma 0.3.2 is a documented learned-force-field option supported by OpenMMForceFields. Its workflow still uses OpenFF molecule handling and PyTorch/DGL dependencies, and the template-generator interface is described as experimental. It was considered, not installed or validated here. Its existence does not establish that a lighter pip-only installation will work in this executor or that its nintedanib parameterization is adequate. A learned charge model would be a declared methodological change requiring its own inspection, not a way to bypass missing charges. [2,4,5]

The existing `work/binding_analysis/run_md_after_structure_qc.py` is a prospective starting script, not a ready-validated production workflow. It requires repaired, hydrogen-complete protein and a bound-pose SDF in the same coordinate frame. Before deployment it still needs tested parameter generation, explicit restart/error handling, improved trajectory storage, and independent preparation review. Its `--structure-qc-approved` flag is a record of a review, not a replacement for that review. No production invocation was made.

## Minimum useful initial campaign and stopping rule

There is no universal number of nanoseconds that proves convergence. For a publishable **exploratory local-pose persistence analysis**, budget an initial **three independently initialized 100 ns replicas per target–ligand complex**, after equilibration. This is a project planning recommendation, not a validated minimum or a guarantee of binding stability. Independent velocities and seeds are necessary; independent equilibrated starts and alternative credible poses/loop models are useful where structural uncertainty is substantial. Separate simulations are required for materially different protonation or starting-pose hypotheses. [6]

| Planned scope | Production aggregate |
|---|---:|
| One selected nintedanib pose/microstate in each of two targets; 3 × 100 ns each | 600 ns |
| Above plus matched apo receptor simulations; 3 × 100 ns for each receptor | 1,200 ns |
| Above plus a well-defined positive-control target–ligand system; 3 × 100 ns | 1,500 ns |

For a kinase-ligand preparation control, the PDGFRA–imatinib crystallographic construct can be considered only after its own construct/gap repair. A different control or an experimentally determined homolog validates the pipeline within that system; it does not validate the other target's predicted binding pose. Adding each alternative pose or microstate at 3 × 100 ns adds 300 ns per receptor. Screening all 16 compounds with three replicas per target would be 9.6 µs before apo controls, and is outside the justified first-stage scope.

Prospective settings: minimization with justified solute restraints, gradual thermalization to 310 K, approximately 0.2 ns restrained NVT followed by 1–2 ns staged NPT at 1 bar and restraint release; explicit TIP3P water, at least 1.2 nm solute padding, neutralization plus 150 mM salt; PME with 1.0 nm real-space cutoff, 2 fs integration with hydrogen constraints. Those durations are starting settings: extend equilibration until density, temperature, energy and relevant structural observables settle. They do not waive structural review or provide evidence of association.

Assess backbone and pocket-aligned ligand heavy-atom RMSD, displacement from the pocket, ligand torsions, hinge/contact occupancy with chemical definitions, pocket hydration and residue fluctuations. Report every independent replica and pose/microstate condition. Use blocks longer than observed correlation times and quantify between-replica uncertainty. Compare early/late segments and replica conformational distributions; extend or change sampling when they disagree materially. Stable RMSD alone does not prove affinity, and absence of dissociation in 100 ns does not establish a favorable binding free energy. Conventional trajectories cannot calibrate the screen to experimental Kd or determine DILI prevention. [6]

## Credible compute and storage plan

A GPU workstation or an existing institutional HPC allocation with a compatible OpenMM CUDA/OpenCL installation is the practical route. A planning request of one modern supported GPU, 4 host CPU cores, 16–32 GiB host RAM and at least 50 GiB free scratch is reasonable for these roughly 60,000-atom systems; exact needs and throughput must be measured on the final prepared complex. This is a resource request, not a claim that access exists. No paid service, external job or deployment was initiated. The available Proto audit identified no deployed conventional MD runner; Boltz and BioEmu outputs are not substitutes for a conventional explicit-solvent MD campaign.

First run the final prepared system briefly on the chosen device, checking energy/force validity and measuring sustained ns/day after compilation and equilibration. At a measured aggregate production rate R ns/day, wall time is aggregate ns/R, plus setup, equilibration, analysis and queue overhead. If the verified device sustains 50, 100 or 200 ns/day, a 600 ns campaign requires 12, 6 or 3 production device-days respectively. These are **conditional arithmetic scenarios**, not measured GPU forecasts. Parallel replicas reduce elapsed time only when separately allocated resources provide that throughput.

At the observed shared one-thread rate, 100 ns would take about 170 days and 600 ns about 1,018 days. These figures explain the decision against local production; they are not estimates for an idle eight-core machine. Linear eightfold scaling was not assumed or tested.

For the measured 57,629-atom system, saving uncompressed float32 coordinates every 10 ps would use about 41.5 GB for 600 ns, before file overhead, checkpoints and analysis. Saving solute-only coordinates every 10 ps and whole-system coordinates every 100 ps lowers that estimate to about 7.9 GB, before overhead, using roughly 5,200 solute atoms. Use compatible topology files, periodic-box metadata, solvent snapshots for hydration analysis, frequent restart checkpoints and a full final state. Check analysis requirements before reducing saved water frames. Save input structures, exact atom mappings, charge/parameter files, simulation XML, software lock, seeds and preparation manifests with checksums.

## Sources

1. Open Force Field Initiative. Installation documentation. https://docs.openforcefield.org/en/latest/install.html (accessed 2026-10-02).
2. OpenMMForceFields maintainers. Supported force fields, template generators and software requirements. https://github.com/openmm/openmmforcefields/blob/main/README.md (accessed 2026-10-02).
3. OpenMM developers. Getting Started. https://docs.openmm.org/latest/userguide/application/01_getting_started.html (accessed 2026-10-02).
4. Espaloma developers. Software and installation. https://github.com/choderalab/espaloma (accessed 2026-10-02).
5. Takaba K et al. Machine-learned molecular mechanics force field for the simulation of protein–ligand systems and beyond. https://arxiv.org/abs/2307.07085. The model is a methodological alternative, not one used in this work.
6. Grossfield A et al. Best Practices for Quantification of Uncertainty and Sampling Quality in Molecular Simulations. Living J Comput Mol Sci. 2018;1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067; https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067.
