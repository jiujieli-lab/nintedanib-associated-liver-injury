"""Audit complete raw docking outputs and create a compact, reproducible run inventory."""
import csv
import hashlib
import json
import math
from pathlib import Path

from rdkit import Chem

P = Path(__file__).resolve().parent
library_path = P.parent / 'candidate_library' / 'Candidate_Screening_Library.csv'
library = list(csv.DictReader(library_path.open()))
seeds = [104729, 130363, 155921]
groups = {
    'primary': [(target, ligand['name'], seed) for target in ['PDGFRA', 'FLT4']
                for ligand in library for seed in seeds],
    'apo_sensitivity': [('PDGFRA_APO_8PQJ', 'nintedanib', seed) for seed in seeds],
    'separate_preflight_control': [('PDGFRA', 'Imatinib_redocking', 104729)],
}
inventory = []
for group, jobs in groups.items():
    for target, name, seed in jobs:
        folder = P / 'runs' / f'{target}_{name}_{seed}'
        result = json.loads((folder / 'result.json').read_text())
        poses = list(Chem.SDMolSupplier(str(folder / 'poses.sdf'), removeHs=True))
        assert poses and all(pose is not None for pose in poses), str(folder)
        scores = [float(pose.GetProp('vina_affinity_kcal_mol')) for pose in poses]
        assert all(math.isfinite(score) for score in scores), str(folder)
        assert len(poses) == len(result['pose_energies']), str(folder)
        assert scores == [row[0] for row in result['pose_energies']], str(folder)
        assert scores[0] == result['best_vina_score_kcal_mol'], str(folder)
        assert result['ligand_MMFF_converged'], str(folder)
        assert result['seed'] == seed and result['exhaustiveness'] == 32, str(folder)
        assert result['target'] == target and result['name'] == name, str(folder)
        inventory.append({'group': group, 'target': target, 'name': name, 'seed': seed,
                          'retained_pose_count': len(poses),
                          'best_vina_score_kcal_mol': scores[0],
                          'output_dir': str(folder.relative_to(P))})
with (P / 'Final_Docking_Run_Inventory.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(inventory[0]))
    writer.writeheader()
    writer.writerows(inventory)
summary = {
    'searches_by_group': {group: sum(row['group'] == group for row in inventory) for group in groups},
    'poses_by_group': {group: sum(row['retained_pose_count'] for row in inventory if row['group'] == group) for group in groups},
    'all_pose_SDFs_parse': True,
    'all_pose_scores_finite_and_match_raw_results': True,
    'all_ligand_MMFF_minimizations_converged': True,
    'primary_and_sensitivity_fixed_seeds': seeds,
    'library_sha256': hashlib.sha256(library_path.read_bytes()).hexdigest(),
    'interpretation': 'Empirical rigid-receptor docking; stochastic searches are not biological replicates or MD. Sensitivity results do not replace the primary campaign.',
}
(P / 'Final_Docking_Execution_QC.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
