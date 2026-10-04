"""Add absent heavy atoms only, without reconstructing missing loops or moving deposited atoms."""
import json,os
os.environ['OPENMM_CPU_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
from pdbfixer import PDBFixer
from openmm import app,unit
P=Path(__file__).resolve().parent
source=P/'PDGFRA_8PQJ_receptor.pdb'
fixer=PDBFixer(filename=str(source))
before={(a.residue.chain.id,a.residue.id,a.name):p.value_in_unit(unit.angstrom) for a,p in zip(fixer.topology.atoms(),fixer.positions)}
fixer.findMissingResidues();fixer.missingResidues={}
fixer.findMissingAtoms()
missing={r.chain.id+':'+r.id:[a.name for a in atoms] for r,atoms in fixer.missingAtoms.items()}
terminals={r.chain.id+':'+r.id:list(atoms) for r,atoms in fixer.missingTerminals.items()}
fixer.addMissingAtoms(seed=104729)
after={(a.residue.chain.id,a.residue.id,a.name):p.value_in_unit(unit.angstrom) for a,p in zip(fixer.topology.atoms(),fixer.positions)}
displacements={':'.join(k):float(np.linalg.norm(np.asarray(after[k])-v)) for k,v in before.items()}
maximum=max(displacements.values())
if maximum>1e-4:raise ValueError(f'Deposited atoms moved by {maximum} A')
with (P/'8PQJ_fixed_heavy.pdb').open('w') as f:app.PDBFile.writeFile(fixer.topology,fixer.positions,f,keepIds=True)
report={'software':'PDBFixer1.12.0/OpenMM8.4.0','seed':104729,'missing_residue_reconstruction':False,'missing_heavy_atoms':missing,'missing_terminal_atoms':terminals,'initial_atoms':len(before),'repaired_atoms':len(after),'max_deposited_atom_displacement_A':maximum,'added_atoms':list(':'.join(k) for k in after if k not in before)}
(P/'8PQJ_heavy_atom_repair_QC.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
