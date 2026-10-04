"""CPU feasibility only: ligand-free solvated FLT4 scaffold, not target-binding MD."""
import json,time,pathlib,os,hashlib,importlib.metadata as md
import numpy as np
import openmm as mm
from openmm import app,unit
from pdbfixer import PDBFixer
ROOT=pathlib.Path(__file__).parent
src=ROOT.parent/'docking/FLT4_AF_receptor.pdb'
out=ROOT/'benchmark';out.mkdir(exist_ok=True)
manifest={'purpose':'Runtime/preparation feasibility only; ligand-free predicted scaffold; not evidence of nintedanib binding or stability','status':'started','versions':{x:md.version(x) for x in ['openmm','pdbfixer','numpy']},'platform':'CPU','threads':1,'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'source':str(src),'shared_executor_docking_campaign_running':True}
def save(): (out/'benchmark_manifest.json').write_text(json.dumps(manifest,indent=2))
save();t=time.monotonic()
fixer=PDBFixer(filename=str(src));fixer.findMissingResidues();fixer.findMissingAtoms()
manifest['detected_missing_residues']={str(k):v for k,v in fixer.missingResidues.items()}
manifest['detected_missing_atoms']={f'{r.chain.id}{r.id}{r.name}':[a.name for a in v] for r,v in fixer.missingAtoms.items()}
manifest['terminal_additions']={f'{r.chain.id}{r.id}{r.name}':v for r,v in fixer.missingTerminals.items()}
if fixer.missingResidues:raise ValueError('Internal residue reconstruction is outside benchmark scope')
fixer.addMissingAtoms(seed=20261002)
ff=app.ForceField('amber14/protein.ff14SB.xml','amber14/tip3p.xml')
model=app.Modeller(fixer.topology,fixer.positions)
model.addHydrogens(ff,pH=7.4)
positions=np.asarray(model.positions.value_in_unit(unit.nanometer))
box=positions.max(axis=0)-positions.min(axis=0)+2.4
manifest['protein_atoms_including_H']=model.topology.getNumAtoms()
model.addSolvent(ff,model='tip3p',boxSize=mm.Vec3(*box)*unit.nanometer,ionicStrength=.15*unit.molar,neutralize=True)
manifest['box_lengths_nm']=box.tolist();manifest['solvated_atoms']=model.topology.getNumAtoms();manifest['solvated_residues']=model.topology.getNumResidues();save()
system=ff.createSystem(model.topology,nonbondedMethod=app.PME,nonbondedCutoff=1.0*unit.nanometer,constraints=app.HBonds,rigidWater=True,ewaldErrorTolerance=5e-4)
integ=mm.LangevinMiddleIntegrator(310*unit.kelvin,1/unit.picosecond,.002*unit.picoseconds);integ.setRandomNumberSeed(20261002)
sim=app.Simulation(model.topology,system,integ,mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
sim.context.setPositions(model.positions)
manifest['setup_seconds']=time.monotonic()-t;save();print('SETUP',manifest['solvated_atoms'],manifest['setup_seconds'],flush=True)
t=time.monotonic();sim.minimizeEnergy(maxIterations=100)
manifest['limited_minimization_seconds']=time.monotonic()-t;save();print('MINIMIZED',flush=True)
sim.context.setVelocitiesToTemperature(310*unit.kelvin,20261002)
sim.step(20)
timings=[];t0=time.monotonic()
for repeat in range(3):
 t=time.monotonic();sim.step(100);elapsed=time.monotonic()-t
 state=sim.context.getState(getEnergy=True)
 epot=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
 if not np.isfinite(epot):raise ValueError('Nonfinite potential energy')
 timings.append({'steps':100,'simulated_ps':.2,'wall_seconds':elapsed,'ns_per_day':.0002*86400/elapsed,'potential_energy_kJ_per_mol':epot})
 manifest['timings']=timings;save();print('TIMED',timings[-1],flush=True)
 if time.monotonic()-t0>45:break
manifest['status']='completed_feasibility_benchmark';manifest['total_integrated_ps']=sim.currentStep*.002;manifest['median_ns_per_day']=float(np.median([x['ns_per_day'] for x in timings]));manifest['interpretation']='One-thread result measured under shared docking load, on protein-only scaffold. Not a production equilibration, target-ligand simulation, or convergence demonstration. No ligand force-field parameters were generated.'
save()
with (out/'prepared_protein_only_benchmark.pdb').open('w') as h:app.PDBFile.writeFile(model.topology,model.positions,h)
