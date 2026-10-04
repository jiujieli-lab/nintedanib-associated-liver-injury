"""Actual AutoDock Vina rigid-receptor docking, with retained inputs and outputs."""
import argparse,csv,json,time
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolAlign
from meeko import MoleculePreparation,PDBQTWriterLegacy,PDBQTMolecule,RDKitMolCreate
from vina import Vina

P=Path(__file__).resolve().parent

def dock(name,smiles,target,seed=104729,exhaustiveness=32,cpu=4,reference=None):
    rec=json.loads((P/'receptor_manifest.json').read_text())[target]
    out=P/'runs'/f'{target}_{name}_{seed}';out.mkdir(parents=True,exist_ok=True)
    if (out/'result.json').exists():return json.loads((out/'result.json').read_text())
    m=Chem.MolFromSmiles(smiles)
    if m is None:raise ValueError('Invalid SMILES')
    canonical=Chem.MolToSmiles(m,True);m=Chem.AddHs(m)
    params=AllChem.ETKDGv3();params.randomSeed=seed
    if AllChem.EmbedMolecule(m,params)!=0:raise ValueError('3D embedding failed')
    converged=AllChem.MMFFOptimizeMolecule(m,maxIters=2000)
    with Chem.SDWriter(str(out/'input_ligand.sdf')) as w:w.write(m)
    setups=MoleculePreparation().prepare(m)
    if len(setups)!=1:raise ValueError('Unexpected multiple ligand setups')
    pdbqt,success,errors=PDBQTWriterLegacy.write_string(setups[0])
    if not success:raise ValueError(errors)
    (out/'input_ligand.pdbqt').write_text(pdbqt)
    config={'target':target,'name':name,'canonical_smiles':canonical,'formal_charge':Chem.GetFormalCharge(m),'seed':seed,'exhaustiveness':exhaustiveness,'cpu':cpu,'n_poses':10,'energy_range':5,'scoring_function':'vina','ligand_MMFF_converged':converged==0,'receptor':rec['receptor'],'box_center':rec['box_center'],'box_size':rec['box_size']}
    (out/'config.json').write_text(json.dumps(config,indent=2))
    start=time.monotonic();v=Vina(sf_name='vina',cpu=cpu,seed=seed,verbosity=1)
    v.set_receptor(str(P/rec['receptor']));v.set_ligand_from_string(pdbqt)
    v.compute_vina_maps(center=rec['box_center'],box_size=rec['box_size'])
    v.dock(exhaustiveness=exhaustiveness,n_poses=10,min_rmsd=1)
    v.write_poses(str(out/'poses.pdbqt'),n_poses=10,energy_range=5,overwrite=True)
    energies=v.energies(n_poses=10,energy_range=5).tolist()
    pq=PDBQTMolecule((out/'poses.pdbqt').read_text(),skip_typing=True)
    mols=RDKitMolCreate.from_pdbqt_mol(pq)
    if len(mols)!=1 or mols[0] is None:raise ValueError('Pose conversion failed')
    pose=mols[0]
    rmsds=[]
    ref=Chem.SDMolSupplier(str(P/reference),removeHs=True)[0] if reference else None
    heavy=Chem.RemoveHs(pose)
    with Chem.SDWriter(str(out/'poses.sdf')) as w:
        for i in range(pose.GetNumConformers()):
            pose.SetProp('vina_affinity_kcal_mol',str(energies[i][0]));pose.SetProp('pose_rank',str(i+1));w.write(pose,confId=i)
            if ref is not None:rmsds.append(float(rdMolAlign.CalcRMS(heavy,ref,prbId=i,refId=0,maxMatches=100000)))
    result={**config,'elapsed_seconds':time.monotonic()-start,'pose_energies':energies,'best_vina_score_kcal_mol':energies[0][0],'reference_RMSD_A_receptor_frame':rmsds,'output_dir':str(out.relative_to(P))}
    (out/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--control',action='store_true');a.add_argument('--library');a.add_argument('--only-name');a.add_argument('--target',choices=['PDGFRA','FLT4','PDGFRA_APO_8PQJ']);a.add_argument('--seed',type=int,default=104729);a.add_argument('--exhaustiveness',type=int,default=32);a.add_argument('--cpu',type=int,default=4)
    args=a.parse_args()
    if args.control:
        rec=json.loads((P/'receptor_manifest.json').read_text())['PDGFRA'];dock('Imatinib_redocking',rec['reference_smiles'],'PDGFRA',args.seed,args.exhaustiveness,args.cpu,rec['reference_ligand'])
    elif args.library:
        rows=list(csv.DictReader(open(args.library)))
        for row in rows:
            if args.only_name is None or row['name']==args.only_name:dock(row['name'],row['smiles'],args.target,args.seed,args.exhaustiveness,args.cpu)
    else:a.error('Specify --control or --library')
