"""Summarize only actual completed Vina runs; no synthetic missing results."""
import csv,json,itertools
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
P=Path(__file__).resolve().parent
library=list(csv.DictReader((P.parent/'candidate_library'/'Candidate_Screening_Library.csv').open()))
allrows=[];summary=[];contacts=[];pose_qc=[]
for target in ['PDGFRA','FLT4']:
    for lig in library:
        completed=[]
        for seed in [104729,130363,155921]:
            f=P/'runs'/f"{target}_{lig['name']}_{seed}"/'result.json'
            if f.exists():
                d=json.loads(f.read_text());allrows.append({**lig,'target':target,'seed':seed,'vina_score_kcal_mol':d['best_vina_score_kcal_mol'],'elapsed_seconds':d['elapsed_seconds'],'output_dir':d['output_dir']});completed.append(d)
        if not completed:continue
        scores=[d['best_vina_score_kcal_mol'] for d in completed]
        rmsds=[]
        for a,b in itertools.combinations(completed,2):
            ma=Chem.SDMolSupplier(str(P/a['output_dir']/'poses.sdf'),removeHs=True)[0];mb=Chem.SDMolSupplier(str(P/b['output_dir']/'poses.sdf'),removeHs=True)[0]
            rmsds.append(rdMolAlign.CalcRMS(ma,mb,maxMatches=100000))
        summary.append({**lig,'target':target,'completed_searches':len(scores),'median_vina_score_kcal_mol':float(np.median(scores)),'min_vina_score_kcal_mol':min(scores),'max_vina_score_kcal_mol':max(scores),'max_pairwise_top_pose_RMSD_A':max(rmsds) if rmsds else '', 'receptor_type':'experimental_WT_6JOL' if target=='PDGFRA' else 'predicted_AlphaFold_v6'})
        if lig['name']=='nintedanib':
            recpath=P/('PDGFRA_6JOL_receptor.pdb' if target=='PDGFRA' else 'FLT4_AF_receptor.pdb')
            rec=[l for l in recpath.read_text().splitlines() if l.startswith('ATOM') and l[76:78].strip()!='H']
            recxyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in rec])
            for d in completed:
                gatekeeper=674 if target=='PDGFRA' else 927
                hinge_cys=677 if target=='PDGFRA' else 930
                for rank,pose in enumerate(Chem.SDMolSupplier(str(P/d['output_dir']/'poses.sdf'),removeHs=True),1):
                    pxyz=pose.GetConformer().GetPositions()
                    distances={}
                    for n in [gatekeeper,gatekeeper+1,gatekeeper+2,hinge_cys]:
                        coords=recxyz[np.array([int(l[22:26])==n for l in rec])]
                        distances[n]=float(np.linalg.norm(pxyz[:,None,:]-coords[None,:,:],axis=2).min())
                    pose_qc.append({'target':target,'seed':d['seed'],'pose_rank':rank,'vina_score_kcal_mol':pose.GetProp('vina_affinity_kcal_mol'),'gatekeeper_residue':gatekeeper,'minimum_gatekeeper_distance_A':distances[gatekeeper],'hinge_cysteine_residue':hinge_cys,'minimum_hinge_cysteine_distance_A':distances[hinge_cys],'minimum_hinge_region_distance_A':min(distances[gatekeeper+1],distances[gatekeeper+2],distances[hinge_cys])})
                mol=Chem.SDMolSupplier(str(P/d['output_dir']/'poses.sdf'),removeHs=True)[0]
                xyz=mol.GetConformer().GetPositions();dist=np.linalg.norm(recxyz[:,None,:]-xyz[None,:,:],axis=2)
                unique={}
                for i,line in enumerate(rec):
                    md=float(dist[i].min());key=(line[21],int(line[22:26]),line[17:20])
                    if md<4 and (key not in unique or md<unique[key]):unique[key]=md
                for (chain,resnum,resname),mind in unique.items():contacts.append({'target':target,'seed':d['seed'],'chain':chain,'residue_number':resnum,'residue_name':resname,'minimum_heavy_atom_distance_A':round(mind,3),'criterion':'distance<4 A; not hydrogen-bond classification'})
                if d['seed']==104729:
                    ligpdb=Chem.MolToPDBBlock(mol).splitlines();out=[]
                    for line in ligpdb:
                        if line.startswith(('ATOM','HETATM')):out.append('HETATM'+line[6:17]+'LIG B'+f'{1:4d}'+line[26:])
                    (P/f'{target}_nintedanib_best_docked_complex.pdb').write_text('\n'.join(rec+['TER']+out+['END'])+'\n')
for filename,rows in [('Docking_All_Completed_Runs.csv',allrows),('Docking_Summary.csv',summary),('Nintedanib_Docking_Contacts.csv',contacts),('Nintedanib_Pose_QC_All_Searches.csv',pose_qc)]:
    if rows:
        with (P/filename).open('w',newline='') as out:
            w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps({'completed_runs':len(allrows),'expected_runs':96,'compound_target_summaries':len(summary)},indent=2))
