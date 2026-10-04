import csv,itertools,json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
P=Path(__file__).resolve().parent
rows=[];poses=[];summary=[]
for label,target,recfile in [('Primary_6JOL','PDGFRA','PDGFRA_6JOL_receptor.pdb'),('Sensitivity_8PQJ','PDGFRA_APO_8PQJ','PDGFRA_8PQJ_receptor.pdb')]:
    lines=[l for l in (P/recfile).read_text().splitlines() if l.startswith('ATOM') and l[76:78].strip()!='H']
    coords={n:np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lines if int(l[22:26])==n]) for n in [674,675,676,677]}
    completed=[]
    for seed in [104729,130363,155921]:
        folder=P/'runs'/f'{target}_nintedanib_{seed}'
        if not (folder/'result.json').exists():continue
        d=json.loads((folder/'result.json').read_text());ms=list(Chem.SDMolSupplier(str(folder/'poses.sdf'),removeHs=True));completed.append((d,ms[0]))
        for rank,m in enumerate(ms,1):
            xyz=m.GetConformer().GetPositions();distances={n:float(np.linalg.norm(xyz[:,None,:]-c[None,:,:],axis=2).min()) for n,c in coords.items()}
            v={'receptor':label,'target':target,'seed':seed,'pose_rank':rank,'vina_score_kcal_mol':float(m.GetProp('vina_affinity_kcal_mol')),'minimum_Thr674_gatekeeper_distance_A':distances[674],'minimum_hinge_675_677_distance_A':min(distances[675],distances[676],distances[677]),'minimum_Cys677_distance_A':distances[677]}
            poses.append(v)
            if rank==1:rows.append(v)
        if label=='Sensitivity_8PQJ' and seed==104729:
            mol=ms[0];lig=[]
            for l in Chem.MolToPDBBlock(mol).splitlines():
                if l.startswith(('ATOM','HETATM')):lig.append('HETATM'+l[6:17]+'LIG B'+f'{1:4d}'+l[26:])
            (P/'PDGFRA_8PQJ_nintedanib_best_docked_complex.pdb').write_text('\n'.join(lines+['TER']+lig+['END'])+'\n')
    scores=[d['best_vina_score_kcal_mol'] for d,m in completed]
    pairrms=[float(rdMolAlign.CalcRMS(a[1],b[1],maxMatches=100000)) for a,b in itertools.combinations(completed,2)]
    if scores:summary.append({'receptor':label,'n_searches_completed':len(scores),'median_vina_score_kcal_mol':float(np.median(scores)),'minimum_score':min(scores),'maximum_score':max(scores),'max_top_pose_pairwise_RMSD_A':max(pairrms) if pairrms else None})
for name,data in [('PDGFRA_Nintedanib_Receptor_Sensitivity.csv',rows),('PDGFRA_Nintedanib_Sensitivity_All_Poses.csv',poses)]:
    if data:
        with (P/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
(P/'PDGFRA_8PQJ_Sensitivity_Summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
