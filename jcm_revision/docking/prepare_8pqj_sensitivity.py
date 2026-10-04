import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from Bio.SeqUtils import seq1
from prepare_receptors import P,ca_records,fit,pdbqt_receptor

native=json.loads((P.parent.parent/'binding_analysis'/'PDGFRA_UniProt.json').read_text())['sequence']['value']
input_pdb=P/'8PQJ_fixed_heavy.pdb'
if not input_pdb.exists():raise FileNotFoundError('First run repair_8pqj_heavy_atoms.py; original 8PQJ has incomplete side chains.')
lines=[l for l in input_pdb.read_text().splitlines() if l.startswith('ATOM') and l[21]=='A']
observed={int(l[22:26]):seq1(l[17:20]) for l in lines if l[12:16].strip()=='CA'}
errors=[(n,aa,native[n-1]) for n,aa in observed.items() if aa!=native[n-1]]
if errors:raise ValueError(errors)
assert observed[674]=='T' and observed[842]=='D'
source,_=ca_records(P/'PDGFRA_6JOL_receptor.pdb');target,_=ca_records(P/'8PQJ.pdb')
src={r.id[1]:r for r in source};tar={r.id[1]:r for r in target};common=sorted(src.keys()&tar.keys())
x=np.array([src[n]['CA'].coord for n in common]);y=np.array([tar[n]['CA'].coord for n in common]);mask=np.ones(len(common),bool)
for _ in range(8):
 r,t=fit(x[mask],y[mask]);dist=np.linalg.norm(x@r+t-y,axis=1);new=dist<max(2.5,np.quantile(dist,0.75))
 if np.array_equal(new,mask):break
 mask=new
r,t=fit(x[mask],y[mask]);dist=np.linalg.norm(x@r+t-y,axis=1)
lig=Chem.SDMolSupplier(str(P/'6JOL_STI_reference.sdf'),removeHs=True)[0];xyz=lig.GetConformer().GetPositions()@r+t
for i,c in enumerate(xyz):lig.GetConformer().SetAtomPosition(i,c.tolist())
with Chem.SDWriter(str(P/'8PQJ_transferred_STI_box_reference.sdf')) as w:w.write(lig)
center=(xyz.min(0)+xyz.max(0))/2;size=np.maximum(xyz.max(0)-xyz.min(0)+12,22)
count=pdbqt_receptor('PDGFRA_8PQJ_receptor','\n'.join(lines)+'\nEND\n')
record={'receptor':'PDGFRA_8PQJ_receptor.pdbqt','source':'8PQJ WT apo PDGFRA, 1.82 A','box_template':'6JOL imatinib transferred by matched native-residue C-alpha core alignment; box only','box_center':center.tolist(),'box_size':size.tolist(),'residue_count':count,'sequence_mismatches':errors,'native_T674':observed[674],'native_D842':observed[842],'alignment_matched_CA':len(common),'alignment_core_CA':int(mask.sum()),'alignment_core_RMSD_A':float(np.sqrt(np.mean(dist[mask]**2))),'transform_rotation':r.tolist(),'transform_translation':t.tolist(),'core_residue_numbers':[n for n,ok in zip(common,mask) if ok],'purpose':'Separate prespecified receptor-conformation sensitivity; never replaces primary 6JOL panel'}
manifest=json.loads((P/'receptor_manifest.json').read_text());manifest['PDGFRA_APO_8PQJ']=record;(P/'receptor_manifest.json').write_text(json.dumps(manifest,indent=2));(P/'PDGFRA_8PQJ_preparation_QC.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
