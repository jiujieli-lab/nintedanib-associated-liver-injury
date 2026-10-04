"""Prepare rigid receptors and reference-ligand boxes. No docking results generated here."""
import json
from pathlib import Path
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import Polymer, PDBQTWriterLegacy

P=Path(__file__).resolve().parent

def ligand_from_pdb(pid, code):
    lines=[l for l in (P/f'{pid}.pdb').read_text().splitlines() if l.startswith('HETATM') and l[17:20]==code]
    block='\n'.join(lines)+'\nEND\n'
    (P/f'{pid}_{code}_reference.pdb').write_text(block)
    smiles=json.loads((P/f'{code}_chemcomp.json').read_text())['rcsb_chem_comp_descriptor']['SMILES_stereo']
    mol=Chem.MolFromPDBBlock(block,removeHs=True)
    if mol is None: raise ValueError('Ligand parse failed')
    mol=AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles),mol)
    with Chem.SDWriter(str(P/f'{pid}_{code}_reference.sdf')) as w:w.write(mol)
    return mol,smiles

def pdbqt_receptor(name, text):
    (P/f'{name}.pdb').write_text(text)
    polymer=Polymer.from_pdb_string(text,allow_bad_res=False,default_altloc='A')
    rigid,flex=PDBQTWriterLegacy.write_from_polymer(polymer)
    assert not flex
    (P/f'{name}.pdbqt').write_text(rigid)
    return len(polymer.monomers)

def ca_records(path):
    chain=PDBParser(QUIET=True).get_structure('s',str(path))[0]['A']
    allowed={'CME':'C','PTR':'Y'}
    rs=[r for r in chain if 'CA' in r and (r.id[0]==' ' or r.resname in allowed)]
    return rs,''.join(allowed.get(r.resname,seq1(r.resname)) for r in rs)

def fit(x,y):
    xc=x.mean(0);yc=y.mean(0)
    u,s,vt=np.linalg.svd((x-xc).T@(y-yc));r=u@vt
    if np.linalg.det(r)<0:u[:,-1]*=-1;r=u@vt
    return r,yc-xc@r

def main():
    report={}
    imatinib,sti_smiles=ligand_from_pdb('6JOL','STI')
    xyz=imatinib.GetConformer().GetPositions()
    center=(xyz.max(0)+xyz.min(0))/2;size=np.maximum(xyz.max(0)-xyz.min(0)+12,22)
    lines=[l for l in (P/'6JOL.pdb').read_text().splitlines() if l.startswith('ATOM') and l[21]=='A']
    n=pdbqt_receptor('PDGFRA_6JOL_receptor','\n'.join(lines)+'\nEND\n')
    report['PDGFRA']={'receptor':'PDGFRA_6JOL_receptor.pdbqt','source':'6JOL, WT PDGFRA imatinib complex, 1.90 A','reference_ligand':'6JOL_STI_reference.sdf','reference_smiles':sti_smiles,'box_center':center.tolist(),'box_size':size.tolist(),'residue_count':n,'caveat':'Experimental construct deletes kinase insert and contains unresolved loops; no WT-to-mutant substitution; all deposited protein heavy atoms retained, waters and ligand removed.'}

    xin,xin_smiles=ligand_from_pdb('3C7Q','XIN')
    af=P.parent/'proto'/'FLT4_AF_kinase_domain.pdb'
    ref,seqx=ca_records(P/'3C7Q.pdb');tar,seqy=ca_records(af)
    aligner=Align.PairwiseAligner();aligner.substitution_matrix=substitution_matrices.load('BLOSUM62');aligner.open_gap_score=-10;aligner.extend_gap_score=-0.5
    al=aligner.align(seqx,seqy)[0];(P/'FLT4_to_VEGFR2_alignment.txt').write_text(str(al))
    pairs=[(i,j) for (a,b),(c,d) in zip(*al.aligned) for i,j in zip(range(a,b),range(c,d))]
    x=np.array([ref[i]['CA'].coord for i,j in pairs]);y=np.array([tar[j]['CA'].coord for i,j in pairs]);mask=np.ones(len(pairs),bool)
    for _ in range(8):
        r,t=fit(x[mask],y[mask]);dist=np.linalg.norm(x@r+t-y,axis=1);new=dist<max(2.5,np.quantile(dist,0.75))
        if np.array_equal(mask,new):break
        mask=new
    r,t=fit(x[mask],y[mask]);dist=np.linalg.norm(x@r+t-y,axis=1)
    newmol=Chem.Mol(xin);conf=newmol.GetConformer()
    transfer=xin.GetConformer().GetPositions()@r+t
    for i,v in enumerate(transfer):conf.SetAtomPosition(i,v.tolist())
    with Chem.SDWriter(str(P/'FLT4_transferred_XIN_box_reference.sdf')) as w:w.write(newmol)
    center=(transfer.max(0)+transfer.min(0))/2;size=np.maximum(transfer.max(0)-transfer.min(0)+12,22)
    n=pdbqt_receptor('FLT4_AF_receptor',af.read_text())
    residues=[]
    for res in tar:
        heavy=np.array([a.coord for a in res if a.element!='H'])
        if np.min(np.linalg.norm(heavy[:,None,:]-transfer[None,:,:],axis=2))<6:residues.append((int(res.id[1]),res.resname,float(res['CA'].bfactor)))
    report['FLT4']={'receptor':'FLT4_AF_receptor.pdbqt','source':'AlphaFold AF-P35916-F1 v6; kinase domain 845-1173','box_template':'VEGFR2 3C7Q XIN; homolog ligand transferred only to define ATP-site box, not a predicted FLT4 ligand pose','box_center':center.tolist(),'box_size':size.tolist(),'residue_count':n,'alignment_matched_CA':len(pairs),'alignment_core_CA':int(mask.sum()),'alignment_core_RMSD_A':float(np.sqrt(np.mean(dist[mask]**2))),'pocket_residue_plddt':residues,'pocket_mean_plddt':float(np.mean([v[2] for v in residues])),'caveat':'Predicted rigid receptor; no experimental FLT4 kinase structure or self-redocking validation available; native pLDDT is confidence, not temperature factor.'}
    (P/'receptor_manifest.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

if __name__=='__main__':main()
