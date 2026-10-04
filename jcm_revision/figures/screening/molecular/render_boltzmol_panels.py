#!/usr/bin/env python3
"""Render the supplied pocket-cropped BoltzMol mmCIF coordinates exactly.

No predictions are run, no coordinates are fitted, and no structure is completed.
Ligand bonds come from mmCIF chem_comp_bond; identities are verified to metadata
and the submitted library. Canonical residue numbers are verified against the
archived full-length UniProt-indexed AlphaFold files and submitted target input.
"""
from pathlib import Path
import csv,json,hashlib
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.interpolate import make_interp_spline
from PIL import Image
import render_molecular_panels as r

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[2]
SPECS=[
 {'target':'PDGFRA','ligand':'imatinib','pres':'pres_El5PkTZ0qk9YvZ8V8XIZ','offset':592,
  'canonical_pdb':'AF-P16234-F1-model_v6.pdb','crop_count':199,'full_count':362,
  'selected':[644,677,836],'roll':20,'flip':False,'asset':'Panel_E_BoltzMol_PDGFRA_imatinib'},
 {'target':'FLT4','ligand':'tivozanib','pres':'pres_D49fhaQle5QZqdqqr6um','offset':844,
  'canonical_pdb':'AF-P35916-F1-model_v6.pdb','crop_count':195,'full_count':329,
  'selected':[896,930,1055],'roll':-20,'flip':True,'asset':'Panel_F_BoltzMol_FLT4_tivozanib'}]

def parse_complex(spec):
    run_dir=BASE/'boltz-experiments'/f'jcm-{spec["target"].lower()}-focused16-v1'
    result=run_dir/'results'/spec['pres']
    path=result/'files/result'/f'{spec["pres"]}_predicted.cif'
    d=MMCIF2Dict(str(path));metadata=json.loads((result/'metadata.json').read_text())
    run=json.loads((run_dir/'run.json').read_text())
    submitted=run['input']['target']['entities'][0]['value']
    assert len(submitted)==spec['full_count']
    canonical={}
    canonical_path=BASE/'proto'/spec['canonical_pdb']
    for line in canonical_path.read_text().splitlines():
        if line.startswith('ATOM') and line[12:16].strip()=='CA':canonical[int(line[22:26])]=line[17:20]
    assert all(seq1(canonical[i+spec['offset']+1])==aa for i,aa in enumerate(submitted))
    assert ''.join(seq1(a) for a in d['_entity_poly_seq.mon_id'])==submitted
    n=len(d['_atom_site.id']);pos=np.array([[float(d[f'_atom_site.Cartn_{k}'][i]) for k in 'xyz'] for i in range(n)])
    ligidx=[i for i in range(n) if d['_atom_site.label_comp_id'][i]=='LIG1']
    pdb=[];protein_positions=[];res=set();ligcoord={d['_atom_site.label_atom_id'][i]:pos[i] for i in ligidx}
    assert set(d['_atom_site.B_iso_or_equiv'])=={'100.000'}
    for i in range(n):
        if d['_atom_site.group_PDB'][i]!='ATOM':continue
        local=int(d['_atom_site.label_seq_id'][i]);rn=local+spec['offset'];name=d['_atom_site.label_comp_id'][i]
        assert canonical[rn]==name and seq1(name)==submitted[local-1]
        res.add(rn);atom=d['_atom_site.label_atom_id'][i];elem=d['_atom_site.type_symbol'][i]
        x,y,z=pos[i];serial=len(pdb)+1
        atomfield=(' '+atom.ljust(3)) if len(elem)==1 else atom.ljust(4)
        pdb.append(f'ATOM  {serial:5d} {atomfield} {name:>3} A{rn:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00100.00          {elem:>2}  ')
        protein_positions.append(pos[i])
    assert len(res)==spec['crop_count']
    prot=Chem.MolFromPDBBlock('\n'.join(pdb)+'\nEND\n',sanitize=False,removeHs=True)
    # The PDB adapter supplies only residue/bond metadata; restore full CIF precision.
    assert prot.GetNumAtoms()==len(protein_positions)
    for i,position in enumerate(protein_positions):prot.GetConformer().SetAtomPosition(i,position)
    # Build the ligand explicitly from the deposited chemical component tables.
    rw=Chem.RWMol();idx={};atom_names=[]
    for i,name in enumerate(d['_chem_comp_atom.atom_id']):
        if d['_chem_comp_atom.comp_id'][i]!='LIG1':continue
        if d['_chem_comp_atom.type_symbol'][i]=='H':continue
        atom=Chem.Atom(d['_chem_comp_atom.type_symbol'][i].title())
        charge=d['_chem_comp_atom.charge'][i];atom.SetFormalCharge(0 if charge in ['.','?'] else int(charge))
        atom.SetIsAromatic(d['_chem_comp_atom.pdbx_aromatic_flag'][i]=='Y')
        atom.SetProp('cif_atom_name',name);idx[name]=rw.AddAtom(atom);atom_names.append(name)
    types={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for i,name in enumerate(d['_chem_comp_bond.atom_id_1']):
        if d['_chem_comp_bond.comp_id'][i]!='LIG1':continue
        name2=d['_chem_comp_bond.atom_id_2'][i]
        if name in idx and name2 in idx:rw.AddBond(idx[name],idx[name2],types[d['_chem_comp_bond.value_order'][i]])
    lig=rw.GetMol();Chem.SanitizeMol(lig)
    conf=Chem.Conformer(lig.GetNumAtoms())
    for name,i in idx.items():conf.SetAtomPosition(i,ligcoord[name])
    lig.AddConformer(conf)
    expected=Chem.MolFromSmiles(metadata['smiles']);assert Chem.MolToSmiles(lig)==Chem.MolToSmiles(expected)
    library=BASE/'candidate_library/Candidate_Screening_Library.csv'
    row=next(row for row in csv.DictReader(library.open()) if row['name']==spec['ligand'])
    assert Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles']))==Chem.MolToSmiles(expected)
    assert metadata['external_id']==spec['ligand']
    contacts=[]
    for rn in sorted(res):
        c=r.min_contact(prot,lig,rn)
        c['ligand_atom_name']=atom_names[c['ligand_atom_index_zero_based']]
        if c['distance_A']<4:contacts.append(c)
    independent=BASE/'md/Boltz_Priority_Pair_Contacts.csv'
    independent_rows=[rr for rr in csv.DictReader(independent.open()) if rr['target']==spec['target'] and rr['compound']==spec['ligand']]
    assert len(independent_rows)==len(contacts)
    for c in contacts:
        rr=next(rr for rr in independent_rows if int(rr['native_residue'])==c['protein']['residue_number'])
        assert rr['protein_atom']==c['protein']['atom_name'] and rr['ligand_atom']==c['ligand_atom_name']
        assert abs(float(rr['minimum_distance_A'])-c['distance_A'])<.000051
    source_paths=[path,result/'metadata.json',run_dir/'run.json',canonical_path,library,independent]
    qc={'target':spec['target'],'ligand':spec['ligand'],'prediction_id':spec['pres'],
        'ligand_smiles':Chem.MolToSmiles(lig),'ligand_formula':rdMolDescriptors.CalcMolFormula(lig),
        'ligand_heavy_atoms':lig.GetNumAtoms(),'ligand_bonds':lig.GetNumBonds(),
        'residue_number_offset':spec['offset'],'protein_coordinate_residues':len(res),'submitted_target_residues':len(submitted),
        'all_coordinate_residues_match_submitted_and_canonical_sequence':True,
        'all_contacts_match_independent_QC_atom_pairs_and_distances':True,
        'B_factor_values':[100.0],'B_factor_used_as_confidence':False,'all_contacts_lt4A':contacts,
        'input_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}}
    return prot,lig,qc

def cropped_backbone(ax,prot,cam,focus,cutoff=9):
    runs=[];run=[];last=None;lastnum=None
    for atom in prot.GetAtoms():
        info=atom.GetPDBResidueInfo()
        if info.GetName().strip()!='CA':continue
        pos=r.xyz(prot)[atom.GetIdx()];rn=info.GetResidueNumber()
        keep=np.min(np.linalg.norm(focus-pos,axis=1))<cutoff
        broken=last is not None and (rn!=lastnum+1 or np.linalg.norm(pos-last)>4.5)
        if not keep or broken:
            if len(run)>1:runs.append(np.array(run))
            run=[]
        if keep:run.append(pos)
        last,lastnum=pos,rn
    if len(run)>1:runs.append(np.array(run))
    for run in runs:
        p=r.projection(run,cam)
        if len(p)>3:p=make_interp_spline(np.arange(len(p)),p,k=3)(np.linspace(0,len(p)-1,len(p)*12))
        ax.plot(p[:,0],p[:,1],color='#E8ECEF',lw=3.2,solid_capstyle='round',zorder=0)
        ax.plot(p[:,0],p[:,1],color='#D5DDE2',lw=.55,solid_capstyle='round',zorder=.1)
    return len(runs)

def render(spec,prot,lig,qc):
    cs=[next(c for c in qc['all_contacts_lt4A'] if c['protein']['residue_number']==n) for n in spec['selected']]
    focus=np.vstack([r.xyz(lig)]+[r.xyz(prot)[r.residue_ids(prot,n)] for n in spec['selected']])
    cam=r.camera(r.xyz(lig),roll=spec['roll'],flip=spec['flip'])
    f,a=r.figure((.01,.145,.98,.70));r.bounds(a,focus,cam,2.5)
    qc['backbone_trace_segments']=cropped_backbone(a,prot,cam,r.xyz(lig),8)
    for n in spec['selected']:r.draw_mol(a,prot,cam,'#8597A3',selection=r.residue_ids(prot,n),lw=.74,zbase=2,atoms=False)
    r.draw_mol(a,lig,cam,r.TEAL,atom_colors=True,lw=1.02,zbase=4)
    # Chlorine is named explicitly; no extra category is inferred from atom color.
    for atom in lig.GetAtoms():
        if atom.GetSymbol()=='Cl':
            p=r.projection(r.xyz(lig)[atom.GetIdx()],cam)
            a.annotate('Cl',p[:2],xytext=(-6,4),ha='right',textcoords='offset points',fontsize=6.8,color=r.TEAL,zorder=8)
    if spec['target']=='PDGFRA':
        placements={644:(.91,.70,'right'),677:(.11,.26,'left'),836:(.94,.45,'right')}
    else:
        placements={896:(.06,.36,'left'),930:(.94,.72,'right'),1055:(.95,.41,'right')}
    for c in cs:
        p=r.contact_line(a,prot,lig,c,cam,'#758C99');n=c['protein']['residue_number'];tx,ty,ha=placements[n]
        a.annotate(f'{c["protein"]["residue"].title()}{n}\n{c["distance_A"]:.3f} Å',xy=p[0,:2],xytext=(tx,ty),
           textcoords='figure fraction',ha=ha,va='center',fontsize=7.0,color=r.INK,zorder=30,
           bbox={'facecolor':'white','edgecolor':'none','alpha':1,'pad':1.0},
           arrowprops={'arrowstyle':'-','color':'#909CA4','lw':.5,'shrinkA':3,'shrinkB':2})
    f.text(.5,.94,'BoltzMol predicted pocket',ha='center',fontsize=7.5,color=r.NAVY)
    f.text(.5,.09,'Pocket crop · selected proximities',ha='center',fontsize=7,color=r.INK)
    f.text(.5,.031,'Dashed: heavy-atom distance <4 Å',ha='center',fontsize=6.9,color='#617078')
    r.save(f,spec['asset'])
    qc['selected_contacts']=cs;qc['camera_center']=cam[0].tolist();qc['camera_rotation']=cam[1].tolist()
    qc['geometry']='Actual mmCIF coordinates; one shared rigid orthographic camera; no ligand fit or coordinate editing'
    qc['semantics']='Predicted cropped pocket; selected heavy-atom proximity contacts only; no hydrogen-bond assignment'
    qc['backbone']='Uniform gray C-alpha trace; cubic interpolation only within consecutive residue runs; no gaps completed'
    qc['figure']={'inches':[2.25,2.4],'dpi':600,'pixels':[1350,1440],'panel_letters':False}
    return qc

def main():
    all_qc=[]
    for spec in SPECS:
        prot,lig,qc=parse_complex(spec);all_qc.append(render(spec,prot,lig,qc))
    (HERE/'BoltzMol_molecular_provenance.json').write_text(json.dumps(all_qc,indent=2)+'\n')
    preview=Image.new('RGB',(2700,1440),'white')
    for i,spec in enumerate(SPECS):preview.paste(Image.open(HERE/(spec['asset']+'.png')).convert('RGB'),(1350*i,0))
    preview.resize((1080,576),Image.Resampling.LANCZOS).save(HERE/'BoltzMol_molecular_preview.png')
    print(json.dumps([{k:q[k] for k in ['target','ligand','ligand_formula','protein_coordinate_residues','selected_contacts']} for q in all_qc],indent=2))

if __name__=='__main__':main()
