#!/usr/bin/env python3
"""Exact-coordinate Figure 7 molecular views; no ligand or receptor fitting.

Uses orthographic camera rotation for display only. SDF bond orders are retained
and aromatic bonds are drawn in a Kekulé representation. Receptor backdrop is
an interpolated C-alpha trace, not a secondary-structure assignment. All displayed
distances are calculated in 3-D before camera projection. No H-bond assignment.
"""
from pathlib import Path
import json, csv, hashlib, os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'.mplconfig'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import patheffects as pe
from scipy.interpolate import make_interp_spline
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2] / 'docking'
TEAL, CORAL, NAVY = '#177E89', '#D66A52', '#244F70'
INK, GRAY = '#263840', '#A9B4BB'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,
 'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5,
 'savefig.facecolor':'white'})
SOURCES = set()

def read_json(rel):
    p=ROOT/rel; SOURCES.add(p); return json.loads(p.read_text())

def ligand(rel, index=0):
    p=ROOT/rel; SOURCES.add(p)
    m=list(Chem.SDMolSupplier(str(p),removeHs=False))[index]
    assert m is not None
    m=Chem.RemoveHs(m)
    return m

def protein(rel):
    p=ROOT/rel; SOURCES.add(p)
    return Chem.MolFromPDBFile(str(p),removeHs=True,sanitize=False)

def xyz(m): return m.GetConformer().GetPositions()

def residue_ids(m,n):
    return [a.GetIdx() for a in m.GetAtoms() if a.GetPDBResidueInfo().GetResidueNumber()==n]

def atom_info(m,i):
    a=m.GetAtomWithIdx(int(i)); r=a.GetPDBResidueInfo()
    return {'chain':r.GetChainId(),'residue':r.GetResidueName().strip(),
      'residue_number':r.GetResidueNumber(),'atom_name':r.GetName().strip(),'index':int(i)}

def min_contact(prot,lig,n):
    ids=residue_ids(prot,n); dist=np.linalg.norm(xyz(prot)[ids,None,:]-xyz(lig)[None,:,:],axis=2)
    i,j=np.unravel_index(dist.argmin(),dist.shape)
    return {'protein':atom_info(prot,ids[i]),'ligand_atom_index_zero_based':int(j),
      'ligand_element':lig.GetAtomWithIdx(int(j)).GetSymbol(),'distance_A':float(dist[i,j])}

def camera(points, roll=0, flip=False):
    center=np.mean(points,axis=0)
    _,_,vh=np.linalg.svd(points-center,full_matrices=False)
    # Put longest direction on vertical axis, then apply one shared camera roll.
    rot=np.array([vh[1],vh[0],vh[2]]).T
    if np.linalg.det(rot)<0: rot[:,2]*=-1
    if flip: rot[:,:2]*=-1
    t=np.deg2rad(roll)
    spin=np.array([[np.cos(t),np.sin(t),0],[-np.sin(t),np.cos(t),0],[0,0,1]])
    return center, rot@spin

def projection(p,cam): return (np.asarray(p)-cam[0])@cam[1]

def figure(area=(.015,.19,.97,.65)):
    f=plt.figure(figsize=(2.25,2.4),dpi=600,facecolor='white')
    a=f.add_axes(area); a.set_aspect('equal'); a.axis('off')
    return f,a

def bounds(ax,points,cam,pad=2.1):
    pp=projection(points,cam); mn=pp[:,:2].min(0)-pad; mx=pp[:,:2].max(0)+pad
    ax.set_xlim(mn[0],mx[0]);ax.set_ylim(mn[1],mx[1])

def backbone(ax,prot,cam,focus,cutoff=9):
    runs=[];run=[]; last=None
    for a in prot.GetAtoms():
        r=a.GetPDBResidueInfo()
        if r.GetName().strip()!='CA': continue
        pos=xyz(prot)[a.GetIdx()]
        keep=np.min(np.linalg.norm(focus-pos,axis=1))<cutoff
        if not keep or (last is not None and np.linalg.norm(pos-last)>4.5):
            if len(run)>1:runs.append(np.array(run))
            run=[]
        if keep:run.append(pos)
        last=pos
    if len(run)>1:runs.append(np.array(run))
    for run in runs:
        p=projection(run,cam)
        if len(p)>3:
            t=np.arange(len(p));p=make_interp_spline(t,p,k=3)(np.linspace(0,len(p)-1,len(p)*12))
        ax.plot(p[:,0],p[:,1],color='#E8ECEF',lw=3.3,solid_capstyle='round',zorder=0)
        ax.plot(p[:,0],p[:,1],color='#D5DDE2',lw=.55,solid_capstyle='round',zorder=.1)

def draw_mol(ax,m,cam,color,atom_colors=False,selection=None,lw=1.0,alpha=1,zbase=3,atoms=True):
    m=Chem.Mol(m)
    if any(b.GetIsAromatic() for b in m.GetBonds()):Chem.Kekulize(m,clearAromaticFlags=True)
    p=projection(xyz(m),cam)
    selected=set(range(m.GetNumAtoms())) if selection is None else set(selection)
    def ac(i):
        return {'N':NAVY,'O':CORAL,'S':'#B69A42'}.get(m.GetAtomWithIdx(i).GetSymbol(),color) if atom_colors else color
    # Render each true bond, then atoms, in approximate camera depth order.
    for b in sorted(m.GetBonds(),key=lambda b: (p[b.GetBeginAtomIdx(),2]+p[b.GetEndAtomIdx(),2])/2):
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        if i not in selected or j not in selected: continue
        a,c=p[i,:2],p[j,:2]; vec=c-a; d=np.linalg.norm(vec)
        if d<1e-6: continue
        normal=np.array([-vec[1],vec[0]])/d
        order=b.GetBondTypeAsDouble()
        offsets=[0] if order<1.9 else ([-.15,.15] if order<2.9 else [-.22,0,.22])
        zz=zbase+(p[i,2]+p[j,2])*.001
        for off in offsets:
            aa=a+normal*off;cc=c+normal*off;mid=(aa+cc)/2
            for p1,p2,col in [(aa,mid,ac(i)),(mid,cc,ac(j))]:
                ax.plot([p1[0],p2[0]],[p1[1],p2[1]],color=col,lw=lw if order<1.9 else lw*.70,
                        alpha=alpha,solid_capstyle='round',zorder=zz)
    if atoms:
        for i in sorted(selected,key=lambda i:p[i,2]):
            ax.scatter(p[i,0],p[i,1],s=(2.3 if m.GetAtomWithIdx(i).GetSymbol()=='C' else 3.5),
                       c=ac(i),edgecolors='none',alpha=alpha,zorder=zbase+.1+p[i,2]*.001)

def label(ax,text,pos,offset,ha='left',color=INK,size=7):
    return ax.annotate(text,xy=pos,xytext=offset,textcoords='offset points',ha=ha,va='center',
      fontsize=size,color=color,zorder=30,
      bbox={'facecolor':'white','edgecolor':'none','alpha':.92,'pad':.8},
      arrowprops={'arrowstyle':'-','color':'#909CA4','lw':.5,'shrinkA':2,'shrinkB':2})

def contact_line(ax,prot,lig,c,cam,color=GRAY):
    a=xyz(prot)[c['protein']['index']];b=xyz(lig)[c['ligand_atom_index_zero_based']]
    p=projection([a,b],cam)
    ax.plot(p[:,0],p[:,1],color=color,lw=.8,ls=(0,(2,2)),zorder=6)
    return p

def legend(f,items,y=.94):
    handles=[Line2D([0],[0],color=color,lw=1.7,marker='o',markersize=2.0,label=name) for name,color in items]
    f.legend(handles=handles,loc='center',bbox_to_anchor=(.5,y),frameon=False,
       ncol=len(items),fontsize=7.1,handlelength=1.35,handletextpad=.45,columnspacing=.9)

def save(f,name):
    for ext in ['png','pdf','svg']:f.savefig(HERE/f'{name}.{ext}',dpi=600)
    plt.close(f)
    im=Image.open(HERE/f'{name}.png'); assert im.size==(1350,1440),im.size

def main():
    out={};pdg=protein('PDGFRA_6JOL_receptor.pdb');flt=protein('FLT4_AF_receptor.pdb')
    ref=ligand('6JOL_STI_reference.sdf');red=ligand('runs/PDGFRA_Imatinib_redocking_104729/poses.sdf')
    rs=read_json('runs/PDGFRA_Imatinib_redocking_104729/result.json')
    # CalcRMS computes same-frame symmetry-aware RMSD without alignment.
    rmsd=rdMolAlign.CalcRMS(Chem.Mol(red),Chem.Mol(ref),maxMatches=100000)
    assert abs(rmsd-rs['reference_RMSD_A_receptor_frame'][0])<1e-6
    c=camera(xyz(ref),roll=24,flip=False)
    f,a=figure((.02,.19,.96,.66));focus=np.vstack([xyz(ref),xyz(red)])
    bounds(a,focus,c,2.8);backbone(a,pdg,c,focus,8)
    draw_mol(a,ref,c,TEAL,lw=1.04,zbase=3)
    draw_mol(a,red,c,CORAL,lw=.84,zbase=4)
    legend(f,[('Crystal (6JOL)',TEAL),('Redocked',CORAL)])
    f.text(.5,.128,f'Heavy-atom RMSD {rmsd:.3f} Å',ha='center',fontsize=8,color=INK)
    f.text(.5,.064,'Receptor frame · no ligand alignment',ha='center',fontsize=6.8,color='#617078')
    save(f,'Panel_E_redocking')
    out['E']={'crystal':'6JOL_STI_reference.sdf','docked':'runs/PDGFRA_Imatinib_redocking_104729/poses.sdf','pose':1,'seed':104729,
      'rmsd_A_recorded_result':rs['reference_RMSD_A_receptor_frame'][0],
      'rmsd_A_recalculated_from_exported_sdf':rmsd,'camera_center':c[0].tolist(),'camera_rotation':c[1].tolist(),
      'RMSD_method':'symmetry-aware heavy-atom CalcRMS; unchanged receptor frame, no ligand alignment'}

    p1=ligand('runs/PDGFRA_nintedanib_104729/poses.sdf',0)
    p3=ligand('runs/PDGFRA_nintedanib_104729/poses.sdf',2)
    scores=read_json('runs/PDGFRA_nintedanib_104729/result.json')['pose_energies']
    qcs=read_json('PDGFRA_nintedanib_pose_hinge_QC.json')
    cs=[min_contact(pdg,m,674) for m in [p1,p3]]
    assert abs(cs[0]['distance_A']-qcs[0]['minimum_gatekeeper_Thr674_distance_A'])<1e-6
    assert abs(cs[1]['distance_A']-qcs[2]['minimum_gatekeeper_Thr674_distance_A'])<1e-6
    focus=np.vstack([xyz(p1),xyz(p3),xyz(pdg)[residue_ids(pdg,674)]])
    c=camera(focus,roll=-18,flip=True)
    f,a=figure((.01,.19,.98,.62));bounds(a,focus,c,2.8);backbone(a,pdg,c,focus,7)
    draw_mol(a,pdg,c,'#6F808D',selection=residue_ids(pdg,674),lw=.75,zbase=2)
    draw_mol(a,p1,c,TEAL,lw=.93,zbase=3)
    draw_mol(a,p3,c,CORAL,lw=.93,zbase=4)
    for m,co,cc in [(p1,TEAL,cs[0]),(p3,CORAL,cs[1])]:contact_line(a,pdg,m,cc,c,co)
    pp=projection(xyz(pdg)[residue_ids(pdg,674)].mean(0),c)
    label(a,'Thr674',pp[:2],(7,12),size=7)
    legend(f,[('Pose 1',TEAL),('Pose 3',CORAL)],.955)
    f.text(.5,.865,'−8.502 / −8.341 kcal mol⁻¹',ha='center',fontsize=7.2,color=INK)
    f.text(.5,.136,'Thr674 minimum proximity',ha='center',fontsize=7.2,color=INK)
    f.text(.5,.075,'8.43 Å / 3.38 Å',ha='center',fontsize=8,color=INK)
    f.text(.5,.022,'Predicted poses · seed 104729',ha='center',fontsize=7,color='#617078')
    save(f,'Panel_F_PDGFRA_pose_ambiguity')
    out['F']={'source':'runs/PDGFRA_nintedanib_104729/poses.sdf','poses':[1,3],'seed':104729,
      'scores_kcal_mol':[scores[0][0],scores[2][0]],'contacts_to_Thr674':cs,
      'camera_center':c[0].tolist(),'camera_rotation':c[1].tolist(),
      'interpretation':'Two predicted poses retained; neither is proven. Supplied hinge-region minima are to gatekeeper Thr674, not a demonstrated hinge hydrogen bond.'}

    g=ligand('runs/FLT4_nintedanib_104729/poses.sdf',0)
    gs=read_json('runs/FLT4_nintedanib_104729/result.json')
    selected=[896,927,934,1055]
    gc=[min_contact(flt,g,n) for n in selected]
    contactpath=ROOT/'Nintedanib_Docking_Contacts.csv';SOURCES.add(contactpath)
    rows=list(csv.DictReader(contactpath.open()))
    for cc in gc:
        row=next(r for r in rows if r['target']=='FLT4' and int(r['seed'])==104729 and int(r['residue_number'])==cc['protein']['residue_number'])
        assert abs(cc['distance_A']-float(row['minimum_heavy_atom_distance_A']))<.00051
        assert cc['distance_A']<4
    # Only actual nearby selected residues, including their molecular connectivity.
    focus=np.vstack([xyz(g)]+[xyz(flt)[residue_ids(flt,n)] for n in selected])
    c=camera(xyz(g),roll=-10,flip=True)
    f,a=figure((.00,.145,1,.725));bounds(a,focus,c,2.2);backbone(a,flt,c,xyz(g),7)
    for n in selected:draw_mol(a,flt,c,'#8C9AA4',selection=residue_ids(flt,n),lw=.7,zbase=2,atoms=False)
    draw_mol(a,g,c,TEAL,atom_colors=True,lw=1.03,zbase=4)
    # Label offsets are editorial placement only; leader lines end at exact receptor atoms.
    placements={896:(.91,.79,'right'),927:(.07,.65,'left'),934:(.93,.245,'right'),1055:(.95,.49,'right')}
    for cc in gc:
        pp=contact_line(a,flt,g,cc,c,'#718492');n=cc['protein']['residue_number'];tx,ty,ha=placements[n]
        nm=cc['protein']['residue'].title()
        a.annotate(f'{nm}{n}\n{cc["distance_A"]:.2f} Å',xy=pp[0,:2],xytext=(tx,ty),
           textcoords='figure fraction',ha=ha,va='center',fontsize=7.0,color=INK,zorder=30,
           bbox={'facecolor':'white','edgecolor':'none','alpha':1,'pad':1.0},
           arrowprops={'arrowstyle':'-','color':'#909CA4','lw':.5,'shrinkA':3,'shrinkB':2})
    f.text(.5,.94,'AlphaFold receptor · predicted pose',ha='center',fontsize=7.3,color=NAVY)
    f.text(.5,.093,'Selected heavy-atom proximities <4 Å',ha='center',fontsize=6.9,color=INK)
    f.text(.5,.034,'Pose 1 · seed 104729',ha='center',fontsize=7,color='#617078')
    save(f,'Panel_G_FLT4_predicted_pocket')
    out['G']={'receptor':'FLT4_AF_receptor.pdb','receptor_type':'AlphaFold AF-P35916-F1 v6 predicted kinase domain',
      'ligand':'runs/FLT4_nintedanib_104729/poses.sdf','pose':1,'seed':104729,'score_kcal_mol':gs['best_vina_score_kcal_mol'],
      'selected_contacts':gc,'camera_center':c[0].tolist(),'camera_rotation':c[1].tolist(),
      'contact_semantics':'selected minimum heavy-atom distances <4 Å, calculated in 3-D; no hydrogen-bond classification'}
    out['global']={'size_inches':[2.25,2.4],'png_pixels':[1350,1440],'dpi':600,
      'backbone':'C-alpha trace, cubic interpolation only within contiguous local residue runs; no invented secondary structure',
      'bonds':'actual SDF bond types; aromatic bonds shown as Kekulé single/double bonds; hydrogens hidden',
      'geometry':'shared rigid orthographic camera per panel; no coordinate editing or ligand superposition',
      'palette':{'teal':TEAL,'coral':CORAL,'navy':NAVY},
      'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SOURCES)}}
    (HERE/'molecular_panel_provenance.json').write_text(json.dumps(out,indent=2)+'\n')
    imgs=[Image.open(HERE/(nm+'.png')).convert('RGB') for nm in ['Panel_E_redocking','Panel_F_PDGFRA_pose_ambiguity','Panel_G_FLT4_predicted_pocket']]
    sheet=Image.new('RGB',(4050,1440),'white')
    for i,im in enumerate(imgs):sheet.paste(im,(1350*i,0))
    sheet.resize((1620,576),Image.Resampling.LANCZOS).save(HERE/'molecular_panels_preview.png')
    print(json.dumps({'rmsd_A':rmsd,'F_contacts':cs,'G_contacts':gc,'output':str(HERE)},indent=2))

if __name__=='__main__':main()
