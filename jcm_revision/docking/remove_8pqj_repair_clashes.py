"""Rotate only newly rebuilt, distal side-chain atoms to remove steric clashes.

No deposited atom is moved. This is a local geometry repair, not MD or an energy calculation.
"""
from pathlib import Path
import json,numpy as np
P=Path(__file__).resolve().parent
path=P/'8PQJ_fixed_heavy.pdb';lines=path.read_text().splitlines();idx=[i for i,l in enumerate(lines) if l.startswith('ATOM')]
records=[lines[i] for i in idx];xyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in records]);original=xyz.copy()
added=set(json.loads((P/'8PQJ_heavy_atom_repair_QC.json').read_text())['added_atoms'])
keys=[l[21]+':'+l[22:26].strip()+':'+l[12:16].strip() for l in records]
resnums=np.array([int(l[22:26]) for l in records]);report=[]
for cycle in range(3):
 for n in [780,785]:
  moving=np.array([i for i,k in enumerate(keys) if resnums[i]==n and k in added])
  fixed=np.array([i for i in range(len(records)) if resnums[i]!=n])
  ca=xyz[keys.index(f'A:{n}:CA')];cb=xyz[keys.index(f'A:{n}:CB')];axis=cb-ca;axis/=np.linalg.norm(axis)
  v=xyz[moving]-cb;best=None
  for degree in range(0,360,15):
   angle=np.deg2rad(degree);rot=v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+(v@axis)[:,None]*axis[None,:]*(1-np.cos(angle));candidate=rot+cb
   distance=np.linalg.norm(candidate[:,None,:]-xyz[fixed][None,:,:],axis=2)
   score=float(np.square(np.maximum(2.8-distance,0)).sum())
   if best is None or score<best[0]:best=(score,degree,candidate,float(distance.min()))
  xyz[moving]=best[2];report.append({'cycle':cycle,'residue':n,'rotation_degrees':best[1],'clash_penalty':best[0],'min_other_residue_distance_A':best[3]})
for i,j in enumerate(idx):
 l=lines[j];lines[j]=l[:30]+f'{xyz[i,0]:8.3f}{xyz[i,1]:8.3f}{xyz[i,2]:8.3f}'+l[54:]
path.write_text('\n'.join(lines)+'\n')
deposited=[i for i,k in enumerate(keys) if k not in added]
assert np.max(np.linalg.norm(xyz[deposited]-original[deposited],axis=1))==0
(P/'8PQJ_rebuilt_sidechain_clash_repair.json').write_text(json.dumps({'method':'Grid search of chi1 rotations in 15-degree increments for newly rebuilt Asn780 and Asp785 side chains only; minimize sum of squared nonbonded-distance deficits below 2.8 A; no force-field affinity interpretation','deposited_atoms_moved':0,'report':report},indent=2));print(json.dumps(report,indent=2))
