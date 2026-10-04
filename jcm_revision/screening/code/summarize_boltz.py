import csv,json,hashlib
from pathlib import Path
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[2]
compounds=json.loads((ROOT/'candidate_library/candidate_structures_final16.json').read_text())
by_id={c['name'].lower():c for c in compounds}
rows=[]; manifests=[]
for target in ['PDGFRA','FLT4']:
    run_dir=ROOT/f'boltz-experiments/jcm-{target.lower()}-focused16-v1'
    run=json.loads((run_dir/'run.json').read_text())
    assert run['status']=='succeeded', f'{target} has not completed'
    results=[json.loads(l) for l in (run_dir/'results/index.jsonl').read_text().splitlines() if l]
    assert len({r['external_id'] for r in results})==len(results)
    by_result={r['external_id']:r for r in results}
    assert len(results)==run['progress']['num_molecules_screened']
    assert len(results)+run['progress']['rejection_summary']['filtered_count']==len(compounds)
    assert run['progress']['num_molecules_failed']==0
    order={r['external_id']:i+1 for i,r in enumerate(sorted(results,key=lambda r:-r['metrics']['binding_confidence']))}
    for cid,c in by_id.items():
        r=by_result.get(cid)
        row={'compound':c['name'],'target':target,'status':'scored' if r else 'filtered',
             'rank_within_target':order.get(cid),'binding_confidence':None,'optimization_score':None,
             'structure_confidence':None,'complex_plddt':None,'complex_iplddt':None,'iptm':None,'ptm':None,
             'boltz_logD':None,'boltz_permeability_score':None,'boltz_solubility_category':None,
             'filter_reason':'' if r else 'Default provider structural-alert filter; individual rule not returned',
             'input_smiles':c['canonical_smiles'],'returned_smiles':r['smiles'] if r else '',
             'input_inchikey':Chem.MolToInchiKey(Chem.MolFromSmiles(c['canonical_smiles'])),
             'returned_smiles_identity_verified':None,'engine':run['engine'],'engine_version':run['engine_version'],
             'pipeline':run['pipeline'],'pipeline_version':run['pipeline_version'],
             'job_id':run['id'],'result_id':r['id'] if r else '',
             'result_structure_path':'','structure_sha256':''}
        if r:
            assert Chem.MolToInchiKey(Chem.MolFromSmiles(r['smiles']))==row['input_inchikey']
            row['returned_smiles_identity_verified']=True
            for k,v in r['metrics'].items():
                if k in row:row[k]=v
            for k,out in [('lipophilicity','boltz_logD'),('permeability','boltz_permeability_score'),('solubility','boltz_solubility_category')]:
                row[out]=r.get('adme',{}).get(k)
            structure=run_dir/r['paths']['structure'];assert structure.is_file() and structure.stat().st_size>0
            row['result_structure_path']=str(structure.relative_to(ROOT))
            row['structure_sha256']=hashlib.file_digest(structure.open('rb'),'sha256').hexdigest()
        rows.append(row)
    manifests.append({'target':target,'job_id':run['id'],'status':run['status'],'input_n':len(compounds),
        'scored_n':len(results),'filtered_n':len(compounds)-len(results),
        'filtered_compounds':[r['compound'] for r in rows if r['target']==target and r['status']=='filtered'],
        'top5_binding_confidence':[{'compound':by_id[r['external_id']]['name'],'score':r['metrics']['binding_confidence']} for r in sorted(results,key=lambda r:-r['metrics']['binding_confidence'])[:5]],
        'started_at':run['started_at'],'completed_at':run['completed_at']})
out=ROOT/'screening'
with (out/'Boltz_All32_Dispositions.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out/'Boltz_Screening_Summary.json').write_text(json.dumps({'runs':manifests,'scoring_scope':'Target-specific AI model outputs; not biochemical affinities or clinical probabilities','filtered_interpretation':'No score generated; filtering does not demonstrate absence of binding','software':'boltz-api 0.43.0'},indent=2))
print(json.dumps(manifests,indent=2))
