import csv
import json
from pathlib import Path
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
compounds = json.loads((ROOT / 'candidate_library/candidate_structures_final16.json').read_text())

def key(s):
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        raise ValueError('Invalid molecular structure')
    return Chem.MolToInchiKey(mol)

by_key = {key(c['canonical_smiles']): c for c in compounds}
assert len(by_key) == 16
sources = [
    ROOT.parent / 'binding_analysis/Inductive_Bio_Raw_Predictions.json',
    ROOT / 'screening/Inductive_Raw_New_Candidates.json',
]
rows = []
identity = []
for source in sources:
    raw = json.loads(source.read_text())
    for block in raw['content']:
        if block['type'] != 'text':
            continue
        model = json.loads(block['text'])
        for pred in model['predictions']:
            returned = pred['smiles']
            ik = key(returned)
            c = by_key[ik]
            row = {
                'compound': c['name'], 'model_id': model['model_id'],
                'model_version': model['model_version'], 'config_version': model['config_version'],
                'value': pred.get('continuous_prediction'),
                'model_lower_bound': pred.get('continuous_prediction_low'),
                'model_upper_bound': pred.get('continuous_prediction_high'),
                'bounds_coverage': 'Not specified by provider; not treated as confidence intervals',
                'category': pred.get('categorical_prediction'),
                'out_of_domain_flag': pred.get('out_of_domain_flag'),
                'low_confidence_flag': pred.get('low_confidence_flag'),
                'error': pred.get('error'), 'status': pred.get('status'),
                'input_smiles': c['canonical_smiles'], 'returned_smiles': returned,
                'verified_inchikey': ik, 'identity_match': True,
                'source_file': str(source.relative_to(ROOT.parent)),
            }
            rows.append(row)
            identity.append({k: row[k] for k in ['compound','model_id','input_smiles','returned_smiles','verified_inchikey','identity_match']})

assert len(rows) == 48
assert len({(r['compound'], r['model_id']) for r in rows}) == 48
refs = {r['model_id']: r for r in rows if r['compound'] == 'nintedanib'}
for r in rows:
    ref = refs[r['model_id']]
    assert (r['model_version'],r['config_version']) == (ref['model_version'],ref['config_version'])
    r['delta_vs_nintedanib'] = r['value'] - ref['value']

out = ROOT / 'screening'
with (out / 'Inductive_Properties_All16.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out / 'Inductive_Structure_Identity_Verification.json').write_text(json.dumps(identity,indent=2))
summary={
    'n_compounds':16,'n_predictions':48,
    'successful':sum(r['status']=='success' and not r['error'] for r in rows),
    'out_of_domain':sum(bool(r['out_of_domain_flag']) for r in rows),
    'low_confidence':sum(bool(r['low_confidence_flag']) for r in rows),
    'model_versions':{k:v['model_version'] for k,v in refs.items()},
    'all_structures_verified_by_stereospecific_inchikey':True,
    'property_scope':'Predicted logD at pH 7.4 and most acidic/basic pKa only; not clinical ADMET or hepatic safety.',
    'bounds_note':rows[0]['bounds_coverage'],
    'nintedanib_predictions':'Reused prior raw predictions with identical model and configuration versions.',
    'nerandomilast_identity_note':'Provider canonicalized the SMILES; stereospecific InChIKey matches the input R-sulfoxide.'
}
(out/'Inductive_Summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
for r in rows:
    if r['model_id']=='mcp_public_logd':
        print(r['compound'],round(r['value'],3),round(r['delta_vs_nintedanib'],3),r['out_of_domain_flag'])
