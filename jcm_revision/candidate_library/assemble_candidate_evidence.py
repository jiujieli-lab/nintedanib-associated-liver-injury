import pathlib,json,csv,hashlib
ROOT=pathlib.Path(__file__).parent
s=json.loads((ROOT/'candidate_structures_final16.json').read_text())
labels=json.loads((ROOT/'label_review.json').read_text())
by={r['compound'].lower():r for r in labels['compounds']}
acts=list(csv.DictReader((ROOT/'All_Target_Activity_Disposition.csv').open()))
pdgfra={'imatinib','crenolanib','avapritinib','ripretinib'}
out=[]
for x in s:
 n=x['name'];l=by[n.lower()]
 row=dict(x);row['clinical_evidence']=l
 if n=='nintedanib':role='Index drug and biochemical positive reference'
 elif n in pdgfra:role='PDGFRA-directed mechanistic comparator'
 elif n=='MAZ51':role='Nonclinical FLT4 inhibition comparator; off-target caution'
 elif n in ['pirfenidone','nerandomilast']:role='Pulmonary clinical nonkinase reference; not validated receptor nonbinder'
 else:role='VEGFR/multikinase mechanistic comparator'
 row['candidate_role']=role
 row['hepatoprotection_established']=False
 row['clinical_substitute_recommendation']='Not established by this screen; no clinical switching recommendation'
 for gene in ['PDGFRA','FLT4']:
  a=[r for r in acts if r['compound']==n and r['gene_symbol']==gene]
  row[gene+'_evidence']={'raw_records':len(a),'eligible_exact_human_mapped_nonvariant':sum(r['eligible_exact_human_mapped_nonvariant']=='True' for r in a),'note':'No variant described is not proof of wild-type construct. No record does not prove no binding.'}
 out.append(row)
(ROOT/'Candidate_Library_Evidence.json').write_text(json.dumps(out,indent=2))
flat=[]
for x in out:
 l=x['clinical_evidence']
 flat.append({'compound':x['name'],'candidate_role':x['candidate_role'],'chembl_id':x['chembl_id'],'selected_smiles':x['canonical_smiles'],'selected_inchi_key':x['standard_inchi_key'],'structure_url':x['structure_url'],'identity_resolution':x.get('identity_resolution','ChEMBL neutral parent matched exact PubChem InChIKey'),'regulatory_status':l['regulatory_status'],'indications':l['approved_indications_summary'],'target_mechanism':l['target_mechanism'],'hepatic_safety_summary':l['hepatic_safety_summary'],'hepatic_warning_category':l.get('hepatic_warning_category',''),'label_url':l.get('label_url',''),'PDGFRA_raw_records':x['PDGFRA_evidence']['raw_records'],'PDGFRA_exact_nonvariant_records':x['PDGFRA_evidence']['eligible_exact_human_mapped_nonvariant'],'FLT4_raw_records':x['FLT4_evidence']['raw_records'],'FLT4_exact_nonvariant_records':x['FLT4_evidence']['eligible_exact_human_mapped_nonvariant'],'hepatoprotection_established':False,'notes':'Target engagement is mechanistic evidence; docking or MD cannot establish hepatoprotection. Clinical reference compounds are not validated nonbinders.'})
with (ROOT/'Candidate_Library_Evidence.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
# Replace the early structure CSV so downstream readers cannot consume stale nerandomilast stereo.
with (ROOT/'candidate_structures.csv').open('w',newline='') as f:
 fields=['name','chembl_id','canonical_smiles','standard_inchi_key','structure_url','identity_resolution'];w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(s)
hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='Candidate_Library_SHA256.json' and '__pycache__' not in str(p)}
(ROOT/'Candidate_Library_SHA256.json').write_text(json.dumps(hashes,indent=2))
print('Assembled',len(out),'compounds;',len(hashes),'files hashed')
