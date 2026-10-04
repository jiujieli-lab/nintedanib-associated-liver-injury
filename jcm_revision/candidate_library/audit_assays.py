import pathlib,json,csv,re,collections,hashlib,statistics
ROOT=pathlib.Path(__file__).parent; RAW=ROOT/'raw_chembl'
names=json.loads((ROOT/'candidate_structures_final16.json').read_text())
target={'CHEMBL2007':'PDGFRA','CHEMBL1955':'FLT4'}
records=[]; docs={}
for n in names:
 p=RAW/(n['name']+'_target_activities.json')
 j=json.loads(p.read_text())
 assert j['page_meta']['total_count']==len(j['activities']),n['name']
 for a in j['activities']:
  d=a.get('document_chembl_id'); dp=RAW/(str(d)+'.json')
  doc=json.loads(dp.read_text()) if dp.exists() else {}
  docs[d]=doc
  desc=a.get('assay_description') or '';low=desc.lower();title=doc.get('title') or ''
  secondary_docs={'CHEMBL1250538','CHEMBL2177017','CHEMBL2417410','CHEMBL4270503','CHEMBL4304778','CHEMBL5150004','CHEMBL5241055','CHEMBL6103541','CHEMBL5365465','CHEMBL5634097'}
  primary='secondary_review_or_reanalysis' if d in secondary_docs or re.search(r'\breview\b|Oxetanes in Drug Discovery Campaigns',title,re.I) else ('patent' if doc.get('doc_type')=='PATENT' else 'publication_primary_status_not_fully_manually_reviewed')
  if d in ['CHEMBL1908390','CHEMBL1150977','CHEMBL1144455','CHEMBL6153081']:primary='verified_primary_publication'
  if a.get('assay_variant_mutation') or re.search(r'mutant|mutation|exon\s*18|TEL.fused|FIP1L1',desc,re.I):wt='variant_or_fusion'
  elif re.search(r'wild.?type',desc,re.I):wt='explicit_wild_type'
  else:wt='no_variant_described_not_proven_wild_type'
  if 'extracted from Sf9 cells' in desc:context='biochemical_recombinant_protein'; format_note='ChEMBL cell-based BAO conflicts with extracted-protein HTRF description; classified biochemical by description'
  elif a.get('bao_label')=='cell-based format':context='cellular';format_note=''
  elif a.get('bao_label')=='single protein format' or re.search(r'recombinant|cytoplasmic domain|kinase domain|Enzyme Inhibition',desc,re.I):context='biochemical';format_note=''
  else:context='assay_context_uncertain';format_note=''
  reasons=[]
  if a.get('target_organism')!='Homo sapiens':reasons.append('target_not_human')
  if wt=='variant_or_fusion':reasons.append('variant_or_fusion')
  if a.get('data_validity_comment'):reasons.append('validity_flag')
  if a.get('potential_duplicate'):reasons.append('potential_duplicate')
  if a.get('standard_type') not in ['Kd','Ki','IC50']:reasons.append('endpoint_outside_Kd_Ki_IC50')
  if a.get('standard_value') is None or a.get('standard_units')!='nM':reasons.append('no_numeric_nM_value')
  if a.get('standard_relation')!='=':reasons.append('nonexact_relation')
  if primary.startswith('secondary_'):reasons.append('secondary_source')
  clean=not reasons
  row={
   'compound':n['name'],'gene_symbol':target[a['target_chembl_id']],
   **{k:a.get(k) for k in ['activity_id','molecule_chembl_id','parent_molecule_chembl_id','target_chembl_id','target_organism','target_tax_id','standard_type','standard_relation','standard_value','standard_units','pchembl_value','standard_text_value','activity_comment','assay_chembl_id','assay_description','assay_type','assay_variant_accession','assay_variant_mutation','bao_format','bao_label','data_validity_comment','potential_duplicate','document_chembl_id','document_year']},
   'wild_type_status':wt,'target_species_provenance':'explicit_human_in_description' if 'human' in low else ('unknown_origin_in_description' if 'unknown origin' in low else 'ChEMBL_human_target_mapping_only'),
   'interpreted_assay_context':context,'assay_format_note':format_note,'source_type':doc.get('doc_type'),'document_title':title,'document_doi':doc.get('doi'),'document_pubmed_id':doc.get('pubmed_id'),'document_primary_status':primary,
   'eligible_exact_human_mapped_nonvariant':clean,
   'eligible_biochemical_exact_human_mapped_nonvariant':clean and context.startswith('biochemical'),
   'eligible_explicit_human_sensitivity':clean and 'human' in low,
   'disposition':'eligible' if clean else ';'.join(reasons),
   'activity_url':'https://www.ebi.ac.uk/chembl/api/data/activity/'+str(a['activity_id'])+'.json',
   'document_url':'https://doi.org/'+doc['doi'] if doc.get('doi') else 'https://www.ebi.ac.uk/chembl/api/data/document/'+str(d)+'.json'}
  records.append(row)

def save_csv(name,rows,fields=None):
 with (ROOT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),extrasaction='ignore');w.writeheader();w.writerows(rows)

save_csv('All_Target_Activity_Disposition.csv',records)
save_csv('Eligible_Exact_Human_Mapped_Nonvariant_Records.csv',[r for r in records if r['eligible_exact_human_mapped_nonvariant']])
summary=[]
for n in names:
 for t in target.values():
  for ep in ['Kd','Ki','IC50']:
   rr=[r for r in records if r['compound']==n['name'] and r['gene_symbol']==t and r['standard_type']==ep]
   ok=[r for r in rr if r['eligible_exact_human_mapped_nonvariant']]
   bio=[r for r in ok if r['eligible_biochemical_exact_human_mapped_nonvariant']]
   summary.append({'compound':n['name'],'gene_symbol':t,'endpoint':ep,'all_records':len(rr),'eligible_exact_nonvariant_records':len(ok),'biochemical_records':len(bio),'cellular_records':sum(r['interpreted_assay_context']=='cellular' for r in ok),'explicit_human_records':sum(r['eligible_explicit_human_sensitivity'] for r in ok),'explicit_wild_type_records':sum(r['wild_type_status']=='explicit_wild_type' for r in ok),'biochemical_values_nM':';'.join(r['standard_value'] for r in bio),'biochemical_activity_ids':';'.join(str(r['activity_id']) for r in bio),'note':'Empty activity values mean no eligible record under these rules, not proven absence of binding; do not pool IC50/Kd/Ki or infer hepatoprotection.'})
save_csv('Endpoint_Separated_Experimental_Summary.csv',summary)
refs=[]
for d,doc in docs.items():
 refs.append({k:doc.get(k) for k in ['document_chembl_id','title','authors','journal','year','volume','first_page','last_page','doi','pubmed_id','doc_type']})
save_csv('Assay_Source_References.csv',refs)
# Retain the matched Davis panel as a benchmark, explicitly labeling duplicate imported rows.
matched=[r for r in records if r['document_chembl_id']=='CHEMBL1908390' and r['standard_type']=='Kd']
save_csv('Davis2011_Matched_Kd_Benchmark.csv',matched)
counts=collections.Counter(r['disposition'] for r in records)
manifest={'retrieval_date':'2026-10-02','compounds':len(names),'activity_records':len(records),'documents':len(docs),'exact_nonvariant_records':sum(r['eligible_exact_human_mapped_nonvariant'] for r in records),'biochemical_records':sum(r['eligible_biochemical_exact_human_mapped_nonvariant'] for r in records),'disposition_counts':dict(counts),'rules':['human target mapping required; exact relation and quantitative nM Kd/Ki/IC50 analyzed separately','assay variant metadata and description identify mutation/fusion records, which are excluded from nonvariant summaries','no variant described does not prove an experimentally verified wild-type construct','potential duplicates, validity flags, and clearly identified secondary-source reviews excluded from primary exact summaries','all censored, duplicate, variant, nonquantitative and nonstandard endpoint records retained in full disposition','Davis2011 benchmark preserves duplicated earlier-source values for matched-panel use; duplicates must not be counted as independent evidence','Sf9 extracted-protein HTRF descriptions override conflicting cellular BAO classification, with raw BAO retained','no synthetic potency values or artificial nonbinder labels for absent records']}
(ROOT/'Assay_Audit_Manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest,indent=2))
