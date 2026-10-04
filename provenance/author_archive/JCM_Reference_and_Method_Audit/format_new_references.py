from pathlib import Path
import json,re,html,hashlib
B=Path(__file__).parent;M=B/'new_reference_metadata'
def clean(s):return html.unescape(re.sub('<[^>]+>','',s)).strip()
def initials(s):
 out=[]
 for w in re.findall(r'[\wÀ-ž]+(?:-[\wÀ-ž]+)*',s,re.U):
  if '-' in w:out.append('-'.join(x[0].upper()+'.' for x in w.split('-')))
  elif w.isupper():out.append(''.join(x+'.' for x in w))
  else:out.append(w[0].upper()+'.')
 return ''.join(out)
def author(a):return clean(a['family'])+', '+initials(a.get('given','')) if a.get('family') else clean(a.get('name',''))
J={'Journal of Chemical Information and Modeling':'J. Chem. Inf. Model.','Nucleic Acids Research':'Nucleic Acids Res.','Nature':'Nature','The American Journal of Pathology':'Am. J. Pathol.','Angiogenesis':'Angiogenesis','Nature Biotechnology':'Nat. Biotechnol.'}
R=[];V=[]
for n in range(98,108):
 notes=[]
 if n not in [100,104]:
  path=M/f'new_ref_{n}.json';m=json.loads(path.read_text())['message'];a=m['author'];authors='; '.join(author(x) for x in a[:10])+('; et al.' if len(a)>10 else '')
  title=clean(m['title'][0]);journal=J.get(m['container-title'][0],m['container-title'][0]);year=str((m.get('published-print') or m.get('published') or m.get('issued'))['date-parts'][0][0]);vol=m.get('volume','');page=m.get('page',m.get('article-number','')).replace('-','–');doi=m['DOI']
  seg=[{'text':authors+' '+title.rstrip('.')+'. '},{'text':journal,'italic':True},{'text':' '},{'text':year,'bold':True},{'text':', '},{'text':vol,'italic':True},{'text':', '+page+'. https://doi.org/'+doi+'.'}]
  notes.append('DOI, title, complete author order, issue year, volume and pages verified against Crossref registration metadata; up to the first ten authors listed under MDPI ACS style.')
  if n==98:notes.append('Cites method publication 1.2.0; Methods must state actually used Vina 1.2.7.')
  if n==99:notes.append('Final peer-reviewed 2025 paper; do not replace with earlier preprint. Methods version 0.8.0.')
  if n==101:notes.append('2024 issue year; online-first publication 2 November 2023. AlphaFold DB model v6 is an entry version, not algorithm AlphaFold-6.')
  V.append({'number':n,'status':'primary_DOI_metadata_verified','doi':doi,'authors_total':len(a),'authors_listed':min(10,len(a)),'source_url':'https://api.crossref.org/works/'+doi,'metadata_file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 elif n==100:
  path=M/'new_ref_100.json';m=json.loads(path.read_text())['data']['attributes'];a=m['creators'];formatted=[]
  for c in a[:10]:
   # Zenodo creators include handles; do not invent corresponding personal identities.
   name=c['name'];parts=name.split()
   formatted.append(parts[-1]+', '+initials(' '.join(parts[:-1])) if len(parts)>1 else name)
  authors='; '.join(formatted)+'; et al.';title=m['titles'][0]['title'];doi=m['doi'];year=str(m['publicationYear']);ver=m['version']
  seg=[{'text':authors+' '},{'text':title,'italic':True},{'text':'; Version '+ver+'; Zenodo: '},{'text':year,'bold':True},{'text':'. https://doi.org/'+doi+'.'}]
  notes=['Exact software release DOI and version verified against DataCite/Zenodo. Creator names are retained as deposited; the handle sriniker is not expanded into an unverified identity. Do not substitute the all-versions DOI. RDKit runtime version 2026.03.6 corresponds to Release_2026_03_6.']
  V.append({'number':n,'status':'software_release_metadata_verified','doi':doi,'authors_total':len(a),'authors_listed':10,'source_url':'https://api.datacite.org/dois/'+doi,'metadata_file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 else:
  path=M/'new_ref_104_rcsb.json';m=json.loads(path.read_text());authors='; '.join(a['name'] for a in m['audit_author']);title=m['struct']['title'];year=m['rcsb_accession_info']['initial_release_date'][:4];doi='10.2210/pdb6JOL/pdb'
  seg=[{'text':authors+' '},{'text':title,'italic':True},{'text':'; Protein Data Bank Entry 6JOL; Protein Data Bank: '},{'text':year,'bold':True},{'text':'. https://doi.org/'+doi+'.'}]
  notes=['Dataset citation, not a journal article. Author order/title from RCSB audit_author and struct records; year denotes initial public release on 25 March 2020, not deposition in 2019. Associated literature remains To be published.']
  V.append({'number':n,'status':'official_structure_dataset_metadata_verified','doi':doi,'authors_total':3,'authors_listed':3,'source_url':'https://data.rcsb.org/rest/v1/core/entry/6JOL','metadata_file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 R.append({'number':n,'doi':doi,'text':''.join(x['text'] for x in seg),'segments':seg,'notes':notes})
for n,org,title,url,fn,note in [
 (108,'Boltz','Screen Small Molecule Libraries','https://api.boltz.bio/docs/api/guides/small-molecule-library-screen/','new_ref_108_boltz.html','Official documentation for the actual boltzmol pipeline, version 1.0, not a citation to the Boltz-2 research model. Original requested URL redirects to this canonical URL.'),
 (109,'Inductive Bio','Inductive Bio MCP Connector','https://www.inductive.bio/mcp-connector','new_ref_109_inductive.html','Official provider description of the MCP property models and model-discovery/prediction tools. Model identifiers/versions and exact returned values must be cited from archived run responses, not inferred from this overview.')]:
 t=f'{org}. {title}. Available online: {url} (accessed on 2 October 2026).';R.append({'number':n,'doi':None,'text':t,'segments':[{'text':t}],'notes':[note]});p=M/fn
 V.append({'number':n,'status':'official_provider_documentation_verified','source_url':url,'metadata_file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(B/'New_References_JCM_ACS.json').write_text(json.dumps(R,ensure_ascii=False,indent=2))
(B/'New_References_JCM_ACS.txt').write_text('\n\n'.join(f"{x['number']}. {x['text']}" for x in R)+'\n')
(B/'New_Reference_Audit.json').write_text(json.dumps(V,ensure_ascii=False,indent=2))
print((B/'New_References_JCM_ACS.txt').read_text())
