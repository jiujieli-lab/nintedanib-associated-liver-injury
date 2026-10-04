from retrieve_library import ROOT,RAW,get
import json,urllib.parse,concurrent.futures

tasks={
 'nintedanib_pubchem_inchikey':'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/XZXHXSATPCNXJR-ZIADKAODSA-N/property/IsomericSMILES,InChI,InChIKey,MolecularFormula/JSON',
 'MAZ51_pubchem_9839842':'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/9839842/synonyms/JSON',
 'nerandomilast_pubchem_166177189':'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/166177189/synonyms/JSON',
 'MAZ51_target_activities':'https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=CHEMBL597366&target_chembl_id__in=CHEMBL2007,CHEMBL1955&limit=1000'
}
for p in RAW.glob('*_target_activities.json'):
 for a in json.loads(p.read_text())['activities']:
  if a.get('document_chembl_id'):
   did=a['document_chembl_id'];tasks[did]='https://www.ebi.ac.uk/chembl/api/data/document/'+did+'.json'

def one(kv):
 k,u=kv
 p=RAW/(k+'.json')
 if p.exists():return
 try:
  j=get(u);p.write_text(json.dumps(j,indent=2));print(k,'OK',flush=True)
 except Exception as e:print(k,str(e),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(one,tasks.items()))
