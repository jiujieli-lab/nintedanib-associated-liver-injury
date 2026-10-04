import concurrent.futures, urllib.request, urllib.parse, json, pathlib, time, hashlib, csv

ROOT=pathlib.Path(__file__).parent
RAW=ROOT/'raw_chembl'; RAW.mkdir(exist_ok=True)
NAMES=['nintedanib','imatinib','crenolanib','avapritinib','ripretinib','axitinib','tivozanib','lenvatinib','sunitinib','pazopanib','sorafenib','regorafenib','ponatinib','MAZ51','pirfenidone','nerandomilast']

def get(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url,timeout=50) as r: return json.load(r)
        except Exception:
            if attempt==2: raise
            time.sleep(1+attempt)

def one(name):
    url='https://www.ebi.ac.uk/chembl/api/data/molecule.json?'+urllib.parse.urlencode({'pref_name__iexact':name.upper(),'limit':100})
    prior=RAW/(name+'_lookup.json')
    data=json.loads(prior.read_text()) if prior.exists() else get(url)
    (RAW/(name+'_lookup.json')).write_text(json.dumps(data,indent=2))
    rows=data.get('molecules',[])
    if not rows:
        url='https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?'+urllib.parse.urlencode({'q':name,'limit':100})
        p=RAW/(name+'_search.json');data=json.loads(p.read_text()) if p.exists() else get(url); p.write_text(json.dumps(data,indent=2)); rows=data.get('molecules',[])
    exact=[m for m in rows if (m.get('pref_name') or '').lower()==name.lower() or any((v.get('molecule_synonym') or '').lower().replace('-','')==name.lower().replace('-','') for v in m.get('molecule_synonyms',[]))]
    if not exact:
        print(name,'NO EXACT MATCH',[(m['molecule_chembl_id'],m['pref_name']) for m in rows],flush=True)
        return {'name':name,'error':'no exact preferred-name match','lookup_url':url}
    m=exact[0]; parent=m['molecule_hierarchy']['parent_chembl_id']
    if parent !=m['molecule_chembl_id']:m=get('https://www.ebi.ac.uk/chembl/api/data/molecule/'+parent+'.json')
    (RAW/(name+'_parent.json')).write_text(json.dumps(m,indent=2))
    s=m.get('molecule_structures') or {}
    out={'name':name,'chembl_id':m['molecule_chembl_id'],'pref_name':m['pref_name'],'max_phase':m['max_phase'],'canonical_smiles':s.get('canonical_smiles'),'standard_inchi':s.get('standard_inchi'),'standard_inchi_key':s.get('standard_inchi_key'),'structure_url':'https://www.ebi.ac.uk/chembl/api/data/molecule/'+m['molecule_chembl_id']+'.json','retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'source_parent':parent,'source_molecule_properties':m.get('molecule_properties')}
    print(name,m['molecule_chembl_id'],'phase',m['max_phase'],flush=True)
    return out

if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        vals=list(ex.map(one,NAMES))
    (ROOT/'candidate_structures.json').write_text(json.dumps(vals,indent=2))
    with (ROOT/'candidate_structures.csv').open('w') as f:
        fields=['name','chembl_id','pref_name','max_phase','canonical_smiles','standard_inchi','standard_inchi_key','structure_url','retrieved_utc','source_parent','error','lookup_url']
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(vals)
