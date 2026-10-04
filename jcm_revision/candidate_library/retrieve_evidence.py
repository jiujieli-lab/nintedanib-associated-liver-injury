from retrieve_library import ROOT,RAW,get
import concurrent.futures, urllib.parse,json,time,csv
names=json.loads((ROOT/'candidate_structures.json').read_text())

def compound(n):
    name=n['name']
    out={'name':name}
    # Independently resolve the supplied common drug name in PubChem.
    u='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.parse.quote(name)+'/property/IsomericSMILES,CanonicalSMILES,InChI,InChIKey,MolecularFormula,MolecularWeight/JSON'
    try:
        p=get(u); (RAW/(name+'_pubchem.json')).write_text(json.dumps(p,indent=2));out['pubchem']=p;out['pubchem_url']=u
    except Exception as e:out['pubchem_error']=str(e)
    if n.get('chembl_id'):
        u='https://www.ebi.ac.uk/chembl/api/data/activity.json?'+urllib.parse.urlencode({'molecule_chembl_id':n['chembl_id'],'target_chembl_id__in':'CHEMBL2007,CHEMBL1955','limit':1000})
        try:
            a=get(u); (RAW/(name+'_target_activities.json')).write_text(json.dumps(a,indent=2));out['activities_url']=u;out['activities_count']=len(a.get('activities',[]));out['activity_page_meta']=a['page_meta']
        except Exception as e:out['activity_error']=str(e)
    print(name,out.get('activities_count'),out.get('pubchem_error','OK'),flush=True)
    return out

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: out=list(ex.map(compound,names))
(ROOT/'cross_validation_and_activities_manifest.json').write_text(json.dumps(out,indent=2))
