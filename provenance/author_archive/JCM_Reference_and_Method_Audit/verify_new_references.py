from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import quote
import concurrent.futures,json
B=Path(__file__).parent/'new_reference_metadata';B.mkdir(exist_ok=True)
ids=['10.1021/acs.jcim.1c00203','10.1021/acs.jcim.5c02271','10.5281/zenodo.22140358','10.1093/nar/gkad1011','10.1093/nar/gkab1061','10.1038/s41586-021-03819-2','10.2210/pdb6JOL/pdb','10.1016/j.ajpath.2020.06.006','10.1007/s10456-020-09718-w','10.1038/nbt1358']
def get(t):
 n,d=t;u=('https://api.datacite.org/dois/' if d.startswith(('10.5281/','10.2210/')) else 'https://api.crossref.org/works/')+quote(d,safe='');o={'number':n,'doi':d,'source_url':u}
 try:
  r=urlopen(Request(u,headers={'User-Agent':'Scholarly bibliographic verification'}),timeout=45);raw=r.read();p=B/f'new_ref_{n}.json';p.write_bytes(raw);m=json.loads(raw);o.update(status='fetched',file=str(p))
 except Exception as e:o.update(status='error',error=str(e))
 return o
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:r=list(ex.map(get,enumerate(ids,98)))
(B/'index.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
