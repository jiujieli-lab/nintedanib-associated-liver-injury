from urllib.request import Request,urlopen
from urllib.parse import quote
from pathlib import Path
import re,json,concurrent.futures
B=Path(__file__).parent;R=json.loads((B/'reference_verification.json').read_text())
work=[]
for r in R:
 if r['doi']:continue
 g=re.search(r'GSE\d+',r['text']);p=re.search(r'PXD\d+',r['text'])
 if g: u=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g[0]}&targ=self&form=text&view=brief'
 elif p:u='https://www.ebi.ac.uk/pride/ws/archive/v2/projects/'+p[0]
 elif r['number']==95:u='https://ddbj.nig.ac.jp/resource/bioproject/PRJDB12477.json'
 else:continue
 work.append((r['number'],u))
def get(nu):
 n,u=nu
 try:
  res=urlopen(Request(u,headers={'User-Agent':'Scholarly dataset citation verification'}),timeout=40);b=res.read();(B/'metadata'/f'ref_{n:03d}_details.txt').write_bytes(b)
  t=b.decode(errors='replace');return {'number':n,'url':u,'http':res.status,'bytes':len(b),'preview':t[:900]}
 except Exception as e:return {'number':n,'url':u,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:out=list(ex.map(get,work))
(B/'dataset_source_verification.json').write_text(json.dumps(out,ensure_ascii=False,indent=2));print(json.dumps(out,ensure_ascii=False,indent=2))
