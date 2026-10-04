from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import quote
import json,concurrent.futures,time,re,hashlib
B=Path(__file__).parent
refs=json.loads((B/'references_original.json').read_text());(B/'metadata').mkdir(exist_ok=True)
def fetch(r):
 n=r['number']; d=r['doi'];out=dict(r)
 try:
  if d:
   u=('https://api.datacite.org/dois/' if d.startswith('10.5281/') else 'https://api.crossref.org/works/')+quote(d,safe='')
  else:
   u=re.search(r'https?://\S+',r['text']).group().rstrip('.,;')
  p=B/'metadata'/f'ref_{n:03d}.json' if d else B/'metadata'/f'ref_{n:03d}.html'
  if p.exists(): b=p.read_bytes();http=200
  else:
   for attempt in range(3):
    try:
     res=urlopen(Request(u,headers={'User-Agent':'Scholarly bibliographic metadata audit'}),timeout=40);b=res.read();http=res.status;p.write_bytes(b);break
    except Exception:
     if attempt==2:raise
     time.sleep(1+attempt)
  out.update(source_url=u,http_status=http,bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
  if d:
   raw=json.loads(b);m=raw.get('message') or raw.get('data',{}).get('attributes',{})
   if d.startswith('10.5281/'):
    title=m.get('titles',[{}])[0].get('title','');authors=m.get('creators',[]);journal=m.get('publisher','');year=m.get('publicationYear');vol='';page='';mdoi=m.get('doi','')
   else:
    title='; '.join(m.get('title',[]));authors=m.get('author',[]);journal='; '.join(m.get('container-title',[]));year=(m.get('published-print') or m.get('published') or m.get('issued') or {}).get('date-parts',[[None]])[0][0];vol=m.get('volume','');page=m.get('page') or m.get('article-number','');mdoi=m.get('DOI','')
   norm=lambda t: re.sub(r'[^a-z0-9]','',t.lower())
   titleok=norm(title) in norm(r['text']);out.update(metadata_title=title,metadata_authors=authors,metadata_journal=journal,metadata_year=year,metadata_volume=vol,metadata_page=page,metadata_doi=mdoi,title_exact_normalized_match=titleok,doi_match=mdoi.lower()==d.lower(),status='doi_metadata_match' if titleok else 'doi_resolved_title_review')
  else:out['status']='official_url_fetched_requires_content_review'
 except Exception as e:out.update(status='unresolved',error=str(e))
 print(n,out['status'],flush=True)
 return out
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:out=list(ex.map(fetch,refs))
(B/'reference_verification.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
print('WROTE',len(out))
