from pathlib import Path
import json,re
R=Path(__file__).resolve().parents[1];audit=[]
for area in ['screening_archive','author_archive']:
 for p in (R/area).rglob('*.json'):
  if p.stat().st_size>15_000_000:continue
  obj=json.loads(p.read_text());changed=[]
  def clean(x,path=''):
   if isinstance(x,dict):return {k:clean(v,path+'/'+k) for k,v in x.items()}
   if isinstance(x,list):return [clean(v,path+'/'+str(i)) for i,v in enumerate(x)]
   if isinstance(x,str) and x.startswith('http') and re.search(r'[?&](?:AWSAccessKeyId|X-Amz-Signature|Signature)=',x,re.I):
    changed.append(path);return 'TEMPORARY_ASSET_TRANSFER_URL_OMITTED'
   return x
  result=clean(obj)
  if changed:
   if isinstance(result,dict):result['_archive_note']='Temporary signed asset-transfer URLs were omitted. Scientific fields and downloaded coordinate files are unchanged.'
   p.write_text(json.dumps(result,indent=2));audit.append({'file':str(p.relative_to(R)),'transport_url_fields_removed':changed})
(R/'packaging/Transport_Metadata_Audit.json').write_text(json.dumps(audit,indent=2));print({'files':len(audit),'url_fields':sum(len(x['transport_url_fields_removed']) for x in audit)})
