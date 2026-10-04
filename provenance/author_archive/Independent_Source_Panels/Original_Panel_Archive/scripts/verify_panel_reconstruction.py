from pathlib import Path
from PIL import Image,ImageChops,ImageDraw
import csv,json,gc,hashlib
Image.MAX_IMAGE_PIXELS=None
ROOT=Path(__file__).resolve().parents[3]
rows=list(csv.DictReader(open(ROOT/'work/audit/panel_crop_manifest.csv')))
out=[]
for f in dict.fromkeys(r['original_figure'] for r in rows):
 rr=[r for r in rows if r['original_figure']==f]
 with Image.open(ROOT/rr[0]['source_file']) as src:
  src.load();reconstructed=Image.new('RGB',src.size,'white')
  for r in rr:
   xy=tuple(int(r['crop_'+x]) for x in ('left','top','right','bottom'))
   with Image.open(ROOT/r['output_file']) as p:
    p.load();mask=Image.new('L',p.size,255)
    for eb in json.loads(r.get('neighbour_text_exclusions','[]')):
     lx,ty,rx,by=eb;mask.paste(0,(lx-xy[0],ty-xy[1],rx-xy[0],by-xy[1]))
    reconstructed.paste(p,xy[:2],mask)
  bbox=ImageChops.difference(src,reconstructed).getbbox()
  out.append({'figure':f,'all_source_pixels_reconstructed':bbox is None,'difference_bbox':str(bbox),'total_panel_count':sum(r['role']=='panel' for r in rr),'supporting_elements':sum(r['role']!='panel' for r in rr)})
  print(f,bbox,flush=True)
  assert bbox is None,(f,bbox)
  del reconstructed
 gc.collect()
with open(ROOT/'work/audit/panel_source_reconstruction_qa.csv','w',newline='') as s:
 w=csv.DictWriter(s,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
print('Verified lossless reconstruction of all source figure pixels from the delivered panels and shared elements.',flush=True)
