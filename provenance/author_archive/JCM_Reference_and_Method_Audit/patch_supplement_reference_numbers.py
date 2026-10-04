from pathlib import Path
import fitz,re,json,hashlib,sys,collections
B=Path(__file__).parent
src=Path(sys.argv[1]) if len(sys.argv)>1 else Path('deliverables/Supplementary_Materials/Additional_file_1_Supplementary_Information.pdf')
out=Path(sys.argv[2]) if len(sys.argv)>2 else B/'Additional_file_1_Reference_Corrected.pdf'
mapping={int(k):int(v) for k,v in json.loads(Path(sys.argv[3]).read_text()).items()} if len(sys.argv)>3 else {54:55,55:58,56:59,57:60,58:61}
d=fitz.open(src);patches=[]
fontdir=Path('/opt/codex/runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/share/fonts/truetype')
for pi,p in enumerate(d):
 for b in p.get_text('rawdict')['blocks']:
  for line in b.get('lines',[]):
   chars=[c for s in line['spans'] for c in s['chars']];text=''.join(c['c'] for c in chars);matches=[]
   for m in re.finditer(r'\[[0-9,; –-]+\]',text):
    for num in re.finditer(r'\d+',m[0]):
     old=int(num[0]);start=m.start()+num.start();end=m.start()+num.end()
     if old in mapping:matches.append((start,end,old))
   if pi>=91:
    m=re.match(r'\s*(\d+)\.\s+[A-Z]',text)
    if m and int(m[1]) in mapping:matches.append((m.start(1),m.end(1),int(m[1])))
   for start,end,old in matches:
    c=chars[start];ss=next(s for s in line['spans'] if any(x is c for x in s['chars']));rect=fitz.Rect(chars[start]['bbox'])
    for cc in chars[start+1:end]:rect|=fitz.Rect(cc['bbox'])
    f=fontdir/('LiberationSerif-Bold.ttf' if 'Bold' in ss['font'] else 'LiberationSerif-Regular.ttf')
    patches.append({'page':pi+1,'old':str(old),'new':str(mapping[old]),'context':text,'rect':list(rect),'origin':c['origin'],'size':ss['size'],'fontfile':str(f)})
for pg in sorted(set(x['page'] for x in patches)):
 p=d[pg-1]
 for x in (x for x in patches if x['page']==pg):
  rect=fitz.Rect(x['rect']);rect.x0+=.03;rect.x1-=.03;p.add_redact_annot(rect,fill=(1,1,1))
 p.apply_redactions(images=0,graphics=0)
 for x in (x for x in patches if x['page']==pg):
  alias='RefB' if 'Bold.ttf' in x['fontfile'] else 'RefR';p.insert_font(fontname=alias,fontfile=x['fontfile']);p.insert_text(x['origin'],x['new'],fontname=alias,fontsize=x['size'])
d.save(out,garbage=4,deflate=True)
a=fitz.open(src);z=fitz.open(out);changed=sorted(set(x['page'] for x in patches));assert len(a)==len(z)==96
for i in range(len(z)):
 if i+1 not in changed:assert a[i].get_text()==z[i].get_text()
 else:
  exp=collections.Counter(''.join(a[i].get_text().split()))
  for x in patches:
   if x['page']==i+1:exp.subtract(x['old']);exp.update(x['new'])
  assert +exp==collections.Counter(''.join(z[i].get_text().split())),i+1
qa=B/'supplement_reference_qa';qa.mkdir(exist_ok=True)
for pg in changed:
 p=z[pg-1];p.get_pixmap(matrix=fitz.Matrix(1.3,1.3),alpha=False).save(qa/f'page_{pg:03d}.png')
log={'source':str(src),'output':str(out),'map':mapping,'pages':len(z),'changed_pages':changed,'patch_count':len(patches),'patches':patches,'unchanged_pages_exact_text':len(z)-len(changed),'changed_pages_character_inventory_verified':True,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
(B/'supplement_reference_patch_log.json').write_text(json.dumps(log,indent=2))
print(json.dumps({k:v for k,v in log.items() if k!='patches'},indent=2))
