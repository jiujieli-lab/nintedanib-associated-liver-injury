from pathlib import Path
import fitz,json,re,hashlib
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];WORK=ROOT.parent
src=WORK/'source/转化医学投稿20260922/Additional_Files/Additional_file_1_Supplementary_Information.pdf'
out=ROOT/'submission/Supplementary_Materials/Additional_file_1_Supplementary_Information.pdf';out.parent.mkdir(exist_ok=True,parents=True)
doc=fitz.open(src)
patches=[p for p in json.loads((WORK/'supplementary/patch_log.json').read_text())['patches'] if p['kind'] in ('line','title')]
for p in patches:
 if p['page']==1:
  if p['new'].startswith('Figures S1'):p['new']='Figures S1–S14 and tables are in files 2–3; the source crosswalk is on page 93.'
  if p['new'].startswith('Data and code'):p['new']='Additional files 4–6 detail targeted analyses; files 7–9 describe compound screening.'
  if p['new'].startswith('binding and virtual'):p['new']='The prospective hepatic perturbation protocol has not been executed.'
 if p['page']==83 and p['new'].startswith('Reproducibility'):p['new']='The author reproducibility archive contains machine-readable inputs and outputs, including'
# Canonical main citation numbers; no duplicated bibliography in this supplement.
m0=json.loads((WORK/'build/reference_number_map.json').read_text());m1=json.loads((ROOT/'manuscript/Final_Reference_Number_Map.json').read_text())
fontdir=Path('/opt/codex/runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/share/fonts/truetype')
pat=re.compile(r'\[([0-9,;\s–-]+)\]')
def numbers(s):
 out=[]
 for part in re.split(r'[,;]\s*',s):
  pair=re.split(r'[–-]',part.strip());out.extend(range(int(pair[0]),int(pair[1])+1) if len(pair)==2 else [int(pair[0])])
 return out
for pg,page in enumerate(doc[:92],1):
 for b in page.get_text('rawdict')['blocks']:
  for line in b.get('lines',[]):
   chars=[(c,s) for s in line['spans'] for c in s['chars']];txt=''.join(c['c'] for c,s in chars)
   matches=[(m,'citation') for m in pat.finditer(txt)]
   if pg in (48,49,69):matches.extend((m,'figure-label') for m in re.finditer(r'Figure 4(?=[A-Z\s).,]|$)',txt))
   for m,kind in matches:
    if kind=='citation':
     vals=sorted(set(m1[str(m0[str(x)])] for x in numbers(m[1])));new='['+', '.join(map(str,vals))+']'
    else:new='Fig. S12'
    old=m.group()
    if old==new:continue
    cc=chars[m.start():m.end()];r=fitz.Rect(cc[0][0]['bbox'])
    for c,s in cc[1:]:r|=fitz.Rect(c['bbox'])
    c,s=cc[0];fontfile=fontdir/('LiberationSerif-Bold.ttf' if 'Bold' in s['font'] else 'LiberationSerif-Regular.ttf')
    font=fitz.Font(fontfile=str(fontfile));size=s['size'];natural=font.text_length(new,fontsize=size)
    # Fit the reference within exactly the original span, preserving adjacent source text.
    size=min(size,size*r.width/natural)
    patches.append({'page':pg,'kind':kind,'old':old,'new':new,'rect':list(r),'origin':list(c['origin']),'size':size,'fontfile':str(fontfile),'context':txt})
# Remove duplicated old-format bibliography, retaining all results and the complete S10 legend.
for pg in sorted(set(p['page'] for p in patches)|{92}):
 page=doc[pg-1];edits=[p for p in patches if p['page']==pg]
 for p in edits:
  r=fitz.Rect(p['rect']);r.x0+=.015;r.x1-=.015;page.add_redact_annot(r,fill=(1,1,1))
 if pg==92:page.add_redact_annot(fitz.Rect(38,468,550,790),fill=(1,1,1))
 page.apply_redactions(images=0,graphics=0)
 for p in edits:
  alias='JCMSerifB' if 'Bold.ttf' in p['fontfile'] else 'JCMSerif';page.insert_font(fontname=alias,fontfile=p['fontfile'])
  if p['kind']=='title':
   for t,xy in zip(p['lines'],p['origins']):page.insert_text(xy,t,fontname=alias,fontsize=p['size'])
  else:page.insert_text(p['origin'],p['new'],fontname=alias,fontsize=p['size'])
 if pg==92:
  page.insert_font(fontname='JCMSerif',fontfile=str(fontdir/'LiberationSerif-Regular.ttf'))
  page.insert_text((72.1,485),'References are numbered as in the main manuscript.',fontname='JCMSerif',fontsize=12)
doc.delete_pages(92,95)
cross=fitz.open(ROOT/'supplementary_information/Source_Figure_Crosswalk.pdf');assert len(cross)==1;doc.insert_pdf(cross)
doc.set_metadata({**doc.metadata,'title':'Supplementary information: Clinical proteomics-guided network integration prioritizes PDGFRA and FLT4 for investigation of nintedanib-associated liver injury'})
doc.save(out,garbage=4,deflate=True)
check=fitz.open(out);orig=fitz.open(src);changed=sorted(set(p['page'] for p in patches)|{92});unchanged=0
for i in range(91):
 if i+1 not in changed:
  assert check[i].get_text()==orig[i].get_text();unchanged+=1
 else:
  expected=Counter(''.join(orig[i].get_text().split()))
  for p in [x for x in patches if x['page']==i+1]:expected.subtract(''.join(p['old'].split()));expected.update(''.join(p['new'].split()))
  actual=Counter(''.join(check[i].get_text().split()))
  assert +expected==actual,(i+1,+(expected-actual),+(actual-expected))
# All pre-bibliography S10 content on the last retained page survives.
oldlast=orig[91].get_text().split('References')[0].splitlines()[1:]
newlast=check[91].get_text()
for line in oldlast:
 if line.strip():assert line.strip() in newlast,line
log={'source':str(src),'output':str(out),'pages':len(check),'changed_pages':changed,'unchanged_pages_exact_text_match':unchanged,'all_result_text_preserved':True,'changed_pages_exact_character_inventory':True,'reference_map_original_to_final':{k:m1[str(v)] for k,v in m0.items()},'patches':patches,'bibliography':'Duplicate source-style list removed; citations point to full ACS list in main manuscript. No Results removed.','sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
(ROOT/'supplementary_information/Final_Supplement1_Patch_Audit.json').write_text(json.dumps(log,indent=2,ensure_ascii=False))
r=ROOT/'supplementary_information/render';r.mkdir(exist_ok=True)
for pg in changed+[93]:check[pg-1].get_pixmap(matrix=fitz.Matrix(1.2,1.2),alpha=False).save(r/f'page_{pg:03}.png')
print(json.dumps({k:v for k,v in log.items() if k not in ('patches','reference_map_original_to_final')},indent=2))
