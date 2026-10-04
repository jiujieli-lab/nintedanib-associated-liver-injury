from pathlib import Path
import hashlib,json,re,gc,numpy as np
from PIL import Image,ImageDraw
import fitz

BASE=Path('/workspace/scratch/c5634a2d588f');REV=BASE/'work/jcm_revision';QA=REV/'figures/asset_qa';ART=REV/'artwork_archive'
Image.MAX_IMAGE_PIXELS=None
manifest=json.loads((REV/'Final_Figure_Asset_Map.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()

decode=[]
for entry in manifest['main_figures']+manifest['supplementary_tiffs']:
 p=BASE/entry['output']
 with Image.open(p) as im:
  im.load();assert im.mode=='RGB';assert list(im.size)==entry['size_px']
 assert sha(p)==entry['source_sha256']
 decode.append({'file':entry['output'],'decode_after_copy_process_exit':True,'source_hash_equal':True})
 gc.collect()
 print('Decoded '+entry['final_figure'],flush=True)

pdf=BASE/manifest['supplementary_pdf']['output'];doc=fitz.open(pdf);assert len(doc)==14
assert doc.get_page_labels()==[{'startpage':0,'prefix':'S','firstpagenum':1,'style':'D'}]
page_checks=[]
for i,entry in enumerate(manifest['supplementary_pdf']['page_order']):
 src=fitz.open(BASE/entry['source']);sp=src[entry['source_page']-1]
 clip=fitz.Rect(0,24,sp.rect.width,sp.rect.height+24)
 a=sp.get_pixmap(matrix=fitz.Matrix(1,1),colorspace=fitz.csRGB,alpha=False)
 b=doc[i].get_pixmap(matrix=fitz.Matrix(1,1),colorspace=fitz.csRGB,alpha=False,clip=clip)
 equal=a.width==b.width and a.height==b.height and a.samples==b.samples
 assert a.width==b.width and a.height==b.height
 delta=np.abs(np.frombuffer(a.samples,dtype=np.uint8).astype(np.int16)-np.frombuffer(b.samples,dtype=np.uint8).astype(np.int16))
 form=next(x[0] for x in doc[i].get_xobjects() if x[1]=='fullpage')
 stream_equal=sp.read_contents()==doc.xref_stream(form)
 assert stream_equal,(i+1,'page content stream changed')
 assert float(delta.mean())<0.02,(i+1,float(delta.mean()))
 text=doc[i].get_text()
 number=bool(re.search(r'(?:Figure|Fig\.?)\s*S\s*'+str(i+1)+r'\b',text,re.I))
 assert number,(i+1,'figure label not found')
 out=QA/f'Supplement_Page_{i+1:02}.png';doc[i].get_pixmap(dpi=90).save(out)
 page_checks.append({'page':i+1,'figure':f'S{i+1}','correct_figure_number_in_text':number,'original_page_content_stream_byte_identical':stream_equal,'source_content_render_pixel_identical_at_72dpi':equal,'render_mean_absolute_8bit_difference':float(delta.mean()),'render_max_absolute_8bit_difference':int(delta.max()),'preview':str(out.relative_to(BASE))})
 src.close();print('Verified PDF S'+str(i+1),flush=True)

for start in (1,5,9,13):
 cards=[]
 for n in range(start,min(start+4,15)):
  im=Image.open(QA/f'Supplement_Page_{n:02}.png').convert('RGB');im.thumbnail((500,700))
  card=Image.new('RGB',(520,740),'#dedede');card.paste(im,((520-im.width)//2,30));ImageDraw.Draw(card).text((10,8),f'S{n}',fill='black');cards.append(card)
 sheet=Image.new('RGB',(520*len(cards),740),'white')
 for n,card in enumerate(cards):sheet.paste(card,(520*n,0))
 sheet.save(QA/f'Supplement_Contact_{start}_{min(start+3,14)}.png')

result={'copied_tiffs':decode,'pdf_page_count':14,'pdf_order':[f'S{i}' for i in range(1,15)],'page_checks':page_checks,'no_added_title_pages':True,'figures_preserved_without_scaling':True,'pdf_visually_verified':False}
(QA/'Final_Asset_Verification.json').write_text(json.dumps(result,indent=2)+'\n')
manifest['verification']=result
manifest['workflow_svg_qa']=json.loads((ART/'Workflow_Panels/Workflow_SVG_QA.json').read_text())
manifest['artwork_archive_inventory']=[{'file':str(p.relative_to(REV)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(ART.rglob('*')) if p.is_file()]
(REV/'Final_Figure_Asset_Map.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
print('Final verification complete.',flush=True)
