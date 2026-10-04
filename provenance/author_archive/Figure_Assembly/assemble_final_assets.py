from pathlib import Path
import hashlib, json, os, shutil, struct
import fitz
from PIL import Image

BASE=Path('/workspace/scratch/c5634a2d588f')
REV=BASE/'work/jcm_revision'
FIG=REV/'figures'
SUB=REV/'submission'
ART=REV/'artwork_archive'
QA=FIG/'asset_qa'
for p in (SUB/'Main_Figures', SUB/'Supplementary_Materials', ART/'Main_Artwork', ART/'Supplementary_TIFFs', ART/'Supplementary_Artwork', ART/'Independent_New_Panels', QA): p.mkdir(parents=True,exist_ok=True)
Image.MAX_IMAGE_PIXELS=None
manifest={'main_figures':[],'supplementary_tiffs':[],'artwork_copies':[],'supplementary_pdf':{}}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''): h.update(chunk)
    return h.hexdigest()

def copy(src,dst,kind='artwork_copies'):
    src=Path(src);dst=Path(dst);dst.parent.mkdir(parents=True,exist_ok=True)
    tmp=dst.with_name(dst.name+'.partial');shutil.copyfile(src,tmp);os.replace(tmp,dst)
    s=sha(src);d=sha(dst);assert s==d
    rec={'source':str(src.relative_to(BASE)),'output':str(dst.relative_to(BASE)),'source_sha256':s,'output_sha256':d,'copy_byte_identical':True,'bytes':dst.stat().st_size}
    if dst.suffix.lower() in ('.tif','.tiff'):
        with dst.open('rb') as f: hdr=f.read(8)
        rec['first_ifd_offset']=struct.unpack('<I' if hdr[:2]==b'II' else '>I',hdr[4:8])[0]
        assert rec['first_ifd_offset']>0
        with Image.open(dst) as im: rec.update(size_px=list(im.size),mode=im.mode,dpi=[float(v) for v in im.info.get('dpi',[])])
    manifest[kind].append(rec);return rec

main={1:FIG/'final/Current_Figure_1.tiff',2:FIG/'repaired_quantitative/Figure_2.tiff',3:FIG/'repaired_quantitative/Figure_3.tiff',4:BASE/'figure4_repair/Figure_4_Repaired_1200dpi.tiff',5:FIG/'final/Current_Figure_6.tiff',6:FIG/'final/Current_Figure_7.tiff',7:FIG/'screening/final/Figure_7_Compound_Comparison.tiff'}
for n,src in main.items():
    rec=copy(src,SUB/f'Main_Figures/Figure_{n}.tiff','main_figures');rec['final_figure']=f'Figure {n}'
    for ext in ('.png','.pdf'):
        companion=src.with_suffix(ext)
        if companion.exists(): copy(companion,ART/f'Main_Artwork/Figure_{n}{ext}')
for ext in ('.png','.pdf','.svg'):
    copy(BASE/f'work/binding_analysis/results/Figure_4_Binding_Evidence{ext}',ART/f'Main_Artwork/Figure_4{ext}')

supp={**{n:FIG/f'repaired_quantitative/Figure_S{n}.tiff' for n in range(1,12)},12:FIG/'final/Current_Figure_5.tiff',13:FIG/'screening/final/Figure_S13_Docking_and_Evidence.tiff',14:FIG/'screening/final/Figure_S14_Properties_and_Receptor_Sensitivity.tiff'}
for n,src in supp.items():
    rec=copy(src,ART/f'Supplementary_TIFFs/Figure_S{n}.tiff','supplementary_tiffs');rec['final_figure']=f'Figure S{n}'
    if src.with_suffix('.png').exists():copy(src.with_suffix('.png'),ART/f'Supplementary_Artwork/Figure_S{n}.png')

original=BASE/'work/source/转化医学投稿20260922/Additional_Files/Additional_file_2_Supplementary_Figures.pdf'
original_doc=fitz.open(original);assert len(original_doc)==10
merged=fitz.open()
page_records=[]
for i,p in enumerate(original_doc):
    single=fitz.open();single.insert_pdf(original_doc,from_page=i,to_page=i)
    single.save(ART/f'Supplementary_Artwork/Figure_S{i+1}.pdf',garbage=3,deflate=True);single.close()
    r=p.rect;out=merged.new_page(width=r.width,height=r.height+24)
    out.show_pdf_page(fitz.Rect(0,24,r.width,r.height+24),original_doc,i,keep_proportion=False)
    out.insert_text((8,15),f'Supplementary Figure S{i+1}',fontname='hebo',fontsize=10,color=(0.12,0.18,0.21))
    page_records.append({'page':i+1,'figure':f'S{i+1}','source':str(original.relative_to(BASE)),'source_page':i+1,'source_sha256':sha(original),'original_page_preserved':True,'content_transform':'translation only: +24 pt y; no scaling or cropping','label_source':'added 24 pt external top margin; source labels retained'})
extras={11:BASE/'deliverables/Supplementary_Figures/Figure_S11_Vector.pdf',12:FIG/'final/Current_Figure_5.pdf',13:FIG/'screening/final/Figure_S13_Docking_and_Evidence.pdf',14:FIG/'screening/final/Figure_S14_Properties_and_Receptor_Sensitivity.pdf'}
for n,src in extras.items():
    d=fitz.open(src);assert len(d)==1
    r=d[0].rect;margin=24
    out=merged.new_page(width=r.width,height=r.height+margin)
    out.show_pdf_page(fitz.Rect(0,margin,r.width,r.height+margin),d,0,keep_proportion=False)
    out.insert_text((8,15),f'Supplementary Figure S{n}',fontname='hebo',fontsize=10,color=(0.12,0.18,0.21))
    page_records.append({'page':n,'figure':f'S{n}','source':str(src.relative_to(BASE)),'source_page':1,'source_sha256':sha(src),'original_page_preserved':True,'content_transform':'translation only: +24 pt y; no scaling or cropping','label_source':'added 24 pt external top margin'})
    copy(src,ART/f'Supplementary_Artwork/Figure_S{n}.pdf')
    d.close()
merged.set_toc([[1,f'Supplementary Figure S{n}',n] for n in range(1,15)])
merged.set_page_labels([{'startpage':0,'prefix':'S','style':'D','firstpagenum':1}])
merged.set_metadata({'title':'Supplementary Figures S1–S14','author':'','subject':'Final supplementary figure sequence','keywords':'S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14'})
output=SUB/'Supplementary_Materials/Additional_file_2_Supplementary_Figures.pdf'
tmp=output.with_suffix('.partial.pdf');merged.save(tmp,garbage=3,deflate=True);merged.close();os.replace(tmp,output)
manifest['supplementary_pdf']={'output':str(output.relative_to(BASE)),'sha256':sha(output),'bytes':output.stat().st_size,'page_count':14,'page_order':page_records,'no_added_title_pages':True}

# Independent new assets, retaining their native data and coordinate representations.
for src in (BASE/'work/binding_analysis/results/Individual_Panels').iterdir():
    if src.is_file() and src.suffix.lower() in ('.png','.pdf','.svg'):copy(src,ART/'Independent_New_Panels/Figure_4'/src.name)
for src in (BASE/'deliverables/Independent_Panels/Figure_S11').iterdir():
    if src.is_file() and src.suffix.lower() in ('.png','.pdf','.svg'):copy(src,ART/'Independent_New_Panels/Figure_S11'/src.name)
for folder,dest in [('fixed_panels','Screening_Fixed_Evidence'),('molecular','Screening_Molecular')]:
    for src in (FIG/'screening'/folder).iterdir():
        if src.is_file() and src.suffix.lower() in ('.png','.pdf','.svg','.json','.py'):copy(src,ART/'Independent_New_Panels'/dest/src.name)
for src in (FIG/'screening').iterdir():
    if src.is_file() and (src.suffix=='.py' or src.name.endswith(('_Legend.txt','_Manifest.json','_QA.json'))):copy(src,ART/'Provenance_and_Code/Screening'/src.name)
for name in ('build_workflow_figures.py','composite_validation.json','export_reopen_validation.json','FIGURE_DELIVERY_NOTES.md'):
    copy(FIG/name,ART/'Provenance_and_Code/Workflows'/name)
copy(BASE/'figure4_repair/validation.json',ART/'Provenance_and_Code/Figure_4_Repair_Validation.json')
copy(FIG/'repaired_quantitative/Original_Figure_Export_Manifest.json',ART/'Provenance_and_Code/Original_Figure_Export_Manifest.json')
copy(FIG/'repaired_quantitative/Separate_Process_Decode_QA.json',ART/'Provenance_and_Code/Original_Figure_Decode_QA.json')

manifest['main_figure_mapping']={'Current_Figure_1':'Figure_1','Original_Figure_2':'Figure_2','Original_Figure_3':'Figure_3','Figure_4_Repaired_1200dpi':'Figure_4','Current_Figure_6':'Figure_5','Current_Figure_7':'Figure_6','New_Compound_Comparison':'Figure_7','Current_Figure_5':'Figure_S12'}
(REV/'Final_Figure_Asset_Map.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'main_tiffs':len(manifest['main_figures']),'supplementary_tiffs':len(manifest['supplementary_tiffs']),'supplementary_pdf_pages':len(page_records),'main_tiff_bytes':sum(x['bytes'] for x in manifest['main_figures']),'supplementary_pdf_bytes':output.stat().st_size,'artwork_copies':len(manifest['artwork_copies'])}))
