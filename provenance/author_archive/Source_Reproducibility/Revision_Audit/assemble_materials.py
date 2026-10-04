from pathlib import Path
import shutil,zipfile,json,hashlib,io,csv
from PIL import Image
from pypdf import PdfReader,PdfWriter
ROOT=Path(__file__).resolve().parents[2];SRC=ROOT/'work/source/转化医学投稿20260922';OUT=ROOT/'deliverables';SUP=OUT/'Supplementary_Materials'
Image.MAX_IMAGE_PIXELS=None
checks=[]
def lossless_tiff(src,dst):
 with Image.open(src) as im:
  dpi=im.info.get('dpi',(300,300));before=hashlib.sha256(im.tobytes()).hexdigest();size=im.size;mode=im.mode
  im.save(dst,compression='tiff_adobe_deflate',dpi=dpi)
 with Image.open(dst) as im:
  after=hashlib.sha256(im.tobytes()).hexdigest();assert before==after and im.size==size and im.mode==mode
 checks.append({'file':str(dst.relative_to(OUT)),'pixels':size,'dpi':[float(d) for d in dpi],'pixel_sha256':before,'source_pixels_unchanged':True})
for n in [1,2,3,5,6,7]:
 shutil.copy2(ROOT/f'work/panels/reassembled/Figure_{n}_Reassembled_Original_Panels.tiff',OUT/f'Main_Figures/Figure_{n}.tiff')
lossless_tiff(ROOT/'work/binding_analysis/results/Figure_4_Binding_Evidence_1200dpi.tiff',OUT/'Main_Figures/Figure_4.tiff')
shutil.copy2(ROOT/'work/binding_analysis/results/Figure_4_Binding_Evidence.pdf',OUT/'Main_Figures/Figure_4_Vector.pdf')
lossless_tiff(SRC/'Graphical_Abstract.tiff',OUT/'Graphical_Abstract.tiff')
figdir=OUT/'Supplementary_Figures';figdir.mkdir(exist_ok=True)
for n in range(1,11):lossless_tiff(SRC/f'supplementary Figures_S1-S10/Figure_S{n}_1500dpi.tiff',figdir/f'Figure_S{n}.tiff')
lossless_tiff(ROOT/'work/new_virtual/Figure_S11_Virtual_Perturbation_1200dpi.tiff',figdir/'Figure_S11.tiff')
shutil.copy2(ROOT/'work/supplementary/Additional_file_1_Supplementary_Information_Revised.pdf',SUP/'Additional_file_1_Supplementary_Information.pdf')
shutil.copy2(SRC/'Additional_Files/Additional_file_3_Supplementary_Tables.xlsx',SUP/'Additional_file_3_Supplementary_Tables.xlsx')
shutil.copy2(ROOT/'work/tables/Additional_file_5_Targeted_Analysis_Tables.xlsx',SUP/'Additional_file_5_Targeted_Analysis_Tables.xlsx')
w=PdfWriter();w.append(SRC/'Additional_Files/Additional_file_2_Supplementary_Figures.pdf');w.append(ROOT/'work/new_virtual/Figure_S11_Virtual_Perturbation.pdf')
w.add_metadata({'/Title':'Supplementary Figures S1–S11','/Author':'Shaoxin Huang','/Subject':'Original supplementary figures and targeted network intervention results'})
for n in range(1,12):w.add_outline_item(f'Supplementary Figure S{n}',n-1)
with (SUP/'Additional_file_2_Supplementary_Figures.pdf').open('wb') as f:w.write(f)
assert len(PdfReader(SUP/'Additional_file_2_Supplementary_Figures.pdf').pages)==11
old=PdfReader(SRC/'Additional_Files/Additional_file_2_Supplementary_Figures.pdf');new=PdfReader(SUP/'Additional_file_2_Supplementary_Figures.pdf')
assert all(old.pages[i].extract_text()==new.pages[i].extract_text() for i in range(10))
# Five original data/code packages retained byte-identically in a separate author archive.
rdir=OUT/'Reproducibility/Original_Source_Packages';rdir.mkdir(parents=True,exist_ok=True)
for f in sorted((SRC/'Additional_Files').glob('*.zip')):shutil.copy2(f,rdir/f.name)
# Original panel archive plus all fifteen new independent panels.
shutil.copy2(ROOT/'work/Original_Panel_Archive.zip',OUT/'Independent_Panels/Original_165_Panels.zip')
for name,base in [('Figure_4',ROOT/'work/binding_analysis/results/Individual_Panels'),('Figure_S11',ROOT/'work/new_virtual/Individual_Panels')]:
 dest=OUT/'Independent_Panels'/name;dest.mkdir(exist_ok=True)
 for p in sorted(base.glob('*.png')):shutil.copy2(p,dest/p.name)
# Reproducible targeted archive with exactly the source inputs expected by executable scripts.
stage=ROOT/'work/build/targeted_archive';stage.mkdir(exist_ok=True)
for sub in ['new_virtual','binding_analysis']:
 src=ROOT/'work'/sub;dst=stage/sub;dst.mkdir(exist_ok=True)
 for p in src.iterdir():
  if p.is_file() and p.suffix in ['.py','.csv','.gz','.json','.txt'] and not p.name.startswith(('boltz_','Boltz_','run_md_','analyze_md_','manuscript_','Kd_')) and p.name not in ['run.log']:
   shutil.copy2(p,dst/p.name)
 if sub=='binding_analysis':
  (dst/'results').mkdir(exist_ok=True)
  for p in (src/'results').glob('*'):
   if p.is_file() and p.suffix in ['.csv','.json']:shutil.copy2(p,dst/'results'/p.name)
# Exclude planned MD from completed-analysis code archive; add actual inputs.
paths=[
'data/Additional_file_5_Liver_and_Lung_Single_Cell_Data/Source_Data/Liver_Single_Cell/hepatocyte_network_matrix.csv.gz',
'data/Additional_file_5_Liver_and_Lung_Single_Cell_Data/Source_Data/Liver_Single_Cell/endothelial_network_matrix.csv.gz',
'data/Additional_file_5_Liver_and_Lung_Single_Cell_Data/Source_Data/Liver_Single_Cell/macrophage_network_matrix.csv.gz',
'data/Additional_file_5_Liver_and_Lung_Single_Cell_Data/Source_Data/Liver_Single_Cell/cell_metadata.csv',
'data/Additional_file_4_Clinical_Evidence_and_Analysis_Code/Analysis/network_diffusion/network_clinical_seed_weights.csv',
'data/Additional_file_4_Clinical_Evidence_and_Analysis_Code/Analysis/network_diffusion/network_specification_scores.csv']
for rel in paths:
 dst=stage/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/'work'/rel,dst)
(stage/'README.txt').write_text('Completed targeted analyses and exact input data\n\nUnzip while preserving all relative folders. From this directory run:\npython new_virtual/run_target_pair_perturbation.py\npython new_virtual/make_figure_s11.py\npython new_virtual/validate_extension.py\npython binding_analysis/analyze_binding_evidence.py\n\nThe Python environment needs numpy, scipy, pandas, scikit-learn, threadpoolctl, matplotlib and Pillow. Actual package versions appear in ENVIRONMENT.json. Calculations use the supplied frozen inputs. Prediction responses from Inductive Bio are archived real service results; rerunning the local binding script does not silently call the service. No conventional docking, Boltz prediction or MD trajectory is included in these completed results.\n\nSupplementary Tables N18 and N19 are new_virtual/gene_level_perturbations.csv.gz and new_virtual/matched_null_metrics.csv.gz.\n')
import importlib.metadata
versions={k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','scikit-learn','threadpoolctl','matplotlib','Pillow']};(stage/'ENVIRONMENT.json').write_text(json.dumps(versions,indent=2))
files=[{'file':str(p.relative_to(stage)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(stage.rglob('*')) if p.is_file() and p.name!='MANIFEST.json']
(stage/'MANIFEST.json').write_text(json.dumps(files,indent=2))
zipout=SUP/'Additional_file_6_Targeted_Analysis_Code.zip'
with zipfile.ZipFile(zipout,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(stage.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(stage))
with zipfile.ZipFile(zipout) as z:assert z.testzip() is None
# Concrete future inputs kept outside completed results.
pending=OUT/'Pending_Structure_and_MD'
for name in ['Boltz_PDGFRA_Nintedanib_Payload.json','Boltz_FLT4_Nintedanib_Payload.json','Boltz_Input_Manifest.json','PDGFRA_Kinase_Domain.fasta','FLT4_Kinase_Domain.fasta','PDGFRA_UniProt.json','FLT4_UniProt.json','Nintedanib_ChEMBL.json','run_md_after_structure_qc.py','analyze_md_trajectory.py','Non_Docking_Structure_and_MD_Execution_Plan.md']:
 shutil.copy2(ROOT/'work/binding_analysis'/name,pending/name)
# Keep document assembly and verification source with local-relative inputs documented.
audit=OUT/'Reproducibility/Revision_Audit';audit.mkdir(exist_ok=True)
for p in [ROOT/'work/build/manuscript_change_log.json',ROOT/'work/build/reference_number_map.json',ROOT/'work/build/original_results_retention.json',ROOT/'work/build/main_build_metrics.json',ROOT/'work/audit/panel_inventory.csv',ROOT/'work/audit/panel_crop_manifest.csv',ROOT/'work/audit/panel_source_reconstruction_qa.csv',ROOT/'work/tables/workbook_qa.json',ROOT/'work/tables/table_manifest.json',ROOT/'work/supplementary/patch_log.json']:
 shutil.copy2(p,audit/p.name)
for p in (ROOT/'work/build').glob('*.py'):shutil.copy2(p,audit/p.name)
shutil.copy2(ROOT/'work/tables/build_targeted_tables.mjs',audit/'build_targeted_tables.mjs')
(audit/'lossless_image_export_checks.json').write_text(json.dumps(checks,indent=2))
print(json.dumps({'new_code_zip_MB':zipout.stat().st_size/1e6,'supplementary_pdf_pages':11,'lossless_tiff_exports':len(checks),'main_figures':7,'supplementary_figures':11},indent=2))
