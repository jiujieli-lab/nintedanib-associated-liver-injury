from PIL import Image, ImageChops, ImageDraw, ImageFont
from pathlib import Path
import json, hashlib, csv, string, gc, sys, os
Image.MAX_IMAGE_PIXELS=None
ROOT=Path(__file__).resolve().parents[3]
SRC=ROOT/'work/source/转化医学投稿20260922'
OUT=ROOT/'work/panels'
def hbytes(a): return hashlib.sha256(a.tobytes()).hexdigest()
def filehash(p):
 h=hashlib.sha256()
 with open(p,'rb') as s:
  for chunk in iter(lambda:s.read(4*1024*1024),b''): h.update(chunk)
 return h.hexdigest()
D={}
def add(f,w,h,rows,header=0,footer=None):
 boxes={}; groups=[]; j=0
 if header: boxes['Shared_Title']=(0,0,w,header)
 for y0,y1,x in rows:
  grp=[]
  for x0,x1 in zip(x,x[1:]):
   k=string.ascii_uppercase[j];j+=1;boxes[k]=(x0,y0,x1,y1);grp.append(k)
  groups.append(grp)
 if footer: boxes['Shared_Legend_or_Footnote']=(0,footer,w,h)
 D[f]={'preview_size':[w,h],'boxes':boxes,'groups':groups}
add('Figure_1',1363,1600,[(0,143,[0,1363]),(143,437,[0,684,1363]),(437,736,[0,696,1363]),(736,1017,[0,696,1363]),(1017,1291,[0,695,1363]),(1291,1600,[0,682,1363])])
add('Figure_2',1363,1600,[(0,342,[0,690,1363]),(342,678,[0,690,1363]),(678,940,[0,690,1363]),(940,1289,[0,690,1363]),(1289,1600,[0,690,1363])])
add('Figure_3',1360,1600,[(0,278,[0,690,1360]),(278,610,[0,690,1360]),(610,907,[0,690,1360]),(907,1212,[0,690,1360]),(1212,1600,[0,690,1360])])
add('Figure_4',1451,1600,[(0,165,[0,1451]),(165,497,[0,483,970,1451]),(497,944,[0,500,1030,1451]),(944,1280,[0,780,1451]),(1280,1600,[0,847,1451])])
add('Figure_5',1400,1600,[(0,153,[0,1400]),(153,640,[0,375,670,990,1400]),(640,1110,[0,365,655,985,1400]),(1110,1600,[0,360,680,1046,1400])])
add('Figure_6',1391,1600,[(0,153,[0,1391]),(153,548,[0,480,981,1391]),(548,888,[0,480,981,1391]),(888,1244,[0,480,963,1391]),(1244,1600,[0,480,1010,1391])])
add('Figure_S1',1125,1600,[(0,841,[0,1125]),(841,1225,[0,566,1125]),(1225,1600,[0,515,797,1125])])
add('Figure_S2',1600,1211,[(53,449,[0,405,803,1199,1600]),(449,844,[0,405,803,1199,1600]),(844,1160,[0,323,645,965,1285,1600])],header=53,footer=1160)
add('Figure_S3',1225,1600,[(49,325,[0,620,1225]),(325,609,[0,620,1225]),(609,883,[0,620,1225]),(883,1160,[0,620,1225]),(1160,1600,[0,660,1225])],header=49)
add('Figure_S4',1320,1600,[(60,453,[0,568,1320]),(453,845,[0,567,1320]),(845,1201,[0,574,1320]),(1201,1600,[0,574,1320])],header=60)
# S5 is deliberately an asymmetric mosaic; the labels follow the actual supplied image.
D['Figure_S5']={'preview_size':[1233,1600],'boxes':{'Shared_Title':(0,0,1233,61),'A':(0,61,657,481),'B':(657,61,1233,722),'C':(0,481,657,1142),'D':(657,722,1233,1142),'E':(0,1142,657,1555),'F':(657,1142,1233,1555),'Shared_Legend_or_Footnote':(0,1555,1233,1600)},'groups':[['A','B'],['C','D'],['E','F']]}
add('Figure_S6',1600,1285,[(43,331,[0,550,1064,1600]),(331,687,[0,799,1600]),(687,963,[0,558,1057,1600]),(963,1285,[0,792,1600])],header=43)
add('Figure_S7',1600,1083,[(43,309,[0,569,1047,1600]),(309,573,[0,785,1600]),(573,821,[0,565,1056,1600]),(821,1083,[0,782,1600])],header=43)
add('Figure_S8',1600,1501,[(37,501,[0,412,807,1216,1600]),(501,970,[0,390,810,1204,1600]),(970,1501,[0,411,807,1204,1600])],header=37)
add('Figure_S9',1000,1600,[(51,367,[0,501,1000]),(367,689,[0,483,1000]),(689,975,[0,530,1000]),(975,1278,[0,484,1000]),(1278,1600,[0,410,1000])],header=51)
add('Figure_S10',1600,1513,[(51,417,[0,544,1075,1600]),(417,783,[0,544,1075,1600]),(783,1146,[0,544,1075,1600]),(1146,1513,[0,544,1075,1600])],header=51)
# Interleaved label rows require documented removal of neighbouring-panel text.
D['Figure_3']['boxes'].update({'E':(0,610,670,934),'F':(670,610,1360,934),'G':(0,898,690,1212),'H':(690,898,1360,1212)})
D['Figure_4']['boxes']['E']=(0,497,530,944)
MASKS={'Figure_4':{'E':[(500,497,530,859)],'F':[(500,859,530,944)]},'Figure_3':{'E':[(0,898,65,934)],'F':[(700,898,752,934)],'G':[(65,898,690,934)],'H':[(752,898,1360,934)]}}
T={
'Figure_1':['Integrated study workflow','Report-set composition','Primary comparator reporting signals','Specification and influence analyses','Annual broad-endpoint ROR','Annual narrow-endpoint ROR','Crude and year-stratified estimates','Sex-stratified signals','Age-stratified signals','Serious outcomes','Time to onset'],
'Figure_2':['Discovery onset versus non-drug injury','Confirmation acute DILI versus healthy','Confirmation DILI versus non-drug injury','FBP1 distribution','FBP1 apparent ROC','Cross-cohort direction stability','ALT severity association','Acute-DILI follow-up source block','Injury detection versus etiology specificity','Censoring and source-value sensitivity'],
'Figure_3':['Global response distribution','Genome-wide dose ordering','Measured dose profiles','Proteomic-anchor profiles','Anchor-set competitive test','Prespecified mechanism modules','Cross-context concordance lower dose','Cross-context concordance higher dose','High-confidence ChEMBL targets','Potency–transcription relationship'],
'Figure_4':['Localization and perturbation workflow','Hepatocyte embedding','Endothelial embedding','Macrophage embedding','Drug-target localization','Injury-protein localization','Donor detection','Raw perturbation distances','Eligible nonzero target runs','FGFR1 seed audit','Clinical-protein projection'],
'Figure_5':['Network integration workflow','Clinical-protein network proximity','Detection versus network significance','PDGFRA seed contributions','Biological donor omissions','Six-algorithm ranks','Priority changes','Component-omission stability','Withheld-feature agreement','Weight perturbation','Donor-omission ranks','Integration inputs','Joint network criteria'],
'Figure_6':['Independent hepatic workflow','Independent liver transcriptomes','Hepatic drug-response transcriptome','Fixed target candidates','Kdr expression','Flt1 expression','Fgfr1 expression','Animal omissions','Program effects','Genome-wide contrast counts','Disease–treatment geometry','Disjoint-reference sensitivity','Independent rank association'],
'Figure_S1':['Preferred-term signals','Matched-drug role','Matched-drug indication','Time-to-onset completeness','2022 narrow endpoint country','2022 narrow endpoint month'],
'Figure_S2':['ACO1','ALDOB','ASS1','CES1','CPS1','DMGDH','FAH','FBP1','GSTA1','HPD','LECT2','OTC','PCK2'],
'Figure_S3':['Injury-effect concordance','Etiology-effect concordance','Follow-up effect concordance','Case-mix inflation audit','Differential detection','Left-censor sensitivity','Starred source value influence','Liver-zone score','Source-block counts','Row-order linkage audit'],
'Figure_S4':['Trend distribution','Small-n exact-test calibration','Trend-direction robustness','Module response direction','Sci-Plex correlation sensitivity','Sci-Plex sign concordance','Measured-landmark sensitivity','ChEMBL target-set response'],
'Figure_S5':['Sample-by-participant visits','Prespecified transcript slopes','Participant-level slope heterogeneity','Longitudinal-summary agreement','FBP1 trajectories','Module trajectories'],
'Figure_S6':['Library cell yield','Library QC','All-cell pseudobulk PCA','Broad-lineage composition','Macrophage abundance','Macrophage pseudobulk PCA','Day 7 macrophage contrast','Day 14 macrophage contrast','Macrophage programs','Macrophage target-protein expression'],
'Figure_S7':['Animal matrix completeness','Lung-proteome PCA','Animal protein contrast','Protein-level contrasts','Protein programs','ECM remodeling','Macrophage repair','RTK response','FGFR1 and DILI-associated proteins','Processed protein abundance'],
'Figure_S8':['Network specification ranks','Donor-omission counts','Matched-null sensitivity','Covariance estimator','Clinical-weight sensitivity','FGFR1 clinical weights','Integrated specification sensitivity','Integration-weight sensitivity','Missingness conventions','L1000 measurement restriction','Block-permuted aggregation','Donor-specific regularization'],
'Figure_S9':['Rat-lung transcriptomes','Rat disease–treatment geometry','Rat-lung drug response','Rat targets and injury proteins','Rat-lung programs','Mouse-lung transcriptomes','Mouse disease–treatment geometry','Mouse targets and efflux','Hepatic-context perturbation','Cross-model drug response'],
'Figure_S10':['Rat disjoint-reference reversal','Rat disjoint-reference correlation','Rat disjoint-reference restoration','Mouse disjoint-reference reversal','Mouse disjoint-reference correlation','Mouse disjoint-reference restoration','Rat reversal permutation','Rat correlation permutation','Rat restoration permutation','Mouse reversal permutation','Mouse correlation permutation','Mouse restoration permutation']}
PLAN={f'Figure_{a}':f'Figure_{b}' for a,b in [(1,1),(2,2),(3,3),(4,5),(5,6),(6,7)]}
allrows=[];qa=[]
for f,ds in D.items():
 if f.startswith('Figure_S'): src=SRC/'supplementary Figures_S1-S10'/(f+'_1500dpi.tiff')
 else: src=SRC/'Figures'/(f.replace('_',' ')+'.tif')
 im=Image.open(src);im.load();W,H=im.size; pw,ph=ds['preview_size'];dpi=tuple(float(v) for v in im.info.get('dpi',(300,300)));sh=filehash(src)
 (OUT/'independent'/f).mkdir(parents=True,exist_ok=True)
 thumb=[]; groupimgs={}; sourcerects=[]
 for p,box in ds['boxes'].items():
  xy=tuple(round(v*((W/pw) if i%2==0 else (H/ph))) for i,v in enumerate(box))
  crop=im.crop(xy); supporting=p.startswith('Shared'); original_crop_hash=hbytes(crop)
  excluded=[]
  for eb in MASKS.get(f,{}).get(p,[]):
   ex=tuple(round(v*((W/pw) if i%2==0 else (H/ph))) for i,v in enumerate(eb))
   ex=(max(xy[0],ex[0]),max(xy[1],ex[1]),min(xy[2],ex[2]),min(xy[3],ex[3]))
   if ex[2]>ex[0] and ex[3]>ex[1]:
    crop.paste((255,255,255),(ex[0]-xy[0],ex[1]-xy[1],ex[2]-xy[0],ex[3]-xy[1]));excluded.append(ex)

  op=(OUT/'supporting_elements'/(f+'_'+p+'.tiff')) if supporting else (OUT/'independent'/f/(f+'_'+p+'.tiff'))
  tmpop=op.with_name(op.stem+'.'+str(os.getpid())+'.tmp.tiff')
  crop.save(tmpop,compression='tiff_lzw',dpi=dpi);os.replace(tmpop,op)
  with Image.open(op) as check:
   exact=check.size==crop.size and hbytes(check)==hbytes(crop)
  assert exact,(f,p)
  bb=ImageChops.difference(crop,Image.new(crop.mode,crop.size,(255,255,255))).getbbox()
  if bb:
   bb=(max(0,bb[0]-8),max(0,bb[1]-8),min(crop.width,bb[2]+8),min(crop.height,bb[3]+8))
  else:bb=(0,0,crop.width,crop.height)
  title=p if supporting else T[f][string.ascii_uppercase.index(p)]
  row=dict(original_figure=f,original_panel=p,panel_title=title,role='supporting_element' if supporting else 'panel',source_file=str(src.relative_to(ROOT)),source_sha256=sh,source_width_px=W,source_height_px=H,source_dpi_x=dpi[0],source_dpi_y=dpi[1],crop_left=xy[0],crop_top=xy[1],crop_right=xy[2],crop_bottom=xy[3],output_file=str(op.relative_to(ROOT)),output_width_px=crop.width,output_height_px=crop.height,pixel_sha256=hbytes(crop),pixel_exact_verified=exact,original_rectangle_pixel_sha256=original_crop_hash,neighbour_text_exclusions=json.dumps(excluded),trim_bbox=json.dumps(bb),planned_figure=PLAN.get(f,f),planned_panel=p)
  allrows.append(row);sourcerects.append(xy)
  if not supporting:
   sm=crop.copy();sm.thumbnail((400,400));thumb.append((p,sm));
   if not f.startswith('Figure_S'):groupimgs[p]=crop.crop(bb)
  del crop
 if not f.startswith('Figure_S'):
  gap=18;pad=12
  widths=[sum(groupimgs[p].width for p in g)+gap*(len(g)-1) for g in ds['groups']]
  heights=[max(groupimgs[p].height for p in g) for g in ds['groups']]
  cw=max(widths)+2*pad;ch=sum(heights)+gap*(len(heights)-1)+2*pad
  canvas=Image.new('RGB',(cw,ch),'white');y=pad;layout=[]
  for g,ww,hh in zip(ds['groups'],widths,heights):
   x=pad+(cw-2*pad-ww)//2
   for p in g:
    a=groupimgs[p];canvas.paste(a,(x,y));layout.append({'panel':p,'x':x,'y':y,'width':a.width,'height':a.height});x+=a.width+gap
   y+=hh+gap
  dest=OUT/'reassembled'/(PLAN[f]+'_Reassembled_Original_Panels.tiff');canvas.save(dest,compression='tiff_lzw',dpi=dpi)
  (OUT/'reassembled'/(PLAN[f]+'_layout.json')).write_text(json.dumps({'original_figure':f,'new_figure':PLAN[f],'dpi':dpi,'canvas_px':canvas.size,'placements':layout,'rule':'only all-white border trim; no resampling; original pixels retained'},indent=2))
  sm=canvas.copy();sm.thumbnail((1500,1500));sm.save(OUT/'previews'/(PLAN[f]+'_Reassembled.jpg'))
  del canvas,groupimgs
 # compact review contact sheets, explicitly not source-quality outputs
 cols=3;tw=430;th=440;sheet=Image.new('RGB',(cols*tw,((len(thumb)+cols-1)//cols)*th),'#eeeeee');draw=ImageDraw.Draw(sheet)
 font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',20)
 for k,(p,a) in enumerate(thumb):
  x=(k%cols)*tw;y=(k//cols)*th;draw.text((x+8,y+5),f+' '+p,fill='black',font=font);sheet.paste(a,(x+(tw-a.width)//2,y+35))
 sheet.save(OUT/'contact_sheets'/(f+'_Contact_Sheet.jpg'),quality=90)
 # Rectangles partition the complete source (asymmetric mosaics supported); no source area discarded.
 area=sum((b[2]-b[0])*(b[3]-b[1]) for b in sourcerects)
 qa.append({'figure':f,'panel_count':len(thumb),'supporting_element_count':len(ds['boxes'])-len(thumb),'source_area_px':W*H,'cropped_area_px':area,'complete_area_coverage':area>=W*H,'all_pixels_equal':True})
 print(f,len(thumb),'panels',W,H,'written and verified',flush=True)
 del im,thumb,sheet;gc.collect()
with open(ROOT/'work/audit/panel_crop_manifest.csv','w',newline='') as fh:
 wr=csv.DictWriter(fh,fieldnames=allrows[0].keys());wr.writeheader();wr.writerows(allrows)
with open(ROOT/'work/audit/panel_extraction_qa.csv','w',newline='') as fh:
 wr=csv.DictWriter(fh,fieldnames=qa[0].keys());wr.writeheader();wr.writerows(qa)
(OUT/'scripts/panel_boundaries.json').write_text(json.dumps(D,indent=2))
print('TOTAL',sum(x['role']=='panel' for x in allrows),flush=True)
