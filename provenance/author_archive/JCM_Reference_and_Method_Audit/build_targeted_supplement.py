from pathlib import Path
import re,json,hashlib
from docx import Document
from docx.shared import Pt,Inches,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
B=Path(__file__).resolve().parent
ROOT=B.parent
src=Path('/workspace/scratch/c5634a2d588f/deliverables/Supplementary_Materials/Additional_file_4_Targeted_Analyses.docx')
source=Document(src);P={i:p.text for i,p in enumerate(source.paragraphs)}
mapfile=ROOT/'manuscript/Final_Reference_Number_Map.json';mp={int(k):v for k,v in json.loads(mapfile.read_text()).items()}
refs={x['number']:x for x in json.loads((ROOT/'manuscript/Final_References.json').read_text())}
local_to_old={1:54,2:56,3:57}
# Confirm local bibliography identities before applying the main bibliography map.
for i,doi in [(47,'10.1038/nbt.1990'),(48,'10.1016/j.cell.2026.07.052'),(49,'10.1016/j.cell.2026.08.025')]:assert doi in P[i]
# Existing supplement references are local; moved Results already use main-source numbering.
def citations(t,local=False):
 def rep(m):
  nums=[]
  for s in m[1].split(','):
   if re.fullmatch(r'\s*\d+\s*',s):nums.append(int(s))
   else:
    a,z=re.split('[–-]',s);nums.extend(range(int(a),int(z)+1))
  vals=[mp[local_to_old[n] if local else n] for n in nums]
  return '['+', '.join(map(str,vals))+']'
 return re.sub(r'\[([\d, –-]+)\]',rep,t)
def clean(t):
 pairs=[('reanalysed the released ChEMBL','reanalysed the ChEMBL'),('The original exact-relation','The specified exact-relation'),('the frozen integrated ranking','the fixed integrated ranking'),('missing kinetic values were not imputed and no conventional molecular docking was performed','missing kinetic values were not imputed'),('The frozen maximum','The primary maximum'),('released ChEMBL activity records and a fresh ChEMBL molecule record','ChEMBL activity records and the ChEMBL molecule record'),(' after live model discovery',''),('We extended the inherited donor-balanced networks with an explicitly model-conditional intervention analysis','We evaluated donor-balanced networks using a model-conditional intervention analysis'),('inherited gene universe','fixed gene universe'),('the inherited 10% expression threshold','the specified 10% expression threshold'),('These extension settings','These settings'),('The extension therefore defines','This analysis defines'),('while retaining the original interpretation of','and supports the interpretation of'),('the earlier favorable percentiles','the percentile-normalized scores'),('The earlier favorable percentiles','The percentile-normalized scores'),('None of the five revised methods','None of the five integration methods'),('We retained similarly informative negative findings from more detailed pulmonary measurements.','Detailed pulmonary measurements characterized the cellular and protein context of treatment response.'),('required to reproduce the additions','required to reproduce the analyses'),('three original normalized compartment matrices','three normalized compartment matrices'),('each completed analysis','each analysis'),('The threshold was retained rather than relaxed to recover a preferred molecule.','The specified threshold was applied to all candidates.')]
 for a,b in pairs:t=t.replace(a,b)
 return t
textfile=B/'content_base/Supplementary_Results_Integration.txt'
blocks=[s.strip() for s in textfile.read_text().split('\n\n') if s.strip()]
headings=['Liver-cell localization and numerical calibration','Network decomposition and sensitivity analyses','Single and joint PDGFRA and FLT4 network interventions','Pulmonary cellular and protein context','Supplementary Figure S12. Liver localization and numerical auditing calibrate virtual perturbation evidence']
groups={};key=None
for t in blocks:
 if t in headings:key=t;groups[key]=[]
 else:assert key;groups[key].append(t)
assert [len(groups[k]) for k in headings]==[3,5,1,1,1]
d=Document()
sec=d.sections[0];sec.page_width=Inches(8.5);sec.page_height=Inches(11)
sec.top_margin=Inches(.75);sec.bottom_margin=Inches(.75);sec.left_margin=Inches(.8);sec.right_margin=Inches(.8)
for sty in ['Normal','Title','Heading 1','Heading 2']:
 s=d.styles[sty];s.font.name='Times New Roman';s.font.color.rgb=RGBColor(0,0,0)
 s.paragraph_format.space_after=Pt(6)
 if sty=='Normal':s.font.size=Pt(11);s.paragraph_format.line_spacing=1.12
 if sty=='Title':s.font.size=Pt(16);s.font.bold=True
 if sty=='Heading 1':s.font.size=Pt(12);s.font.bold=True;s.paragraph_format.space_before=Pt(10)
 if sty=='Heading 2':s.font.size=Pt(11);s.font.bold=True;s.paragraph_format.space_before=Pt(8)
# Remove theme font overrides and the default Title border.
for element in d.styles.element.xpath('.//w:rFonts'):
 for attr in list(element.attrib):
  if 'theme' in attr.lower():del element.attrib[attr]
for element in d.styles.element.xpath('.//w:pBdr'):
 element.getparent().remove(element)
pat=re.compile(r'10([−-]\d+)|(?<![A-Za-z])([KI])([di])\b|\b(IC)(50)\b|\b(log)(2|10)\b|\b(pK)(a)\b')
def para(t,style='Normal',source_id=None):
 p=d.add_paragraph(style=style);pos=0
 for m in pat.finditer(t):
  p.add_run(t[pos:m.start()])
  if m[1]:p.add_run('10');p.add_run(m[1]).font.superscript=True
  else:
   g=m.groups();a=next(i for i in (1,3,5,7) if g[i] is not None);p.add_run(g[a]);p.add_run(g[a+1]).font.subscript=True
  pos=m.end()
 p.add_run(t[pos:]);p.paragraph_format.widow_control=True
 if source_id is not None:disposition.append({'source':source_id,'output_paragraph':len(d.paragraphs)-1,'output_text':t})
 return p
num=0;disposition=[]
def heading(t):
 global num
 num+=1;para(f'S{num}. {t}','Heading 1')
def old(i):return para(clean(citations(P[i],True)),source_id=f'existing:{i}')
def moved(k):
 for j,t in enumerate(groups[k]):para(clean(citations(t)),source_id=f'moved:{k}:{j+1}')
para('Targeted biochemical and network validation','Title')
para(P[1])
para('This supplement describes endpoint-resolved biochemical evidence, molecular-property predictions and network interventions for PDGFRA and FLT4, together with cellular localization and sensitivity analyses. Figure 4 presents biochemical binding records; Supplementary Figures S11 and S12 present target-pair interventions and numerical calibration. Tables N1–N17 are supplied in Additional file 5, and complete gene-level and matched-null outputs N18–N19 accompany Additional file 6. Reference numbers follow the main manuscript bibliography.')
heading('Assay stratified nintedanib target engagement');old(4)
heading('Biochemical engagement and assay dependent selectivity');old(6);old(7);old(12)
heading('Molecular property prediction methods');old(9)
heading('Molecular property prediction results');old(11)
heading('Single and dual target network intervention methods');old(14);old(15);old(16)
heading('Target pair intervention results');moved(headings[2])
para('The intervention magnitudes are changes in network scores, not predicted reductions in injury. The results provide no additional evidence for biological synergy, receptor-to-receptor regulation or causal mediation of liver injury.',source_id='existing:18 unique qualifications; quantitative content retained in moved paragraph')
old(19);old(20)
heading('Liver cell localization and numerical calibration');moved(headings[0])
heading('Network decomposition and sensitivity analyses');moved(headings[1])
heading('Pulmonary cellular and protein context');moved(headings[3])
heading('Supplementary figure legends')
para(P[22],'Heading 2');old(23)
para(headings[4],'Heading 2');moved(headings[4])
heading('Supplementary table and source data guide')
for i in range(25,44):old(i)
heading('Reproducibility');old(45)
heading('References')
for n in sorted(mp[x] for x in local_to_old.values()):
 r=refs[n];p=d.add_paragraph();p.add_run(f'{n}. ')
 for seg in r['segments']:
  z=p.add_run(seg['text']);z.bold=seg.get('bold',False);z.italic=seg.get('italic',False)
# Page field in the footer.
f=sec.footer.paragraphs[0];f.alignment=2;f.add_run('Page ')
fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');f._p.append(fld)
out=ROOT/'submission/Supplementary_Materials/Additional_file_4_Targeted_Analyses.docx';out.parent.mkdir(parents=True,exist_ok=True);d.save(out)
# Validate complete moved paragraphs after citation/neutral editorial normalization.
outtext='\n'.join(p.text for p in d.paragraphs)
for h in headings:
 for t in groups[h]:assert clean(citations(t)) in outtext,(h,t[:50])
assert 'References are numbered independently' not in outtext
assert '[1]' not in outtext and '[2]' not in outtext and '[3]' not in outtext
for x in ['[73]','[75]','[76]','Supplementary Figure S11','Supplementary Figure S12','95731','70696','4.07×10−20','9.71×10−17']:assert x in outtext,x
report={'source_file':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'output':str(out),'final_reference_map':str(mapfile),'final_reference_map_sha256':hashlib.sha256(mapfile.read_bytes()).hexdigest(),'local_reference_mapping':{k:mp[v] for k,v in local_to_old.items()},'moved_results_paragraphs':10,'moved_S12_legend_paragraphs':1,'all_moved_text_present_after_documented_editorial_and_reference_normalization':True,'source_paragraph_18':'Quantitative content merged with complete moved target-pair paragraph; unique causal/score interpretation retained separately.','source_content_preserved':True,'words':len(outtext.split()),'section_count':num,'render_status':'pending'}
(B/'Targeted_Supplement_Build_Report.json').write_text(json.dumps(report,indent=2))
(B/'Targeted_Supplement_Disposition.json').write_text(json.dumps(disposition,ensure_ascii=False,indent=2))
(B/'Targeted_Supplement_Text.txt').write_text(outtext)
print(json.dumps(report,indent=2))
