from pathlib import Path
import json,re,html,csv
B=Path(__file__).parent
R=json.loads((B/'reference_verification.json').read_text())
J={
'Acta Pharmacologica Sinica':'Acta Pharmacol. Sin.','Advances in Therapy':'Adv. Ther.','American Journal of Physiology-Gastrointestinal and Liver Physiology':'Am. J. Physiol. Gastrointest. Liver Physiol.','American Journal of Physiology-Lung Cellular and Molecular Physiology':'Am. J. Physiol. Lung Cell. Mol. Physiol.','American Journal of Respiratory Cell and Molecular Biology':'Am. J. Respir. Cell Mol. Biol.','American Journal of Respiratory and Critical Care Medicine':'Am. J. Respir. Crit. Care Med.','Annals of Internal Medicine':'Ann. Intern. Med.','Basic & Clinical Pharmacology & Toxicology':'Basic Clin. Pharmacol. Toxicol.','BioMed Research International':'Biomed Res. Int.','Bioinformatics':'Bioinformatics','British Journal of Pharmacology':'Br. J. Pharmacol.','Cancer Research':'Cancer Res.','Cell':'Cell','Clinical Pharmacokinetics':'Clin. Pharmacokinet.','Communications Biology':'Commun. Biol.','Drug Metabolism and Disposition':'Drug Metab. Dispos.','Drug Metabolism and Pharmacokinetics':'Drug Metab. Pharmacokinet.','Drug Safety':'Drug Saf.','ERJ Open Research':'ERJ Open Res.','European Journal of Pharmaceutics and Biopharmaceutics':'Eur. J. Pharm. Biopharm.','European Respiratory Journal':'Eur. Respir. J.','Frontiers in Pharmacology':'Front. Pharmacol.','Gastroenterology':'Gastroenterology','Gut':'Gut','Hepatology':'Hepatology','International Immunopharmacology':'Int. Immunopharmacol.','JCI Insight':'JCI Insight','Journal of Aerosol Medicine and Pulmonary Drug Delivery':'J. Aerosol Med. Pulm. Drug Deliv.','Journal of Hepatology':'J. Hepatol.','Journal of Medicinal Chemistry':'J. Med. Chem.','Journal of Multivariate Analysis':'J. Multivar. Anal.','Journal of the Royal Statistical Society Series B: Statistical Methodology':'J. R. Stat. Soc. Ser. B Stat. Methodol.','Nature':'Nature','Nature Biotechnology':'Nat. Biotechnol.','Nature Communications':'Nat. Commun.','Nature Human Behaviour':'Nat. Hum. Behav.','Nature Metabolism':'Nat. Metab.','Nature Methods':'Nat. Methods.','Nature Reviews Disease Primers':'Nat. Rev. Dis. Primers','Nature Reviews Drug Discovery':'Nat. Rev. Drug Discov.','New England Journal of Medicine':'N. Engl. J. Med.','Nucleic Acids Research':'Nucleic Acids Res.','PLoS Computational Biology':'PLoS Comput. Biol.','Patterns':'Patterns','Pharmaceuticals':'Pharmaceuticals','Science':'Science','Science Advances':'Sci. Adv.','Scientific Reports':'Sci. Rep.','Statistical Applications in Genetics and Molecular Biology':'Stat. Appl. Genet. Mol. Biol.','The American Journal of Human Genetics':'Am. J. Hum. Genet.','The Annals of Statistics':'Ann. Stat.','The Journal of Clinical Pharmacology':'J. Clin. Pharmacol.','The Lancet Respiratory Medicine':'Lancet Respir. Med.','Translational Gastroenterology and Hepatology':'Transl. Gastroenterol. Hepatol.'}
def clean(s):return html.unescape(re.sub('<[^>]+>','',s)).strip()
def initials(s):
 words=re.findall(r'[\wÀ-ž]+(?:-[\wÀ-ž]+)*',s,re.U);out=[]
 for w in words:
  if '-' in w:out.append('-'.join(p[0].upper()+'.' for p in w.split('-')))
  elif w.isupper():out.append(''.join(x+'.' for x in w))
  else:out.append(w[0].upper()+'.')
 return ''.join(out)
def auth(a):
 if a.get('family'):return clean(a['family'])+', '+initials(a.get('given',''))
 return clean(a.get('name',''))
O=[]
for r in R:
 n=r['number'];notes=[]
 if r['doi'] and n!=93:
  A=[dict(x) for x in r['metadata_authors']]
  if n==75:A[4].update(family='Carcione',given='Claudia')
  authors='; '.join(auth(a) for a in A[:10])+('; et al.' if len(A)>10 else '')
  if n==17:authors='European Association for the Study of the Liver.';notes.append('Preserved corporate author from source citation; Crossref additionally lists guideline contributors.')
  if n==75:notes.append('Crossref inverts the fifth author given/family fields. Source Carcione C retained and full identity Claudia Carcione confirmed on the official Ri.MED staff profile https://www.fondazionerimed.eu/bioresearcher/claudia-carcione/?lang=en and ORCID0000-0002-0305-1781.')
  title=clean(r['metadata_title']);journal=J.get(clean(r['metadata_journal']),clean(r['metadata_journal']));year=str(r['metadata_year']);vol=r['metadata_volume'];page=r['metadata_page']
  if n==43:page='Article32';notes.append('Pagination completed from PubMed PMID16646851; Crossref omits article number.')
  if n==52:page='1–26';notes.append('Crossref omits pagination; original1–26 retained, supported by publisher/primary bibliographic citation.')
  if n==20 and page=='9-9':page='9'
  page=page.replace('-','–')
  text=f'{authors} {title.rstrip(".")}. {journal} {year}, {vol}, {page}. https://doi.org/{r["metadata_doi"]}.'
  seg=[{'text':authors+' '+title.rstrip('.')+'. '},{'text':journal,'italic':True},{'text':' '},{'text':year,'bold':True},{'text':', '},{'text':vol,'italic':True},{'text':', '+page+'. https://doi.org/'+r['metadata_doi']+'.'}]
 elif n==93:
  authors='Lauer, D.; Gote-Schniering, J.';text=authors+' Radioproteomics Stratifies Molecular Response to Antifibrotic Treatment in Pulmonary Fibrosis; Version 1; Zenodo: 2024. https://doi.org/10.5281/zenodo.11395642.';seg=[{'text':text}]
 else:
  t=r['text'];u=re.search(r'https?://\S+',t)[0].rstrip('.,;');prefix=t[:t.index(u)].strip();prefix=re.sub(r'\s*Available from:$|\s*Available online:$','',prefix)
  prefix=re.sub(r'\s*\[Internet\]','',prefix);prefix=re.sub(r'\s*\[cited [^\]]+\]','',prefix)
  # Preserve frozen access dates rather than recasting as new retrieval dates.
  if n==84:access='19 September 2026'
  elif n in [8,9,10,14]:access='5 September 2026'
  else:access='4 September 2026'
  # Keep original corporate title and remove bibliographic publication-place clauses for websites.
  if n in [8,9]:prefix=t.split(' [Internet]')[0]+'.'
  if n==10:prefix='European Medicines Agency. Jascayd: EPAR—Medicine Overview.'
  if n==82:u='https://maayanlab.cloud/Harmonizome/resource/Sci-Plex';access='2 October 2026'
  if n==97:u='https://string-db.org/help/api/';access='2 October 2026'
  if n==14:prefix='Boehringer Ingelheim Pharmaceuticals, Inc. OFEV (Nintedanib) Capsules, for Oral Use: Full Prescribing Information.'
  text=prefix.rstrip('.;')+'. Available online: '+u+' (accessed on '+access+').';seg=[{'text':text}]
 O.append({'number':n,'doi':r['doi'],'text':text,'segments':seg,'notes':notes})
(B/'References_JCM_ACS.json').write_text(json.dumps(O,ensure_ascii=False,indent=2))
(B/'References_JCM_ACS.txt').write_text('\n\n'.join(f"{x['number']}. {x['text']}" for x in O))
print('Formatted',len(O))
