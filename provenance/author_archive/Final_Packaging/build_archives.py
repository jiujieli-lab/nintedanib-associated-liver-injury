from pathlib import Path
import shutil,json,csv,hashlib,zipfile,re
R=Path(__file__).resolve().parents[1];W=R.parent;BASE=W.parent
S=R/'submission/Supplementary_Materials';S.mkdir(exist_ok=True,parents=True)
for n,name in [(3,'Supplementary_Tables.xlsx'),(5,'Targeted_Analysis_Tables.xlsx'),(6,'Targeted_Analysis_Code.zip')]:
 fn=f'Additional_file_{n}_{name}';shutil.copy2(BASE/'deliverables/Supplementary_Materials'/fn,S/fn)
shutil.copy2(R/'tables/output/Additional_file_8_Screening_Tables.xlsx',S/'Additional_file_8_Screening_Tables.xlsx')
P=R/'screening_archive';P.mkdir(exist_ok=True)
def copytree(src,dst,omit=()):
 for f in src.rglob('*'):
  if not f.is_file() or any(x in omit or x.startswith('.') or x=='__pycache__' for x in f.relative_to(src).parts):continue
  out=dst/f.relative_to(src);out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,out)
for d in ['candidate_library','screening','boltz-experiments','docking']:
 copytree(R/d,P/'jcm_revision'/d,('env','manuscript_audit'))
proto=P/'jcm_revision/proto';proto.mkdir(exist_ok=True)
for name in ['AF-P35916-F1-model_v6.pdb','FLT4_Proto_Boltz2_Affinity_NOT_SUBMITTED.json','FLT4_alphafold_fetch_result.json','PDGFRA_alphafold_fetch_result.json','AlphaFold_domain_confidence.json','FLT4_AF_kinase_domain.pdb','PDGFRA_AF_kinase_domain.pdb','AF-P16234-F1-model_v6.pdb','Proto_execution_audit.md','workspace_status_latest.json']:
 shutil.copy2(R/'proto'/name,proto/name)
qa=P/'jcm_revision/md';qa.mkdir(exist_ok=True)
for f in (R/'md').glob('Boltz*'):
 if f.is_file():shutil.copy2(f,qa/f.name)
raw=P/'binding_analysis';raw.mkdir(exist_ok=True)
shutil.copy2(W/'binding_analysis/Inductive_Bio_Raw_Predictions.json',raw/'Inductive_Bio_Raw_Predictions.json')
idx=P/'source_indices';idx.mkdir(exist_ok=True)
shutil.copy2(R/'proto/manuscript_audit/Additional3_Final_Figure_Index_Map.csv',idx/'Additional3_Final_Figure_Index_Map.csv')
shutil.copy2(R/'supplementary_information/Source_Figure_Crosswalk.pdf',idx/'Source_Figure_Crosswalk.pdf')
# Quantitative figure source code, data and coordinate checks; rendered artwork is supplied separately.
for f in (R/'figures/screening').rglob('*'):
 if f.is_file() and f.suffix in ('.py','.json','.csv','.txt','.md') and not any(x.startswith('.') or x=='__pycache__' for x in f.parts):
  q=P/'jcm_revision/figures/screening'/f.relative_to(R/'figures/screening');q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,q)
for f in (R/'tables').glob('*.mjs'):
 q=P/'jcm_revision/tables'/f.name;q.parent.mkdir(exist_ok=True);shutil.copy2(f,q)
for name in ['Table_Manifest.json','Workbook_QA.json']:
 shutil.copy2(R/'tables/output'/name,P/'jcm_revision/tables'/name)
(P/'README.txt').write_text('''FOCUSED COMPOUND SCREENING: SOURCE DATA AND CODE

This archive contains completed computational outputs for 16 compounds against PDGFRA and FLT4, with all scored, filtered and retained results. There are 96 primary Vina searches, three separately specified apo-PDGFRA sensitivity searches and one separate crystallographic redocking control. The three seeds are computational searches, not biological replicates. No target–ligand production molecular-dynamics trajectories are included or claimed.

Directory layout
source_indices: complete89-row mapping of analysis-source figure IDs to final manuscript panels.
jcm_revision/candidate_library: exact compound identities, label sources, all191 experimental records, exclusions, endpoint summaries and the source-matched Kd panel.
jcm_revision/boltz-experiments: both completed provider jobs, all24 exported pocket models and raw metrics. Filtered candidates have no returned model score. Exported coordinates are cropped pocket structures.
jcm_revision/screening: complete32-input disposition, 48 molecular-property predictions, exact provider responses and summarization scripts.
jcm_revision/docking: receptor structures, preparation, grids, molecular inputs, all search logs/poses and terminal validation. Python environment binaries are omitted; resolved software versions are preserved in methods and manifests.
jcm_revision/proto: actual AlphaFold retrieval and tool-execution provenance. No unexecuted Proto prediction is represented as a result.
jcm_revision/md: Boltz coordinate/contact quality checks only; the directory name does not indicate an MD simulation.
jcm_revision/figures/screening: plotting scripts and underlying numerical/coordinate records for Figure7 and FiguresS13–S14.
binding_analysis/Inductive_Bio_Raw_Predictions.json: the matched-version raw nintedanib predictions used with the15 further compounds.

Reproduction
Run `python verify_archive.py` from this directory for an offline integrity and count check (Python standard library only). All SHA-256 entries must match. Analysis scripts use paths relative to their own locations unless an archived runtime path is explicitly retained as provenance. Summarization requires RDKit2026.03.6 and the standard numerical stack. Docking used AutoDockVina1.2.7, Meeko0.8.0 and RDKit2026.03.6. See docking/Docking_Methods_and_Limits.md and Additionalfile7 for exact preparation, seeds, grid and exhaustiveness. Rerunning provider inference requires the respective accounts and may incur new charges; no such rerun is necessary to inspect these archived results. Do not execute the campaign launcher merely to review results.

Interpretation
Binding confidence is a target-specific model output, not a measured dissociation constant. The returned_smiles_identity_verified field checks returned molecular metadata against input stereospecific InChIKeys; it is not an independent check of coordinate chirality. Contact distances do not prove hydrogen bonds. Rigid docking scores do not measure affinity, clinical benefit or hepatic safety. Source record inequalities and duplicate flags are retained.
''')
(P/'verify_archive.py').write_text('''from pathlib import Path
import hashlib,csv,json
r=Path(__file__).resolve().parent
for line in (r/'SHA256SUMS.txt').read_text().splitlines():
 digest,name=line.split('  ',1)
 assert hashlib.sha256((r/name).read_bytes()).hexdigest()==digest,name
counts={'screening/Boltz_All32_Dispositions.csv':32,'screening/Inductive_Properties_All16.csv':48,'docking/Docking_All_Completed_Runs.csv':96,'docking/Docking_Summary.csv':32,'candidate_library/All_Target_Activity_Disposition.csv':191,'candidate_library/Davis2011_Matched_Kd_Benchmark.csv':12}
for name,n in counts.items():
 assert len(list(csv.DictReader((r/'jcm_revision'/name).open())))==n,(name,n)
for target in ('pdgfra','flt4'):
 d=json.loads((r/f'jcm_revision/boltz-experiments/jcm-{target}-focused16-v1/run.json').read_text())
 assert d['status']=='succeeded' and d['progress']['num_molecules_screened']==12
print('PASS: all file hashes and primary analysis counts match.')
''')
lines=[]
for f in sorted(P.rglob('*')):
 if f.is_file() and f.name!='SHA256SUMS.txt':lines.append(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+str(f.relative_to(P)))
(P/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n')
zpath=S/'Additional_file_9_Screening_Code_and_Data.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(P.rglob('*')):
  if f.is_file():z.write(f,f.relative_to(P))
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
print(json.dumps({'additional9_files':len(lines)+1,'size_MB':zpath.stat().st_size/1e6},indent=2))
