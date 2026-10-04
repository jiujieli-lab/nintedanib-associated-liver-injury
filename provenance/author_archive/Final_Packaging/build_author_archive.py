from pathlib import Path
import shutil,json,hashlib
R=Path(__file__).resolve().parents[1];BASE=R.parent.parent;A=R/'author_archive';A.mkdir(exist_ok=True)
def tree(src,dst,extensions=None):
 for p in src.rglob('*'):
  if not p.is_file() or any(x.startswith('.') or x=='__pycache__' for x in p.relative_to(src).parts):continue
  if extensions and p.suffix.lower() not in extensions:continue
  q=dst/p.relative_to(src);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
tree(BASE/'deliverables/Reproducibility',A/'Source_Reproducibility')
tree(BASE/'deliverables/Independent_Panels',A/'Independent_Source_Panels')
tree(R/'journal_audit',A/'JCM_Reference_and_Method_Audit',extensions={'.json','.py','.txt','.md','.tsv','.html'})
tree(R/'manuscript',A/'Manuscript_Assembly',extensions={'.json','.py','.txt','.md'})
tree(R/'supplementary_information',A/'Supplementary_Information_Audit',extensions={'.json','.py'})
tree(R/'proto/manuscript_audit',A/'Independent_Review',extensions={'.json','.md','.csv'})
tree(R/'tables',A/'Table_Assembly',extensions={'.json','.mjs','.ndjson'})
for f in (R/'figures').glob('*'):
 if f.is_file() and f.suffix in ('.json','.md','.py'):
  q=A/'Figure_Assembly'/f.name;q.parent.mkdir(exist_ok=True);shutil.copy2(f,q)
tree(R/'figures/manuscript_qa',A/'Figure_Assembly/Manuscript_QA',extensions={'.json','.md'})
for f in (R/'figures/screening').glob('*'):
 if f.is_file() and f.suffix in ('.json','.md','.py','.txt'):
  q=A/'Figure_Assembly/Screening_Figures'/f.name;q.parent.mkdir(exist_ok=True);shutil.copy2(f,q)
for name in ['MD_Feasibility_and_Prospective_Protocol.md','MD_Resource_Estimate.json','Resource_and_Installation_Audit.json','MD_Audit_SHA256.json','benchmark_cpu.py']:
 p=R/'md'/name;q=A/'MD_Feasibility_Only'/name;q.parent.mkdir(exist_ok=True);shutil.copy2(p,q)
tree(R/'md/benchmark',A/'MD_Feasibility_Only/benchmark')
tree(R/'proto',A/'Plugin_Execution_Provenance',extensions={'.json','.md'})
q=A/'Figure4_Repair_Audit';q.mkdir(exist_ok=True);shutil.copy2(BASE/'figure4_repair/validation.json',q/'validation.json')
tree(R/'packaging',A/'Final_Packaging',extensions={'.py','.json','.txt'})
(A/'README.txt').write_text('''AUTHOR REPRODUCIBILITY ARCHIVE

This author-facing archive preserves the original source-analysis packages, the165 independent source panels, targeted-analysis panels, reference verification, paragraph-disposition records, assembly code and completed quality checks. It supplements the journal-facing JCM_Submission_Package.zip and the separate high-resolution artwork archive. It is not a substitute for the numbered supplementary files.

Current journal figure numbering is defined by the manuscript and the source crosswalk on page93 of Additionalfile1. Raw source-package filenames and source-analysis identifiers are deliberately preserved. They may differ from final figure labels. The exact89-row workbook crosswalk is in Additionalfile9/source_indices and in this archive's Independent_Review folder.

The MD_Feasibility_Only folder contains a ligand-free, sub-picosecond CPU feasibility benchmark and a prospective production protocol. No target–ligand production trajectory, convergence result, binding free energy or drug-stability conclusion is claimed. Receptor reconstruction, validated ligand parameterization and an appropriate computing backend are still required for the proposed production simulations.

The submitted scientific text must be the version approved by all authors. The author checklist identifies unresolved author-controlled attestations and the need for an accessible repository location if large source archives cannot be uploaded with the journal files. Never upload this entire author archive as part of the120MB portal set without checking the combined size.

Archived audit documents describe their own analysis or assembly dates; final manuscript, supplementary files and final file manifest take precedence for publication numbering and delivery status. No software environment caches, API credentials or account tokens are included.
''')
print(A)
