from pathlib import Path
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
