"""Wait for primary campaign completion, then execute three separately reported sensitivity searches."""
import concurrent.futures,json,os,subprocess,sys,time
from pathlib import Path
P=Path(__file__).resolve().parent
while True:
    path=P/'campaign_progress.json'
    try:
        data=json.loads(path.read_text())
        if data['completed']==data['total']:
            if any(r['returncode']!=0 for r in data['results']):raise RuntimeError('Resolve failed primary jobs first')
            break
    except FileNotFoundError:pass
    except json.JSONDecodeError:pass
    time.sleep(15)
def run(seed):
    out=P/'logs'/f'PDGFRA_APO_8PQJ_nintedanib_{seed}.log'
    cmd=[sys.executable,str(P/'run_docking.py'),'--library',str(P.parent/'candidate_library'/'Candidate_Screening_Library.csv'),'--only-name','nintedanib','--target','PDGFRA_APO_8PQJ','--seed',str(seed),'--cpu','2','--exhaustiveness','32']
    with out.open('w') as log:
        r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'1'})
    print(json.dumps({'seed':seed,'returncode':r.returncode}),flush=True)
    return {'seed':seed,'returncode':r.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as e:results=list(e.map(run,[104729,130363,155921]))
(P/'8PQJ_sensitivity_execution_status.json').write_text(json.dumps(results,indent=2))
