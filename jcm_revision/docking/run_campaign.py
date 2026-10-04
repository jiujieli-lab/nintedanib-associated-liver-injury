"""Bounded local CPU campaign; first seed across all compounds before repeats."""
import concurrent.futures,csv,json,os,subprocess,sys,time
from pathlib import Path
P=Path(__file__).resolve().parent
LIB=P.parent/'candidate_library'/'Candidate_Screening_Library.csv'
SEEDS=[104729,130363,155921]
rows=list(csv.DictReader(LIB.open()))
jobs=[(seed,row['name'],target) for seed in SEEDS for row in rows for target in ['PDGFRA','FLT4']]
(P/'logs').mkdir(exist_ok=True)
config={'n_compounds':len(rows),'targets':['PDGFRA','FLT4'],'seeds':SEEDS,'exhaustiveness':32,'concurrent_processes':4,'threads_per_process':2,'n_jobs':len(jobs),'interpretation':'Three stochastic searches per compound/target; these are not biological or MD replicates.'}
(P/'campaign_config.json').write_text(json.dumps(config,indent=2))
def run(job):
    seed,name,target=job;log=P/'logs'/f'{target}_{name}_{seed}.log'
    command=[sys.executable,str(P/'run_docking.py'),'--library',str(LIB),'--only-name',name,'--target',target,'--seed',str(seed),'--exhaustiveness','32','--cpu','2']
    with log.open('w') as out:
        result=subprocess.run(command,stdout=out,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'1'})
    return {'target':target,'name':name,'seed':seed,'returncode':result.returncode,'log':str(log.relative_to(P))}
done=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for result in pool.map(run,jobs):
        done.append(result);(P/'campaign_progress.json').write_text(json.dumps({'completed':len(done),'total':len(jobs),'results':done},indent=2));print(json.dumps(result),flush=True)
print('COMPLETE',flush=True)
