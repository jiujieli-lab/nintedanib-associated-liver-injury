#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from PIL import Image
OUT=Path(__file__).resolve().parent
NET=OUT.parent/'data/Additional_file_4_Clinical_Evidence_and_Analysis_Code/Analysis/network_diffusion'
g=pd.read_csv(OUT/'gene_level_perturbations.csv.gz');o=pd.read_csv(NET/'network_specification_scores.csv');o=o[o.is_primary_spec&o.weight_profile.eq('discovery_DO_vs_HV')]
x=g[g.intervention.eq('sham')].merge(o,left_on=['compartment','gene'],right_on=['compartment','gene_symbol'])
baseline_error=float((x.baseline_score-x.diffusion_score).abs().max());assert baseline_error<1e-12 and len(x)==50
m=pd.read_csv(OUT/'perturbation_metrics.csv');d=pd.read_csv(OUT/'donor_target_expression.csv')
r=m[m.omitted_donor.ne('none')&m.intervention.isin(['PDGFRA','FLT4','PDGFRA+FLT4'])].groupby(['compartment','intervention']).clinical_weighted_relative_change.agg(['min','max']);r.to_csv(OUT/'donor_omission_ranges.csv')
assert (d.groupby('compartment').n_cells.sum()==np.array([1688,7002,2384])).all()  # each compartment is represented once per target
assert m[m.intervention.isin(['sham','restored_graph'])].max_abs_gene_change.eq(0).all()
cal=pd.read_csv(OUT/'primary_calibration.csv');assert len(cal)==9 and cal.q_bh_9.ge(.05).all()
s=json.loads((OUT/'validation_summary.json').read_text());s['baseline_reproduction_against_archived_primary_graph_max_abs_difference']=baseline_error;s['n_baseline_target_compartment_scores_reproduced']=len(x)
Image.MAX_IMAGE_PIXELS=None
with Image.open(OUT/'Figure_S11_Virtual_Perturbation_1200dpi.tiff') as im:
 assert im.size==(10080,10080) and im.info['dpi']==(1200.,1200.)
 s['figure_size_pixels']=im.size;s['figure_dpi']=tuple(map(float,im.info['dpi']))
s['individual_panels']=len(list((OUT/'Individual_Panels').glob('*.png')));assert s['individual_panels']==9
(OUT/'validation_summary.json').write_text(json.dumps(s,indent=2))
rows=[]
for p in sorted(OUT.rglob('*')):
 if p.is_file() and p.name not in ['output_manifest.csv','run.log']:
  rows.append({'file':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(rows).to_csv(OUT/'output_manifest.csv',index=False)
print(json.dumps(s,indent=2))
