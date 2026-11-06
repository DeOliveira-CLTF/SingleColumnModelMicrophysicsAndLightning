from scm_lightning.config import load_config
"""Aggregate already grouped flashes with frozen area and explicit coverage."""
from pathlib import Path
import argparse
import json
import sys
import pandas as pd
from scm_lightning.observations import aggregate
from scm_lightning.provenance import ROOT,write_json,sha256

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',default=str(ROOT/'configs/cases/argentina_20181213.yaml'));args=ap.parse_args()
    c=load_config(args.config)
    if c['area'] is None or c['observation_coverage'] is None:
        raise ValueError('Fixed area and explicit coverage required')
    out=ROOT/c.get('paths',{}).get('results','results/argentina_20181213')/'observations'
    f=pd.read_csv(out/'flashes.csv.gz')
    f['start']=pd.to_datetime(f.start,utc=True,format='mixed')
    coverage=c['observation_coverage']
    if isinstance(coverage,str) and coverage=='full_case':
        coverage=[True]*len(pd.date_range(c['start'],c['end'],freq=f"{c['interval_minutes']}min",inclusive='left'))
    result=aggregate(f,c['start'],c['end'],c['area'],coverage,c['interval_minutes'])
    result.to_csv(out/'counts_15min.csv',index=False)
    aggregate(f,c['start'],c['end'],c['area'],[v for v in coverage for _ in range(3)],5).to_csv(out/'counts_5min.csv',index=False)
    components=[result.assign(type='total')]
    for kind in ['IC','CG']:
        components.append(aggregate(f[f.kind==kind],c['start'],c['end'],c['area'],coverage,c['interval_minutes']).assign(type=kind))
    pd.concat(components).to_csv(out/'counts_by_type_15min.csv',index=False)
    m=json.loads((out/'manifest.json').read_text())
    m.update(flag_mapping_confirmed=c['flag_mapping_confirmed'],classification_status='confirmed_by_user',
        flag_confirmation=c.get('flag_confirmation','same confirmed export convention'),coverage_status='valid_for_full_case_confirmed_by_user',
        coverage_confirmation=c['coverage_confirmation'],area_wsen=c['area'],
        area_selection=c['area_selection'],count_support='first detected event time and location; half-open UTC intervals',
        flashes_sha256=sha256(out/'flashes.csv.gz'),counts_sha256=sha256(out/'counts_15min.csv'))
    write_json(out/'manifest.json',m)
    print(json.dumps({'case_flashes_in_fixed_area':int(result['count'].sum()),'bins':len(result)}))

if __name__=='__main__':main()
