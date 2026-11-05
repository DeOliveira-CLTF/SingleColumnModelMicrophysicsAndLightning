"""Aggregate already grouped flashes with frozen area and explicit coverage."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
from leona.lightning import aggregate
from leona.provenance import ROOT,write_json,sha256

def main():
    c=json.loads((ROOT/'configs/leona.json').read_text())
    if c['area'] is None or c['observation_coverage'] is None:
        raise ValueError('Fixed area and explicit coverage required')
    out=ROOT/'products/leona/observations'
    f=pd.read_csv(out/'flashes.csv.gz')
    f['start']=pd.to_datetime(f.start,utc=True,format='mixed')
    result=aggregate(f,c['start'],c['end'],c['area'],c['observation_coverage'],c['interval_minutes'])
    result.to_csv(out/'counts_15min.csv',index=False)
    components=[result.assign(type='total')]
    for kind in ['IC','CG']:
        components.append(aggregate(f[f.kind==kind],c['start'],c['end'],c['area'],c['observation_coverage'],c['interval_minutes']).assign(type=kind))
    pd.concat(components).to_csv(out/'counts_by_type_15min.csv',index=False)
    m=json.loads((out/'manifest.json').read_text())
    m.update(flag_mapping_confirmed=c['flag_mapping_confirmed'],classification_status='confirmed_by_user',
        flag_confirmation=c['flag_confirmation'],coverage_status='valid_for_full_case_confirmed_by_user',
        coverage_confirmation=c['coverage_confirmation'],area_wsen=c['area'],
        area_selection=c['area_selection'],count_support='first detected event time and location; half-open UTC intervals',
        flashes_sha256=sha256(out/'flashes.csv.gz'),counts_sha256=sha256(out/'counts_15min.csv'))
    write_json(out/'manifest.json',m)
    print(json.dumps({'case_flashes_in_fixed_area':int(result['count'].sum()),'bins':len(result)}))

if __name__=='__main__':main()
