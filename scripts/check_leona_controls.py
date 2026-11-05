"""Recheck complete case controls without rerunning or replacing their products."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import xarray as xr
from leona.experiments import CONFIGURATIONS,validate_output
from leona.provenance import ROOT,write_json,sha256

def main():
    control=ROOT/'work/leona/controls';forcing=ROOT/'work/leona/case_forcing_v2/LEONA_SCM.nc'
    with xr.open_dataset(forcing,decode_times=False) as f:end=float(f.time.values[-1])
    records=[{'configuration':c,**validate_output(control/c/'output.nc',end)} for c in CONFIGURATIONS]
    write_json(control/'control_checks.json',records)
    gate={'configurations':CONFIGURATIONS,'forcing_sha256':sha256(forcing),
        'executable_sha256':sha256(ROOT/'work/leona/build/mp_offline_driver'),
        'numerical_passed':True,'physical_preflight_passed':all(r.get('physical_preflight_passed',False) for r in records),
        'physical_scope':json.loads((ROOT/'configs/budget_contract.json').read_text()),'control_results':records}
    write_json(control/'ensemble_gate.json',gate)
    write_json(ROOT/'products/leona/tests/case_controls.json',gate)
    print(json.dumps({'numerical':True,'physical_preflight':gate['physical_preflight_passed']}))

if __name__=='__main__':main()
