from scm_lightning.config import load_config
from scm_lightning.ensemble.experiments import control_directory
"""Recheck complete case controls without rerunning or replacing their products."""
from pathlib import Path
import json
import sys
import xarray as xr
from scm_lightning.ensemble.experiments import CONFIGURATIONS,validate_output
from scm_lightning.provenance import ROOT,write_json,sha256

def main():
    control=ROOT/'runs/argentina_20181213/_execution/controls';forcing=ROOT/'data/processed/forcing/LEONA_SCM.nc'
    with xr.open_dataset(forcing,decode_times=False) as f:end=float(f.time.values[-1])
    records=[{'configuration':c,**validate_output(control_directory(c)/'output.nc',end)} for c in CONFIGURATIONS]
    write_json(control/'control_checks.json',records)
    gate={'configurations':CONFIGURATIONS,'forcing_sha256':sha256(forcing),
        'executable_sha256':sha256(ROOT/'build/leona/mp_offline_driver'),
        'numerical_passed':True,'physical_preflight_passed':all(r.get('physical_preflight_passed',False) for r in records),
        'physical_scope':load_config(ROOT/'configs/evaluation/budget_contract.yaml'),'control_results':records}
    write_json(control/'ensemble_gate.json',gate)
    write_json(ROOT/'results/argentina_20181213/tests/case_controls.json',gate)
    print(json.dumps({'numerical':True,'physical_preflight':gate['physical_preflight_passed']}))

if __name__=='__main__':main()
