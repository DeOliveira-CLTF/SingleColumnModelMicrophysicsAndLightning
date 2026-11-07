"""Versioned MYNN surface/plume/coupling audit; separate deterministic runs.

Build isolated binaries with ``make BUILD=build/mynn_audit_vN``. Stages:
run (original instrumentation), corrected (opt-in proven boundary fixes),
reference-v2 and no-shadow-v3 (identity checks), predictors (uncalibrated),
analyze (CSV/three PNGs), finalize (JSON and <=15 terminal lines).
Existing run identities deliberately reject changed binaries or forcings.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import threading
import numpy as np
import pandas as pd
import xarray as xr
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.ensemble.experiments import CONFIGURATIONS,control_directory,namelist,execute,validate_output

RUNS=ROOT/'runs/argentina_20181213/mynn_audit_v1'
OUT=ROOT/'results/argentina_20181213/mynn_audit_v1'
FORCING=ROOT/'data/processed/forcing/LEONA_SCM.nc'
EXE=ROOT/'build/mynn_audit_v1/mp_offline_driver'
LOCK=threading.Lock()
SURFACE='tsk qsfc hfx lh qfx ust rmol_sf pblh z0 mavail dx previous_pblh maxwidth maxmf plume_top wstar psfc surface_theta lowest_rv dt'.split()
GATES='fltv superadiabatic maxwidth_pre minwidth dx_cap pbl_cap cloud_cap flux_cap Psig_w Psig_shcu max_resolved_w nup flux_limiter ktop fltv_entry wind_area_factor'.split()


def protected():
    paths=[]
    for c in CONFIGURATIONS:paths.append(control_directory(c)/'output.nc')
    paths+=list((ROOT/'results/argentina_20191028/calibration').rglob('*'))
    paths+=list((ROOT/'results/argentina_20181213/article').rglob('*'))
    return [{'path':str(p),'sha256':sha256(p)} for p in paths if p.is_file()]


def run(phase='original'):
    OUT.mkdir(parents=True,exist_ok=True)
    protection=OUT/'protected_inputs.json'
    if not protection.exists():write_json(protection,protected())
    configurations=CONFIGURATIONS if phase in ['original','reference_v2'] else ['TH_CF','WS_CF']
    exe=EXE if phase=='original' else ROOT/('build/mynn_audit_v3/mp_offline_driver' if phase=='no_shadow_v3' else 'build/mynn_audit_v2/mp_offline_driver')
    if phase=='alternative':
        raise ValueError('Optional alternative deferred: original boundary audit found errors; validate corrected boundary first.')
    def one(cfg):
        folder=RUNS/phase/cfg;folder.mkdir(parents=True,exist_ok=True);output=folder/'output.nc'
        text=namelist(cfg,FORCING,output).replace('write_water_budget=.true.,',
            'write_water_budget=.true.,write_vertical_diagnostics=.true.,write_mynn_audit=.true.,')
        if phase in ['no_shadow','no_shadow_v3']:text=text.replace('mynn_edmf_type=0,','mynn_edmf_type=0,mynn_audit_no_shadow=.true.,')
        if phase=='corrected':text=text.replace('mynn_edmf_type=0,','mynn_edmf_type=0,correct_mynn_boundary_units=.true.,')
        if phase=='alternative':text=text.replace('mynn_edmf_type=0,','mynn_edmf_type=1,')
        identity={'configuration':cfg,'phase':phase,'executable_sha256':sha256(exe),'forcing_sha256':sha256(FORCING),'namelist':text}
        if output.exists():
            if json.loads((folder/'identity.json').read_text())!=identity:raise ValueError('Changed audit inputs: choose fresh version')
        else:
            (folder/'namelist.mp_offline').write_text(text);write_json(folder/'identity.json',identity);execute(exe,folder)
        with LOCK:
            result=validate_output(output,172800)
            baseline=control_directory(cfg)/'output.nc'
            with xr.open_dataset(output,decode_times=False) as a,xr.open_dataset(baseline,decode_times=False) as b:
                if phase in ['original','reference_v2','no_shadow','no_shadow_v3']:
                    assert set(a.variables)==set(b.variables)
                    for name in a.variables:np.testing.assert_array_equal(a[name],b[name],err_msg=cfg+'/'+name)
                    result['all_native_variables_bit_identical']=True
                else:result['all_native_variables_bit_identical']=False
        result.update(configuration=cfg,phase=phase,output=str(output),baseline_sha256=sha256(baseline))
        print(phase+' '+cfg+': '+('bit-identical' if result['all_native_variables_bit_identical'] else 'experimental'),flush=True)
        return result
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,configurations))
    write_json(OUT/(phase+'_run_tests.json'),{'runs':results,'phase':phase,'no_ensemble_reruns':True})
    return results


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['run','no-shadow','no-shadow-v3','reference-v2','corrected','alternative','analyze','predictors','finalize']);a=p.parse_args()
    if a.stage=='analyze':
        from scm_lightning.evaluation.mynn_audit_products import analyze
        analyze()
    elif a.stage=='predictors':
        from scm_lightning.evaluation.mynn_audit_products import corrected_predictors
        corrected_predictors()
    elif a.stage=='finalize':
        from scm_lightning.evaluation.mynn_audit_products import finalize
        finalize()
    else:run({'run':'original','no-shadow':'no_shadow','no-shadow-v3':'no_shadow_v3','reference-v2':'reference_v2','corrected':'corrected','alternative':'alternative'}[a.stage])


if __name__=='__main__':main()
