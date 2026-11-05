"""Historical forcing smoke test and main-driver instrumentation invariance.

These are NOT LEONA case controls. Physical trajectory compared bit for bit
with process diagnostics disabled/enabled, for the six configurations.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import xarray as xr
from leona.experiments import CONFIGURATIONS,namelist,validate_output
from leona.provenance import ROOT,sha256,write_json

def linux(path):
    p=Path(path).resolve().as_posix()
    return '/mnt/'+p[0].lower()+p[2:] if len(p)>1 and p[1]==':' else p

def execute(executable,folder):
    with (folder/'execution.log').open('w') as log:
        if os.name=='nt':
            command='cd '+shlex.quote(linux(folder))+' && exec '+shlex.quote(linux(executable))
            subprocess.run(['wsl','-d','Ubuntu-24.04','--','bash','-lc',command],stdout=log,stderr=subprocess.STDOUT,check=True)
        else:
            subprocess.run([str(executable)],cwd=folder,stdout=log,stderr=subprocess.STDOUT,check=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='work/leona/driver_checks');ap.add_argument('--resume',action='store_true');args=ap.parse_args()
    output=ROOT/args.output;output.mkdir(parents=True,exist_ok=args.resume)
    exe=ROOT/'work/leona/build/mp_offline_driver'
    forcing=ROOT/'data_goamazon/GOAMAZON_SINGLE_SCM_driver.nc'
    records=[]
    for configuration in CONFIGURATIONS:
        physical=None
        for diagnostics in [False,True]:
            folder=output/(configuration+('_hooks' if diagnostics else '_plain'));folder.mkdir(exist_ok=args.resume)
            nml=namelist(configuration,forcing,folder/'output.nc')
            if not diagnostics:nml=nml.replace('write_water_budget=.true.','write_water_budget=.false.')
            if os.name=='nt':
                nml=nml.replace(forcing.resolve().as_posix(),linux(forcing)).replace((folder/'output.nc').resolve().as_posix(),linux(folder/'output.nc'))
            nml=nml.replace('legacy_temperature_tendency=.false.,',f'legacy_temperature_tendency=.false.,write_process_diagnostics={".true." if diagnostics else ".false."},')
            nmlfile=folder/'namelist.mp_offline'
            if nmlfile.exists() and nmlfile.read_text()!=nml:raise ValueError('Different resumed test namelist')
            nmlfile.write_text(nml)
            if not (folder/'output.nc').exists():execute(exe,folder)
            with xr.open_dataset(forcing,decode_times=False) as source:
                first=float(source.time.values[0]);last=float(source.time.values[-1])
                end=first+30*np.ceil((last-first)/30) # preserve historical loop's terminal partial forcing interval
            result=validate_output(folder/'output.nc',end)
            with xr.open_dataset(folder/'output.nc',decode_times=False) as ds:
                arrays={name:ds[name].values.copy() for name in ds if not name.startswith('flash_rate')}
                if physical is None: physical=arrays
                else:
                    for name in arrays:np.testing.assert_array_equal(arrays[name],physical[name],err_msg=configuration+'/'+name)
            records.append({'configuration':configuration,'diagnostics':diagnostics,**result})
        records[-1]['instrumentation_bit_identical']=True
    write_json(ROOT/'products/leona/tests/driver_checks.json',{'scope':'historical GOAMAZON forcing; not LEONA controls','executable_sha256':sha256(exe),'forcing_sha256':sha256(forcing),'runs':records})
    print('12 historical-forcing runs passed; six instrumentation comparisons bit identical')

if __name__=='__main__':main()
