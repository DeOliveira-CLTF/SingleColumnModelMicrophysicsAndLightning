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
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='work/leona/driver_checks');args=ap.parse_args()
    output=ROOT/args.output;output.mkdir(parents=True,exist_ok=False)
    exe=ROOT/'work/leona/build/mp_offline_driver'
    forcing=ROOT/'data_goamazon/GOAMAZON_SINGLE_SCM_driver.nc'
    records=[]
    for configuration in CONFIGURATIONS:
        physical=None
        for diagnostics in [False,True]:
            folder=output/(configuration+('_hooks' if diagnostics else '_plain'));folder.mkdir()
            nml=namelist(configuration,forcing,folder/'output.nc')
            if os.name=='nt':
                nml=nml.replace(forcing.resolve().as_posix(),linux(forcing)).replace((folder/'output.nc').resolve().as_posix(),linux(folder/'output.nc'))
            nml=nml.replace('legacy_temperature_tendency=.false.,',f'legacy_temperature_tendency=.false.,write_process_diagnostics={".true." if diagnostics else ".false."},')
            (folder/'namelist.mp_offline').write_text(nml)
            execute(exe,folder)
            result=validate_output(folder/'output.nc')
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
