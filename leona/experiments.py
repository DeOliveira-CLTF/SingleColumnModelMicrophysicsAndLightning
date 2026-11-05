"""Driver-principal runs only. Controls gate paired synthetic ensemble runs."""
from pathlib import Path
import os
import subprocess
import shlex
import numpy as np
import xarray as xr
from .provenance import write_json,sha256

CONFIGURATIONS = [f'{s}_{c}' for s in ['TH','WS'] for c in ['C0','CL','CF']]

def model_path(path):
    p=Path(path).resolve().as_posix()
    if os.name=='nt' and len(p)>1 and p[1]==':':
        return '/mnt/'+p[0].lower()+p[2:]
    return p

def execute(executable,folder):
    with (Path(folder)/'execution.log').open('w') as log:
        if os.name=='nt':
            command='cd '+shlex.quote(model_path(folder))+' && exec '+shlex.quote(model_path(executable))
            subprocess.run(['wsl','-d','Ubuntu-24.04','--','bash','-lc',command],stdout=log,stderr=subprocess.STDOUT,check=True)
        else:
            subprocess.run([str(executable)],cwd=folder,stdout=log,stderr=subprocess.STDOUT,check=True)

def namelist(configuration, forcing, output):
    scheme,coupling = configuration.split('_')
    if configuration not in CONFIGURATIONS:
        raise ValueError(configuration)
    return f'''&mp_offline
 path_input="{model_path(forcing)}",
 path_output="{model_path(output)}",
 nz=60,ztop=20000.,dz0=50.,dt_microp=30.,
 use_mynn=.true.,mynn_edmf=1,mynn_edmf_type=0,mynn_uedmf_mp=0,
 apply_rt_adv=.false.,legacy_temperature_tendency=.false.,
 couple_mynn_condensate={'.true.' if coupling != 'C0' else '.false.'},
 mynn_condensate_local_only={'.true.' if coupling == 'CL' else '.false.'},
 microp_scheme="{'thompson' if scheme == 'TH' else 'wsm6'}",
/
'''

def validate_output(path):
    with xr.open_dataset(path,decode_times=False) as ds:
        required={'theta','pressure','qv','qc','qr','qi','qs','qg','rainnc','time'}
        if not required<=set(ds.variables):
            raise ValueError(f'Missing output variables: {required-set(ds.variables)}')
        for name in ['theta','temperature','pressure','qv','qc','qr','qi','qs','qg','rainnc']:
            if name not in ds:
                continue
            a = ds[name].values
            if not np.isfinite(a).all():
                raise ValueError(f'Nonfinite {name}')
            if name in ['qv','qc','qr','qi','qs','qg'] and np.min(a)<-1e-12:
                raise ValueError(f'Negative mixing ratio {name}')
        if np.any(np.diff(ds.rainnc)<-1e-7):
            raise ValueError('Missing/reset total precipitation')
        if 'time' not in ds or np.any(np.diff(ds.time)<=0):
            raise ValueError('Nonmonotonic output time')
        return {'status':'passed_basic_numerical_checks','steps':ds.sizes['time'],
            'sha256':sha256(path),'physical_budget_gate':'pending_external_budget_review'}

def run_controls(executable, forcing, directory):
    executable = Path(executable).resolve()
    directory = Path(directory).resolve()
    results=[]
    for configuration in CONFIGURATIONS:
        folder = directory/configuration
        folder.mkdir(parents=True,exist_ok=True)
        output = folder/'output.nc'
        if output.exists():
            raise FileExistsError('Never overwrite prior run; choose a fresh directory')
        nml = namelist(configuration,forcing,output)
        (folder/'namelist.mp_offline').write_text(nml)
        execute(executable,folder)
        results.append({'configuration':configuration,**validate_output(output)})
    write_json(directory/'control_checks.json',results)
    return results

def perturb_scm(control, destination, delta):
    """Perturb only initial T/log(rv); retain ALL time-varying forcing arrays."""
    from .contracts import exner
    with xr.open_dataset(control,decode_times=False) as source:
        f = source.load()
    if Path(destination).exists():
        raise FileExistsError(destination)
    # Derived thermodynamic aliases must remain consistent; use dephy conversion
    # after changing initial theta/rv to avoid stale ta/qv in the resulting input.
    Case,_ = __import__('leona.forcing',fromlist=['load_dephy']).load_dephy()
    case = Case('LEONA_ensemble'); case.read(str(control))
    t = f.ta.values.squeeze()+delta[0]
    rv = f.rv.values.squeeze()*np.exp(delta[1])
    coordinate = case.variables['ua']
    for name in ['ta','theta','thetal','qv','qt','rv','rt','ql','qi','hur']:
        if name in case.variables:
            case.remove_variable(name)
    coords={'lev':coordinate.level,'height':case.variables['zh'],'pressure':case.variables['pa']}
    case.add_init_temp(t,**coords)
    case.add_init_rv(rv,**coords)
    # Recompute aliases using the official Case methods, preserving prescribed
    # pa/zh rather than invoking the hydrostatic completion path again.
    for name,method in [('theta','compute_theta'),('thetal','compute_thetal')]:
        case.add_init_variable(name,getattr(case,method)(pressure=case.variables['pa']),**coords)
    for name,method in [('rt','compute_rt'),('qv','compute_qv'),('qt','compute_qt')]:
        case.add_init_variable(name,getattr(case,method)(),**coords)
    # Condensed dry-air mixing ratios stay fixed; corresponding moist-air
    # specific contents follow the new total dry-air mixing ratio.
    denom=1+case.variables['rt'].data
    for name,source in [('ql','rl'),('qi','ri')]:
        case.add_init_variable(name,case.variables[source].data/denom,**coords)
    case.add_init_variable('hur',case.compute_hur(),**coords)
    # Case.read creates per-variable axis IDs. This is still an SCM file:
    # retain one common lev dimension for ALL initial and forcing profiles.
    from .forcing import common_grid
    common_grid(case,f.lev.values)
    case.write(str(destination))
    with xr.open_dataset(destination,decode_times=False) as changed:
        for name in f.data_vars:
            if 'time' in f[name].dims:
                np.testing.assert_array_equal(f[name].values,changed[name].values)
