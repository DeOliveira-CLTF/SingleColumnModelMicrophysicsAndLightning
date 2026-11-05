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
 apply_rt_adv=.false.,legacy_temperature_tendency=.false.,write_water_budget=.true.,
 couple_mynn_condensate={'.true.' if coupling != 'C0' else '.false.'},
 mynn_condensate_local_only={'.true.' if coupling == 'CL' else '.false.'},
 microp_scheme="{'thompson' if scheme == 'TH' else 'wsm6'}",
/
'''

def validate_output(path,expected_end=None):
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
        if expected_end is not None and not np.isclose(float(ds.time.values[-1]),expected_end,atol=1e-5,rtol=0):
            raise ValueError('Simulation does not cover complete forcing period')
        temperature=ds.theta.values*(ds.pressure.values/1.e5)**(2/7)
        if np.any((temperature<100)|(temperature>400)) or np.any(ds.pressure.values<=0):
            raise ValueError('Unphysical temperature or pressure')
        if np.any(np.diff(ds.pressure.values,axis=-1)>=0):
            raise ValueError('Pressure must decrease with height')
        result={'status':'passed_basic_numerical_checks','steps':ds.sizes['time'],
            'sha256':sha256(path),'physical_budget_gate':'not_instrumented'}
    budget=Path(str(path)+'.water_budget.csv')
    if budget.exists():
        from .io import numeric_csv
        b=numeric_csv(budget)
        finite=np.isfinite(b.values).all()
        residual=float(b.residual_kg_m2.abs().max())
        min_vapor=float(b.min_qv_before_mp.min())
        passed=bool(finite and min_vapor>=-1e-12)
        result.update(physical_budget_gate='diagnostic_available' if finite else 'diagnostic_invalid',
            physical_preflight_passed=passed,max_entry_density_water_residual_kg_m2=residual,
            entry_density_residual_within_1e_minus7=bool(residual<=1e-7),
            min_entry_vapor_kg_kg=min_vapor,budget_sha256=sha256(budget),
            budget_scope='ENTRY native-density weighting; Thompson recomputes rho internally (module_mp_thompson.F); this residual is NOT an exact variable-density closure contract. Prescribed forcing and MYNN are separate.')
    return result

def run_controls(executable, forcing, directory):
    executable = Path(executable).resolve()
    directory = Path(directory).resolve()
    results=[]
    with xr.open_dataset(forcing,decode_times=False) as source:
        expected_end=float(source.time.values[-1])
    for configuration in CONFIGURATIONS:
        folder = directory/configuration
        folder.mkdir(parents=True,exist_ok=True)
        output = folder/'output.nc'
        if output.exists():
            raise FileExistsError('Never overwrite prior run; choose a fresh directory')
        nml = namelist(configuration,forcing,output)
        nml=nml.replace('legacy_temperature_tendency=.false.,','legacy_temperature_tendency=.false.,write_process_diagnostics=.true.,')
        (folder/'namelist.mp_offline').write_text(nml)
        execute(executable,folder)
        results.append({'configuration':configuration,**validate_output(output,expected_end)})
    write_json(directory/'control_checks.json',results)
    gate={'configurations':CONFIGURATIONS,'forcing_sha256':sha256(forcing),'numerical_passed':True,
        'executable_sha256':sha256(executable),
        'physical_preflight_passed':all(r.get('physical_preflight_passed',False) for r in results),
        'physical_scope':'finite/positive species; pressure/temperature bounds; entry vapor nonnegative; finite entry-density water diagnostics. No exact mass/energy closure claim under prescribed pressure.',
        'control_results':results}
    write_json(directory/'ensemble_gate.json',gate)
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
    if not np.isfinite(t).all() or not np.isfinite(rv).all() or np.any(t<=0) or np.any(rv<=0):
        raise ValueError('Invalid realized T/log(rv) perturbation; no clipping or artificial moistening')
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
