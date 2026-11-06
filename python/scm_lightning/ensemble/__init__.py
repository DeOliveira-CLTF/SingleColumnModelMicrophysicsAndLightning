from scm_lightning.provenance import resolve_path
"""Paired synthetic initial perturbations from joint ERA5 EDA anomalies."""
import numpy as np

def joint_perturbations(temperature, rv, members=20, seed=20181213):
    """EDA arrays (10,z) -> synthetic (member,2,z), not ten reanalysis members.

    Low-rank Gaussian with sample joint T/log(r) covariance; no wind perturbation.
    Realized draws must be saved by caller, not reconstructed from seed alone.
    """
    t,r = np.asarray(temperature,float),np.asarray(rv,float)
    if t.shape != r.shape or t.ndim != 2 or t.shape[0] != 10:
        raise ValueError('Expected ten collocated EDA profiles (10,z)')
    if not np.isfinite(t).all() or not np.isfinite(r).all() or np.any(r<=0):
        raise ValueError('EDA profiles must be finite with positive rv')
    anomalies = np.concatenate([t,np.log(r)],axis=1)
    anomalies -= anomalies.mean(axis=0)
    rng = np.random.default_rng(seed)
    weights = rng.normal(size=(members,10))/np.sqrt(9)
    samples = weights@anomalies
    return samples.reshape(members,2,t.shape[1]), weights

def apply_initial(t,rv,perturbation):
    return np.asarray(t)+perturbation[0], np.asarray(rv)*np.exp(perturbation[1])

def collocate(config,control,pl_path,sl_path,output):
    """Ten regional EDA members -> common initial T/log(r) profiles.

    Surface fields follow the same 2m-to-0m policy as the deterministic forcing.
    No initial wind draw is generated; pressure remains prescribed by control.
    """
    import json
    from pathlib import Path
    import pandas as pd
    import xarray as xr
    from scm_lightning.forcing import normalize
    from scm_lightning.forcing.spatial import subset
    from scm_lightning.provenance import sha256
    if Path(output).exists():raise FileExistsError(output)
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    if config['area'] is None or config['area_selection']['status']!='fixed_before_lightning_evaluation':
        raise ValueError('Meteorologically fixed area required for EDA collocation')
    pl,sl=normalize(xr.load_dataset(pl_path)),normalize(xr.load_dataset(sl_path))
    if pl.sizes.get('number')!=10 or sl.sizes.get('number')!=10:
        raise ValueError('Exactly ten EDA reanalysis members required')
    stamp=pd.Timestamp(config['start']).tz_localize(None)
    pl=pl.sel(time=stamp);sl=sl.sel(time=stamp)
    w,s,e,n=config['area']
    pl,pw=subset(pl,config['area'])
    sl,sw=subset(sl,config['area'])
    if pl.sizes['latitude']<1 or pl.sizes['longitude']<1:
        raise ValueError('Area has no EDA grid cell')
    with xr.open_dataset(control,decode_times=False) as f:
        z=f.lev.values.copy()
    mean=lambda x:x.weighted(pw).mean(('latitude','longitude'))
    # Exclude subsurface pressure levels before averaging, for each member.
    sl_on_pl=sl[['sp']].interp(latitude=pl.latitude,longitude=pl.longitude)
    above=(pl.level*100 < sl_on_pl.sp).all(('latitude','longitude'))
    temperatures=[];mixing_ratios=[]
    for member in pl.number.values:
        p=pl.sel(number=member);surface=sl.sel(number=member)
        valid=above.sel(number=member)
        p=p.isel(level=np.flatnonzero(valid.values))
        origin=float(mean(surface.z)/9.80665)
        hz=mean(p.z).values/9.80665-origin
        ix=np.argsort(hz);hz=hz[ix]
        t=mean(p.t).values[ix]
        r=mean(p.q/(1-p.q)).values[ix]
        ps=surface.sp
        es=611.2*np.exp(17.67*(surface.d2m-273.15)/(surface.d2m-29.65))
        rs=float(mean((287.04/461.5)*es/(ps-es)))
        ts=float(mean(surface.t2m))
        if len(hz)<2 or hz[0]<=0 or hz[-1]<z[-1] or np.any(r<=0) or rs<=0:
            raise ValueError('EDA terrain/top/moisture support incomplete; no extrapolation')
        temperatures.append(np.interp(z,np.r_[0,hz],np.r_[ts,t]))
        mixing_ratios.append(np.exp(np.interp(z,np.r_[0,hz],np.log(np.r_[rs,r]))))
    if not np.isfinite(temperatures).all() or not np.isfinite(mixing_ratios).all():
        raise ValueError('Nonfinite EDA collocation')
    np.savez_compressed(output,temperature=np.asarray(temperatures),rv=np.asarray(mixing_ratios),z=z,
        forcing_sha256=sha256(control),eda_pl_sha256=sha256(pl_path),eda_sl_sha256=sha256(sl_path),
        area_json=json.dumps(config['area']),initial_time_utc=config['start'],
        interpolation='native pressure levels above terrain; T and log(r) linear in AGL height; 2m surface mapped to 0m')
    return Path(output)

def stability(statistic_by_member, sizes=(5,10,15,20)):
    """Numerical convergence diagnostics only; no automatic stopping claim."""
    a = np.asarray(statistic_by_member)
    return [{'members':n,'mean':float(a[:n].mean()),'std':float(a[:n].std(ddof=1))} for n in sizes if n<=len(a)]

def prepare(config, control, eda_profiles, control_checks, output, inherit=None):
    """Validated, collocated EDA NPZ must contain temperature,rv (10,z).

    Generate 20 paired inputs + the control shared across the six configurations.
    Requires all six case controls with numerical and physical checks passed.
    """
    import json
    from pathlib import Path
    from scm_lightning.ensemble.experiments import CONFIGURATIONS,perturb_scm,namelist
    from scm_lightning.provenance import write_json,sha256
    gate=json.loads(Path(control_checks).read_text())
    if set(gate.get('configurations',[])) != set(CONFIGURATIONS) or not gate.get('numerical_passed') or not gate.get('physical_preflight_passed'):
        raise ValueError('Six LEONA controls must pass numerical and physical gates')
    if gate.get('forcing_sha256') != sha256(control):
        raise ValueError('Controls were not run with this forcing')
    if config['area'] is None:
        raise ValueError('Fixed meteorological support required')
    data=np.load(eda_profiles,allow_pickle=False)
    if 'forcing_sha256' not in data or str(data['forcing_sha256']) != sha256(control):
        raise ValueError('EDA collocation must be tied to the control forcing checksum')
    import xarray as xr
    with xr.open_dataset(control,decode_times=False) as f:
        if 'z' not in data or not np.array_equal(data['z'],f.lev.values):
            raise ValueError('EDA vertical grid differs from control')
    n,seed=config['ensemble']['perturbed_members'],config['ensemble']['seed']
    delta,weights=joint_perturbations(data['temperature'],data['rv'],n,seed)
    prior=None;prior_runs={};prior_files={};accepted=[];old_n=0
    if inherit is not None:
        prior_path=Path(inherit);prior=json.loads(prior_path.read_text())
        old_n=prior['synthetic_members']
        if n<=old_n or prior['seed']!=seed:raise ValueError('Extension must increase members with same seed')
        for key in ['area','start','end','reference_configuration','split']:
            if prior['config'][key]!=config[key]:raise ValueError('Extension changes the fixed experiment design')
        old_draws=np.load(prior_path.parent/'realized_perturbations.npz')
        np.testing.assert_allclose(delta[:old_n],old_draws['delta'],rtol=1e-12,atol=1e-14)
        np.testing.assert_array_equal(weights[:old_n],old_draws['weights'])
        delta[:old_n]=old_draws['delta'] # exact saved prefix, independent of BLAS tiling
        prior_runs={(r['configuration'],r['member']):r for r in prior['runs']}
        prior_files={r['member']:resolve_path(r['forcing']) for r in prior['runs']}
        accepted=[sha256(prior_path),*prior.get('accepted_prior_manifest_sha256',[])]
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/'realized_perturbations.npz',delta=delta,weights=weights,
        eda_temperature=data['temperature'],eda_rv=data['rv'],seed=seed)
    files=[Path(control).resolve()]
    for i in range(n):
        if i+1<=old_n:
            files.append(prior_files[i+1]);continue
        target=out/f'member_{i+1:02}.nc'
        perturb_scm(control,target,delta[i]);files.append(target.resolve())
    runs=[]
    for configuration in CONFIGURATIONS:
        for member,forcing in enumerate(files):
            if (configuration,member) in prior_runs:
                runs.append(prior_runs[(configuration,member)]);continue
            folder=out/configuration/f'{member:02}';folder.mkdir(parents=True)
            (folder/'namelist.mp_offline').write_text(namelist(configuration,forcing,folder/'output.nc'))
            runs.append({'configuration':configuration,'member':member,'forcing':str(forcing),'forcing_sha256':sha256(forcing),'folder':str(folder.resolve())})
    write_json(out/'manifest.json',{'case':config['case'],'seed':seed,'distribution':'low-rank Gaussian joint T/log(rv) from sample EDA covariance',
        'synthetic_members':n,'reanalysis_members_used_for_covariance':10,'control_member':0,'paired_across_configurations':True,
        'pressure_balance':'prescribed pressure; no exact hydrostatic balance','wind_perturbations':False,
        'controls_sha256':sha256(control_checks),'config':config,'runs':runs,
        'control_executable_sha256':gate.get('executable_sha256'),
        'accepted_prior_manifest_sha256':accepted,'inherited_members':old_n,
        'realized_perturbations_sha256':sha256(out/'realized_perturbations.npz')})
    return runs

def run_prepared(executable,manifest):
    """Run the principal driver from a checked preparation manifest; resumable."""
    import json
    from pathlib import Path
    from scm_lightning.ensemble.experiments import execute,validate_output
    from scm_lightning.provenance import sha256,write_json
    from scm_lightning.download import NETCDF_LOCK
    from concurrent.futures import ThreadPoolExecutor,as_completed
    import xarray as xr
    path=Path(manifest)
    plan=json.loads(path.read_text())
    identity={'manifest_sha256':sha256(path),'executable_sha256':sha256(executable)}
    if plan.get('control_executable_sha256')!=identity['executable_sha256']:
        raise ValueError('Ensemble executable differs from the checked controls')
    records=[]
    def one(run):
        folder=resolve_path(run['folder']);forcing=resolve_path(run['forcing'])
        if not folder.is_absolute():
            from scm_lightning.provenance import ROOT
            folder=ROOT/folder
        # Manifest paths are workspace-local. Inputs and namelists are immutable
        # over a resumed run, recorded before the first process launch.
        if sha256(forcing)!=run['forcing_sha256']:
            raise ValueError('Prepared forcing changed')
        output=folder/'output.nc';tag=folder/'run_identity.json'
        expected=dict(identity,forcing_sha256=run['forcing_sha256'],namelist_sha256=sha256(folder/'namelist.mp_offline'))
        if output.exists():
            recorded=json.loads(tag.read_text()) if tag.exists() else {}
            accepted={identity['manifest_sha256'],*plan.get('accepted_prior_manifest_sha256',[])}
            if recorded.get('manifest_sha256') not in accepted or any(recorded.get(k)!=v for k,v in expected.items() if k!='manifest_sha256'):
                raise ValueError('Existing run has a different/unrecorded identity')
        else:
            write_json(tag,expected)
            execute(Path(executable).resolve(),folder)
        with NETCDF_LOCK:
            with xr.open_dataset(forcing,decode_times=False) as source:expected_end=float(source.time.values[-1])
            checked=validate_output(output,expected_end)
        return dict(run,**checked)
    # Independent processes have separate physics-module state. NetCDF checks
    # run under a lock because the local C library is not assumed thread-safe.
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(one,run) for run in plan['runs']]
        for future in as_completed(futures):
            records.append(future.result())
            records.sort(key=lambda r:(r['configuration'],r['member']))
            write_json(path.parent/'run_checks.json',records)
    return records
