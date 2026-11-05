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

def stability(statistic_by_member, sizes=(5,10,15,20)):
    """Numerical convergence diagnostics only; no automatic stopping claim."""
    a = np.asarray(statistic_by_member)
    return [{'members':n,'mean':float(a[:n].mean()),'std':float(a[:n].std(ddof=1))} for n in sizes if n<=len(a)]

def prepare(config, control, eda_profiles, control_checks, output):
    """Validated, collocated EDA NPZ must contain temperature,rv (10,z).

    Generate 20 paired inputs + the control shared across the six configurations.
    Requires all six case controls with numerical and physical checks passed.
    """
    import json
    from pathlib import Path
    from .experiments import CONFIGURATIONS,perturb_scm,namelist
    from .provenance import write_json,sha256
    gate=json.loads(Path(control_checks).read_text())
    if set(gate.get('configurations',[])) != set(CONFIGURATIONS) or not gate.get('numerical_passed') or not gate.get('physical_budget_passed'):
        raise ValueError('Six LEONA controls must pass numerical and physical gates')
    if gate.get('forcing_sha256') != sha256(control):
        raise ValueError('Controls were not run with this forcing')
    if config['area'] is None:
        raise ValueError('Fixed meteorological support required')
    data=np.load(eda_profiles,allow_pickle=False)
    if 'forcing_sha256' not in data or str(data['forcing_sha256']) != sha256(control):
        raise ValueError('EDA collocation must be tied to the control forcing checksum')
    n,seed=config['ensemble']['perturbed_members'],config['ensemble']['seed']
    delta,weights=joint_perturbations(data['temperature'],data['rv'],n,seed)
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/'realized_perturbations.npz',delta=delta,weights=weights,
        eda_temperature=data['temperature'],eda_rv=data['rv'],seed=seed)
    files=[Path(control).resolve()]
    for i in range(n):
        target=out/f'member_{i+1:02}.nc'
        perturb_scm(control,target,delta[i]);files.append(target.resolve())
    runs=[]
    for configuration in CONFIGURATIONS:
        for member,forcing in enumerate(files):
            folder=out/configuration/f'{member:02}';folder.mkdir(parents=True)
            (folder/'namelist.mp_offline').write_text(namelist(configuration,forcing,folder/'output.nc'))
            runs.append({'configuration':configuration,'member':member,'forcing':str(forcing),'forcing_sha256':sha256(forcing),'folder':str(folder)})
    write_json(out/'manifest.json',{'case':config['case'],'seed':seed,'distribution':'low-rank Gaussian joint T/log(rv) from sample EDA covariance',
        'synthetic_members':n,'reanalysis_members_used_for_covariance':10,'control_member':0,'paired_across_configurations':True,
        'pressure_balance':'prescribed pressure; no exact hydrostatic balance','wind_perturbations':False,
        'controls_sha256':sha256(control_checks),'config':config,'runs':runs,
        'realized_perturbations_sha256':sha256(out/'realized_perturbations.npz')})
    return runs
