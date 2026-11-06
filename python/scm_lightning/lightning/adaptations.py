from scm_lightning.config import load_config
"""Exploratory 1D diagnostics; never aliases of the original CRM formulations.

EDMF flux proxies below assume that grid-column hydrometeors represent the
active area-mean plume. They are NOT covariant plume species fluxes.
"""
import argparse
import json
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import xarray as xr
from scm_lightning.provenance import ROOT, sha256, write_json

VERSION = 'SCM_vertical_diagnostics_v1'
ASSUMPTION = ('Exploratory: active area-mean EDMF plume velocity combined with '
              'grid-column hydrometeors assumed representative of those plumes; '
              'no resolved plume hydrometeor covariance or entrainment correction.')


def validate_units(dataset):
    expected={'w_LS':'m s-1','w_LS_interface':'m s-1','w_EDMF':'m s-1',
              'mf_w':'kg m-2 s-1','mf_a':'1','rho_MYNN':'kg m-3',
              'rhoz':'kg m-3','z_interface':'m','layer_thickness':'m',
              'temperature':'K','pressure':'Pa'}
    for name, unit in expected.items():
        if dataset[name].attrs.get('units') != unit:
            raise ValueError('Unexpected diagnostic units: '+name)


def interface_density(rho, dz):
    """Exactly bl_mynn_subroutines.F90 rhoz(k), excluding the unused lid."""
    rho, dz = np.asarray(rho), np.asarray(dz)
    if rho.shape[-1] != len(dz) or np.any(dz <= 0) or not np.isfinite(rho).all() or np.any(rho <= 0):
        raise ValueError('Invalid MYNN density or layer geometry')
    return (rho[..., :-1]*dz[1:]+rho[..., 1:]*dz[:-1])/(dz[1:]+dz[:-1])


def recover_velocity(mf_w, rhoz, mf_a):
    """mf_w/rhoz, NOT mf_w/(rhoz*mf_a); inactive plumes stay missing."""
    mf_w, rhoz, mf_a = np.broadcast_arrays(mf_w, rhoz, mf_a)
    if not (np.isfinite(mf_w).all() and np.isfinite(rhoz).all() and np.isfinite(mf_a).all()):
        raise ValueError('Missing native plume diagnostics')
    if np.any(rhoz <= 0) or np.any((mf_a < 0) | (mf_a > 1)):
        raise ValueError('Invalid interface density or area fraction')
    return np.divide(mf_w, rhoz, out=np.full(mf_w.shape, np.nan), where=mf_a > 0)


def crossing(temperature, z, target):
    """Lowest downward isotherm crossing. No invented/extrapolated crossing."""
    t = np.asarray(temperature)
    mask = (t[..., :-1] >= target) & (t[..., 1:] <= target) & (t[..., :-1] > t[..., 1:])
    index = np.argmax(mask, axis=-1)
    valid = mask.any(axis=-1)
    lo = np.take_along_axis(t, index[..., None], axis=-1)[..., 0]
    hi = np.take_along_axis(t, (index+1)[..., None], axis=-1)[..., 0]
    fraction = np.divide(lo-target, lo-hi, out=np.zeros_like(lo), where=lo != hi)
    height = np.asarray(z)[index]+fraction*(np.asarray(z)[index+1]-np.asarray(z)[index])
    return np.where(valid, height, np.nan)


def at_height(values, z, height):
    values = np.asarray(values)
    return np.asarray([np.interp(h, z, row, left=np.nan, right=np.nan)
                       if np.isfinite(h) else np.nan for row, h in zip(values, height)])


def surface_cape(p, temperature, rv, layer_selection='bottom_top'):
    """MetPy surface parcel, liquid pseudoadiabatic CAPE (explicit adaptation).

    Parcel starts in the lowest simulated cell, not an invented surface level.
    qv is dry-air mixing ratio, so vapor pressure is p*r/(epsilon+r).
    This diagnostic is NOT asserted to reproduce Romps' isentropic CAPE.
    """
    from metpy.calc import dewpoint, surface_based_cape_cin, parcel_profile_with_lcl, cape_cin
    from metpy.units import units
    p, temperature, rv = map(np.asarray, (p, temperature, rv))
    if not np.isfinite([p, temperature, rv]).all() or np.any(p <= 0) or np.any(np.diff(p) >= 0) or np.any(rv < 0):
        raise ValueError('Invalid thermodynamic profile for CAPE')
    if rv[0] == 0:
        return np.nan  # A zero-vapor parcel has no finite dewpoint for this provider.
    vapor_pressure = p*rv/(0.6219569100577033+rv)
    # At exactly dry upper levels, drop only those dry samples rather than fabricate humidity.
    valid = vapor_pressure > 0
    if valid.sum() < 4:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        td = dewpoint(vapor_pressure[valid]*units.Pa)
        if layer_selection=='bottom_top':
            cape, _ = surface_based_cape_cin(p[valid]*units.Pa, temperature[valid]*units.K, td)
        elif layer_selection=='most_cape':
            pressure,environment,dewpoints,parcel=parcel_profile_with_lcl(p[valid]*units.Pa,temperature[valid]*units.K,td)
            cape, _ = cape_cin(pressure,environment,dewpoints,parcel,which_lfc='most_cape',which_el='most_cape')
        else:
            raise ValueError('Unknown a priori CAPE layer-selection policy')
    value=float(cape.to('J/kg').magnitude)
    # Never turn a negative signed integral into invented positive energy.
    return value if layer_selection=='bottom_top' or value>=0 else np.nan


def lpi_1d(temperature, dz, qc, qr, qi, qs, qg, w):
    """Yair epsilon, volume -> depth mean, no horizontal-neighborhood test."""
    def ratio(a, b):
        return np.divide(a, b, out=np.zeros_like(a), where=b > 0)
    ice = qg*(ratio(np.sqrt(qs*qg), qs+qg)+ratio(np.sqrt(qi*qg), qi+qg))
    liquid = qc+qr
    epsilon = ratio(2*np.sqrt(ice*liquid), ice+liquid)
    mask = (temperature >= 253.15) & (temperature <= 273.15)
    depth = np.sum(mask*dz, axis=-1)
    return np.divide(np.sum(np.where(mask, epsilon*w*w*dz, 0), axis=-1), depth,
                     out=np.full(depth.shape, np.nan), where=depth > 0)


def diagnose_file(vertical_path, native_path, origin, minutes=15):
    """Step-begin snapshots × actual step precipitation; sum complete 15-min bins."""
    with xr.open_dataset(vertical_path, decode_times=False) as source:
        d = source.load()
    validate_units(d)
    with xr.open_dataset(native_path, decode_times=False) as source:
        native = source[['time', 'rainnc']].load()
    z = d.lev.values
    edges = np.zeros(len(z)+1)
    for k, center in enumerate(z):
        edges[k+1] = 2*center-edges[k]
    dz = d.layer_thickness.values
    np.testing.assert_allclose(dz, np.diff(edges), atol=1e-10, rtol=0)
    np.testing.assert_allclose(d.z_interface, edges[1:-1], atol=1e-10, rtol=0)
    expected_rho = interface_density(d.rho_MYNN.values, dz)
    np.testing.assert_allclose(d.rhoz.values, expected_rho, rtol=2e-14, atol=0)
    velocity = recover_velocity(d.mf_w.values, d.rhoz.values, d.mf_a.values)
    np.testing.assert_array_equal(velocity, d.w_EDMF.values)
    np.testing.assert_array_equal(d.plume_active.values > 0, d.mf_a.values > 0)
    step = np.diff(native.time.values)
    np.testing.assert_array_equal(d.time.values, native.time.values[:-1])
    if not np.all(step > 0):
        raise ValueError('Invalid model time axis')
    precipitation = np.diff(native.rainnc.values)/step
    if np.any(precipitation < -1e-10):
        raise ValueError('Precipitation reset')
    precipitation = np.maximum(precipitation, 0)
    temperature, pressure, rv = d.temperature.values, d.pressure.values, d.qv.values
    rho_d = pressure/(287.*temperature*(1+rv/(287./461.6)))
    qc, qr, qi, qs, qg = (d[q].values for q in ['qc', 'qr', 'qi', 'qs', 'qg'])
    vii = np.sum(rho_d*(qi+qs+qg)*dz, axis=-1)
    raw = {'McCaul_F2_1D_adaptado': .20*vii}
    reach = pd.DataFrame({'time_s': d.time.values})
    for target, label in [(258.15, 'minus15C'), (260., '260K')]:
        height = crossing(temperature, z, target)
        area = at_height(d.mf_a.values, d.z_interface.values, height)
        # Both bracketing interfaces must contain active plumes, preventing a
        # fictitious plume crossing an inactive gap through interpolation.
        bracket = np.clip(np.searchsorted(d.z_interface.values, height), 1, len(z)-2)
        active = d.mf_a.values > 0
        rows = np.arange(len(height))
        reached = np.isfinite(height) & (height >= d.z_interface.values[0]) & (height <= d.z_interface.values[-1])
        reached &= active[rows, bracket-1] & active[rows, bracket]
        reach[label+'_height_m'] = height
        reach[label+'_active_plume'] = reached
        reach[label+'_mf_a_interpolated'] = area
    for representation in ['LS', 'EDMF']:
        if representation == 'LS':
            w = d.w_LS.values
        else:
            # Only for depth-mean / maximum diagnostics: inactive interfaces
            # mean no mass-flux ascent, while the exported plume velocity stays NaN.
            w = np.asarray([np.interp(z, d.z_interface, np.nan_to_num(row, nan=0.), left=0., right=0.) for row in velocity])
        wmax = np.max(w if representation == 'LS' else np.nan_to_num(velocity, nan=0.), axis=-1)
        suffix = '_'+representation+'_1D_adaptado'
        raw['PR92_wmax'+suffix] = 5.e-6*np.maximum(0, wmax)**4.5
        raw['PW10'+suffix] = precipitation*(wmax > 10.)
        for target, method in [(260., 'IFluxT'), (258.15, 'McCaul_F1')]:
            height = reach['260K_height_m' if target == 260 else 'minus15C_height_m'].values
            if representation == 'LS':
                speed = at_height(w, z, height)
            else:
                speed = at_height(np.nan_to_num(velocity, nan=0.), d.z_interface.values, height)
                active = reach['260K_active_plume' if target == 260 else 'minus15C_active_plume'].values
                speed = np.where(active, speed, np.where(np.isfinite(height), 0., np.nan))
            # Interpolate the two factors separately to the same isotherm.
            if method == 'IFluxT':
                value = at_height(rho_d*qi, z, height)*speed*(speed > 1.)
            else:
                value = .042*at_height(qg*1000., z, height)*np.maximum(0, speed)
            raw[method+suffix] = value
        raw['McCaul_F3'+suffix] = .95*raw['McCaul_F1'+suffix]+.05*raw['McCaul_F2_1D_adaptado']
        raw['LPI'+suffix] = lpi_1d(temperature, dz, qc, qr, qi, qs, qg, w)
    seconds = minutes*60
    groups = (d.time.values//seconds).astype(int)
    records = []
    for group in np.unique(groups):
        mask = groups == group
        if not np.isclose(step[mask].sum(), seconds):
            raise ValueError('Incomplete diagnostic bin')
        timestamp = pd.Timestamp(origin)+pd.Timedelta(seconds=int(group*seconds))
        first = np.flatnonzero(mask)[0]
        cape = surface_cape(pressure[first], temperature[first], rv[first])
        records.append({'time': timestamp, 'method': 'CAPExP_SBCAPE15_1D_adaptado',
                        'predictor': cape*np.sum(precipitation[mask]*step[mask])})
        for method, values in raw.items():
            value = np.sum(values[mask]*step[mask]) if np.isfinite(values[mask]).all() else np.nan
            records.append({'time': timestamp, 'method': method, 'predictor': value})
    predictions = pd.DataFrame(records)
    predictions['status'] = np.where(np.isfinite(predictions.predictor), 'available', 'unavailable')
    return predictions, reach, d


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-root', default='runs/argentina_20181213/_diagnostic_audit_v2')
    args = parser.parse_args()
    c = load_config(ROOT/'configs/cases/argentina_20181213.yaml')
    out = ROOT/'results/argentina_20181213/vertical_diagnostics'
    out.mkdir(parents=True, exist_ok=True)
    summary, all_predictions = [], []
    for folder in sorted((ROOT/args.audit_root).iterdir()):
        if not folder.is_dir():
            continue
        native = folder/'output.nc'
        vertical = folder/'output.nc.vertical.nc'
        result, reach, dataset = diagnose_file(vertical, native, c['start'], c['interval_minutes'])
        destination = out/folder.name
        destination.mkdir(exist_ok=True)
        result['configuration'] = folder.name
        result['member'] = 0
        result.to_csv(destination/'predictors_15min.csv', index=False)
        reach.to_csv(destination/'isotherm_reach.csv.gz', index=False)
        all_predictions.append(result)
        summary.append({'configuration': folder.name, 'mynn_edmf_type': int(dataset.attrs['mynn_edmf_type']),
            'mynn_uedmf_mp': int(dataset.attrs['mynn_uedmf_mp']), 'samples': len(reach),
            'active_260K_samples': int(reach['260K_active_plume'].sum()),
            'active_minus15C_samples': int(reach['minus15C_active_plume'].sum()),
            'graupel_plume_flux_available': False, 'species_channel_order': 'qr,qi,qs,qg',
            'reason': 'Main MYNN type=0 mode=0 does not compute uedmf_sp_flux',
            'units_density_velocity_masks_verified': True,
            'vertical_source': str(vertical.relative_to(ROOT)), 'vertical_sha256': sha256(vertical)})
    if len(summary) != 6:
        raise ValueError('Expected all six configurations')
    pd.concat(all_predictions).to_csv(out/'predictors_15min.csv.gz', index=False)
    pd.DataFrame(summary).to_csv(out/'availability_and_reach.csv', index=False)
    write_json(out/'manifest.json', {'version': VERSION, 'runs': summary, 'EDMF_hydrometeor_assumption': ASSUMPTION,
        'velocity_sum': False, 'physics_changed': False, 'ensemble_scope': 'Six deterministic controls; not a rerun of the 40-member ensemble',
        'origin_utc':c['start'], 'vertical_height_reference':'m AGL; internal top-of-layer interfaces; surface and lid omitted',
        'CAPE': 'MetPy liquid pseudoadiabatic surface parcel from lowest simulated cell; diagnosed every 15 minutes, held constant in bin. Adaptation of Romps CAPExP.',
        'time_integration': 'Step-begin snapshots times step duration; P from exact native rainnc increments. No snownc addition.',
        'isotherm': 'Lowest downward column-temperature crossing; both neighboring interfaces active for EDMF reach.',
        'EDMF_original_flux_status': 'Unavailable without collocated plume species. Exploratory column-representative versions have distinct method IDs.',
        'LPI': '1D depth mean with center-cell 0 to -20 C mask; no horizontal neighborhood criterion; explicitly adapted.'})
    print('Saved six-configuration exploratory predictors and native isotherm-reach tables.')


def evaluate_main():
    """Optional, separately labeled calibration of the six exploratory controls."""
    from scm_lightning.evaluation import fit_scale
    if load_config(ROOT/'configs/cases/argentina_20181213.yaml').get('role')=='evaluation_only':
        raise ValueError('No December 2018 calibration, including exploratory versions; retain raw LS/EDMF diagnostics')
    parser=argparse.ArgumentParser(description='Exploratory controls only; does not modify the main ensemble calibration')
    parser.add_argument('--output', default='results/argentina_20181213/vertical_diagnostics/control_evaluation')
    args=parser.parse_args()
    folder=ROOT/args.output
    if folder.exists():
        raise FileExistsError('Choose a fresh exploratory evaluation directory')
    c=load_config(ROOT/'configs/cases/argentina_20181213.yaml')
    root=ROOT/'results/argentina_20181213'
    predictors=pd.read_csv(root/'vertical_diagnostics/predictors_15min.csv.gz')
    observed=pd.read_csv(root/'observations/counts_15min.csv')
    table=predictors.merge(observed[['time','count']],on='time',validate='many_to_one')
    table['time']=pd.to_datetime(table.time,utc=True)
    coefficients=[]
    for (method,configuration), rows in table.groupby(['method','configuration']):
        development=rows[rows.time<pd.Timestamp(c['split']['development_end'])]
        scale,status=fit_scale(development.predictor.values,development['count'].values)
        n=int((np.isfinite(development.predictor)&np.isfinite(development['count'])).sum())
        coefficients.append({'method':method,'configuration':configuration,'scale':scale,'status':status,'n_development':n})
    coefficient=pd.DataFrame(coefficients)
    held=table[table.time>=pd.Timestamp(c['split']['evaluation_start'])].merge(coefficient,on=['method','configuration'])
    ref=coefficient[coefficient.configuration==c['reference_configuration']][['method','scale']].rename(columns={'scale':'reference_scale'})
    held=held.merge(ref,on='method',validate='many_to_one')
    held['prediction_A']=held.predictor*held.scale
    held['prediction_B']=held.predictor*held.reference_scale
    # A zero known predictor remains a known zero forecast if its scale is unidentifiable.
    for name in ['prediction_A','prediction_B']:
        held.loc[held.predictor==0,name]=0.
    metrics=[]
    for (method,configuration), rows in held.groupby(['method','configuration']):
        for comparison in ['A','B']:
            pred=rows['prediction_'+comparison].values;obs=rows['count'].values
            valid=np.isfinite(pred)&np.isfinite(obs)
            metrics.append({'method':method,'configuration':configuration,'comparison':comparison,'n':int(valid.sum()),
                'rmse':float(np.sqrt(np.mean((pred[valid]-obs[valid])**2))) if valid.any() else np.nan,
                'mae':float(np.mean(np.abs(pred[valid]-obs[valid]))) if valid.any() else np.nan})
    folder.mkdir(parents=True)
    coefficient.to_csv(folder/'coefficients.csv',index=False)
    held.to_csv(folder/'heldout_predictions.csv',index=False)
    pd.DataFrame(metrics).to_csv(folder/'metrics.csv',index=False)
    write_json(folder/'manifest.json',{'scope':'Exploratory deterministic controls only; separate from 40-member primary ensemble',
        'coefficient_fit':'Development control trajectories only; no member-specific fits',
        'split':c['split'],'reference':c['reference_configuration'],
        'EDMF_hydrometeor_assumption':ASSUMPTION,'prediction_spread':'Not estimated from six physical controls'})


if __name__ == '__main__':
    main()
