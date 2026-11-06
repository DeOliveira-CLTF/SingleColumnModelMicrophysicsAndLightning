"""Predeclared meteorology-only pulse design and 18 deterministic main-driver runs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from scm_lightning.config import load_config
from scm_lightning.download import NETCDF_LOCK
from scm_lightning.forcing import normalize
from scm_lightning.forcing.spatial import subset
from scm_lightning.lightning.adaptations import crossing
from scm_lightning.ensemble.experiments import CONFIGURATIONS, execute, namelist, control_directory, validate_output
from scm_lightning.provenance import ROOT, git, sha256, write_json

RUNS=ROOT/'runs/argentina_20181213/updraft_sensitivity'
OUT=ROOT/'results/argentina_20181213/updraft_sensitivity'
CONFIG=ROOT/'configs/experiments/updraft_sensitivity.yaml'
FORCING=ROOT/'data/processed/forcing/LEONA_SCM.nc'


def pulse(z, time, start, duration, center, sigma, amplitude):
    if duration<=0 or sigma<=0 or amplitude<0:
        raise ValueError('Invalid pulse parameters')
    tau=np.asarray(time)-start
    shape=np.where((tau>=0)&(tau<=duration),.5*(1-np.cos(2*3.14159265358979*tau/duration)),0.)
    return amplitude*shape[...,None]*np.exp(-.5*((np.asarray(z)-center)/sigma)**2)


def design():
    """No lightning observations, simulated outcomes or calibration products read."""
    if CONFIG.exists():
        raise FileExistsError('Frozen pulse design exists; do not retune amplitudes')
    case=load_config(ROOT/'configs/cases/argentina_20181213.yaml')
    pressure_source=ROOT/'data/intermediate/era5/reanalysis_pl.nc'
    surface_source=ROOT/'data/intermediate/era5/reanalysis_sl.nc'
    ds=normalize(xr.load_dataset(pressure_source)).sel(time=slice(case['start'][:19],case['end'][:19]))
    sl=normalize(xr.load_dataset(surface_source))
    ds,weights=subset(ds,case['area'])
    sl,sw=subset(sl,case['area'])
    terrain=float((sl.z/9.80665).weighted(sw).mean(('latitude','longitude','time')))
    pressure=ds.level*100.
    rv=ds.q/(1-ds.q)
    rho=pressure/(287.04*ds.t*(1+rv/(287.04/461.5)))*(1+rv)
    meteorology=xr.Dataset({'temperature':ds.t,'height':ds.z/9.80665-terrain,'w':-ds.w/(rho*9.80665)})
    meteorology=meteorology.weighted(weights).mean(('latitude','longitude'))
    order=np.argsort(meteorology.height.mean('time').values)
    t=meteorology.temperature.values[:,order];z=meteorology.height.values[:,order];w=meteorology.w.values[:,order]
    def heights(target):
        return np.array([crossing(row[None,:],zz,target)[0] for row,zz in zip(t,z)])
    z15,z5,z20=heights(258.15),heights(268.15),heights(253.15)
    center=float(np.nanmedian(z15));sigma=float(.5*np.nanmedian(z20-z5))
    if not (center>0 and sigma>0 and center-3*sigma>0):
        raise ValueError('Mixed-phase pulse envelope intersects ground; review meteorological support')
    gaussian=np.exp(-.5*((z-center)/sigma)**2)
    mixed=(t>=253.15)&(t<=268.15)
    reference=float(np.quantile(w[mixed & (w>0)],.95))
    score=np.sum(np.maximum(w,0)*gaussian,axis=1)/np.sum(gaussian,axis=1)
    times=(pd.DatetimeIndex(meteorology.time.values)-pd.Timestamp(case['start']).tz_localize(None)).total_seconds().values
    eligible=(times>=1800)&(times<=172800-1800)
    peak=int(np.flatnonzero(eligible)[np.argmax(score[eligible])])
    left=right=peak
    while left>0 and score[left-1]>=.5*score[peak]:left-=1
    while right<len(score)-1 and score[right+1]>=.5*score[peak]:right+=1
    width=float(times[right]-times[left]+3600.)
    duration=float(np.clip(width,3600.,10800.))
    duration=min(duration,2*times[peak],2*(172800-times[peak]))
    start=float(times[peak]-duration/2)
    with xr.open_dataset(ROOT/'runs/argentina_20181213/thompson_C0/control/output.nc',decode_times=False) as control:
        grid=control.lev.values.copy()
    edge=np.zeros(len(grid)+1)
    for k,c in enumerate(grid):edge[k+1]=2*c-edge[k]
    dz=np.diff(edge)
    g=np.exp(-.5*((grid-center)/sigma)**2)
    cap=float(.20*np.min(dz/np.maximum(g,1.e-300))/30.)
    strong=min(reference,cap)
    config={'case':case['case'],'design_status':'frozen_before_runs_and_any_ENTLN_access',
        'regimes':{'ERA5':0.,'weak':.5*strong,'strong':strong},
        'pulse':{'start_s':start,'duration_s':duration,'center_m_AGL':center,'sigma_m':sigma},
        'design_basis':{'amplitude':'0, 0.5, 1.0 times 95th percentile of positive ERA5 geometric w in -5 to -20 C; capped by incremental Courant<=0.20',
            'reference_positive_w_p95_m_s':reference,'incremental_courant_cap_amplitude_m_s':cap,
            'timing':'Maximum hourly Gaussian-weighted positive ERA5 ascent in fixed column; centered pulse duration from contiguous half-maximum episode, bounded to 1-3 h',
            'duration_unbounded_s':width,'height':'Median ERA5 -15 C height AGL; sigma=half median thickness from -5 to -20 C',
            'interpretation_scope':'Prescribed advective sensitivity, not a dynamically balanced reconstruction or resolved convective plume',
            'centered_Euler':'Courant is an exposure diagnostic, not a stability proof for centered forward Euler; numerical and clipping checks remain mandatory'},
        'forcing_sha256':sha256(FORCING),'meteorology_sources':{str(p.relative_to(ROOT)):sha256(p) for p in [pressure_source,surface_source]},
        'historical_pulse_source_sha256':sha256(ROOT/'archive/goamazon/src/updraft_forcing.F90'),
        'initial_conditions':'Exact unmodified full-case forcing; no saturation adjustment',
        'mynn':'type=0, explicit plume MP=0, w2=0; pulse is not supplied to MYNN as vertical velocity',
        'continuity':'Horizontal advection and pressure remain prescribed ERA5. Added vertical motion is not balanced by a horizontal mass-divergence adjustment; diagnose this inconsistency explicitly.',
        'primary_ensemble':'40-member T/humidity ensemble remains unchanged; only 18 deterministic runs',
        'ENTLN_accessed':False}
    CONFIG.write_text(yaml.safe_dump(config,sort_keys=False),encoding='utf-8')
    OUT.mkdir(parents=True,exist_ok=True)
    write_json(OUT/'design_frozen.json',{'created_utc':pd.Timestamp.now(tz='UTC').isoformat(),'configuration':config,'configuration_sha256':sha256(CONFIG)})
    meteorology.to_netcdf(OUT/'design_meteorology.nc')
    pd.DataFrame({'time':meteorology.time.values,'ascent_score_m_s':score,'z_minus15C_m_AGL':z15,
                  'z_minus5C_m_AGL':z5,'z_minus20C_m_AGL':z20}).to_csv(OUT/'design_basis.csv',index=False)
    print(json.dumps(config['regimes']),config['pulse'])


def run():
    cfg=load_config(CONFIG);frozen=json.loads((OUT/'design_frozen.json').read_text())
    if sha256(CONFIG)!=frozen['configuration_sha256'] or sha256(FORCING)!=cfg['forcing_sha256']:
        raise ValueError('Frozen design or complete forcing changed')
    executable=ROOT/'build/updraft_sensitivity/mp_offline_driver'
    baseline=git('show','26e13bc:src/mp_offline_driver.F90')
    current=(ROOT/'src/mp_offline_driver.F90').read_text(encoding='utf-8')
    pattern=r'  subroutine apply_forcing_and_advance\(.*?  end subroutine apply_forcing_and_advance'
    assert re.search(pattern,baseline,re.S).group()==re.search(pattern,current,re.S).group()
    original_checks=ROOT/'runs/argentina_20181213/_execution/ensemble40/run_checks.json'
    ensemble_before=sha256(original_checks)
    def one(item):
        regime,configuration=item
        folder=RUNS/regime/configuration;folder.mkdir(parents=True,exist_ok=True)
        output=folder/'output.nc';p=cfg['pulse'];amplitude=cfg['regimes'][regime]
        text=namelist(configuration,FORCING,output).replace('write_water_budget=.true.,','write_water_budget=.true.,write_vertical_diagnostics=.true.,')
        text=text.replace('/\n',f' pulse_amplitude={amplitude:.17g},pulse_start={p["start_s"]:.17g},pulse_duration={p["duration_s"]:.17g},\n pulse_center={p["center_m_AGL"]:.17g},pulse_sigma={p["sigma_m"]:.17g},\n/\n')
        identity={'executable_sha256':sha256(executable),'forcing_sha256':sha256(FORCING),'design_sha256':sha256(CONFIG),'namelist':text}
        tag=folder/'identity.json'
        if output.exists():
            if json.loads(tag.read_text())!=identity:raise ValueError('Existing run identity differs; retain it and use another experiment')
            with NETCDF_LOCK:validate_output(output,172800.)
        else:
            (folder/'namelist.mp_offline').write_text(text,encoding='utf-8');write_json(tag,identity)
            execute(executable,folder)
        with NETCDF_LOCK:
            checks=validate_output(output,172800.)
            with xr.open_dataset(output,decode_times=False) as native, xr.open_dataset(str(output)+'.vertical.nc',decode_times=False) as diag:
                assert diag.sizes['time']==5760 and native.sizes['time']==5761
                assert diag.attrs['mynn_edmf_type']==0 and diag.attrs['mynn_uedmf_mp']==0
                expected=pulse(diag.lev.values,diag.time.values,p['start_s'],p['duration_s'],p['center_m_AGL'],p['sigma_m'],amplitude)
                np.testing.assert_allclose(diag.w_perturbation,expected,atol=2e-15,rtol=2e-13)
                np.testing.assert_array_equal(diag.w_effective.values,diag.w_LS.values+diag.w_perturbation.values)
                if regime=='ERA5':
                    with xr.open_dataset(control_directory(configuration)/'output.nc',decode_times=False) as old:
                        for variable in old.variables:np.testing.assert_array_equal(native[variable].values,old[variable].values)
                max_courant=float((abs(diag.w_effective)*30/diag.layer_thickness).max())
        result={'regime':regime,'configuration':configuration,'folder':str(folder.relative_to(ROOT)),**checks,
                'max_layer_courant':max_courant,'control_bit_identical':regime=='ERA5','vertical_sha256':sha256(str(output)+'.vertical.nc')}
        write_json(folder/'checks.json',result)
        print('Finished',regime,configuration,flush=True)
        return result
    items=[(regime,configuration) for regime in cfg['regimes'] for configuration in CONFIGURATIONS]
    with ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(one,items))
    assert sha256(original_checks)==ensemble_before
    write_json(OUT/'run_checks.json',{'runs':results,'design_sha256':sha256(CONFIG),'integrator_unchanged':True,
        'six_ERA5_controls_bit_identical':True,'main_ensemble_checks_sha256_unchanged':ensemble_before,
        'ENTLN_accessed':False,'run_count':len(results)})


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['design','run','postprocess'])
    args=parser.parse_args()
    if args.action=='design':design()
    elif args.action=='run':run()
    else:
        from scm_lightning.evaluation.updraft_sensitivity import main as postprocess
        postprocess()


if __name__=='__main__':main()
