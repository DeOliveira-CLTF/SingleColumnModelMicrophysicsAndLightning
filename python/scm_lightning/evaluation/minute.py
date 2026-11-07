"""Native 30s output -> real one-minute bins; independent 2019-only calibration.

Existing 15min publication and frozen v1 coefficients remain archived. This
module does not split 15min counts or interpolate 15min lightning forecasts.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import os
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import xarray as xr
from scm_lightning.provenance import ROOT,sha256,write_json,resolve_path,git
from scm_lightning.config import load_config
from scm_lightning.calibration.independent import METHODS,normalize_predictors,fit,freeze,load_frozen,apply
from scm_lightning.evaluation.independent import control_predictors,evaluate_predictions,metric_values,tables
from scm_lightning.forcing.contracts import area_km2
from scm_lightning.observations import aggregate

OUT=ROOT/'results/argentina_20181213/article'
PREVIOUS=ROOT/'results/argentina_20181213/legacy/article_independent_15min_20261106'
FROZEN=ROOT/'results/argentina_20191028/calibration/v2_1min'
CAL=ROOT/'configs/cases/argentina_20191028.yaml'
CASE=ROOT/'configs/cases/argentina_20181213.yaml'
CAPE_CACHE=ROOT/'data/intermediate/cape_exact/most_cape.sqlite3'


def archive_previous():
    if PREVIOUS.exists():return
    if not OUT.resolve().is_relative_to((ROOT/'results/argentina_20181213').resolve()):raise ValueError('Unsafe archive source')
    files=[{'path':str(p.relative_to(OUT)),'sha256':sha256(p)} for p in OUT.rglob('*') if p.is_file()]
    shutil.copytree(OUT,PREVIOUS)
    for record in files:
        if sha256(PREVIOUS/record['path'])!=record['sha256']:raise IOError('15min archive verification failed')
    write_json(PREVIOUS/'archive_identity.json',{'original_interval_minutes':15,'files':files,'native_simulations_changed':False})
    for record in files:(OUT/record['path']).unlink()


def config(path):
    c=load_config(path);c['interval_minutes']=1
    return c


def observations(c):
    folder=ROOT/'results'/c['event_id']/'observations'
    target=folder/'counts_1min.csv'
    if target.exists():return
    frame=pd.read_csv(folder/'flashes.csv.gz',usecols=['start','latitude','longitude','kind'])
    frame['start']=pd.to_datetime(frame.start,utc=True,format='mixed')
    original=load_config(CAL if c['event_id']=='argentina_20191028' else CASE)
    n=len(pd.date_range(c['start'],c['end'],freq='min',inclusive='left'))
    coverage=original['observation_coverage']
    if coverage=='full_case' if isinstance(coverage,str) else False:covered=[True]*n
    else:
        multiplier=original['interval_minutes']
        covered=np.repeat(np.asarray(coverage,bool),multiplier).tolist()
    result=aggregate(frame,c['start'],c['end'],c['area'],covered,1)
    result.to_csv(target,index=False)
    parts=[]
    for kind in ['IC','CG']:
        parts.append(aggregate(frame[frame.kind==kind],c['start'],c['end'],c['area'],covered,1).assign(type=kind))
    pd.concat(parts).to_csv(folder/'counts_by_type_1min.csv',index=False)
    write_json(folder/'minute_manifest.json',{'source_grouped_flashes_sha256':sha256(folder/'flashes.csv.gz'),
        'interval_minutes':1,'bins':n,'counts_sha256':sha256(target),'grouping':c['grouping'],
        'flag_mapping':c['flag_mapping'],'area_wsen':c['area'],'uniform_split_of_15min_counts':False})


def member(run):
    c=config(CASE);folder=resolve_path(run['folder']);source=folder/'output.nc'
    if 'updraft_sensitivity' in source.parts:raise ValueError('Deterministic pulse scenarios are not initial-condition ensemble members')
    target=OUT/'data/members'/run['configuration']/f'{run["member"]:02}'
    tag=target/'minute_identity.json'
    if tag.exists():
        identity=json.loads(tag.read_text())
        manifest=json.loads((target/'manifest.json').read_text())
        if identity['source_sha256']==run['sha256'] and identity['interval_minutes']==1 and manifest['forcing_sha256']==sha256(ROOT/'data/processed/forcing/LEONA_SCM.nc'):return
        raise ValueError('Mismatched minute cache')
    if sha256(source)!=run['sha256']:raise ValueError('Changed native simulation')
    old=PREVIOUS/'data/members'/run['configuration']/f'{run["member"]:02}'/'thermodynamics_15min.csv'
    control_predictors(source,ROOT/'data/processed/forcing/LEONA_SCM.nc',c,target,
        configuration=run['configuration'],member_id=run['member'],export_native=False,cape_cache=old,cape_memo_path=CAPE_CACHE)
    write_json(tag,{'source_sha256':run['sha256'],'interval_minutes':1,'native_snapshots_seconds':30,
        'CAPE':'native profile at minute start; same most_cape pseudoadiabatic definition as 15min version',
        'uniform_disaggregation':False,'code_sha256':sha256(__file__)})


def seed_exact_cache(runs):
    from scm_lightning.lightning.cape_cache import CapeCache
    memo=CapeCache(CAPE_CACHE)
    for run in runs:
        folder=resolve_path(run['folder'])
        source=OUT/'data/members'/run['configuration']/f'{run["member"]:02}'/'thermodynamics_1min.csv'
        cadence=2
        if not source.exists() or not (source.parent/'minute_identity.json').exists():
            source=PREVIOUS/'data/members'/run['configuration']/f'{run["member"]:02}'/'thermodynamics_15min.csv';cadence=30
        frame=pd.read_csv(source)
        with xr.open_dataset(folder/'output.nc',decode_times=False) as d:
            p=d.pressure.values[::cadence][:-1];t=d.theta.values[::cadence][:-1]*(p/1e5)**(2/7);r=d.qv.values[::cadence][:-1]
        if len(frame)!=len(p):raise ValueError('CAPE cache seed is not aligned to native snapshots')
        for pressure,temperature,rv,value in zip(p,t,r,frame.cape_J_kg):memo.seed(pressure,temperature,rv,value)
    memo.close()


def prepare(workers=12):
    archive_previous();OUT.mkdir(parents=True,exist_ok=True)
    write_json(OUT/'minute_mode.json',{'interval_minutes':1,'version':'native_30s_to_1min_v2','previous_results':str(PREVIOUS)})
    for path in [CASE,CAL]:observations(config(path))
    c=config(CASE);root=ROOT/c['ensemble']['active_directory']
    plan=json.loads((root/'manifest.json').read_text())['runs']
    checks={(r['configuration'],r['member']):r for r in json.loads((root/'run_checks.json').read_text())}
    runs=[dict(r,sha256=checks[(r['configuration'],r['member'])]['sha256']) for r in plan]
    if len(runs)!=246:raise ValueError('All six existing 40-member configurations + controls required')
    seed_exact_cache(runs)
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs=[pool.submit(member,r) for r in runs]
        for i,job in enumerate(as_completed(jobs),1):
            job.result();write_json(OUT/'minute_progress.json',{'completed':i,'expected':246,'SCM_reruns':0,'interval_minutes':1})
            if i%6==0:print(f'Native one-minute diagnostics: {i}/246',flush=True)


def calibrate():
    if (FROZEN/'frozen.json').exists():load_frozen(FROZEN);return
    c=config(CAL);run=ROOT/'runs/argentina_20191028/thompson_C0/control'
    identity=json.loads((run/'run_identity.json').read_text())
    if identity['configuration']!='TH_C0' or identity['member']!=0 or not identity['result'].get('physical_preflight_passed') or sha256(run/'output.nc')!=identity['result']['sha256']:
        raise ValueError('Validated deterministic 2019 TH_C0 required')
    target=ROOT/'results/argentina_20191028/diagnostics_1min'
    old=ROOT/'results/argentina_20191028/diagnostics/thermodynamics_15min.csv'
    p,p5=control_predictors(run/'output.nc',ROOT/c['paths']['forcing']/'LEONA_SCM.nc',c,target,cape_cache=old)
    obs=ROOT/'results/argentina_20191028/observations'
    coef,baseline=fit(c,p,pd.read_csv(obs/'counts_1min.csv'),p5,pd.read_csv(obs/'counts_5min.csv'))
    FROZEN.mkdir(parents=True,exist_ok=True);write_json(FROZEN/'configuration_used.json',c)
    freeze(FROZEN,c,coef,baseline,[FROZEN/'configuration_used.json',CAL,run/'run_identity.json',run/'output.nc',
        target/'manifest.json',obs/'minute_manifest.json',obs/'counts_1min.csv',obs/'counts_5min.csv'])


def evaluate():
    c=config(CASE);coef,meta=load_frozen(FROZEN)
    paths=sorted((OUT/'data/members').glob('*/*/predictors_1min.csv'))
    if len(paths)!=246:raise ValueError('Complete native one-minute ensemble required')
    p=normalize_predictors(pd.concat([pd.read_csv(f) for f in paths],ignore_index=True),config(CAL))
    p=apply(p,coef,area_km2(c['area']))
    o=pd.read_csv(ROOT/'results/argentina_20181213/observations/counts_1min.csv')
    constant=meta['constant_density_flashes_km2_per_interval']*area_km2(c['area'])
    metric,stats,mask=evaluate_predictions(p,o,constant,minutes=1)
    data=OUT/'data';data.mkdir(parents=True,exist_ok=True)
    p.to_csv(data/'independent_member_predictions.csv.gz',index=False);stats.to_csv(data/'ensemble_statistics.csv',index=False)
    o.to_csv(data/'observations_1min.csv',index=False);mask.to_csv(data/'common_valid_intervals.csv',index=False)
    diagnostics=[];control_thermodynamics=[]
    for path in sorted((OUT/'data/members').glob('*/*/thermodynamics_1min.csv')):
        f=pd.read_csv(path);configuration=path.parent.parent.name;member_id=int(path.parent.name)
        diagnostics.append({'configuration':configuration,'member':member_id,'interval_minutes':1,'n':len(f),
            'n_missing_cape':int(f.cape_J_kg.isna().sum()),'n_positive_cape':int((f.cape_J_kg>0).sum()),
            'n_positive_precipitation':int((f.precipitation_amount_kg_m2>0).sum()),
            'n_positive_product':int((f.CAPExP_J_m2>0).sum()),'max_cape_J_kg':f.cape_J_kg.max(),
            'total_precipitation_mm':f.precipitation_amount_kg_m2.sum(),'coefficient_applied':False})
        if member_id==0:control_thermodynamics.append(f.assign(configuration=configuration,member=member_id))
    pd.DataFrame(diagnostics).to_csv(data/'CAPE_precipitation_audit_1min.csv',index=False)
    pd.concat(control_thermodynamics,ignore_index=True).to_csv(data/'CAPE_precipitation_controls_1min.csv.gz',index=False)
    cube=p[['time','method','configuration','member','predictor','prediction_density','prediction']].copy()
    cube.time=pd.to_datetime(cube.time,utc=True).dt.tz_localize(None)
    ds=cube.set_index(['time','method','configuration','member']).to_xarray()
    ds['calibration_coefficient']=('method',coef.set_index('method').scale.reindex(ds.method.values).values)
    ds['observation_count']=('time',o.assign(time=pd.to_datetime(o.time,utc=True).dt.tz_localize(None)).set_index('time')['count'].reindex(ds.time.values).values)
    ds['prediction_count']=ds.prediction.copy()
    ds.prediction_count.attrs.update(units='1',long_name='forecast flashes in the one-minute interval')
    ds.prediction.attrs.update(units='min-1',long_name='total IC+CG flashes per one-minute interval; numerically flashes/min')
    ds.prediction_density.attrs.update(units='km-2 min-1',interval_seconds=60)
    ds.observation_count.attrs.update(units='1',long_name='observed total IC+CG flashes during the one-minute interval')
    ds.predictor.attrs['units']='method-dependent native proxy time integral over 60 seconds; see Table_02'
    ds['time_bounds']=(('time','bounds'),np.column_stack([ds.time.values,ds.time.values+np.timedelta64(60,'s')]))
    ds.time.attrs.update(bounds='time_bounds',long_name='start of the UTC one-minute interval')
    ds.attrs.update(interval_seconds=60,native_snapshot_seconds=30,calibration_event='argentina_20191028',
        evaluation_event='argentina_20181213',calibration_configuration='TH_C0',calibration_member=0,
        coefficient_sha256=meta['coefficient_sha256'],uniform_split_of_15min_counts='false',
        evaluation_area_km2=area_km2(c['area']),calibration_area_km2=meta['area_km2'],
        evaluation_area_wsen=json.dumps(c['area']),calibration_area_wsen=json.dumps(meta['area_wsen']),
        rate_convention='one-minute interval count / 1 minute; cumulative flashes integrate rate over these bounds')
    for key in ['method','configuration']:ds[key]=xr.Variable((key,),np.asarray(ds[key].values,dtype='S40'))
    time_encoding={'units':'minutes since '+pd.Timestamp(c['start']).strftime('%Y-%m-%d %H:%M:%S'),
        'calendar':'proleptic_gregorian'}
    for name in ['time','time_bounds']:ds[name].encoding.update(time_encoding)
    pending=data/'independent_predictions_1min.pending.nc';ds.to_netcdf(pending);pending.replace(data/'independent_predictions_1min.nc')
    tables(c,coef)
    metric.to_csv(OUT/'tables/Table_04_independent_evaluation.csv',index=False)
    metric[metric.configuration.isin(['TH_C0','WS_C0'])].to_csv(OUT/'tables/Table_04_reference_C0.csv',index=False)
    contrasts=[]
    for method in METHODS:
        f=metric[metric.method==method].set_index('configuration')
        pairs=[('WS_'+mode,'TH_'+mode) for mode in ['C0','CL','CF']]+[(s+'_'+m,s+'_C0') for s in ['TH','WS'] for m in ['CL','CF']]
        for a,b in pairs:
            for key in ['rmse','mae','skill_zero','skill_constant','bias','accumulated_error']:
                contrasts.append({'method':method,'configuration_a':a,'configuration_b':b,'contrast':a+' minus '+b,'metric':key,'difference':f.loc[a,key]-f.loc[b,key]})
    pd.DataFrame(contrasts).to_csv(OUT/'tables/Table_05_microphysics_sensitivity.csv',index=False)
    conventions=json.loads((OUT/'tables/metric_conventions.json').read_text())
    for key in ['rmse','mae','bias']:conventions[key]='flashes/min; '+('prediction minus observation' if key=='bias' else 'one-minute bins')
    conventions['interval_minutes']=1;write_json(OUT/'tables/metric_conventions.json',conventions)
    write_json(OUT/'independent_manifest.json',{'status':'evaluated','interval_minutes':1,'calibration_event':'argentina_20191028',
        'calibration_directory':str(FROZEN),'coefficient_sha256':meta['coefficient_sha256'],
        'constant_baseline_flashes_per_min':constant,'common_valid_intervals':int(mask.common_valid.sum()),
        'physical_2018_runs_repeated':0,'physical_2019_runs_repeated':0,'uniform_disaggregation':False,
        'calibration_1min':'2019 TH_C0 only; McCaul additional peak mapping retains published 5min peak support'})


def check():
    coef,meta=load_frozen(FROZEN)
    p=pd.read_csv(OUT/'data/independent_member_predictions.csv.gz');o=pd.read_csv(OUT/'data/observations_1min.csv')
    constant=meta['constant_density_flashes_km2_per_interval']*area_km2(config(CASE)['area'])
    actual,stats,mask=evaluate_predictions(p,o,constant,minutes=1)
    pd.testing.assert_frame_equal(actual,pd.read_csv(OUT/'tables/Table_04_independent_evaluation.csv'),check_dtype=False,rtol=1e-12)
    assert len(coef)==3 and (p.groupby('method').scale.nunique()==1).all()
    tags=list((OUT/'data/members').glob('*/*/minute_identity.json'))
    assert len(tags)==246
    forcing_hash=sha256(ROOT/'data/processed/forcing/LEONA_SCM.nc')
    for tag in tags:
        row=json.loads(tag.read_text());identity=json.loads((tag.parent/'manifest.json').read_text())
        assert row['source_sha256']==identity['simulation_sha256'] and identity['native_output_unchanged']
        assert identity['forcing_sha256']==forcing_hash
    assert not git('diff','--','src')
    write_json(OUT/'audit/one_minute_validation.json',{'status':'passed','interval_seconds':60,'source_step_seconds':30,
        'members_per_configuration':40,'separate_control':True,'configurations':6,'metric_rows':18,
        'common_bins':int(mask.common_valid.sum()),'coefficient_sha256':meta['coefficient_sha256'],
        'uniform_15min_disaggregation':False,'native_physics_changed':False})
    write_json(OUT/'pending_dependencies.json',{'blocking':[],
        'initial_boundary_sounding_substitution':'not applied: provided soundings outside fixed boxes and/or after initialization; ERA5 retained',
        'independent_precipitation_validation':'Figure 6 uses forcing ERA5 as explicitly requested, not independent observations'})
    print('One-minute independent audit passed.',flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','calibrate','evaluate','check']);ap.add_argument('--workers',type=int,default=12);a=ap.parse_args()
    active=OUT/'independent_manifest.json'
    if active.exists() and json.loads(active.read_text()).get('version')==3:
        raise ValueError('Historical v2 driver cannot overwrite v3; use scripts/evaluate.py independent or the preserved v2 archive')
    if a.stage=='prepare':prepare(a.workers)
    else:{'calibrate':calibrate,'evaluate':evaluate,'check':check}[a.stage]()


if __name__=='__main__':main()
