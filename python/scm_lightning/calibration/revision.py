"""Six immutable 2019-only proportional OLS coefficients; no December fitting."""
import importlib.metadata
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scm_lightning.calibration.independent import normalize_predictors
from scm_lightning.forcing.contracts import area_km2
from scm_lightning.provenance import sha256,write_json

MU_METHOD='MUCAPExP_300hPa_SCM'
METHODS=[MU_METHOD,'IxG','McCaulF3_LS_1D_adapted']
CONSTANTS={'k1':.042,'k2':.20,'r1':.95,'r2':.05}


def normalize(frame,config):
    f=frame.copy();f['method']=f.method.replace({MU_METHOD:'CAPExP'})
    f=normalize_predictors(f,config);f['method']=f.method.replace({'CAPExP':MU_METHOD})
    return f


def clock(frame,config):
    f=frame.copy();f['time']=pd.to_datetime(f.time,utc=True,format='mixed')
    if ((f.time<pd.Timestamp(config['start']))|(f.time>=pd.Timestamp(config['end']))).any():
        raise ValueError('Only October 2019 intervals may enter calibration')
    return f


def fit_controls(config,predictors,observations,five_predictors,five_observations):
    if config.get('event_id')!='argentina_20191028' or config.get('role')!='calibration_only':raise ValueError('2019 calibration event only')
    if pd.Timestamp(config['start'])!=pd.Timestamp('2019-10-28T00:00Z') or pd.Timestamp(config['end'])!=pd.Timestamp('2019-10-30T00:00Z'):
        raise ValueError('Only the fixed October 28-29 event may calibrate')
    if config.get('area') is None or config['area_selection']['status']!='fixed_before_lightning_evaluation':raise ValueError('Fixed weather-only area required')
    if config['calibration']['published_McCaul_coefficients']!=CONSTANTS:raise ValueError('Published constants changed')
    for f in [predictors,five_predictors]:
        if set(f.configuration)!={'TH_C0','WS_C0'} or set(f.member)!={0}:raise ValueError('Only the two deterministic C0 controls may calibrate')
    p,p5=(normalize(clock(f,config),config) for f in [predictors,five_predictors])
    o,o5=(clock(f,config) for f in [observations,five_observations])
    for f in [o,o5]:
        if f.time.duplicated().any():raise ValueError('Duplicate observed intervals')
        if 'covered' not in f:raise ValueError('Explicit observation coverage required')
    area=area_km2(config['area']);rows=[]
    for scheme in ['TH','WS']:
        for method in METHODS:
            minutes=5 if method==METHODS[2] else 1
            x=p5 if minutes==5 else p;obs=o5 if minutes==5 else o
            x=x[(x.configuration==scheme+'_C0')&(x.method==method)]
            f=x.merge(obs[['time','count','covered']],on='time',validate='one_to_one')
            expected=pd.date_range(config['start'],config['end'],freq=f'{minutes}min',inclusive='left')
            if len(f)!=len(expected) or set(f.time)!=set(expected):raise ValueError('Complete declared calibration clock required')
            # Exactly a=sum(X*N)/sum(X^2), in whole fixed-box interval counts.
            values=f.predictor.to_numpy(float)*area;target=f['count'].to_numpy(float)
            valid=np.isfinite(values)&np.isfinite(target)&f.covered.to_numpy(bool)
            if np.any(values[valid]<0) or np.any(target[valid]<0):raise ValueError('Nonnegative predictor and target required')
            denominator=float(np.dot(values[valid],values[valid]));numerator=float(np.dot(values[valid],target[valid]))
            status='fitted' if denominator>0 else ('identically_zero_unrecoverable' if valid.any() else 'no_valid_calibration_data')
            scale=max(0.,numerator/denominator) if denominator>0 else np.nan
            rows.append({'method':method,'microphysics':scheme,'configuration':scheme+'_C0','member':0,
                'scale':scale,'status':status,'fit_method':'nonnegative_proportional_OLS_origin_all_valid_intervals',
                'fit_interval_minutes':minutes,'n_intervals':len(f),'n_valid':int(valid.sum()),
                'n_missing_or_uncovered':int((~valid).sum()),'n_positive_predictor':int((values[valid]>0).sum()),
                'n_observed_active':int((target[valid]>0).sum()),'OLS_numerator':numerator,'OLS_denominator':denominator,
                'fit_support':'whole fixed-box interval count equivalent: normalized diagnostic density times calibration area',
                'coefficient_units':{MU_METHOD:'flashes J-1','IxG':'flashes km-2 (kg2 s m-5)-1',METHODS[2]:'observed flashes / nominal-cell flash equivalent'}[method],
                'target_units':f'flashes per {minutes}min','calibration_area_km2':area,
                'start':config['start'],'end_exclusive':config['end'],'published_constants_changed':False})
    selected=o[o.covered & o['count'].notna()]
    baseline=float(selected['count'].mean()/area) if len(selected) else np.nan
    return pd.DataFrame(rows),baseline


def freeze(destination,config,coefficients,baseline,sources):
    destination=Path(destination)
    if (destination/'frozen.json').exists():raise FileExistsError('Frozen calibration cannot be overwritten')
    expected={(m,s) for m in METHODS for s in ['TH','WS']}
    if len(coefficients)!=6 or set(zip(coefficients.method,coefficients.microphysics))!=expected:raise ValueError('Exactly six identifiable slots required')
    if set(coefficients.configuration)!={'TH_C0','WS_C0'} or set(coefficients.member)!={0}:raise ValueError('Invalid calibration controls')
    if any('20181213' in str(p) for p in sources):raise ValueError('December data cannot enter coefficient freeze')
    destination.mkdir(parents=True,exist_ok=True)
    table=destination/'coefficients.csv';coefficients.to_csv(table,index=False)
    snapshot=destination/'calibration_implementation_at_freeze.py';snapshot.write_bytes(Path(__file__).read_bytes())
    write_json(destination/'frozen.json',{'schema_version':3,'calibration_event':'argentina_20191028',
        'controls':['TH_C0','WS_C0'],'member':0,'coefficient_count':6,'interval_minutes':1,
        'start':config['start'],'end_exclusive':config['end'],'area_wsen':config['area'],'area_km2':area_km2(config['area']),
        'constant_density_flashes_km2_per_interval':baseline,'coefficient_sha256':sha256(table),
        'implementation_sha256':sha256(snapshot),'evaluation_observations_used':False,
        'published_constants':CONSTANTS,'McCaul_fit_minutes':5,'McCaul_prediction_minutes':1,
        'CAPE_definition':'official MetPy most_unstable_cape_cin; lowest 300hPa max theta-e; not Romps isentropic or SBCAPE',
        'shared_across':'each microphysics: C0 CL CF; 40 members plus separate controls; no per-member fit',
        'all_identifiable':bool((coefficients.status=='fitted').all()),
        'versions':{name:importlib.metadata.version(name) for name in ['MetPy','numpy','pandas','xarray','scipy']},
        'sources':[{'path':str(p),'sha256':sha256(p)} for p in sources]})


def load_frozen(destination):
    destination=Path(destination);meta=json.loads((destination/'frozen.json').read_text())
    if meta['calibration_event']!='argentina_20191028' or meta['controls']!=['TH_C0','WS_C0'] or meta['evaluation_observations_used']:
        raise ValueError('Invalid six-coefficient provenance')
    if sha256(destination/'coefficients.csv')!=meta['coefficient_sha256']:raise ValueError('Coefficient hash mismatch')
    if sha256(destination/'calibration_implementation_at_freeze.py')!=meta['implementation_sha256']:raise ValueError('Implementation snapshot changed')
    return pd.read_csv(destination/'coefficients.csv'),meta


def apply(predictors,coefficients,area,comparison='microphysics_specific'):
    expected={(m,s) for m in METHODS for s in ['TH','WS']}
    if len(coefficients)!=6 or coefficients.duplicated(['method','microphysics']).any() or set(zip(coefficients.method,coefficients.microphysics))!=expected:
        raise ValueError('Exactly six unique method/microphysics coefficients required')
    f=predictors.copy();f['microphysics']=f.configuration.str[:2]
    if not set(f.configuration)<={s+'_'+c for s in ['TH','WS'] for c in ['C0','CL','CF']}:raise ValueError('Unknown physical configuration')
    if comparison=='microphysics_specific':
        factors=coefficients[['method','microphysics','scale','status']]
    elif comparison=='common_TH_same_MUCAPE_OLS':
        shared=coefficients[coefficients.microphysics=='TH'][['method','scale','status']]
        factors=pd.concat([shared.assign(microphysics=s) for s in ['TH','WS']],ignore_index=True)
    else:raise ValueError('Unknown independently labeled comparison')
    f=f.merge(factors.rename(columns={'status':'coefficient_status'}),on=['method','microphysics'],validate='many_to_one')
    f['prediction_density']=f.predictor*f.scale
    zero=(f.predictor==0)&(f.coefficient_status=='identically_zero_unrecoverable')
    f.loc[zero,'prediction_density']=0.
    f['prediction']=f.prediction_density*area;f['comparison']=comparison
    return f
