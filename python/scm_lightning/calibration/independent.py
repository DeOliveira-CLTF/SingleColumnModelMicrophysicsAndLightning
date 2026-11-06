"""2019-only deterministic calibration; immutable coefficients transferred as densities.

No input, model choice or coefficient may be selected using evaluation observations.
McCaul's original constants are retained. The extra single-event peak mapping is
explicitly an SCM extension, not the published multi-storm F1/F2 regressions.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scm_lightning.evaluation import fit_scale
from scm_lightning.forcing.contracts import area_km2
from scm_lightning.provenance import sha256, write_json

METHODS = ['CAPExP', 'IxG', 'McCaulF3_LS_1D_adapted']
UNITS = {'CAPExP':'J km-2', 'IxG':'kg2 s m-5', 'McCaulF3_LS_1D_adapted':'published nominal-cell flashes km-2 per bin'}


def nominal_mccaul_cell_area(config):
    grid = config['calibration']
    step = grid['published_McCaul_grid_spacing_degrees']
    lat = grid['published_McCaul_reference_latitude_degrees']
    return area_km2([-step/2,lat-step/2,step/2,lat+step/2])


def normalize_predictors(frame, config):
    """Native time integrals -> explicit transferable column/density support."""
    f = frame.copy()
    f['method'] = f.method.replace({'McCaulF3':'McCaulF3_LS_1D_adapted'})
    if not set(f.method) <= set(METHODS):
        raise ValueError('Only the three principal methods enter this calibration')
    f.loc[f.method=='CAPExP','predictor'] *= 1.e6  # J/m2 -> J/km2
    f.loc[f.method=='McCaulF3_LS_1D_adapted','predictor'] /= 300*nominal_mccaul_cell_area(config)
    f['predictor_units'] = f.method.map(UNITS)
    return f


def validate_context(config, predictors):
    if config.get('role') != 'calibration_only' or config.get('event_id') != 'argentina_20191028':
        raise ValueError('Calibration accepts only the declared independent October 2019 event')
    if config['calibration']['configuration'] != 'TH_C0' or config['calibration']['member'] != 0:
        raise ValueError('Only the predeclared deterministic Thompson C0 is allowed')
    if set(predictors.configuration) != {'TH_C0'} or set(predictors.member) != {0}:
        raise ValueError('Perturbed members and other physical configurations cannot calibrate')
    times = pd.to_datetime(predictors.time,utc=True)
    if ((times < pd.Timestamp(config['start'])) | (times >= pd.Timestamp(config['end']))).any():
        raise ValueError('Calibration predictors outside the independent event')
    if config['area'] is None or config['area_selection']['status'] != 'fixed_before_lightning_evaluation':
        raise ValueError('Weather-selected fixed area is required')


def fit(config, predictors, observations, five_predictors, five_observations):
    validate_context(config,predictors)
    validate_context(config,five_predictors)
    for f in [observations, five_observations]:
        t=pd.to_datetime(f.time,utc=True)
        if ((t<pd.Timestamp(config['start']))|(t>=pd.Timestamp(config['end']))).any():
            raise ValueError('Evaluation observations cannot enter calibration or its constant baseline')
        if f.time.duplicated().any():raise ValueError('Duplicate observation intervals')
        if not f.covered.all():raise ValueError('Complete explicitly confirmed calibration coverage required')
    p=normalize_predictors(predictors,config)
    p5=normalize_predictors(five_predictors,config)
    area=area_km2(config['area'])
    rows=[]
    for method in METHODS:
        minutes=5 if method=='McCaulF3_LS_1D_adapted' else 15
        x=p5 if minutes==5 else p
        obs=five_observations if minutes==5 else observations
        f=x[x.method==method].merge(obs[['time','count']],on='time',validate='one_to_one')
        expected=len(pd.date_range(config['start'],config['end'],freq=f'{minutes}min',inclusive='left'))
        if len(f)!=expected:raise ValueError('Calibration must retain the complete meteorological period')
        target=f['count'].to_numpy(float)/area
        values=f.predictor.to_numpy(float)
        valid=np.isfinite(values)&np.isfinite(target)
        if np.any(values[valid]<0) or np.any(target[valid]<0):raise ValueError('Nonnegative physical inputs required')
        if not valid.any():
            scale,status=np.nan,'no_valid_calibration_data'
        elif minutes==5:
            maximum=values[valid].max()
            scale,status=(float(target[valid].max()/maximum),'fitted') if maximum>0 else (np.nan,'identically_zero_unrecoverable')
        else:
            scale,status=fit_scale(values,target)
        rows.append({'method':method,'scale':scale,'status':status,
            'fit_method':'single_event_peak_density_ratio_5min_SCM_extension' if minutes==5 else 'proportional_OLS_origin_SCM_recalibration',
            'fit_interval_minutes':minutes,'configuration':'TH_C0','member':0,
            'n_intervals':len(f),'n_valid':int(valid.sum()),'n_positive_predictor':int((values[valid]>0).sum()),
            'n_observed_active':int((target[valid]>0).sum()),'predictor_min':float(values[valid].min()) if valid.any() else np.nan,
            'predictor_q50':float(np.median(values[valid])) if valid.any() else np.nan,
            'predictor_q95':float(np.quantile(values[valid],.95)) if valid.any() else np.nan,
            'predictor_max':float(values[valid].max()) if valid.any() else np.nan,
            'predictor_units':UNITS[method],'target_units':f'flashes km-2 per {minutes}min',
            'calibration_area_km2':area,'start':config['start'],'end_exclusive':config['end'],
            'published_constants_changed':False})
    # Independent density baseline transfers to any evaluation area without using its lightning.
    baseline=float(observations['count'].mean()/area)
    return pd.DataFrame(rows),baseline


def freeze(destination, config, coefficients, baseline, sources):
    destination=Path(destination)
    destination.mkdir(parents=True,exist_ok=True)
    manifest=destination/'frozen.json'
    if manifest.exists():raise FileExistsError('Frozen calibration is immutable; choose an explicit new version')
    if len(coefficients)!=3 or set(coefficients.method)!=set(METHODS):raise ValueError('Exactly three coefficients required')
    if set(coefficients.configuration)!={'TH_C0'} or set(coefficients.member)!={0}:raise ValueError('Only deterministic TH_C0 coefficients can be frozen')
    table=destination/'coefficients.csv'
    coefficients.to_csv(table,index=False)
    (destination/'calibration_implementation_at_freeze.py').write_bytes(Path(__file__).read_bytes())
    write_json(manifest,{'schema_version':1,'calibration_event':'argentina_20191028','configuration':'TH_C0','member':0,
        'start':config['start'],'end_exclusive':config['end'],'area_wsen':config['area'],'area_km2':area_km2(config['area']),
        'constant_density_flashes_km2_per_15min':baseline,'coefficient_sha256':sha256(table),
        'all_identifiable':bool((coefficients.status=='fitted').all()),'published_constants':config['calibration']['published_McCaul_coefficients'],
        'mccaul_nominal_cell_area_km2':nominal_mccaul_cell_area(config),
        'mccaul_support_limit':config['calibration']['published_McCaul_area_convention'],
        'literature_limit':config['calibration']['literature_limit'],'grouping':config['grouping'],'flag_mapping':config['flag_mapping'],
        'selection_sha256':config['area_selection']['decision_sha256'],
        'sources':[{'path':str(p),'sha256':sha256(p)} for p in sources],
        'code_sha256':sha256(__file__),'evaluation_observations_used':False})


def load_frozen(directory):
    directory=Path(directory)
    meta=json.loads((directory/'frozen.json').read_text(encoding='utf-8'))
    if meta['calibration_event']!='argentina_20191028' or meta['configuration']!='TH_C0' or meta['member']!=0 or meta['evaluation_observations_used']:
        raise ValueError('Invalid calibration provenance')
    if sha256(directory/'coefficients.csv')!=meta['coefficient_sha256']:raise ValueError('Frozen coefficient checksum mismatch')
    snapshot=directory/'calibration_implementation_at_freeze.py'
    if snapshot.exists() and sha256(snapshot)!=meta['code_sha256']:raise ValueError('Frozen calibration implementation checksum mismatch')
    return pd.read_csv(directory/'coefficients.csv'),meta


def apply(predictors, coefficients, area):
    if coefficients.method.duplicated().any() or set(coefficients.method)!=set(METHODS):
        raise ValueError('One shared coefficient per method is required')
    f=predictors.merge(coefficients[['method','scale','status']].rename(columns={'status':'coefficient_status'}),on='method',validate='many_to_one')
    f['prediction_density']=f.predictor*f.scale
    zero=(f.predictor==0)&(f.coefficient_status=='identically_zero_unrecoverable')
    f.loc[zero,'prediction_density']=0.
    f['prediction']=f.prediction_density*area
    return f
