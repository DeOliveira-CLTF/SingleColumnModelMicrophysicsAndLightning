"""Version 3: MUCAPE and microphysics-specific independent-event calibration."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import xarray as xr
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.config import load_config

CAL=ROOT/'configs/cases/argentina_20191028.yaml'
CASE=ROOT/'configs/cases/argentina_20181213.yaml'
FROZEN=ROOT/'results/argentina_20191028/calibration/v3_mucape_microphysics_OLS'
OUT=ROOT/'results/argentina_20181213/article'
PREVIOUS=ROOT/'results/argentina_20181213/legacy/article_v2_1min_most_cape'
MU_METHOD='MUCAPExP_300hPa_SCM'


def diagnostics_root(event):return ROOT/'results'/event/'diagnostics_v3'


def plan():
    from scm_lightning.provenance import resolve_path
    c=load_config(CASE);root=ROOT/c['ensemble']['active_directory']
    native=json.loads((root/'manifest.json').read_text())['runs']
    checks={(r['configuration'],r['member']):r for r in json.loads((root/'run_checks.json').read_text())}
    runs=[dict(event='argentina_20181213',folder=str(resolve_path(r['folder'])),configuration=r['configuration'],
        member=r['member'],sha256=checks[(r['configuration'],r['member'])]['sha256']) for r in native]
    if len(runs)!=246:raise ValueError('Six controls and 240 existing perturbed members required')
    runs.sort(key=lambda r:(r['member']!=0 or r['configuration'] not in ['TH_C0','WS_C0'],r['configuration'],r['member']))
    for scheme,name in [('TH','thompson'),('WS','wsm6')]:
        folder=ROOT/'runs/argentina_20191028'/(name+'_C0')/'control'
        identity=json.loads((folder/'run_identity.json').read_text())
        if not identity['result'].get('physical_preflight_passed'):raise ValueError('Invalid calibration control')
        runs.append(dict(event='argentina_20191028',folder=str(folder),configuration=scheme+'_C0',member=0,sha256=identity['result']['sha256']))
    return runs


def diagnose_run(run):
    from scm_lightning.lightning.thermodynamics import diagnose,DEFINITION,DEPTH_PA
    from scm_lightning.evaluation.independent import control_predictors
    source=Path(run['folder'])/'output.nc';event=run['event'];cfg=run['configuration'];member=run['member']
    if 'updraft_sensitivity' in source.parts:raise ValueError('Pulse scenarios cannot enter this ensemble')
    destination=diagnostics_root(event)/cfg/f'{member:02}';destination.mkdir(parents=True,exist_ok=True)
    algorithm=sha256(ROOT/'python/scm_lightning/lightning/thermodynamics.py')
    tag=destination/'identity.json'
    if tag.exists():
        previous=json.loads(tag.read_text())
        if previous['source_sha256']==run['sha256'] and previous['algorithm_sha256']==algorithm:return
        raise ValueError('MUCAPE cache identity mismatch')
    if sha256(source)!=run['sha256']:raise ValueError('Native source identity changed')
    c=load_config(CAL if event=='argentina_20191028' else CASE);c['interval_minutes']=1
    if event=='argentina_20181213':
        old=(PREVIOUS if PREVIOUS.exists() else OUT)/'data/members'/cfg/f'{member:02}'
        record=json.loads((old/'minute_identity.json').read_text())
        if record['source_sha256']!=run['sha256']:raise ValueError('Ice integrals do not match native output')
        ice=pd.read_csv(old/'predictors_1min.csv');five=pd.read_csv(old/'predictors_5min.csv')
    elif cfg=='TH_C0':
        old=ROOT/'results/argentina_20191028/diagnostics_1min'
        if json.loads((old/'manifest.json').read_text())['simulation_sha256']!=run['sha256']:raise ValueError('2019 ice source changed')
        ice=pd.read_csv(old/'predictors_1min.csv');five=pd.read_csv(old/'predictors_5min.csv')
    else:
        old=destination/'native_ice_integrals'
        ice,five=control_predictors(source,ROOT/c['paths']['forcing']/'LEONA_SCM.nc',c,old,
            configuration=cfg,export_native=True,compute_cape=False)
    with xr.open_dataset(source,decode_times=False) as original:
        d=original[['time','lev','pressure','theta','qv','rainnc']].load()
    duration=(pd.Timestamp(c['end'])-pd.Timestamp(c['start'])).total_seconds()
    if not (d.time.values[0]==0 and d.time.values[-1]==duration and np.all(np.diff(d.time.values)==30)):
        raise ValueError('Complete native 30s profile sequence required')
    p=d.pressure.values;t=d.theta.values*(p/1e5)**(2/7);rain=d.rainnc.values
    if np.any(np.diff(rain)<0):raise ValueError('Reset accumulated precipitation')
    compare=event=='argentina_20181213' and cfg in ['TH_C0','WS_C0'] and member==0
    rows=[]
    for i in range(0,len(d.time)-1,2):
        row=diagnose(p[i],t[i],d.qv.values[i],d.lev.values,comparison=compare)
        amount=float(rain[i+2]-rain[i]);row.update(time=str(pd.Timestamp(c['start'])+pd.Timedelta(seconds=float(d.time.values[i]))),
            configuration=cfg,member=member,precipitation_amount_kg_m2=amount,CAPExP_J_m2=row['cape_J_kg']*amount)
        if compare:
            for name in ['SBCAPE','previous_most_cape']:
                row[name+'_P_J_m2']=row.get(name+'_J_kg',np.nan)*amount
        rows.append(row)
    therm=pd.DataFrame(rows);therm.to_csv(destination/'thermodynamics_mucape_1min.csv',index=False)
    cap=therm[['time','CAPExP_J_m2','configuration','member']].rename(columns={'CAPExP_J_m2':'predictor'})
    cap['method']=MU_METHOD;cap['status']=therm.status.values
    ice=ice[ice.method!='CAPExP'];result=pd.concat([ice,cap],ignore_index=True)
    result.to_csv(destination/'predictors_1min.csv',index=False);five.to_csv(destination/'predictors_5min.csv',index=False)
    if sha256(source)!=run['sha256']:raise IOError('Native simulation changed during diagnosis')
    forcing=ROOT/(c['paths']['forcing'] if event=='argentina_20191028' else 'data/processed/forcing')/'LEONA_SCM.nc'
    write_json(tag,dict(run,source_sha256=run['sha256'],algorithm_sha256=algorithm,definition=DEFINITION,
        search_depth_Pa=DEPTH_PA,forcing_sha256=sha256(forcing),comparison=compare,native_output_unchanged=True,
        ice_integral_source=str(old),ice_predictor_sha256=sha256(old/'predictors_1min.csv'),
        statuses=therm.status.value_counts().to_dict(),time_step_seconds=60,native_snapshot_seconds=30))


def prepare(workers=12):
    runs=plan()
    environment={name:importlib.metadata.version(name) for name in ['MetPy','numpy','pandas','xarray','scipy']}
    identity=diagnostics_root('argentina_20181213')/'preparation_versions.json'
    if identity.exists() and json.loads(identity.read_text())!=environment:raise ValueError('Changed diagnostic environment; use a fresh version')
    write_json(identity,environment)
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(diagnose_run,r) for r in runs]
        for n,future in enumerate(as_completed(futures),1):
            future.result();write_json(diagnostics_root('argentina_20181213')/'progress.json',{'completed':n,'expected':len(runs),'SCM_2018_reruns':0})
            if n%6==0 or n==len(runs):print(f'Native MUCAPE diagnostics: {n}/{len(runs)}',flush=True)


def calibrate():
    from scm_lightning.calibration.revision import fit_controls,freeze,load_frozen
    if (FROZEN/'frozen.json').exists():load_frozen(FROZEN);return
    c=load_config(CAL);c['interval_minutes']=1
    roots=[diagnostics_root('argentina_20191028')/cfg/'00' for cfg in ['TH_C0','WS_C0']]
    for folder in roots:
        identity=json.loads((folder/'identity.json').read_text())
        if identity['event']!='argentina_20191028' or identity['member']!=0:raise ValueError('Invalid calibration diagnostics')
        if sha256(Path(identity['folder'])/'output.nc')!=identity['source_sha256']:raise ValueError('Changed calibration simulation')
    p=pd.concat([pd.read_csv(r/'predictors_1min.csv') for r in roots],ignore_index=True)
    p5=pd.concat([pd.read_csv(r/'predictors_5min.csv') for r in roots],ignore_index=True)
    observations=ROOT/'results/argentina_20191028/observations'
    coef,baseline=fit_controls(c,p,pd.read_csv(observations/'counts_1min.csv'),p5,pd.read_csv(observations/'counts_5min.csv'))
    sources=[CAL,ROOT/'configs/evaluation/validation_v3_mucape.yaml',observations/'counts_1min.csv',observations/'counts_5min.csv',
        observations/'flashes.csv.gz',ROOT/c['paths']['forcing']/'LEONA_SCM.nc']
    sources+=sorted((ROOT/c['paths']['era5']).glob('*.nc'))
    for root in roots:
        identity=json.loads((root/'identity.json').read_text());native=Path(identity['folder'])
        sources.extend([root/'identity.json',root/'predictors_1min.csv',root/'predictors_5min.csv',root/'thermodynamics_mucape_1min.csv',
            native/'output.nc',native/'run_identity.json',native/'namelist.mp_offline'])
    freeze(FROZEN,c,coef,baseline,sources)
    print('Frozen six 2019-only coefficients.',flush=True)


def prepare_calibration():
    """Prioritize the two controls; no 2019 ensemble is constructed."""
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
    runs=[r for r in plan() if r['event']=='argentina_20191028']
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(diagnose_run,r) for r in runs]):future.result()
    print('Both 2019 MUCAPE controls prepared.',flush=True)


def preserve_previous():
    if PREVIOUS.exists():return
    allowed=(ROOT/'results/argentina_20181213').resolve()
    if not OUT.resolve().is_relative_to(allowed) or not PREVIOUS.resolve().is_relative_to(allowed):raise ValueError('Archive outside event')
    records=[{'path':str(f.relative_to(OUT)),'sha256':sha256(f)} for f in OUT.rglob('*') if f.is_file()]
    shutil.copytree(OUT,PREVIOUS)
    for record in records:
        if sha256(PREVIOUS/record['path'])!=record['sha256']:raise IOError('Previous methodological archive failed')
    write_json(PREVIOUS/'v2_archive_identity.json',{'files':records,'version':'v2_1min_most_cape_TH_common_peak_mapping',
        'native_outputs_modified':False})


def compare_controls():
    """Profile changes and parcel/LFC/EL switches; never read lightning here."""
    from metpy.calc import equivalent_potential_temperature,dewpoint,vapor_pressure
    from metpy.units import units
    from metpy.calc import parcel_profile_with_lcl
    from scm_lightning.lightning.parcel_levels import integration_limits
    folder=ROOT/'results/argentina_20181213/thermodynamic_revision';folder.mkdir(parents=True,exist_ok=True)
    cache=folder/'comparison_contract.json'
    helper_hash=sha256(ROOT/'python/scm_lightning/lightning/parcel_levels.py')
    implementation_hash=sha256(__file__)
    if cache.exists():
        previous=json.loads(cache.read_text())
        if previous.get('status')=='passed' and previous.get('comparison_code_sha256')==implementation_hash and previous.get('integration_limit_code_sha256')==helper_hash:
            if all(sha256(r['thermodynamic_source'])==r['thermodynamic_source_sha256'] for r in previous['controls']):
                print('Validated control CAPE comparison already current.',flush=True);return
    summaries=[]
    for cfg in ['TH_C0','WS_C0']:
        f=pd.read_csv(diagnostics_root('argentina_20181213')/cfg/'00/thermodynamics_mucape_1min.csv')
        if 'previous_most_cape_LFC_Pa' not in f:raise ValueError('Required control definition comparison is incomplete')
        source=Path(json.loads((diagnostics_root('argentina_20181213')/cfg/'00/identity.json').read_text())['folder'])/'output.nc'
        with xr.open_dataset(source,decode_times=False) as d:
            p=d.pressure.values[::2][:-1];t=d.theta.values[::2][:-1]*(p/1e5)**(2/7);rv=d.qv.values[::2][:-1];z=d.lev.values
        if len(f)!=len(p):raise ValueError('Unaligned control thermodynamic clock')
        # Raw diagnostic cache also retains uncorrected-temperature roots.
        # Reconstruct the actual virtual-temperature roots used inside cape_cin.
        for label in ['SBCAPE','previous_most_cape','MUCAPE']:
            for kind in ['LFC','EL']:
                f[label+'_'+kind+'_uncorrected_temperature_Pa']=f[label+'_'+kind+'_Pa']
        reconstructed=[]
        for i in range(len(f)):
            record={}
            if f.parcel_index.notna().iloc[i]:
                pressure=p[i]*units.Pa;temperature=t[i]*units.K
                dewpoints=dewpoint(vapor_pressure(pressure,rv[i]*units.dimensionless))
                prepared=parcel_profile_with_lcl(pressure,temperature,dewpoints)
                mu_index=int(f.parcel_index.iloc[i])
                for label,start,selection in [('SBCAPE',0,'bottom'),('previous_most_cape',0,'most_cape'),('MUCAPE',mu_index,'bottom')]:
                    values=integration_limits(pressure[start:],temperature[start:],dewpoints[start:],
                        which_lfc=selection,which_el='most_cape' if selection=='most_cape' else 'top',prepared=prepared if start==0 else None)
                    record.update({label+'_'+key:value for key,value in values.items()})
            reconstructed.append(record)
        actual=pd.DataFrame(reconstructed)
        for name in actual:f[name]=actual[name]
        for label,column in [('SBCAPE','SBCAPE_J_kg'),('previous_most_cape','previous_most_cape_J_kg'),('MUCAPE','raw_cape_J_kg')]:
            valid=np.isfinite(f[column])&np.isfinite(f[label+'_signed_integral_J_kg'])
            np.testing.assert_allclose(f.loc[valid,label+'_signed_integral_J_kg'],f.loc[valid,column],rtol=1e-10,atol=1e-7)
        f['MUCAPE_top_limited']=(f.raw_cape_J_kg>0)&f.MUCAPE_EL_Pa.isna()
        f['max_abs_delta_T_K']=np.r_[np.nan,np.max(abs(np.diff(t,axis=0)),axis=1)]
        f['max_abs_delta_rv_kg_kg']=np.r_[np.nan,np.max(abs(np.diff(rv,axis=0)),axis=1)]
        with np.errstate(divide='ignore',invalid='ignore'):
            td=dewpoint(vapor_pressure(p*units.Pa,rv*units.dimensionless))
            theta_e=equivalent_potential_temperature(p*units.Pa,t*units.K,td).to('K').magnitude
        nearest_top=np.argmin(abs(p-(p[:,0]-30000)[:,None]),axis=1)
        in_layer=np.arange(p.shape[1])[None,:]<=nearest_top[:,None]
        ranked=np.sort(np.where(in_layer & np.isfinite(theta_e),theta_e,-np.inf),axis=1)
        f['theta_e_leading_candidate_gap_K']=ranked[:,-1]-ranked[:,-2]
        valid=f.parcel_index.notna().to_numpy();idx=f.parcel_index.fillna(0).astype(int).to_numpy()
        f['selected_theta_e_K']=np.where(valid,theta_e[np.arange(len(f)),idx],np.nan)
        f['parcel_at_nearest_native_search_top']=valid & (idx==nearest_top)
        f['parcel_index_changed']=f.parcel_index.diff().abs().gt(0)
        small=(f.max_abs_delta_T_K<.2)&(f.max_abs_delta_rv_kg_kg<.0002)
        f['small_profile_change']=small
        for label,column in [('SBCAPE','SBCAPE_J_kg'),('previous_most_cape','previous_most_cape_J_kg'),('MUCAPE','cape_J_kg')]:
            change=f[column].diff();f[label+'_delta_J_kg']=change
            jump=(change.abs()>500)|((change.abs()>100)&(change.abs()/f[column].shift().abs().clip(lower=100)>.5))
            level_switch=f[label+'_LFC_Pa'].diff().abs().gt(5000)|f[label+'_EL_Pa'].diff().abs().gt(5000)
            f[label+'_large_change_flag']=jump;f[label+'_LFC_or_EL_switch_flag']=level_switch
            f[label+'_selection_change_small_profile_delta']=jump & small & (level_switch|f.parcel_index_changed)
            f[label+'_large_change_with_larger_profile_delta']=jump & ~small
        f.to_csv(folder/f'CAPE_profile_transition_diagnostics_{cfg}.csv.gz',index=False)
        clock=pd.to_datetime(f.time,utc=True)
        window=(clock>=pd.Timestamp('2018-12-13T20:00Z'))&(clock<pd.Timestamp('2018-12-15T00:00Z'))
        f[window].to_csv(folder/f'CAPE_21UTC_and_Dec14_transition_diagnostics_{cfg}.csv',index=False)
        thermodynamic_source=diagnostics_root('argentina_20181213')/cfg/'00/thermodynamics_mucape_1min.csv'
        summaries.append({'configuration':cfg,'native_sha256':sha256(source),'intervals':len(f),'statuses':f.status.value_counts().to_dict(),
            'thermodynamic_source':str(thermodynamic_source),'thermodynamic_source_sha256':sha256(thermodynamic_source),
            'positive_top_limited_MUCAPE_intervals':int(f.get('MUCAPE_top_limited',pd.Series(False,index=f.index)).fillna(False).sum()),
            'selection_change_small_profile_delta_counts':{name:int(f[name+'_selection_change_small_profile_delta'].sum()) for name in ['SBCAPE','previous_most_cape','MUCAPE']}})
    write_json(folder/'comparison_contract.json',{'status':'passed','controls':summaries,'diagnostic_only':True,'lightning_observations_read':False,
        'primary_definition':'official MetPy MUCAPE, lowest 300hPa max theta-e; chosen independently of December lightning',
        'search_depth_preflight':'100/200/300hPa searches over predeclared profile snapshots; 300hPa spans approximately lower 3km and admits elevated inflow',
        'MetPy_native_boundary':'get_layer(interpolate=False) uses nearest native layer-top level; selected pressure/height retained explicitly',
        'SBCAPE_surface':'lowest simulated cell (25m AGL); official surface_based_cape_cin; no extrapolated ground sample',
        'transition_flags_not_causal_proof':True,'flag_thresholds':{'large_absolute_CAPE_change_J_kg':500,'relative_change_fraction':.5,
            'minimum_absolute_change_for_relative_flag_J_kg':100,'relative_denominator_reference_J_kg':100,
            'small_max_delta_T_K':.2,'small_max_delta_rv_kg_kg':.0002,'LFC_EL_switch_Pa':5000},
        'incomplete_CIN_Cape_policy':'invalid/dry/missing/insufficient-depth profiles remain NaN with status; valid physical zero retained',
        'upper_boundary_policy':'official MetPy integrates to EL or native top, no extrapolation; positive open-EL controls flagged top-limited',
        'LFC_EL':'actual virtual-temperature integration limits reconstructed from installed MetPy; uncorrected-temperature roots retained separately',
        'signed_area_audit':'positive plus negative buoyancy contributions reconstruct provider CAPE; neither is substituted as a predictor',
        'integration_limit_code_sha256':helper_hash,'comparison_code_sha256':implementation_hash})
    print('Saved control CAPE/CIN, parcel/LFC/EL and profile-transition diagnostics.',flush=True)


def contrasts(metric,comparison):
    from scm_lightning.calibration.revision import METHODS
    rows=[]
    for method in METHODS:
        f=metric[metric.method==method].set_index('configuration')
        pairs=[('WS_'+mode,'TH_'+mode) for mode in ['C0','CL','CF']]+[(s+'_'+mode,s+'_C0') for s in ['TH','WS'] for mode in ['CL','CF']]
        for a,b in pairs:
            for key in ['rmse','mae','skill_zero','skill_constant','bias','accumulated_error','pearson','spearman','relative_peak_error','peak_time_error_minutes']:
                rows.append({'comparison':comparison,'method':method,'configuration_a':a,'configuration_b':b,
                    'contrast':a+' minus '+b,'metric':key,'difference':f.loc[a,key]-f.loc[b,key]})
    return pd.DataFrame(rows)


def evaluate():
    from scm_lightning.calibration.revision import METHODS,normalize,apply,load_frozen
    from scm_lightning.calibration.independent import load_frozen as load_old
    from scm_lightning.evaluation.independent import evaluate_predictions
    from scm_lightning.forcing.contracts import area_km2
    contract=ROOT/'results/argentina_20181213/thermodynamic_revision/comparison_contract.json'
    if not contract.exists() or json.loads(contract.read_text()).get('status')!='passed':raise ValueError('Validate mandatory control CAPE comparison before promoting v3')
    coef,meta=load_frozen(FROZEN);c=load_config(CASE)
    paths=sorted(diagnostics_root('argentina_20181213').glob('*/*/predictors_1min.csv'))
    if len(paths)!=246:raise ValueError('Full native MUCAPE ensemble required before replacing results')
    preserve_previous()
    raw=pd.concat([pd.read_csv(p) for p in paths],ignore_index=True)
    p=normalize(raw,load_config(CAL));obs=pd.read_csv(ROOT/'results/argentina_20181213/observations/counts_1min.csv')
    area=area_km2(c['area']);baseline=meta['constant_density_flashes_km2_per_interval']*area
    previous_mask=None;contrast_tables=[]
    for comparison,base in [('microphysics_specific',OUT),('common_TH_same_MUCAPE_OLS',OUT/'supplementary/common_TH_same_MUCAPE_OLS')]:
        prediction=apply(p,coef,area,comparison)
        metric,statistics,mask=evaluate_predictions(prediction,obs,baseline,minutes=1,methods=METHODS)
        if previous_mask is not None:pd.testing.assert_frame_equal(mask,previous_mask)
        previous_mask=mask
        for folder in ['data','tables']:(base/folder).mkdir(parents=True,exist_ok=True)
        prediction.to_csv(base/'data/independent_member_predictions.csv.gz',index=False)
        statistics.to_csv(base/'data/ensemble_statistics.csv',index=False)
        obs.to_csv(base/'data/observations_1min.csv',index=False);mask.to_csv(base/'data/common_valid_intervals.csv',index=False)
        metric['comparison']=comparison;metric.to_csv(base/'tables/Table_04_independent_evaluation.csv',index=False)
        metric[metric.configuration.isin(['TH_C0','WS_C0'])].to_csv(base/'tables/Table_04_reference_C0.csv',index=False)
        contrast_tables.append(contrasts(metric,comparison))
        cube=prediction[['time','method','configuration','member','predictor','prediction_density','prediction']].copy()
        cube.time=pd.to_datetime(cube.time,utc=True).dt.tz_localize(None)
        ds=cube.set_index(['time','method','configuration','member']).to_xarray()
        status=prediction[prediction.method==MU_METHOD][['time','configuration','member','status']].copy()
        status['time']=pd.to_datetime(status.time,utc=True).dt.tz_localize(None)
        codes={'valid_positive':0,'valid_zero':1,'invalid_profile':2,'incomplete_profile':3,'undefined_dewpoint':4,'invalid_calculation':5}
        status['MUCAPE_status']=status.status.map(codes).astype('i1')
        ds['MUCAPE_status']=status.set_index(['time','configuration','member']).MUCAPE_status.to_xarray()
        ds.MUCAPE_status.attrs.update(flag_values=np.arange(6,dtype='i1'),flag_meanings=' '.join(codes))
        factors=prediction.groupby(['method','configuration']).scale.first().unstack('configuration')
        ds['calibration_coefficient']=(('method','configuration'),factors.reindex(index=ds.method.values,columns=ds.configuration.values).values)
        ds['observation_count']=('time',obs['count'].values)
        ds['common_valid']=('time',mask.common_valid.to_numpy(dtype='i1'))
        ds['time_bounds']=(('time','bounds'),np.column_stack([ds.time.values,ds.time.values+np.timedelta64(60,'s')]))
        ds.prediction.attrs.update(units='min-1',long_name='IC+CG flashes in the one-minute interval / one minute')
        ds.prediction_density.attrs.update(units='km-2 min-1')
        ds.predictor.attrs['units']='method-dependent normalized interval integral; see Table_02_lightning_parameterizations.csv'
        ds.calibration_coefficient.attrs['units']='method-dependent; see Table_03_calibration_coefficients.csv'
        ds.observation_count.attrs.update(units='1',long_name='IC+CG flashes observed during one minute')
        ds.time.attrs.update(bounds='time_bounds',long_name='UTC interval start')
        ds.attrs.update(interval_seconds=60,native_snapshot_seconds=30,calibration_event='argentina_20191028',evaluation_event='argentina_20181213',
            comparison=comparison,coefficient_sha256=meta['coefficient_sha256'],CAPE='MetPy MUCAPE, lowest 300 hPa, not SBCAPE or Romps isentropic CAPE',
            McCaul_fit='OLS all valid five-minute intervals; native 60s prediction integrals /300s /nominal published analysis-cell area',
            evaluation_area_wsen=json.dumps(c['area']),evaluation_area_km2=area,calibration_area_wsen=json.dumps(meta['area_wsen']),
            calibration_area_km2=meta['area_km2'],MetPy_version=meta['versions']['MetPy'],McCaul_velocity='prescribed LS only; EDMF never added',
            ensemble='members 1..40 separate by microphysics and coupling; member0 control excluded from mean/std')
        for name in ['method','configuration']:ds[name]=xr.Variable((name,),np.asarray(ds[name].values,dtype='S40'))
        for name in ['time','time_bounds']:ds[name].encoding.update(units='minutes since 2018-12-13 00:00:00',calendar='proleptic_gregorian')
        target=base/'data/independent_predictions_1min.nc';temp=target.with_suffix('.pending.nc');ds.to_netcdf(temp);temp.replace(target)
        write_json(base/'evaluation_manifest.json',{'version':3,'comparison':comparison,'interval_minutes':1,
            'coefficient_sha256':meta['coefficient_sha256'],'baseline_flashes_per_min':baseline,
            'common_valid_intervals':int(mask.common_valid.sum()),'evaluation_observations_sha256':sha256(ROOT/'results/argentina_20181213/observations/counts_1min.csv'),
            '2018_observations_used_for_calibration':False})
    pd.concat(contrast_tables,ignore_index=True).to_csv(OUT/'tables/Table_05_microphysics_sensitivity.csv',index=False)
    coef.to_csv(OUT/'tables/Table_03_calibration_coefficients.csv',index=False)
    pd.DataFrame([
        {'event':'argentina_20191028','role':'calibration only','physical_configurations':'TH_C0; WS_C0',
            'members':'two deterministic controls, member0; no calibration ensemble','start':'2019-10-28T00:00Z','end_exclusive':'2019-10-30T00:00Z'},
        {'event':'argentina_20181213','role':'independent evaluation only','physical_configurations':'TH/WS x C0/CL/CF',
            'members':'40 perturbed + separate control per configuration','start':c['start'],'end_exclusive':c['end']}]).to_csv(OUT/'tables/Table_01_experimental_design.csv',index=False)
    # Keep the original Thompson-common v2 assessment with its own CAPE and peak mapping.
    historical=OUT/'supplementary/legacy_original_TH_v2';historical.mkdir(parents=True,exist_ok=True)
    for filename in ['Table_04_independent_evaluation.csv','Table_05_microphysics_sensitivity.csv']:
        shutil.copy2(PREVIOUS/'tables'/filename,historical/filename)
    for path in (PREVIOUS/'figures').glob('Figure_0[34]*'):shutil.copy2(path,historical/path.name)
    write_json(historical/'comparison_manifest.json',{'version':'v2_original_TH_common','calibration':'results/argentina_20191028/calibration/v2_1min',
        'CAPE':'previous lowest-cell most_cape LFC/EL','McCaul':'original frozen 5min peak ratio',
        'not_pooled_with_v3_metrics':True,'reason':'preserve exact original common-Thompson calibration comparison'})
    definitions=pd.DataFrame([
        {'method':MU_METHOD,'equation':'MUCAPE * exact next-60s rainnc increment','definition':'official MetPy most_unstable_cape_cin; maximum theta-e in lowest 300hPa',
            'adaptation':'SCM MUCAPE; liquid pseudoadiabatic; not SBCAPE; not exact Romps isentropic CAPE','predictor_units':'J km-2 per 1min'},
        {'method':'IxG','equation':'integral (rho_d*qi)*(rho_d*qg) dz dt','definition':'cloud ice only; original hydrometeors unchanged',
            'adaptation':'whole-column SCM density mapping','predictor_units':'kg2 s m-5 per 1min'},
        {'method':METHODS[2],'equation':'.95*.042*(w_LS_positive*qg[g/kg])[-15C]+.05*.20*VII','definition':'published F1 F2 F3 constants unchanged',
            'adaptation':'LS SCM score; native integral /300s /published nominal 0.009 degree analysis cell at 34.7N; additional 2019 OLS factor',
            'predictor_units':'nominal-cell flash equivalent km-2 per interval'}])
    definitions.to_csv(OUT/'tables/Table_02_lightning_parameterizations.csv',index=False)
    write_json(OUT/'tables/metric_conventions.json',{'interval_minutes':1,'rates':'flashes/min from actual one-minute integrals',
        'comparison':'microphysics_specific main; common_TH_same_MUCAPE_OLS sensitivity; legacy original TH v2 separate',
        'support':'identical common finite primary predictors across all methods/configurations/members and both v3 coefficient comparisons',
        'RMSE_MAE_bias':'flashes/min; prediction minus observation','skills':'1-MSE/MSE_reference; zero and 2019 observed mean density times 2018 area',
        'correlations':'Pearson and Spearman; undefined constant series NaN','peaks':'first maximum; signed UTC timing error minutes; constant series NaN'})
    write_json(OUT/'independent_manifest.json',{'version':3,'status':'evaluated','interval_minutes':1,
        'calibration_directory':str(FROZEN),'coefficient_count':6,'coefficient_sha256':meta['coefficient_sha256'],
        'main_comparison':'microphysics_specific','sensitivity':'common_TH_same_MUCAPE_OLS',
        'legacy_comparison':'supplementary/legacy_original_TH_v2','prior_article_archive':str(PREVIOUS),
        'physical_2018_runs_repeated':0,'physical_2019_runs_added':['WS_C0'],'new_ensemble_members':0})
    print('Evaluated separate six-coefficient and common-Thompson v3 comparisons.',flush=True)


def check():
    from scm_lightning.calibration.revision import METHODS,apply,load_frozen
    from scm_lightning.evaluation.independent import evaluate_predictions
    from scm_lightning.forcing.contracts import area_km2
    from scm_lightning.provenance import git
    coef,meta=load_frozen(FROZEN)
    if not meta['all_identifiable']:raise ValueError('Unidentifiable coefficients remain; do not claim complete calibration')
    for source in meta['sources']:
        if '20181213' in source['path'] or sha256(source['path'])!=source['sha256']:raise ValueError('Frozen calibration provenance changed')
    area=area_km2(load_config(CASE)['area']);baseline=meta['constant_density_flashes_km2_per_interval']*area
    masks=[]
    for base,comparison in [(OUT,'microphysics_specific'),(OUT/'supplementary/common_TH_same_MUCAPE_OLS','common_TH_same_MUCAPE_OLS')]:
        p=pd.read_csv(base/'data/independent_member_predictions.csv.gz');o=pd.read_csv(base/'data/observations_1min.csv')
        actual,stats,mask=evaluate_predictions(p,o,baseline,minutes=1,methods=METHODS)
        actual['comparison']=comparison
        pd.testing.assert_frame_equal(actual,pd.read_csv(base/'tables/Table_04_independent_evaluation.csv'),check_dtype=False,rtol=1e-12,atol=1e-12)
        recalculated=apply(p.drop(columns=['scale','coefficient_status','prediction','prediction_density','comparison']),coef,area,comparison)
        np.testing.assert_allclose(recalculated.prediction,p.prediction,rtol=1e-12,equal_nan=True)
        if comparison=='microphysics_specific' and not p.groupby(['method','microphysics']).scale.nunique().eq(1).all():raise ValueError('Member/coupling-specific coefficients detected')
        if comparison!='microphysics_specific' and not p.groupby('method').scale.nunique().eq(1).all():raise ValueError('Shared-physics comparison has separate coefficients')
        masks.append(mask)
    pd.testing.assert_frame_equal(*masks)
    native=[]
    for run in plan():
        identity=json.loads((diagnostics_root(run['event'])/run['configuration']/f'{run["member"]:02}'/'identity.json').read_text())
        if identity['source_sha256']!=run['sha256'] or sha256(Path(run['folder'])/'output.nc')!=run['sha256']:raise ValueError('Native physical output changed')
        if not identity['native_output_unchanged']:raise ValueError('Diagnostic output identity failed')
        native.append({'event':run['event'],'configuration':run['configuration'],'member':run['member'],'sha256':run['sha256']})
    if git('diff','--','src'):raise ValueError('Fortran working tree changed')
    for ext in ['png','pdf']:
        name='figures/Figure_05_microphysical_response.'+ext
        if sha256(OUT/name)!=sha256(PREVIOUS/name):raise ValueError('Original microphysical Figure 5 changed')
    meteorology=[]
    for event in ['argentina_20181213','argentina_20191028']:
        folder=ROOT/'results'/event/'meteorological_validation';m=json.loads((folder/'MERGE_validation_manifest.json').read_text())
        if m['area_wsen']!=load_config(ROOT/'configs/cases'/(event+'.yaml'))['area'] or m['common_valid_intervals']!=48:raise ValueError('Incompatible precipitation support')
        for source in m['MERGE_sources']+m['ERA5_sources']+m['SCM_sources']:
            if sha256(source['path'])!=source['sha256']:raise ValueError('Meteorological source changed')
        meteorology.append({'event':event,'manifest_sha256':sha256(folder/'MERGE_validation_manifest.json')})
    write_json(OUT/'audit/v3_validation.json',{'status':'passed','version':3,'six_coefficients_2019_only':True,
        'comparison_modes':['microphysics_specific','common_TH_same_MUCAPE_OLS'],'coefficient_sha256':meta['coefficient_sha256'],
        'native_physical_outputs_checked':native,'original_2018_outputs_unchanged':246,'2019_new_controls':1,
        'Fortran_unchanged':True,'Figure_05_unchanged':True,'LS_EDMF_summed':False,'native_seconds':30,'evaluation_seconds':60,
        'common_valid_intervals':int(masks[0].common_valid.sum()),'meteorology':meteorology,
        'versions':meta['versions'],'MERGE_algorithm_version_limit':'historical archive GRIB does not identify algorithm release; collection v1, actual files and hashes recorded',
        'code':[{'path':str(p.relative_to(ROOT)),'sha256':sha256(p)} for p in [Path(__file__),ROOT/'python/scm_lightning/lightning/thermodynamics.py',
            ROOT/'python/scm_lightning/calibration/revision.py',ROOT/'python/scm_lightning/evaluation/merge.py',ROOT/'python/scm_lightning/plotting/revision.py',
            ROOT/'python/scm_lightning/lightning/parcel_levels.py']]})
    write_json(OUT/'pending_dependencies.json',{'blocking':[],
        'MERGE_historical_algorithm_version':'not encoded in legacy GRIB; STAC collection v1, descriptors, source hashes and actual acquisition preserved',
        'CAPE_upper_boundary':'MetPy uses EL or native profile top without extrapolation; open-EL control samples flagged',
        'initial_boundary_sounding_substitution':'ERA5 retained; distant/post-initialization soundings remain environmental observations'})
    print('v3 audit passed: six coefficients, two independent comparisons, 248 unchanged native files, MERGE support valid.',flush=True)


def run_control():
    """Only the missing WSM6 control; exact Thompson input/numerics/executable."""
    from scm_lightning.ensemble.experiments import execute,validate_output,model_path
    reference=ROOT/'runs/argentina_20191028/thompson_C0/control'
    folder=ROOT/'runs/argentina_20191028/wsm6_C0/control'
    c=load_config(CAL);forcing=ROOT/c['paths']['forcing']/'LEONA_SCM.nc'
    identity=json.loads((reference/'run_identity.json').read_text())
    executable=ROOT/'build/updraft_sensitivity/mp_offline_driver'
    if sha256(executable)!=identity['executable_sha256'] or sha256(forcing)!=identity['forcing_sha256']:
        raise ValueError('2019 reference executable/forcing identity changed')
    if sha256(reference/'output.nc')!=identity['result']['sha256']:raise ValueError('Changed Thompson reference')
    duration=(pd.Timestamp(c['end'])-pd.Timestamp(c['start'])).total_seconds()
    output=folder/'output.nc';folder.mkdir(parents=True,exist_ok=True)
    nml=(reference/'namelist.mp_offline').read_text()
    expected=nml.replace('microp_scheme="thompson"','microp_scheme="wsm6"').replace(model_path(reference/'output.nc'),model_path(output))
    if expected==nml or 'microp_scheme="wsm6"' not in expected:raise ValueError('Cannot construct matched control')
    if output.exists():
        if (folder/'namelist.mp_offline').read_text()!=expected:raise ValueError('Existing WSM6 input differs')
        result=validate_output(output,duration)
    else:
        (folder/'namelist.mp_offline').write_text(expected,encoding='utf-8')
        execute(executable,folder);result=validate_output(output,duration)
    if not result.get('physical_preflight_passed'):raise ValueError('WSM6 physical preflight failed')
    write_json(folder/'run_identity.json',{'configuration':'WS_C0','member':0,'forcing_sha256':sha256(forcing),
        'executable_sha256':sha256(executable),'reference_TH_C0_sha256':identity['result']['sha256'],
        'only_namelist_changes':['microp_scheme','path_output'],'result':result})
    print('Validated 2019 WSM6 C0 control.',flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['run-control','prepare','prepare-calibration','compare-controls','calibrate','evaluate','check']);p.add_argument('--workers',type=int,default=12);a=p.parse_args()
    if a.stage=='run-control':run_control()
    elif a.stage=='prepare':prepare(a.workers)
    elif a.stage=='calibrate':calibrate()
    elif a.stage=='evaluate':evaluate()
    elif a.stage=='compare-controls':compare_controls()
    elif a.stage=='check':check()
    elif a.stage=='prepare-calibration':prepare_calibration()


if __name__=='__main__':main()
