"""Independent-event products. Never fit with December 2018 observations."""
import argparse
import json
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import spearmanr
from scm_lightning.config import load_config
from scm_lightning.provenance import ROOT, sha256, write_json, resolve_path, git
from scm_lightning.forcing.contracts import area_km2
from scm_lightning.lightning.adaptations import surface_cape, crossing, at_height
from scm_lightning.lightning.postprocess import integrate_bins
from scm_lightning.calibration.independent import METHODS, normalize_predictors, fit, freeze, load_frozen, apply

CALIBRATION_CONFIG=ROOT/'configs/cases/argentina_20191028.yaml'
EVALUATION_CONFIG=ROOT/'configs/cases/argentina_20181213.yaml'
OUT=ROOT/'results/argentina_20181213/article'
FROZEN=ROOT/'results/argentina_20191028/calibration/v1'
LEGACY=ROOT/'results/argentina_20181213/legacy/article_within_event_20261106'


def preserve_legacy():
    """Atomic directory rename, inside the verified results root; no native run moves."""
    if LEGACY.exists():return
    source=OUT.resolve();target=LEGACY.resolve();allowed=(ROOT/'results/argentina_20181213').resolve()
    if not source.is_relative_to(allowed) or not target.is_relative_to(allowed):raise ValueError('Archive outside case root')
    if source.exists():
        files=[{'path':str(p.relative_to(source)),'sha256':sha256(p)} for p in source.rglob('*') if p.is_file()]
        target.parent.mkdir(parents=True,exist_ok=True)
        # Windows viewers may hold a directory handle and prevent its rename.
        # Copy and hash-verify before removing only superseded derived products.
        shutil.copytree(source,target)
        for record in files:
            if sha256(target/record['path'])!=record['sha256']:raise IOError('Legacy archive verification failed')
        for record in files:
            p=source/record['path']
            if not p.resolve().is_relative_to(allowed):raise ValueError('Archive cleanup escaped case root')
            if not record['path'].replace('\\','/').startswith('data/members/'):
                p.unlink()
        write_json(target/'preservation.json',{'files':files,'reason':'superseded within-event calibration, not independent validation',
            'source_directory':str(source),'git_commit':git('rev-parse','HEAD'),'native_simulations_moved':False})


def control_predictors(simulation, forcing, config, destination):
    """Read native snapshots; export pre-calibration fields, components and units."""
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    before=sha256(simulation)
    with xr.open_dataset(simulation,decode_times=False) as original:
        d=original[['theta','pressure','qv','qi','qs','qg','rainnc','time','lev']].load()
    z=d.lev.values;seconds=d.time.values
    duration=(pd.Timestamp(config['end'])-pd.Timestamp(config['start'])).total_seconds()
    if seconds[0]!=0 or seconds[-1]!=duration or not np.all(np.diff(seconds)==30):
        raise ValueError('Complete native 30-second case required')
    edges=np.zeros(len(z)+1)
    for k,center in enumerate(z):edges[k+1]=2*center-edges[k]
    dz=np.diff(edges)
    if (dz<=0).any():raise ValueError('Invalid model layer thickness')
    t=d.theta.values*(d.pressure.values/1.e5)**(2/7)
    rho=d.pressure.values/(287*t*(1+d.qv.values/(287/461.6)))
    with xr.open_dataset(forcing,decode_times=False) as f:
        if f.time.values[0]>0 or f.time.values[-1]<duration:raise ValueError('Incomplete temporal forcing')
        height=np.array([np.interp(z,f.lev.values,row) for row in f.wa.values])
        w=np.column_stack([np.interp(seconds,f.time.values,height[:,k]) for k in range(len(z))])
    h=crossing(t,z,258.15)
    f1=.042*np.maximum(at_height(w,z,h),0)*at_height(1000*d.qg.values,z,h)
    vii=np.sum(rho*(d.qi.values+d.qs.values+d.qg.values)*dz,axis=1)
    f2=.20*vii;f3=.95*f1+.05*f2
    ixg=np.sum(rho**2*d.qi.values*d.qg.values*dz,axis=1)
    clock=pd.Timestamp(config['start'])+pd.to_timedelta(seconds,unit='s')
    raw=pd.concat([pd.DataFrame({'time':clock,'method':name,'value':v}) for name,v in
        [('IxG',ixg),('McCaulF3',f3),('McCaulF1',f1),('McCaulF2',f2)]],ignore_index=True)
    five=integrate_bins(raw,config['start'],config['end'],5)
    bins=integrate_bins(raw,config['start'],config['end'],15)
    cape=[];rain=d.rainnc.values
    if (np.diff(rain)<0).any():raise ValueError('Reset/negative total precipitation increment')
    step=30  # thirty native 30s steps per 15min, following interval
    for i in range(0,len(seconds)-1,step):
        value=surface_cape(d.pressure.values[i],t[i],d.qv.values[i],layer_selection='most_cape')
        amount=float(rain[i+step]-rain[i])
        cape.append({'time':str(clock[i]),'cape_J_kg':value,'precipitation_amount_kg_m2':amount,
                     'CAPExP_J_m2':value*amount,'parcel_height_AGL_m':float(z[0]),
                     'parcel_pressure_Pa':float(d.pressure.values[i,0]),'parcel_temperature_K':float(t[i,0]),
                     'parcel_mixing_ratio_kg_kg':float(d.qv.values[i,0]),'parcel_definition':'lowest simulated cell, liquid pseudoadiabatic MetPy; most_cape LFC/EL'})
    therm=pd.DataFrame(cape);therm.to_csv(destination/'thermodynamics_15min.csv',index=False)
    cap=therm[['time','CAPExP_J_m2']].rename(columns={'CAPExP_J_m2':'predictor'})
    cap['method']='CAPExP';cap['status']=np.where(np.isfinite(cap.predictor),'available','unavailable')
    bins['time']=bins.time.astype(str)
    primary=pd.concat([bins[bins.method.isin(['IxG','McCaulF3'])],cap],ignore_index=True)
    five=five[five.method=='McCaulF3'].copy();five['time']=five.time.astype(str)
    for frame in [primary,five]:frame['configuration']='TH_C0';frame['member']=0
    primary.to_csv(destination/'predictors_15min.csv',index=False)
    five.to_csv(destination/'predictors_5min.csv',index=False)
    bins.to_csv(destination/'components_15min.csv',index=False)
    native=xr.Dataset({'w_LS':(('time','height'),w), 'rho_d':(('time','height'),rho),
        'cloud_ice_concentration':(('time','height'),rho*d.qi.values),'graupel_concentration':(('time','height'),rho*d.qg.values),
        'IxG':('time',ixg),'McCaul_F1_LS':('time',f1),'McCaul_F2':('time',f2), 'McCaul_F3_LS':('time',f3),
        'minus15C_height':('time',h),'rainnc':('time',rain)},coords={'time':seconds,'height':z})
    for name,units in {'w_LS':'m s-1','rho_d':'kg m-3','cloud_ice_concentration':'kg m-3','graupel_concentration':'kg m-3',
        'IxG':'kg2 m-5','minus15C_height':'m','rainnc':'kg m-2','height':'m'}.items():native[name].attrs['units']=units
    for name in ['McCaul_F1_LS','McCaul_F2','McCaul_F3_LS']:native[name].attrs['units']='published analysis-cell flashes per 300 seconds; SCM score'
    native.time.attrs['units']='seconds since '+pd.Timestamp(config['start']).strftime('%Y-%m-%d %H:%M:%S')+' UTC'
    native.attrs.update(velocity_support='prescribed LS only; no EDMF addition',CAPE_definition='liquid pseudoadiabatic surface parcel; not asserted isentropic',
        precipitation='rainnc already includes solid water equivalent; never add snownc',source_sha256=before)
    native.to_netcdf(destination/'diagnostics_native.nc')
    if sha256(simulation)!=before:raise IOError('Diagnostic processing modified native model output')
    write_json(destination/'manifest.json',{'simulation':str(simulation),'simulation_sha256':before,'forcing_sha256':sha256(forcing),
        'configuration':'TH_C0','member':0,'case':config['case'],'origin_UTC':config['start'],'native_output_unchanged':True})
    return primary,five


def metric_values(y,p,constant,minutes=15):
    y,p=np.asarray(y,float),np.asarray(p,float)
    if len(y)!=len(p) or not np.isfinite([y,p]).all():raise ValueError('Metrics require an explicit common valid support')
    names=['rmse','mae','skill_zero','skill_constant','pearson','spearman','bias','relative_peak_error','peak_time_error_minutes','accumulated_error']
    result={name:np.nan for name in names};reasons={}
    result['n']=len(y)
    if not len(y):return dict(result,undefined_reasons=json.dumps({'all':'no common valid intervals'}))
    e=p-y;mse=float(np.mean(e**2));zero=float(np.mean(y**2));fixed=float(np.mean((constant-y)**2))
    result.update(rmse=float(np.sqrt(mse)),mae=float(np.mean(abs(e))),bias=float(np.mean(e)),accumulated_error=float(e.sum()))
    for key,den in [('skill_zero',zero),('skill_constant',fixed)]:
        if den>0:result[key]=1-mse/den
        else:reasons[key]='reference MSE is zero'
    if len(y)>1 and np.ptp(y)>0 and np.ptp(p)>0:
        result['pearson']=float(np.corrcoef(y,p)[0,1]);result['spearman']=float(spearmanr(y,p).statistic)
    else:reasons.update(pearson='constant or insufficient series',spearman='constant or insufficient series')
    if y.max()>0:result['relative_peak_error']=float((p.max()-y.max())/y.max())
    else:reasons['relative_peak_error']='observed peak is zero'
    if np.ptp(p)>0 and np.ptp(y)>0:result['peak_time_error_minutes']=int((p.argmax()-y.argmax())*minutes)
    else:reasons['peak_time_error_minutes']='constant series has no distinct peak time'
    result['undefined_reasons']=json.dumps(reasons,sort_keys=True)
    return result


def evaluate_predictions(predictions,observations,constant):
    """One global time mask across methods, configurations and all 40 members."""
    p=predictions.copy();o=observations.copy()
    from scm_lightning.ensemble.experiments import CONFIGURATIONS
    if set(p.method)!=set(METHODS) or set(p.configuration)!=set(CONFIGURATIONS):raise ValueError('Unexpected primary method or physical configuration')
    p.time=pd.to_datetime(p.time,utc=True);o.time=pd.to_datetime(o.time,utc=True)
    if p.duplicated(['time','method','configuration','member']).any() or o.time.duplicated().any():raise ValueError('Duplicate metric support')
    groups=p.groupby(['time','method','configuration']).member.agg(lambda x:tuple(sorted(x)))
    if any(v!=tuple(range(41)) for v in groups):raise ValueError('Expected separate control and exactly 40 paired perturbed members')
    total=p.groupby('time').agg(n=('prediction','size'),valid=('prediction',lambda x:np.isfinite(x).sum()))
    valid=total[(total.n==3*6*41)&(total.valid==3*6*41)].index
    covered=o[o.covered & o['count'].notna()].time
    valid=valid.intersection(covered).sort_values()
    mask=o[['time','covered']].copy();mask['common_valid']=mask.time.isin(valid)
    mask['nonfinite_predictions']=mask.time.map(total.n-total.valid).fillna(3*6*41).astype(int)
    mask['exclusion_reason']=np.where(mask.common_valid,'',np.where(~mask.covered,'observation_not_covered','incomplete_or_nonfinite_primary_predictors'))
    rows=[];stats=[]
    for (method,configuration),g in p[p.member!=0].groupby(['method','configuration']):
        f=g.groupby('time').prediction.agg(['mean','std','count']).rename(columns={'count':'n_members'})
        f.loc[f.n_members!=40,['mean','std']]=np.nan
        f['method']=method;f['configuration']=configuration;stats.append(f.reset_index())
        selected=f.reindex(valid)
        values=metric_values(o.set_index('time').loc[valid,'count'].values,selected['mean'].values,constant)
        # Peak timing uses actual UTC positions when the shared support contains gaps.
        if np.isfinite(values['peak_time_error_minutes']):
            y=o.set_index('time').loc[valid,'count'].values
            values['peak_time_error_minutes']=float((valid[np.argmax(selected['mean'])]-valid[np.argmax(y)]).total_seconds()/60)
        rows.append({'method':method,'configuration':configuration,'coefficient':float(g.scale.iloc[0]),**values})
    return pd.DataFrame(rows),pd.concat(stats,ignore_index=True),mask


def cape_audit(cache):
    rows=[];all_frames=[]
    for path in sorted(Path(cache).glob('*/*/thermodynamics_15min.csv')):
        f=pd.read_csv(path);cape=f.cape_J_kg;rain=f.precipitation_amount_kg_m2
        if (rain.dropna()<0).any():raise ValueError('Negative total precipitation')
        f['configuration']=path.parent.parent.name;f['member']=int(path.parent.name)
        f['CAPExP_J_m2']=cape*rain;all_frames.append(f)
        rows.append({'configuration':f.configuration.iloc[0],'member':f.member.iloc[0],'n':len(f),
            'n_missing_cape':int(cape.isna().sum()),'n_negative_signed_integral':int((cape<0).sum()),'n_positive_cape':int((cape>0).sum()),'n_positive_precipitation':int((rain>0).sum()),
            'n_positive_product':int((cape*rain>0).sum()),'cape_max_J_kg':cape.max(),'total_precipitation_kg_m2':rain.sum(),
            'CAPExP_max_J_m2':(cape*rain).max(),'CAPE_definition':'MetPy liquid pseudoadiabatic; lowest simulated cell',
            'coefficient_applied':False,'rainnc_includes_snow':True,'clock':'native profile at bin start; next exact 900s rainnc increment'})
    if not rows:raise FileNotFoundError('No native-snapshot thermodynamic diagnostic cache')
    out=OUT/'data';out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'CAPE_precipitation_audit.csv',index=False)
    pd.concat(all_frames,ignore_index=True).to_csv(out/'CAPE_precipitation_uncalibrated_15min.csv.gz',index=False)


def tables(config,calibration=None):
    out=OUT/'tables';out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame([{'event':'argentina_20191028','role':'calibration only','physical_configurations':'TH_C0', 'members':'deterministic control 0',
        'start':'2019-10-28T00:00Z','end_exclusive':'2019-10-30T00:00Z'},
        {'event':'argentina_20181213','role':'independent evaluation only','physical_configurations':'TH/WS x C0/CL/CF','members':'40 perturbed + separate control per configuration',
         'start':config['start'],'end_exclusive':config['end']}]).to_csv(out/'Table_01_experimental_design.csv',index=False)
    pd.DataFrame([
        {'method':METHODS[0],'equation':'CAPE * P','reference':'Romps 2014; 2019 eq1','published_coefficient':'2014 CG eta/E is not transferred as total-flash coefficient',
         'adaptation':'MetPy lowest-cell liquid pseudoadiabatic CAPE; not isentropic','predictor_units':'J km-2 per 15min', 'requirements':'Pa, K, kg/kg dry air, total rainnc increment'},
        {'method':METHODS[1],'equation':'integral (rho_d*qi)*(rho_d*qg) dz','reference':'Romps 2019 eq4','published_coefficient':'none universal',
         'adaptation':'whole SCM column, cloud ice only','predictor_units':'kg2 s m-5 per bin','requirements':'rho_d kg/m3; qi and qg kg/kg dry air; dz m'},
        {'method':METHODS[2],'equation':'0.95*0.042*(w_LS_positive*qg_g/kg)[258.15K] + 0.05*0.20*VII','reference':'McCaul 2009 eq3-5',
         'published_coefficient':'k1=.042;k2=.20;r1=.95;r2=.05; original analysis cell per 300s',
         'adaptation':'prescribed LS; nominal published grid area and 300s conversion explicit; extra 2019-only peak mapping',
         'predictor_units':'nominal-cell flash equivalent km-2 per bin','requirements':'collocated w_LS and qg at -15C; VII kg/m2'}
    ]).to_csv(out/'Table_02_lightning_parameterizations.csv',index=False)
    if calibration is not None:calibration.to_csv(out/'Table_03_calibration_coefficients.csv',index=False)
    write_json(out/'metric_conventions.json',{'rmse':'flashes/15min','mae':'flashes/15min','bias':'prediction minus observation, flashes/15min',
        'relative_peak_error':'(predicted maximum - observed maximum)/observed maximum, dimensionless',
        'peak_time_error_minutes':'predicted minus observed peak UTC; earliest peak on ties; constant series -> NaN',
        'accumulated_error':'predicted minus observed flashes on identical valid bins','correlations':'dimensionless; constant/insufficient -> NaN with reason',
        'skills':'1 - model MSE / reference MSE; negative is worse than reference; zero denominator -> NaN',
        'ensemble':'mean and sample std of members 1..40; member 0 separate','precipitation':'total rainnc, never rainnc+snownc'})


def evaluate_event():
    c=load_config(EVALUATION_CONFIG);cal=load_config(CALIBRATION_CONFIG)
    coefficient,meta=load_frozen(FROZEN)
    if c['grouping']!=meta['grouping'] or c['flag_mapping']!=meta['flag_mapping']:raise ValueError('Inconsistent inter-event flash convention')
    cache=OUT/'data/members'
    paths=sorted(cache.glob('*/*/predictors.csv'))
    if len(paths)!=246:raise ValueError('All six existing configurations, 40 members and separate controls required')
    p=normalize_predictors(pd.concat([pd.read_csv(path) for path in paths],ignore_index=True),cal)
    p=apply(p,coefficient,area_km2(c['area']))
    o=pd.read_csv(ROOT/'results/argentina_20181213/observations/counts_15min.csv')
    constant=meta['constant_density_flashes_km2_per_15min']*area_km2(c['area'])
    metric,stats,mask=evaluate_predictions(p,o,constant)
    data=OUT/'data';data.mkdir(parents=True,exist_ok=True)
    p.to_csv(data/'independent_member_predictions.csv.gz',index=False);stats.to_csv(data/'ensemble_statistics.csv',index=False)
    mask.to_csv(data/'common_valid_intervals.csv',index=False);o.to_csv(data/'observations_15min.csv',index=False)
    cube=p[['time','method','configuration','member','predictor','prediction_density','prediction']].copy()
    cube.time=pd.to_datetime(cube.time,utc=True).dt.tz_localize(None)
    ds=cube.set_index(['time','method','configuration','member']).to_xarray()
    factors=coefficient.set_index('method').scale.reindex(ds.method.values)
    ds['calibration_coefficient']=('method',factors.values)
    ds['observation_count']=('time',o.assign(time=pd.to_datetime(o.time,utc=True).dt.tz_localize(None)).set_index('time')['count'].reindex(ds.time.values).values)
    ds['common_valid']=('time',mask.assign(time=pd.to_datetime(mask.time,utc=True).dt.tz_localize(None)).set_index('time').common_valid.reindex(ds.time.values).values.astype('int8'))
    ds.prediction.attrs.update(units='1',long_name='predicted total IC+CG flashes during following 900-second interval')
    ds.prediction_density.attrs.update(units='km-2',interval_seconds=900)
    ds.predictor.attrs['units']='method-dependent; see Table_02_lightning_parameterizations.csv'
    ds.attrs.update(calibration_event='argentina_20191028',evaluation_event='argentina_20181213',
        calibration_configuration='TH_C0',calibration_member=0,calibration_coefficient_sha256=meta['coefficient_sha256'],
        evaluation_area_km2=area_km2(c['area']),interval_seconds=900,member_zero='separate deterministic control, excluded from ensemble mean/std',
        CAPE_definition='MetPy liquid pseudoadiabatic lowest simulated cell; most_cape LFC/EL; no positive floor',
        McCaul_velocity='prescribed LS, no EDMF addition, SCM adaptation',precipitation='rainnc only, includes frozen water equivalent')
    # Fixed-width ASCII avoids Pandas StringDtype / older xarray codec conflicts.
    for key in ['method','configuration']:
        ds[key]=xr.Variable((key,),np.asarray(ds[key].values,dtype='S40'))
    pending=data/'independent_predictions_15min.pending.nc'
    ds.to_netcdf(pending)
    with xr.open_dataset(pending) as checked:
        if dict(checked.sizes)!=dict(time=192,method=3,configuration=6,member=41):raise ValueError('Incomplete forecast NetCDF export')
    pending.replace(data/'independent_predictions_15min.nc')
    tables(c,coefficient)
    metric.to_csv(OUT/'tables/Table_04_independent_evaluation.csv',index=False)
    metric[metric.configuration.isin(['TH_C0','WS_C0'])].to_csv(OUT/'tables/Table_04_reference_C0.csv',index=False)
    contrasts=[]
    for method in METHODS:
        f=metric[metric.method==method].set_index('configuration')
        pairs=[('WS_'+mode,'TH_'+mode) for mode in ['C0','CL','CF']]
        pairs += [(s+'_'+m,s+'_C0') for s in ['TH','WS'] for m in ['CL','CF']]
        for a,b in pairs:
            for name in ['rmse','mae','skill_zero','skill_constant','bias','accumulated_error']:
                contrasts.append({'method':method,'configuration_a':a,'configuration_b':b,'contrast':a+' minus '+b,'metric':name,'difference':f.loc[a,name]-f.loc[b,name]})
    pd.DataFrame(contrasts).to_csv(OUT/'tables/Table_05_microphysics_sensitivity.csv',index=False)
    write_json(OUT/'independent_manifest.json',{'status':'evaluated','calibration_manifest_sha256':sha256(FROZEN/'frozen.json'),
        'calibration_event':meta['calibration_event'],'evaluation_event':'argentina_20181213','calibration_configuration':'TH_C0','calibration_member':0,
        'coefficients_shared_across_all_members_and_physics':True,'coefficient_sha256':meta['coefficient_sha256'],
        'zero_baseline_flashes_per_15min':0,'constant_baseline_flashes_per_15min':constant,'constant_source':'2019 observed mean density times 2018 fixed area',
        'common_valid_intervals':int(mask.common_valid.sum()),'evaluation_observation_sha256':sha256(ROOT/'results/argentina_20181213/observations/counts_15min.csv'),
        'physical_simulations_rerun':0,'new_ensembles':0,'git_commit':git('rev-parse','HEAD')})


def calibrate():
    c=load_config(CALIBRATION_CONFIG)
    run=ROOT/'runs/argentina_20191028/thompson_C0/control'
    identity=json.loads((run/'run_identity.json').read_text())
    if identity['configuration']!='TH_C0' or identity['member']!=0 or not identity['result'].get('physical_preflight_passed'):
        raise ValueError('Calibration requires the validated deterministic TH_C0 control')
    if identity['result']['sha256']!=sha256(run/'output.nc'):raise ValueError('Calibration native output identity mismatch')
    nml=(run/'namelist.mp_offline').read_text()
    for flag in ['microp_scheme="thompson"','couple_mynn_condensate=.false.','mynn_uedmf_mp=0','apply_rt_adv=.false.','legacy_temperature_tendency=.false.']:
        if flag not in nml:raise ValueError('Unexpected calibration physics: '+flag)
    diagnostics=ROOT/'results/argentina_20191028/diagnostics'
    p,p5=control_predictors(run/'output.nc',ROOT/c['paths']['forcing']/'LEONA_SCM.nc',c,diagnostics)
    obs=ROOT/'results/argentina_20191028/observations'
    o=pd.read_csv(obs/'counts_15min.csv');o5=pd.read_csv(obs/'counts_5min.csv')
    coefficients,baseline=fit(c,p,o,p5,o5)
    freeze(FROZEN,c,coefficients,baseline,[CALIBRATION_CONFIG,run/'output.nc',run/'namelist.mp_offline',
        run/'run_identity.json',diagnostics/'manifest.json',obs/'manifest.json',obs/'counts_15min.csv',obs/'counts_5min.csv'])


def prepare():
    preserve_legacy()
    OUT.mkdir(parents=True,exist_ok=True)
    cape_audit(LEGACY/'data/members');tables(load_config(EVALUATION_CONFIG))
    write_json(OUT/'independent_manifest.json',{'status':'awaiting_2019_calibration','calibration_event':'argentina_20191028',
        'evaluation_event':'argentina_20181213','old_calibration_used':False,'legacy_archive':str(LEGACY),
        'pending':['complete 2019 ERA5 acquisition','weather-only column selection','2019 TH_C0 deterministic run','2019 ENTLN grouped flashes and coefficient freeze'],
        'evaluation_event_used_for_calibration':False})


def run_control():
    from scm_lightning.ensemble.experiments import namelist,execute,validate_output
    c=load_config(CALIBRATION_CONFIG)
    forcing=ROOT/c['paths']['forcing']/'LEONA_SCM.nc'
    folder=ROOT/'runs/argentina_20191028/thompson_C0/control';folder.mkdir(parents=True,exist_ok=True)
    output=folder/'output.nc'
    if output.exists():raise FileExistsError('Never overwrite a native simulation')
    executable=ROOT/'build/updraft_sensitivity/mp_offline_driver'
    (folder/'namelist.mp_offline').write_text(namelist('TH_C0',forcing,output).replace('legacy_temperature_tendency=.false.,',
        'legacy_temperature_tendency=.false.,write_process_diagnostics=.true.,'),encoding='utf-8')
    execute(executable,folder)
    duration=(pd.Timestamp(c['end'])-pd.Timestamp(c['start'])).total_seconds()
    write_json(folder/'run_identity.json',{'configuration':'TH_C0','member':0,'forcing_sha256':sha256(forcing),
        'executable_sha256':sha256(executable),'result':validate_output(output,duration)})


def refresh_member_cape(run):
    folder=resolve_path(run['folder']);source=folder/'output.nc'
    old=LEGACY/'data/members'/run['configuration']/f'{run["member"]:02}'
    target=OUT/'data/members'/run['configuration']/f'{run["member"]:02}'
    target.mkdir(parents=True,exist_ok=True)
    tag=target/'independent_cape_identity.json'
    digest=sha256(source)
    identity={'native_sha256':digest,'CAPE_layer_selection':'most_cape LFC and EL; MetPy liquid pseudoadiabatic',
        'diagnostic_code_sha256':sha256(ROOT/'python/scm_lightning/lightning/adaptations.py')}
    if digest!=json.loads((old/'identity.json').read_text())['simulation_sha256']:raise ValueError('Native output identity differs from audited cache')
    if tag.exists() and json.loads(tag.read_text())==identity:return
    with xr.open_dataset(source,decode_times=False) as ds:
        d=ds[['theta','pressure','qv','rainnc','lev','time']].load()
    legacy=pd.read_csv(old/'thermodynamics_15min.csv')
    therm=[];t=d.theta.values*(d.pressure.values/1.e5)**(2/7)
    for j,i in enumerate(range(0,len(d.time)-1,30)):
        value=surface_cape(d.pressure.values[i],t[i],d.qv.values[i],layer_selection='most_cape')
        amount=float(d.rainnc.values[i+30]-d.rainnc.values[i])
        if amount<0:raise ValueError('Reset total precipitation')
        therm.append({'time':str(pd.Timestamp('2018-12-13T00:00Z')+pd.Timedelta(seconds=float(d.time.values[i]))),
            'cape_J_kg':value,'legacy_bottom_top_signed_integral_J_kg':legacy.cape_J_kg.iloc[j],
            'precipitation_amount_kg_m2':amount,'CAPExP_J_m2':value*amount,
            'parcel_pressure_Pa':float(d.pressure.values[i,0]),'parcel_temperature_K':float(t[i,0]),
            'parcel_mixing_ratio_kg_kg':float(d.qv.values[i,0]),'parcel_height_AGL_m':float(d.lev.values[0]),
            'layer_selection':'most_cape LFC/EL; no CAPE floor or thermodynamic adjustment'})
    f=pd.DataFrame(therm);f.to_csv(target/'thermodynamics_15min.csv',index=False)
    p=pd.read_csv(old/'predictors.csv');p.time=pd.to_datetime(p.time,utc=True).astype(str)
    cap=f[['time','CAPExP_J_m2']].rename(columns={'CAPExP_J_m2':'predictor'}).assign(method='CAPExP',configuration=run['configuration'],member=run['member'])
    cap['status']=np.where(cap.predictor.notna(),'available','unavailable')
    pd.concat([p[p.method!='CAPExP'],cap],ignore_index=True).to_csv(target/'predictors.csv',index=False)
    if sha256(source)!=digest:raise IOError('Postprocessing changed native simulation')
    write_json(tag,identity)


def refresh_cape():
    c=load_config(EVALUATION_CONFIG);base=ROOT/c['ensemble']['active_directory']
    runs=json.loads((base/'manifest.json').read_text())['runs']
    if len(runs)!=246:raise ValueError('Expected unchanged 246 native outputs')
    completed=0
    with ProcessPoolExecutor(max_workers=6) as pool:
        futures=[pool.submit(refresh_member_cape,r) for r in runs]
        for f in as_completed(futures):
            f.result();completed+=1
            write_json(OUT/'CAPE_progress.json',{'completed':completed,'expected':246,'model_reruns':0})
            if completed%12==0:print('CAPE diagnostics',completed,'/246',flush=True)
    cape_audit(OUT/'data/members')


def check_main():
    """Technical delivery audit, without fitting or executing a physical model."""
    coefficient,meta=load_frozen(FROZEN)
    p=pd.read_csv(OUT/'data/independent_member_predictions.csv.gz')
    o=pd.read_csv(OUT/'data/observations_15min.csv')
    c=load_config(EVALUATION_CONFIG)
    constant=meta['constant_density_flashes_km2_per_15min']*area_km2(c['area'])
    actual,stats,mask=evaluate_predictions(p,o,constant)
    stored=pd.read_csv(OUT/'tables/Table_04_independent_evaluation.csv')
    pd.testing.assert_frame_equal(actual,stored,check_dtype=False,rtol=1e-12)
    assert (p.groupby('method').scale.nunique()==1).all() and len(coefficient)==3
    tags=list((OUT/'data/members').glob('*/*/independent_cape_identity.json'))
    assert len(tags)==246
    for path in tags:
        tag=json.loads(path.read_text());old=LEGACY/'data/members'/path.parent.relative_to(OUT/'data/members')/'identity.json'
        assert tag['native_sha256']==json.loads(old.read_text())['simulation_sha256']
    assert not git('diff','--','src')
    out=OUT/'audit';out.mkdir(parents=True,exist_ok=True)
    package=ROOT/'python/scm_lightning'
    write_json(out/'implementation_manifest.json',{'git_commit':git('rev-parse','HEAD'),'working_tree_changes_included':True,
        'files':[{'path':str(q.relative_to(ROOT)),'sha256':sha256(q)} for folder in [package,ROOT/'configs',ROOT/'tests'] for q in folder.rglob('*')
                 if q.is_file() and '__pycache__' not in q.parts and q.suffix in ['.py','.yaml','.json']]})
    sources=[('microphysical_processes',ROOT/'results/argentina_20181213/control_diagnostics'),
             ('six_configuration_process_figures',LEGACY/'supplementary'),
             ('deterministic_weak_strong_18_cases',ROOT/'results/argentina_20181213/updraft_sensitivity'),
             ('LS_EDMF_raw_and_isotherm_activity',ROOT/'results/argentina_20181213/vertical_diagnostics'),
             ('calibration_2019',ROOT/'results/argentina_20191028')]
    write_json(OUT/'supplementary/source_index.json',{'physics_supplements_preserved':True,
        'legacy_electrical_coefficients_not_used':True,'sources':[{'kind':kind,'path':str(path),'exists':path.exists()} for kind,path in sources]})
    write_json(out/'independent_validation.json',{'status':'passed','methods':METHODS,'coefficients':3,
        'coefficient_sha256':meta['coefficient_sha256'],'calibration_only_2019_TH_C0_member0':True,
        'existing_2018_outputs_checked_before_and_after_diagnostics':246,'new_2018_simulations':0,
        'common_intervals':int(mask.common_valid.sum()),'excluded_intervals':int((~mask.common_valid).sum()),
        'metric_rows_reconstructed':len(actual),'original_Fortran_unchanged_in_working_tree':True,
        'constant_baseline_source':'2019 observed mean density only','LS_EDMF_summed':False,
        'subdaily_independent_meteorological_observations':'unavailable; CHIRPS2 daily is CFSv2-disaggregated',
        'event_window_confirmation_sha256':sha256(ROOT/'results/argentina_20191028/support/event_window_confirmation.json')})
    print('Independent delivery audit passed: 3 coefficients, 18 metric rows, 246 unchanged native outputs.',flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','refresh-cape','run-control','calibrate','evaluate','check']);a=ap.parse_args()
    {'prepare':prepare,'refresh-cape':refresh_cape,'run-control':run_control,'calibrate':calibrate,'evaluate':evaluate_event,'check':check_main}[a.stage]()


if __name__=='__main__':main()
