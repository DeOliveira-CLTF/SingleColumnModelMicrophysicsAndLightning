"""Derived-only article data: three methods, existing 40-member paired ensemble.

No model execution, Fortran editing, native NetCDF writing, or scenario pooling.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import spearmanr
from scm_lightning.config import load_config
from scm_lightning.provenance import ROOT, resolve_path, sha256, write_json
from scm_lightning.lightning.adaptations import surface_cape, crossing, at_height
from scm_lightning.lightning.postprocess import integrate_bins
from scm_lightning.evaluation import evaluate, split_times

OUT=ROOT/'results/argentina_20181213/article'
CASE=ROOT/'configs/cases/argentina_20181213.yaml'
METHODS=['CAPExP','IxG','McCaulF3']
VERSION='article_native_snapshot_LS_v1'
_forcing=None


def prescribed_w(time,z):
    global _forcing
    if _forcing is None:
        with xr.open_dataset(ROOT/'data/processed/forcing/LEONA_SCM.nc',decode_times=False) as d:
            _forcing=(d.time.values.copy(),d.lev.values.copy(),d.wa.values.copy())
    ft,fz,fw=_forcing
    # Collocated prescribed field at the saved snapshot clock. This is a
    # diagnostic interpolation, not an extra advective tendency.
    height=np.array([np.interp(z,fz,row) for row in fw])
    return np.column_stack([np.interp(time,ft,height[:,k]) for k in range(len(z))])


def member(run):
    folder=resolve_path(run['folder']);source=folder/'output.nc'
    if 'updraft_sensitivity' in source.parts:raise ValueError('Pulse scenario in primary ensemble')
    target=OUT/'data/members'/run['configuration']/f'{run["member"]:02}'
    target.mkdir(parents=True,exist_ok=True)
    tag=target/'identity.json';csv=target/'predictors.csv'
    identity={'version':VERSION,'simulation_sha256':run['sha256'],
              'forcing_sha256':sha256(ROOT/'data/processed/forcing/LEONA_SCM.nc')}
    if csv.exists() and tag.exists() and json.loads(tag.read_text())==identity:
        return str(csv)
    if sha256(source)!=run['sha256']:raise ValueError('Native simulation identity failed')
    with xr.open_dataset(source,decode_times=False) as original:
        d=original[['theta','pressure','qv','qi','qs','qg','rainnc','time','lev']].load()
    z=d.lev.values;seconds=d.time.values
    np.testing.assert_array_equal(seconds,np.arange(5761)*30.)
    edges=np.zeros(len(z)+1)
    for k,center in enumerate(z):edges[k+1]=2*center-edges[k]
    dz=np.diff(edges)
    t=d.theta.values*(d.pressure.values/1.e5)**(2/7)
    rho=d.pressure.values/(287.*t*(1+d.qv.values/(287./461.6)))
    w=prescribed_w(seconds,z)
    height=crossing(t,z,258.15)
    speed=at_height(w,z,height);graupel=at_height(1000*d.qg.values,z,height)
    vii=np.sum(rho*(d.qi.values+d.qs.values+d.qg.values)*dz,axis=1)
    ixg=np.sum(rho**2*d.qi.values*d.qg.values*dz,axis=1)
    f1=.042*np.maximum(speed,0)*graupel
    f3=.95*f1+.05*.20*vii
    origin=pd.Timestamp('2018-12-13T00:00:00Z')
    clock=origin+pd.to_timedelta(seconds,unit='s')
    native=pd.concat([pd.DataFrame({'time':clock,'method':name,'value':values}) for name,values in [('IxG',ixg),('McCaulF3',f3)]],ignore_index=True)
    bins=integrate_bins(native,origin,origin+pd.Timedelta(days=2),15)
    cape=[]
    rain=d.rainnc.values
    if np.any(np.diff(rain)<-1.e-10):raise ValueError('Precipitation reset')
    for i in range(0,5760,30):
        value=surface_cape(d.pressure.values[i],t[i],d.qv.values[i])
        cape.append({'time':clock[i],'method':'CAPExP','predictor':value*(rain[i+30]-rain[i]),
                     'status':'available' if np.isfinite(value) else 'unavailable','cape_J_kg':value,
                     'precipitation_amount_kg_m2':rain[i+30]-rain[i]})
    thermodynamics=pd.DataFrame(cape)
    thermodynamics.to_csv(target/'thermodynamics_15min.csv',index=False)
    bins=pd.concat([bins,thermodynamics[['time','method','predictor','status']]],ignore_index=True)
    bins['configuration']=run['configuration'];bins['member']=run['member']
    # Independent reconciliation against the already validated IxG pipeline.
    old=ROOT/'results/argentina_20181213/members'/run['configuration']/f'{run["member"]:02}'/'predictors.csv'
    prior=pd.read_csv(old);prior=prior[prior.method=='IxG'].sort_values('time')
    np.testing.assert_allclose(bins[bins.method=='IxG'].sort_values('time').predictor,prior.predictor,rtol=1e-12,atol=1e-18)
    bins.to_csv(csv,index=False)
    write_json(tag,identity)
    return str(csv)


def old_metric_audit():
    root=ROOT/'results/argentina_20181213/evaluation'
    h=pd.read_csv(root/'heldout_predictions.csv');metrics=pd.read_csv(root/'metrics.csv')
    contrasts=pd.read_csv(root/'fixed_scale_skill_contrasts.csv')
    rows=[]
    for (method,configuration),f in h[h.member!=0].groupby(['method','configuration']):
        g=f.groupby('time').agg(observation=('count','first'),A=('prediction_A','mean'),B=('prediction_B','mean'))
        for comparison in ['A','B']:
            valid=g[comparison].notna()&g.observation.notna()
            actual=float(np.sqrt(np.mean((g.loc[valid,comparison]-g.loc[valid,'observation'])**2))) if valid.any() else np.nan
            stored=metrics[(metrics.method==method)&(metrics.configuration==configuration)&(metrics.comparison==comparison)].rmse.iloc[0]
            np.testing.assert_allclose(actual,stored,rtol=1e-12,equal_nan=True)
            rows.append({'method':method,'configuration':configuration,'calibration':comparison,'recomputed_rmse':actual,'stored_rmse':stored})
    for row in contrasts.itertuples():
        np.testing.assert_allclose(row.rmse_change,row.rmse-row.reference_rmse,equal_nan=True)
        if np.isfinite(row.rmse):assert abs(row.rmse_change)<=row.rmse+row.reference_rmse+1e-10
    out=OUT/'audit';out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'legacy_metric_reconciliation.csv',index=False)
    write_json(out/'legacy_audit.json',{'stored_metrics_reproduce_series':True,'stored_B_contrasts_reproduce_B_metrics':True,
        'old_figure_mismatch':'Panel (a) used individually calibrated A RMSE; panel (b) used fixed-reference B RMSE differences. They are different forecast systems, not comparable absolute/delta columns.',
        'repair':'New principal Figure 4 uses A RMSE and A contrasts throughout; fixed-scale B is separately labeled in supplementary outputs.',
        'native_data_changed':False,'legacy_files_overwritten':False})


def statistics(predictions,value):
    f=predictions[predictions.member.between(1,40)]
    return f.groupby(['time','method','configuration'])[value].agg(
        n='count',mean='mean',std='std',median='median',p05=lambda x:x.quantile(.05),
        p95=lambda x:x.quantile(.95),minimum='min',maximum='max').reset_index()


def performance(observed,predicted):
    y,x=np.asarray(observed,float),np.asarray(predicted,float);valid=np.isfinite(y)&np.isfinite(x)
    y,x=y[valid],x[valid]
    if not len(y):return dict(n=0,rmse=np.nan,mae=np.nan,bias=np.nan,pearson=np.nan,spearman=np.nan,peak_time_error_minutes=np.nan,accumulated_error=np.nan)
    residual=x-y
    nonconstant=len(y)>1 and np.ptp(y)>0 and np.ptp(x)>0
    original_index=np.flatnonzero(valid)
    return dict(n=len(y),rmse=float(np.sqrt(np.mean(residual**2))),mae=float(np.mean(abs(residual))),
        bias=float(np.mean(residual)),pearson=float(np.corrcoef(y,x)[0,1]) if nonconstant else np.nan,
        spearman=float(spearmanr(y,x).statistic) if nonconstant else np.nan,
        peak_time_error_minutes=int((original_index[np.argmax(x)]-original_index[np.argmax(y)])*15),
        accumulated_error=float(np.sum(x)-np.sum(y)))


def assemble():
    if load_config(CASE).get('role')=='evaluation_only':
        raise ValueError('2018 is evaluation-only. Use scripts/evaluate.py independent evaluate; within-event calibration is archived.')
    c=load_config(CASE);paths=sorted((OUT/'data/members').glob('*/*/predictors.csv'))
    if len(paths)!=246:raise ValueError('All 40 members + separate control in all six configurations are required')
    p=pd.concat([pd.read_csv(path) for path in paths],ignore_index=True)
    o=pd.read_csv(ROOT/'results/argentina_20181213/observations/counts_15min.csv')
    assert len(o)==192 and o.covered.all() and not o['count'].isna().any()
    assert not p.duplicated(['time','method','configuration','member']).any()
    assert (p.groupby(['time','method','configuration']).member.nunique()==41).all()
    assert set(p.method)==set(METHODS)
    result=evaluate(p,o,c['reference_configuration'],c['split']['development_end'],c['split']['evaluation_start'])
    coefficient=result[0];held=result[1]
    p.time=pd.to_datetime(p.time,utc=True);o.time=pd.to_datetime(o.time,utc=True)
    predictions=p.merge(coefficient[['method','configuration','scale','status']].rename(columns={'status':'fit_status'}),on=['method','configuration'],validate='many_to_one')
    reference=coefficient[coefficient.configuration==c['reference_configuration']][['method','scale']].rename(columns={'scale':'reference_scale'})
    predictions=predictions.merge(reference,on='method',validate='many_to_one')
    predictions['prediction_A']=predictions.predictor*predictions.scale
    predictions['prediction_B']=predictions.predictor*predictions.reference_scale
    predictions['split']=split_times(predictions.time,c['split']['development_end'],c['split']['evaluation_start'])
    if not np.isfinite(predictions[['prediction_A','prediction_B']]).all().all():
        raise ValueError('Unavailable calibrated primary forecast; do not silently complete requested figures')
    stats=statistics(predictions,'prediction_A');statsB=statistics(predictions,'prediction_B')
    assert (stats.n==40).all() and (statsB.n==40).all()
    tables=OUT/'tables';data=OUT/'data'
    tables.mkdir(parents=True,exist_ok=True)
    predictions.to_csv(data/'calibrated_member_predictions.csv.gz',index=False)
    stats.to_csv(data/'ensemble_statistics_A.csv',index=False);statsB.to_csv(data/'ensemble_statistics_B.csv',index=False)
    o.to_csv(data/'observations_15min.csv',index=False)
    coefficient.to_csv(data/'calibration_coefficients.csv',index=False)
    metrics=[]
    for protocol,frame in [('A',stats),('B',statsB)]:
        for (method,configuration),f in frame.groupby(['method','configuration']):
            f=f[f.time>=pd.Timestamp(c['split']['evaluation_start'])].sort_values('time').merge(o[['time','count']],on='time',validate='one_to_one')
            row=coefficient[(coefficient.method==method)&(coefficient.configuration==configuration)].iloc[0]
            scale=row.scale if protocol=='A' else reference[reference.method==method].reference_scale.iloc[0]
            assert len(f)==92
            values=performance(f['count'],f['mean'])
            metrics.append({'method':method,'configuration':configuration,'microphysics':'Thompson' if configuration.startswith('TH') else 'WSM6',
                'coupling':configuration.split('_')[1],'calibration':protocol,'coefficient':scale,**values})
    metric=pd.DataFrame(metrics)
    contrasts=[]
    for protocol in ['A','B']:
        f=metric[metric.calibration==protocol]
        for method in METHODS:
            lookup=f[f.method==method].set_index('configuration').rmse
            pairs=[('WS_'+mode,'TH_'+mode,'microphysics',mode) for mode in ['C0','CL','CF']]
            pairs += [(scheme+'_'+mode,scheme+'_C0','coupling',scheme+' '+mode+'-C0') for scheme in ['TH','WS'] for mode in ['CL','CF']]
            for a,b,kind,label in pairs:
                delta=lookup[a]-lookup[b]
                assert abs(delta)<=lookup[a]+lookup[b]+1e-10
                forecasts=stats if protocol=='A' else statsB
                left=forecasts[(forecasts.method==method)&(forecasts.configuration==a)&(forecasts.time>=pd.Timestamp(c['split']['evaluation_start']))][['time','mean']]
                right=forecasts[(forecasts.method==method)&(forecasts.configuration==b)&(forecasts.time>=pd.Timestamp(c['split']['evaluation_start']))][['time','mean']]
                pair=left.merge(right,on='time',validate='one_to_one',suffixes=('_a','_b'))
                assert len(pair)==92
                distance=float(np.sqrt(np.mean((pair.mean_a-pair.mean_b)**2)))
                assert abs(delta)<=distance+1e-8
                contrasts.append({'calibration':protocol,'method':method,'kind':kind,'contrast':label,'configuration_a':a,'configuration_b':b,
                                  'rmse_a':lookup[a],'rmse_b':lookup[b],'delta_rmse':delta,'forecast_RMS_distance':distance})
    metric[metric.calibration=='A'].to_csv(tables/'Table_03_observational_performance.csv',index=False)
    metric[metric.calibration=='A'].sort_values('rmse').to_csv(tables/'Table_03_observational_performance_ranked.csv',index=False)
    metric[metric.calibration=='B'].to_csv(data/'fixed_scale_B_performance.csv',index=False)
    pd.DataFrame(contrasts).to_csv(data/'rmse_contrasts.csv',index=False)
    # Numerical audits are completed before rendering any principal figure.
    for (method,configuration),f in predictions.groupby(['method','configuration']):
        d=f[(f.member!=0)&(f.split=='development')].groupby('time').predictor.mean()
        obs=o.set_index('time').loc[d.index,'count']
        independent=max(0,float(np.sum(d.values*obs.values)/np.sum(d.values*d.values)))
        a=coefficient[(coefficient.method==method)&(coefficient.configuration==configuration)].scale.iloc[0]
        np.testing.assert_allclose(a,independent,rtol=1e-12)
        row=metric[(metric.method==method)&(metric.configuration==configuration)&(metric.calibration=='A')].iloc[0]
        series=f[(f.member!=0)&(f.split=='evaluation')].pivot(index='time',columns='member',values='prediction_A').sort_index()
        y=o.set_index('time').loc[series.index,'count'].to_numpy()
        direct=np.sqrt(np.sum((series.to_numpy().mean(axis=1)-y)**2)/92)
        np.testing.assert_allclose(direct,row.rmse,rtol=1e-12)
        np.testing.assert_allclose(row.accumulated_error,92*row.bias,atol=1e-8)
    write_json(OUT/'audit/current_statistics.json',{'status':'passed','members_per_distribution':40,'control_excluded':True,
        'methods':METHODS,'calibration_bins':96,'gap_bins':4,'evaluation_bins':92,
        'RMSE_independently_reconstructed':True,'coefficient_reconstructed_from_development_only':True,
        'contrast_triangle_bounds':True,'reverse_triangle_forecast_distance_bound':True,'counts_vs_accumulation_reconciled':True,'extra_area_or_time_factor':False,
        'main_protocol':'A: one development ensemble-mean scale for each method/configuration; same scale on each member',
        'fixed_reference_protocol':'B: TH_C0 scale frozen across physical configurations; supplementary only',
        'display_band':'mean +/- sample std (ddof=1); lower edge clipped to zero only in figures',
        'CAPE_clock':'native simulated snapshot at bin start, lowest-cell MetPy pseudoadiabatic parcel; times exact native output clock',
        'CAPExP_integration':'CAPE_at_bin_start * exact rainnc difference over following 900 seconds; no second 900 factor',
        'McCaulF3_integration':'native-snapshot F3 trapezoid over 900s, raw published cell/5min score * seconds; empirical coefficient maps this once to flash counts',
        'LS_support':'unperturbed full ERA5 wa interpolated to snapshot time and model height; no EDMF or weak/strong addition'})


def main():
    if load_config(CASE).get('role')=='evaluation_only':
        raise ValueError('Use scripts/evaluate.py independent prepare|evaluate. No 2018 fitting is allowed.')
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--workers',type=int,default=6);ap.add_argument('--assemble-only',action='store_true')
    args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    old_metric_audit()
    if not args.assemble_only:
        c=load_config(CASE);base=ROOT/c['ensemble']['active_directory']
        checks=json.loads((base/'run_checks.json').read_text());plan=json.loads((base/'manifest.json').read_text())
        assert len(checks)==246 and all(r['physical_preflight_passed'] for r in checks)
        assert len(plan['runs'])==246
        records=[dict(run,sha256=next(r['sha256'] for r in checks if r['configuration']==run['configuration'] and r['member']==run['member'])) for run in plan['runs']]
        write_json(OUT/'audit/source_identity.json',{'ensemble_manifest_sha256':sha256(base/'manifest.json'),'run_checks_sha256':sha256(base/'run_checks.json'),
            'observations_sha256':sha256(ROOT/'results/argentina_20181213/observations/counts_15min.csv'),
            'configuration_sha256':sha256(CASE),'version':VERSION,'model_reruns':0,'weak_strong_in_primary':False})
        finished=0
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(member,run) for run in records]
            for future in as_completed(futures):
                future.result();finished+=1
                write_json(OUT/'progress.json',{'completed':finished,'expected':246,'model_reruns':0})
                if finished%12==0:print('Derived member products',finished,'/246',flush=True)
    assemble()
    print('Primary three-method analysis and metric audit completed.',flush=True)


if __name__=='__main__':main()
