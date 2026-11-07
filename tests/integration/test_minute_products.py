import json
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from scm_lightning.provenance import ROOT,sha256
from scm_lightning.evaluation.minute import OUT,PREVIOUS,FROZEN,config,CASE
from scm_lightning.calibration.independent import load_frozen,METHODS
from scm_lightning.evaluation.independent import metric_values
from scm_lightning.forcing.contracts import area_km2

# Retain the v2 methodological checks against their preserved products.
_v2=ROOT/'results/argentina_20181213/legacy/article_v2_1min_most_cape'
if _v2.exists():OUT=_v2


def test_real_minute_flash_counts_conserve_totals_and_are_not_uniform_splits():
    for event in ['argentina_20181213','argentina_20191028']:
        root=ROOT/'results'/event/'observations';path=root/'counts_1min.csv'
        if not path.exists():pytest.skip('One-minute observations not processed locally')
        one=pd.read_csv(path);old=pd.read_csv(root/'counts_15min.csv')
        one.time=pd.to_datetime(one.time,utc=True);old.time=pd.to_datetime(old.time,utc=True)
        assert len(one)==2880 and (one['count'].dropna()%1==0).all()
        summed=one.groupby(one.time.dt.floor('15min'))['count'].sum().reindex(old.time)
        np.testing.assert_array_equal(summed,old['count'])
        assert one.groupby(one.time.dt.floor('15min'))['count'].nunique().max()>1


def test_native_ice_proxies_conserve_integrals_across_resolution():
    new=OUT/'data/members/TH_C0/00/predictors_1min.csv'
    if not new.exists():pytest.skip('Native minute predictor not yet available')
    old=pd.read_csv(PREVIOUS/'data/members/TH_C0/00/predictors.csv');new=pd.read_csv(new)
    for f in [old,new]:f.time=pd.to_datetime(f.time,utc=True)
    for method in ['IxG','McCaulF3']:
        a=new[new.method==method];b=old[old.method==method].sort_values('time')
        actual=a.groupby(a.time.dt.floor('15min')).predictor.sum(min_count=15).reindex(b.time)
        np.testing.assert_allclose(actual,b.predictor,rtol=2e-12,atol=1e-17,equal_nan=True)


def test_minute_cape_and_rain_use_actual_native_profiles_and_increments():
    from scm_lightning.lightning.adaptations import surface_cape
    path=OUT/'data/members/TH_C0/00/thermodynamics_1min.csv'
    if not path.exists():pytest.skip('Native minute thermodynamics not yet available')
    frame=pd.read_csv(path)
    identity=json.loads((path.parent/'manifest.json').read_text())
    with xr.open_dataset(identity['simulation'],decode_times=False) as d:
        np.testing.assert_allclose(frame.precipitation_amount_kg_m2,np.diff(d.rainnc.values[::2]),rtol=1e-12,atol=1e-14)
        # Include clocks between the archived 15-minute profiles: no interpolation.
        for minute in [1,17,421,1399,2879]:
            p=d.pressure.values[2*minute];t=d.theta.values[2*minute]*(p/1e5)**(2/7)
            actual=surface_cape(p,t,d.qv.values[2*minute],layer_selection='most_cape')
            assert frame.cape_J_kg.iloc[minute]==pytest.approx(actual,rel=1e-12,nan_ok=True)
    np.testing.assert_allclose(frame.CAPExP_J_m2,frame.cape_J_kg*frame.precipitation_amount_kg_m2,rtol=1e-12,equal_nan=True)


def test_minute_calibration_remains_2019_control_only():
    if not (FROZEN/'frozen.json').exists():pytest.skip('One-minute calibration unavailable')
    coef,m=load_frozen(FROZEN)
    assert len(coef)==3 and m['interval_minutes']==1 and m['configuration']=='TH_C0' and m['member']==0
    assert not m['evaluation_observations_used']
    assert set(coef[coef.method!='McCaulF3_LS_1D_adapted'].n_intervals)=={2880}
    for r in m['sources']:assert '20181213' not in r['path'] and sha256(r['path'])==r['sha256']


def test_minute_ensemble_shared_scales_and_metrics_csv_netcdf():
    path=OUT/'data/independent_predictions_1min.nc'
    if not path.exists():pytest.skip('One-minute ensemble processing not complete')
    p=pd.read_csv(OUT/'data/independent_member_predictions.csv.gz')
    o=pd.read_csv(OUT/'data/observations_1min.csv');stats=pd.read_csv(OUT/'data/ensemble_statistics.csv')
    mask=pd.read_csv(OUT/'data/common_valid_intervals.csv');valid=mask[mask.common_valid].time
    table=pd.read_csv(OUT/'tables/Table_04_independent_evaluation.csv');coef,meta=load_frozen(FROZEN)
    assert len(table)==18 and (p.groupby('method').scale.nunique()==1).all()
    baseline=meta['constant_density_flashes_km2_per_interval']*area_km2(config(CASE)['area'])
    for r in table.itertuples():
        f=stats[(stats.method==r.method)&(stats.configuration==r.configuration)].set_index('time').loc[valid]
        y=o.set_index('time').loc[valid,'count'].values
        direct=metric_values(y,f['mean'].values,baseline,minutes=1)
        for name in ['rmse','mae','skill_zero','skill_constant','bias','pearson','spearman','accumulated_error']:
            assert getattr(r,name)==pytest.approx(direct[name],nan_ok=True,rel=1e-12)
    with xr.open_dataset(path) as d:
        assert d.sizes['time']==2880 and d.sizes['member']==41 and d.attrs['interval_seconds']==60
        assert d.prediction.attrs['units']=='min-1'
        np.testing.assert_array_equal(d.time_bounds.values[:,0],d.time.values)
        assert ((d.time_bounds.values[:,1]-d.time_bounds.values[:,0])/np.timedelta64(1,'s')==60).all()


def test_sounding_priority_separate_stations_and_era5_precipitation_support():
    from scm_lightning.observations.profiles import selected
    if not (ROOT/'data/processed/radiosondes/inventory.json').exists():pytest.skip('Normalized profiles unavailable')
    campaign=selected('argentina_20181213');calibration=selected('argentina_20191028')
    assert len(campaign)==8 and set(r['source'] for r in campaign)=={'CSU_RELAMPAGO_in_projetoCHUVA'}
    assert {r['station'] for r in calibration}=={'83928','83937','87576'}
    assert not any(r['cases']['argentina_20191028']['initial_condition_eligible'] for r in calibration)
    path=ROOT/'results/argentina_20181213/meteorological_validation/precipitation_comparison_ERA5_hourly.csv'
    if path.exists():
        f=pd.read_csv(path);assert len(f)==48*6
        t=pd.to_datetime(f.time,utc=True);e=pd.to_datetime(f.end_time,utc=True)
        assert ((e-t).dt.total_seconds()==3600).all()
        assert f.groupby('configuration').ERA5_mm.sum().nunique()==1
