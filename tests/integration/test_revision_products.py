import json
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from scm_lightning.provenance import ROOT,sha256
from scm_lightning.evaluation.revision import OUT,FROZEN,PREVIOUS,diagnostics_root
from scm_lightning.calibration.revision import load_frozen,METHODS
from scm_lightning.config import load_config


def test_merge_real_complete_48_hour_support_and_unmodified_controls():
    for event in ['argentina_20181213','argentina_20191028']:
        folder=ROOT/'results'/event/'meteorological_validation';path=folder/'MERGE_validation_manifest.json'
        if not path.exists():pytest.skip('Local official MERGE acquisition required')
        m=json.loads(path.read_text());c=load_config(ROOT/'configs/cases'/(event+'.yaml'))
        f=pd.read_csv(folder/'precipitation_MERGE_ERA5_SCM_hourly.csv')
        assert len(f)==48 and f.common_valid.all() and m['area_wsen']==c['area']
        t=pd.to_datetime(f.time,utc=True);end=pd.to_datetime(f.end_time,utc=True)
        assert t.iloc[0]==pd.Timestamp(c['start']) and end.iloc[-1]==pd.Timestamp(c['end'])
        assert ((end-t).dt.total_seconds()==3600).all()
        assert m['no_daily_to_hourly_interpolation'] and not m['ERA5_independent_observation']
        assert len(m['MERGE_sources'])==48
        for source in m['SCM_sources']:assert sha256(source['path'])==source['sha256']
        with xr.open_dataset(folder/'MERGE_fixed_box_hourly.nc') as d:
            assert d.sizes['time']==48 and d.MERGE_mm.attrs['interval_seconds']==3600
            np.testing.assert_allclose(d.area_mean_precipitation,f.MERGE_mm,rtol=1e-12)
            actual=(d.MERGE_mm*d.cell_overlap_area_km2).sum(['latitude','longitude'])/d.cell_overlap_area_km2.sum()
            np.testing.assert_allclose(actual,f.MERGE_mm,rtol=1e-12)
            assert ((d.time_bounds.values[:,1]-d.time_bounds.values[:,0])/np.timedelta64(1,'h')==1).all()


def test_new_2019_wsm6_control_is_matched_to_thompson_and_numerically_valid():
    root=ROOT/'runs/argentina_20191028'
    if not (root/'wsm6_C0/control/run_identity.json').exists():pytest.skip('Local WSM6 control unavailable')
    th=json.loads((root/'thompson_C0/control/run_identity.json').read_text());ws=json.loads((root/'wsm6_C0/control/run_identity.json').read_text())
    assert ws['configuration']=='WS_C0' and ws['member']==0
    assert ws['forcing_sha256']==th['forcing_sha256'] and ws['executable_sha256']==th['executable_sha256']
    assert ws['result']['physical_preflight_passed'] and ws['result']['steps']==5761
    assert sha256(root/'wsm6_C0/control/output.nc')==ws['result']['sha256']


def test_six_frozen_coefficients_are_2019_only_identifiable_OLS():
    if not (FROZEN/'frozen.json').exists():pytest.skip('v3 calibration pending local diagnostics')
    coef,m=load_frozen(FROZEN)
    assert len(coef)==6 and m['all_identifiable'] and not m['evaluation_observations_used']
    assert set(coef.configuration)=={'TH_C0','WS_C0'} and set(coef.member)=={0}
    assert set(coef.fit_method)=={'nonnegative_proportional_OLS_origin_all_valid_intervals'}
    assert coef[coef.method==METHODS[2]].fit_interval_minutes.eq(5).all()
    assert coef[coef.method!=METHODS[2]].fit_interval_minutes.eq(1).all()
    np.testing.assert_allclose(coef.scale,np.maximum(0,coef.OLS_numerator/coef.OLS_denominator))
    for record in m['sources']:
        assert '20181213' not in record['path']
        assert sha256(record['path'])==record['sha256']


def test_v3_comparisons_have_shared_support_and_matching_netcdf_coefficients():
    if not (OUT/'evaluation_manifest.json').exists():pytest.skip('v3 evaluation pending local diagnostics')
    coef,_=load_frozen(FROZEN);masks=[]
    for base,own in [(OUT,True),(OUT/'supplementary/common_TH_same_MUCAPE_OLS',False)]:
        f=pd.read_csv(base/'data/ensemble_statistics.csv');mask=pd.read_csv(base/'data/common_valid_intervals.csv');masks.append(mask)
        assert set(f.method)==set(METHODS) and f.groupby(['method','configuration']).time.size().eq(2880).all()
        with xr.open_dataset(base/'data/independent_predictions_1min.nc') as d:
            assert dict(d.sizes)=={'time':2880,'method':3,'configuration':6,'member':41,'bounds':2}
            for i,method in enumerate(d.method.values.astype(str)):
                for j,cfg in enumerate(d.configuration.values.astype(str)):
                    scheme=cfg[:2] if own else 'TH'
                    expected=coef[(coef.method==method)&(coef.microphysics==scheme)].scale.iloc[0]
                    assert d.calibration_coefficient.values[i,j]==pytest.approx(expected)
            assert ((d.time_bounds.values[:,1]-d.time_bounds.values[:,0])/np.timedelta64(1,'s')==60).all()
            status=d.MUCAPE_status.values
            assert np.isin(status,np.arange(6)).all()
            assert d.attrs['MetPy_version']=='1.7.1'
        table=pd.read_csv(base/'tables/Table_04_independent_evaluation.csv')
        assert len(table)==18 and table.n.nunique()==1 and table.n.iloc[0]==mask.common_valid.sum()
    pd.testing.assert_frame_equal(*masks)
    assert len(list(diagnostics_root('argentina_20181213').glob('*/*/identity.json')))==246
    for ext in ['png','pdf']:
        assert sha256(OUT/f'figures/Figure_05_microphysical_response.{ext}')==sha256(PREVIOUS/f'figures/Figure_05_microphysical_response.{ext}')


def test_control_virtual_roots_signed_areas_and_zero_invalid_distinction():
    folder=ROOT/'results/argentina_20181213/thermodynamic_revision'
    if not (folder/'comparison_contract.json').exists():pytest.skip('CAPE control reconstruction pending')
    contract=json.loads((folder/'comparison_contract.json').read_text())
    assert contract['status']=='passed' and not contract['lightning_observations_read']
    for cfg in ['TH_C0','WS_C0']:
        f=pd.read_csv(folder/f'CAPE_profile_transition_diagnostics_{cfg}.csv.gz')
        assert len(f)==2880
        valid=np.isfinite(f.MUCAPE_signed_integral_J_kg)&np.isfinite(f.raw_cape_J_kg)
        np.testing.assert_allclose(f.loc[valid,'MUCAPE_signed_integral_J_kg'],f.loc[valid,'raw_cape_J_kg'],rtol=1e-10,atol=1e-7)
        negative=f.raw_cape_J_kg<0
        assert negative.any() and f.loc[negative,'cape_J_kg'].isna().all()
        zero=f.status=='valid_zero';assert f.loc[zero,'cape_J_kg'].eq(0).all()
