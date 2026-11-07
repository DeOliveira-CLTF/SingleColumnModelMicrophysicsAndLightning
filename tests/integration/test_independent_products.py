"""Real local products; synthetic leakage and undefined-metric tests are in unit/."""
import json
import numpy as np
import pandas as pd
import pytest
from scm_lightning.evaluation.independent import OUT,FROZEN,LEGACY,metric_values
from scm_lightning.calibration.independent import load_frozen,METHODS
from scm_lightning.provenance import ROOT,sha256
from scm_lightning.forcing.contracts import area_km2
from scm_lightning.config import load_config

# Preserve the original independent-15min checks against their immutable archive.
_archived=ROOT/'results/argentina_20181213/legacy/article_independent_15min_20261106'
if _archived.exists():OUT=_archived


def test_legacy_products_preserved_with_original_checksums():
    if not (LEGACY/'preservation.json').exists():pytest.skip('Local legacy archive unavailable')
    m=json.loads((LEGACY/'preservation.json').read_text())
    assert m['native_simulations_moved'] is False
    for f in m['files']:assert sha256(LEGACY/f['path'])==f['sha256']


def test_no_native_2018_output_changed_during_cape_refresh():
    tags=list((OUT/'data/members').glob('*/*/independent_cape_identity.json'))
    if len(tags)!=246:pytest.skip('Native CAPE refresh in progress')
    for path in tags:
        m=json.loads(path.read_text());relative=path.parent.relative_to(OUT/'data/members')
        old=json.loads((LEGACY/'data/members'/relative/'identity.json').read_text())
        assert m['native_sha256']==old['simulation_sha256']
        assert m['CAPE_layer_selection'].startswith('most_cape')


def test_calibration_provenance_and_constant_baseline_are_2019_only():
    if not (FROZEN/'frozen.json').exists():pytest.skip('2019 control/calibration not yet available')
    coef,meta=load_frozen(FROZEN)
    assert len(coef)==3 and set(coef.method)==set(METHODS)
    assert set(coef.configuration)=={'TH_C0'} and set(coef.member)=={0}
    assert meta['evaluation_observations_used'] is False
    assert all('20181213' not in r['path'] for r in meta['sources'])
    o=pd.read_csv(ROOT/'results/argentina_20191028/observations/counts_15min.csv')
    assert len(o)==192 and o.covered.all()
    assert meta['constant_density_flashes_km2_per_15min']==pytest.approx(o['count'].mean()/meta['area_km2'])
    for r in meta['sources']:assert sha256(r['path'])==r['sha256']


def test_independent_tables_and_figures_match_shared_frozen_coefficients():
    file=OUT/'tables/Table_04_independent_evaluation.csv'
    if not file.exists():pytest.skip('Independent evaluation pending')
    p=pd.read_csv(OUT/'data/independent_member_predictions.csv.gz');stats=pd.read_csv(OUT/'data/ensemble_statistics.csv')
    o=pd.read_csv(OUT/'data/observations_15min.csv');mask=pd.read_csv(OUT/'data/common_valid_intervals.csv')
    metric=pd.read_csv(file);coef,meta=load_frozen(FROZEN)
    assert len(metric)==18 and set(metric.method)==set(METHODS)
    assert (p.groupby('method').scale.nunique()==1).all()
    import xarray as xr
    with xr.open_dataset(OUT/'data/independent_predictions_15min.nc') as d:
        assert dict(d.sizes)==dict(time=192,method=3,configuration=6,member=41)
        assert d.attrs['calibration_event']=='argentina_20191028'
        labels=[v.decode() if isinstance(v,bytes) else str(v) for v in d.method.values]
        for r in coef.itertuples():assert float(d.calibration_coefficient.isel(method=labels.index(r.method)))==pytest.approx(r.scale)
    for r in coef.itertuples():np.testing.assert_allclose(p[p.method==r.method].scale,r.scale,equal_nan=True)
    valid=mask[mask.common_valid].time
    c=load_config(ROOT/'configs/cases/argentina_20181213.yaml')
    constant=meta['constant_density_flashes_km2_per_15min']*area_km2(c['area'])
    for r in metric.itertuples():
        f=stats[(stats.method==r.method)&(stats.configuration==r.configuration)].set_index('time').loc[valid]
        y=o.set_index('time').loc[valid,'count'].values
        direct=metric_values(y,f['mean'].values,constant)
        for key in ['rmse','mae','skill_zero','skill_constant','pearson','spearman','bias','accumulated_error']:
            assert getattr(r,key)==pytest.approx(direct[key],nan_ok=True,rel=1e-12)
        members=p[(p.method==r.method)&(p.configuration==r.configuration)&(p.member!=0)].groupby('time').prediction
        np.testing.assert_allclose(members.mean().reindex(f.index),f['mean'],rtol=1e-12)
        np.testing.assert_allclose(members.std().reindex(f.index),f['std'],rtol=1e-12)
    figure=OUT/'data/Figure_04_plotted_metrics.csv'
    if figure.exists():pd.testing.assert_frame_equal(pd.read_csv(figure),metric)
    figure=OUT/'data/Figure_03_C0_plotted_statistics.csv'
    if figure.exists():
        actual=pd.read_csv(figure).sort_values(['method','configuration','time']).reset_index(drop=True)
        expected=stats[stats.configuration.isin(['TH_C0','WS_C0'])].sort_values(['method','configuration','time']).reset_index(drop=True)
        pd.testing.assert_frame_equal(actual,expected)


def test_daily_external_precipitation_is_not_used_as_subdaily_truth():
    file=ROOT/'results/argentina_20181213/meteorological_validation/manifest.json'
    if not file.exists():pytest.skip('External precipitation unavailable')
    meta=json.loads(file.read_text());assert meta['independent_of_SCM_and_ERA5_forcing']
    assert not meta['timing_validation_15min'] and not meta['independent_daily_gauge_truth']


def test_2019_forcing_and_control_keep_interfaces_units_complete_utc_and_physics():
    import xarray as xr
    forcing=ROOT/'data/processed/forcing/argentina_20191028/LEONA_SCM.nc'
    if not forcing.exists():pytest.skip('2019 forcing unavailable')
    with xr.open_dataset(forcing,decode_times=False) as d:
        np.testing.assert_array_equal(d.time,[*np.arange(0,172801,600.)])
        assert '2019-10-28' in d.time.attrs['units']
        assert d.wa.attrs['units']=='m s-1' and d.pa_forc.attrs['units']=='Pa'
        assert d.tntheta_adv.attrs['units']=='K s-1' and d.tnrv_adv.attrs['units']=='s-1'
        for name in ['wa','pa_forc','tntheta_adv','tnrv_adv']:assert d[name].dims==('time','lev') and np.isfinite(d[name]).all()
        np.testing.assert_allclose(d.ta,d.theta*(d.pa/1e5)**(2/7),rtol=2e-6)
        assert (np.diff(d.pa_forc.values,axis=-1)<0).all()
    root=ROOT/'runs/argentina_20191028/thompson_C0/control'
    if not (root/'run_identity.json').exists():pytest.skip('2019 native run not complete')
    identity=json.loads((root/'run_identity.json').read_text());assert identity['configuration']=='TH_C0' and identity['member']==0
    assert identity['result']['physical_preflight_passed']
    assert sha256(root/'output.nc')==identity['result']['sha256']
    nml=(root/'namelist.mp_offline').read_text()
    assert 'mynn_uedmf_mp=0' in nml and 'couple_mynn_condensate=.false.' in nml


def test_main_exports_resolution_and_legends_have_no_external_fig_legend():
    file=OUT/'figure_manifest.json'
    if not file.exists():pytest.skip('Independent exports pending')
    from PIL import Image
    meta=json.loads(file.read_text());assert meta['legends']=='inside axes'
    paths=list((OUT/'figures').glob('Figure_*.png'))
    assert len(paths)==(4 if meta['missing_figures_3_4'] else 6)
    for path in paths:
        with Image.open(path) as img:assert min(img.info['dpi'])>=299.9
