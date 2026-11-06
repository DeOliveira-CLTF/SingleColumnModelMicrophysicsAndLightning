import copy
import numpy as np
import pandas as pd
import pytest
from scm_lightning.config import load_config
from scm_lightning.provenance import ROOT
from scm_lightning.download import era5_request
from scm_lightning.calibration.independent import METHODS, fit, normalize_predictors, apply, freeze, load_frozen, nominal_mccaul_cell_area
from scm_lightning.evaluation.independent import metric_values, evaluate_predictions
from scm_lightning.forcing.contracts import area_km2


def inputs():
    c=load_config(ROOT/'configs/cases/argentina_20191028.yaml')
    c['area']=[-60,-33,-59,-32];c['area_selection'].update(status='fixed_before_lightning_evaluation',decision_sha256='synthetic_test_only')
    t=pd.date_range(c['start'],c['end'],freq='15min',inclusive='left').astype(str)
    five=pd.date_range(c['start'],c['end'],freq='5min',inclusive='left').astype(str)
    rows=[dict(time=stamp,method=m,configuration='TH_C0',member=0,predictor=float(i+1)) for m in ['CAPExP','IxG','McCaulF3'] for i,stamp in enumerate(t)]
    p=pd.DataFrame(rows)
    p5=pd.DataFrame([dict(time=stamp,method='McCaulF3',configuration='TH_C0',member=0,predictor=float(i+1)) for i,stamp in enumerate(five)])
    o=pd.DataFrame({'time':t,'count':np.arange(1,len(t)+1)*2.,'covered':True})
    o5=pd.DataFrame({'time':five,'count':np.arange(1,len(five)+1)*2.,'covered':True})
    return c,p,o,p5,o5


@pytest.mark.parametrize('mutation',['2018_observation','2018_predictor','member','physics','event','unfrozen_area'])
def test_calibration_rejects_leakage_and_nonreference_inputs(mutation):
    c,p,o,p5,o5=inputs()
    if mutation=='2018_observation':o.loc[0,'time']='2018-12-13T00:00:00Z'
    if mutation=='2018_predictor':p.loc[0,'time']='2018-12-13T00:00:00Z'
    if mutation=='member':p.loc[0,'member']=1
    if mutation=='physics':p.loc[0,'configuration']='WS_C0'
    if mutation=='event':c['event_id']='argentina_20181213'
    if mutation=='unfrozen_area':c['area']=None
    with pytest.raises(ValueError):fit(c,p,o,p5,o5)


def test_three_coefficients_only_and_mccaul_peak_uses_five_minutes():
    c,p,o,p5,o5=inputs();coefficient,baseline=fit(c,p,o,p5,o5)
    assert len(coefficient)==3 and set(coefficient.method)==set(METHODS)
    assert set(coefficient.configuration)=={'TH_C0'} and set(coefficient.member)=={0}
    row=coefficient.set_index('method').loc[METHODS[2]]
    assert row.fit_interval_minutes==5 and row.n_intervals==576
    assert row.scale==pytest.approx(2*300*nominal_mccaul_cell_area(c)/area_km2(c['area']))
    assert baseline==pytest.approx(o['count'].mean()/area_km2(c['area']))


def test_explicit_unit_conversions_and_published_constants_preserved():
    c,p,*_=inputs();saved=p.copy(deep=True);f=normalize_predictors(p,c)
    pd.testing.assert_frame_equal(p,saved)
    assert f[f.method=='CAPExP'].predictor.iloc[0]==1e6
    assert f[f.method=='IxG'].predictor.iloc[0]==1
    assert f[f.method==METHODS[2]].predictor.iloc[0]==pytest.approx(1/(300*nominal_mccaul_cell_area(c)))
    assert c['calibration']['published_McCaul_coefficients']==dict(k1=.042,k2=.20,r1=.95,r2=.05)


def test_unidentifiable_predictor_never_creates_a_coefficient():
    c,p,o,p5,o5=inputs();p.loc[p.method=='CAPExP','predictor']=0
    coef,_=fit(c,p,o,p5,o5);row=coef.set_index('method').loc['CAPExP']
    assert np.isnan(row.scale) and row.status=='identically_zero_unrecoverable'


def test_one_coefficient_shared_across_all_members_configurations_and_area_conversion():
    c,p,o,p5,o5=inputs();coef,_=fit(c,p,o,p5,o5)
    normalized=normalize_predictors(p,c);normalized['status']='available'
    f=pd.concat([normalized.assign(configuration=cfg,member=member) for cfg in ['TH_C0','WS_CF'] for member in [0,1,40]])
    a=apply(f,coef,100.);b=apply(f,coef,200.)
    assert (a.groupby('method').scale.nunique()==1).all()
    np.testing.assert_allclose(b.prediction,2*a.prediction)
    np.testing.assert_allclose(b.prediction_density,a.prediction_density)


def test_frozen_coefficients_have_checksum_and_cannot_be_overwritten(tmp_path):
    c,p,o,p5,o5=inputs();coef,baseline=fit(c,p,o,p5,o5)
    freeze(tmp_path,c,coef,baseline,[])
    actual,meta=load_frozen(tmp_path);assert meta['evaluation_observations_used'] is False
    pd.testing.assert_frame_equal(actual,coef)
    with pytest.raises(FileExistsError):freeze(tmp_path,c,coef,baseline,[])
    (tmp_path/'coefficients.csv').write_text('tampered')
    with pytest.raises(ValueError,match='checksum'):load_frozen(tmp_path)


def test_skill_sign_independent_constant_and_undefined_metrics():
    r=metric_values([1,2],[4,4],3)
    assert r['skill_zero']<0 and r['skill_constant']<0
    assert r['bias']==2.5 and r['accumulated_error']==5
    assert np.isnan(r['pearson']) and 'constant' in r['undefined_reasons']
    z=metric_values([0,0],[0,0],0)
    assert np.isnan(z['skill_zero']) and np.isnan(z['skill_constant'])
    assert np.isnan(z['relative_peak_error']) and np.isnan(z['peak_time_error_minutes'])


def test_global_common_mask_and_control_excluded_from_ensemble():
    t=pd.date_range('2018-12-13',periods=3,freq='15min',tz='UTC')
    p=pd.DataFrame([dict(time=stamp,method=m,configuration=s+'_'+mode,member=k,prediction=(1.e9 if k==0 else float(i+1)),scale=2.)
        for i,stamp in enumerate(t) for m in METHODS for s in ['TH','WS'] for mode in ['C0','CL','CF'] for k in range(41)])
    p.loc[(p.time==t[1])&(p.method=='IxG')&(p.configuration=='WS_CF')&(p.member==1),'prediction']=np.nan
    o=pd.DataFrame({'time':t,'count':[1.,2.,3.],'covered':True})
    metric,stats,mask=evaluate_predictions(p,o,5.)
    assert set(metric.n)=={2} and mask.common_valid.tolist()==[True,False,True]
    assert stats.loc[stats.time==t[0],'mean'].eq(1.).all()
    assert (metric.rmse==0).all()


def test_era5_request_is_configurable_and_preserves_2018_default():
    c,*_=inputs();r=era5_request(c,single=True)
    assert r['year']==['2019'] and r['month']==['10'] and r['day']==['27','28','29','30','31']
    assert 'total_precipitation' in r['variable'] and r['area']==[-25,-70,-45,-45]
    assert era5_request({})['day']==['12','13','14','15']


def test_multilayer_cape_retains_signed_legacy_and_selects_a_buoyant_layer_without_changing_state():
    import json
    from scm_lightning.lightning.adaptations import surface_cape
    f=json.loads((ROOT/'tests/fixtures/cape_multilayer_model_profile.json').read_text())
    p,t,r=map(np.array,[f['pressure_Pa'],f['temperature_K'],f['mixing_ratio_kg_kg']])
    old=[a.copy() for a in [p,t,r]]
    signed=surface_cape(p,t,r);cape=surface_cape(p,t,r,layer_selection='most_cape')
    assert signed==pytest.approx(-110.57973,rel=2e-5)
    assert cape==pytest.approx(213.515192,rel=2e-5)
    for a,b in zip([p,t,r],old):np.testing.assert_array_equal(a,b)
    # A strongly stable, dry surface parcel must remain zero; no positive floor.
    assert surface_cape([100000,90000,80000,70000,60000,50000],
        [280,285,280,275,270,265],[.001]*6,layer_selection='most_cape')==0


def test_daily_total_precipitation_uses_only_cumulative_rainnc_and_rejects_resets():
    from scm_lightning.evaluation.meteorology import daily_amount
    assert daily_amount([0,43200,86400],[2,5,9],0,86400)==7
    with pytest.raises(ValueError):daily_amount([0,86400],[2,1],0,86400)
    with pytest.raises(ValueError):daily_amount([0,43200],[0,1],0,86400)


def test_first_event_total_counts_are_causal_at_end_of_target_window():
    from scm_lightning.observations import group_flashes,aggregate
    t=pd.Timestamp('2019-10-30T00:00Z')
    events=pd.DataFrame({'time':[t-pd.Timedelta('0.2s'),t+pd.Timedelta('0.2s'),t+pd.Timedelta('0.3s')],
        'event_id':['a','b','c'],'latitude':[-32,-32,-32],'longitude':[-60,-60,-50],'kind':['IC','CG','IC']})
    before,_=group_flashes(events.iloc[:1]);after,_=group_flashes(events)
    args=('2019-10-29T23:45Z','2019-10-30T00:00Z',[-61,-33,-59,-31],[True],15)
    pd.testing.assert_series_equal(aggregate(before,*args)['count'],aggregate(after,*args)['count'])
    assert after.iloc[0].kind=='CG'  # Later events may change type/multiplicity, never earlier total anchor count.
