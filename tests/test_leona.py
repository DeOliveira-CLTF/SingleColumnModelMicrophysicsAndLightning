import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from leona.contracts import *
from leona.lightning import group_flashes,read_events,haversine_km,aggregate
from leona.evaluate import evaluate,fit_scale,split_times
from leona.proxies import ice_product,ice_flux,mccaul,diagnose
from leona.ensemble import joint_perturbations

def events(milliseconds, latitudes=None, kinds=None):
    n=len(milliseconds)
    return pd.DataFrame({'time':pd.Timestamp('2018-12-13T23:59:59.9Z')+pd.to_timedelta(milliseconds,unit='ms'),
        'latitude':latitudes if latitudes is not None else [0.]*n,'longitude':[0.]*n,
        'kind':kinds if kinds else ['IC']*n,'event_id':[str(i) for i in range(n)]})

def test_theta_temperature_contract():
    p=np.array([100000.,50000.]); t=np.array([300.,250.])
    np.testing.assert_allclose(theta_from_temperature(t,p)*exner(p),t)
    assert theta_from_temperature(1.,50000.) != 1.

def test_humidity_tendency_chain_rule():
    q=.02; dq=1e-6
    assert q_to_r(q) == pytest.approx(q/(1-q))
    assert (q_to_r(q+dq*1e-3)-q_to_r(q))/1e-3 == pytest.approx(q_tendency_to_r(dq,q),rel=1e-6)

def test_horizontal_advection_has_no_vertical_component():
    a=xr.DataArray(np.broadcast_to(np.arange(3.)[:,None],(3,3)),dims=['latitude','longitude'],coords={'latitude':[0.,1.,2.],'longitude':[0.,1.,2.]})
    np.testing.assert_allclose(horizontal_advection(a,0.,2.),-2*180/np.pi/6371000)

def test_spherical_box_weights_include_partial_edge_cells():
    from leona.spatial import box_weights
    box=[-.75,-.75,.75,.75]
    weights=box_weights([-1.,0.,1.],[-1.,0.,1.],box)
    assert float(weights.sum())==pytest.approx(area_km2(box))
    assert float(weights[1,0])<float(weights[1,1])

def test_total_precipitation_not_snow_added():
    np.testing.assert_allclose(precipitation_rate([0,2,3],[0,2,4]),[1,.5])
    with pytest.raises(ValueError): precipitation_rate([2,1],[0,1])

def test_fortran_three_digit_exponent_is_preserved(tmp_path):
    from leona.io import numeric_csv
    p=tmp_path/'fortran.csv';p.write_text('a,b\n1.2345678901234567-107,2.0E-003\n')
    f=numeric_csv(p)
    assert f.a.iloc[0]==pytest.approx(1.2345678901234567e-107,abs=0)
    assert f.b.iloc[0]==.002

def test_mass_concentration_integral():
    assert ice_product(np.array([2.,3.]),np.array([.1,.1]),np.array([.2,.2]),[10,20]) == pytest.approx(4.4)

def test_isotherm_lowest_crossing_and_missing():
    assert isotherm([270,250,270,240],[0,100,200,300],260,[0,2,4,6]) == 1
    assert np.isnan(isotherm([280,270],[0,100],260,[0,1]))

def test_mccaul_g_per_kg_and_minus15():
    # w=10, qg=.001kg/kg=1g/kg -> F1=.42 at -15C.
    assert mccaul([270,250],[0,100],[.001,.001],[10,10],2) == pytest.approx(.95*.42+.05*.4)

def test_ifluxt_cloud_ice_260_threshold():
    assert ice_flux([270,250],[0,100],np.ones(2),np.ones(2)*.001,np.ones(2)*2) == pytest.approx(.002)
    assert ice_flux([270,250],[0,100],np.ones(2),np.ones(2)*.001,np.ones(2)) == 0

def test_no_mutation_and_unavailable_not_zero():
    a=np.ones(2)*.001; before=a.copy()
    rows=diagnose({'rho_d':np.ones(2),'qi':a,'qg':a,'dz':[10,10]}, {'mass_basis':'dry_air'})
    np.testing.assert_array_equal(a,before)
    assert next(r for r in rows if r['method']=='PW10')['status']=='unavailable'

def test_anchor_no_unlimited_chain_mixed_midnight():
    f,a=group_flashes(events([0,600,1200],kinds=['IC','CG','IC']))
    assert list(f.multiplicity)==[2,1]
    assert list(f.kind)==['CG','IC']
    assert len(a)==3 and a.event_id.nunique()==3
    assert f.iloc[0].end.day==14

def test_temporal_inclusive_boundary():
    f,_=group_flashes(events([0,700,700.000001]))
    assert list(f.multiplicity)==[2,1]

def test_space_anchored_not_neighbor_chained():
    f,_=group_flashes(events([0,100,200],latitudes=[0, .08,.16]))
    assert list(f.multiplicity)==[2,1]

def test_spatial_boundary_and_tie_break():
    angle=10/6371*180/np.pi
    assert haversine_km(0,0,angle,0)==pytest.approx(10)
    f,_=group_flashes(events([0,100],latitudes=[0,angle]))
    assert len(f)==1
    f,_=group_flashes(events([0,100,200],latitudes=[0,.15,.075]))
    assert list(f.multiplicity)==[2,1]

def test_azambuja_kept_separate():
    f,_=group_flashes(events([0,600,1200]),'azambuja2017_v1')
    assert list(f.multiplicity)==[3]

def test_csv_nanoseconds_duplicate_provenance(tmp_path):
    f=pd.DataFrame({'year':[2018]*2,'month':[12]*2,'day':[14]*2,'hour':[0]*2,'minute':[0]*2,'second':[0]*2,
       'nanosec':[123456789]*2,'latitude':[-30]*2,'longitude':[-60]*2,'corrente_pico_kA':[1]*2,'ic_cg_flag':[1]*2})
    path=tmp_path/'events.csv'; f.to_csv(path,index=False)
    e,d=read_events([path],{'0':'CG','1':'IC'})
    assert e.time.iloc[0].value%1_000_000_000==123456789
    assert len(e)==1 and len(d)==1 and d.source_row.iloc[0]==3
    f['nanosec']=f.nanosec.astype(float);f.loc[0,'nanosec']=1.5;f.to_csv(path,index=False)
    with pytest.raises(ValueError,match='nanosecond'):read_events([path],{'0':'CG','1':'IC'})

def test_spatial_index_agrees_with_exhaustive_anchor_rule():
    from collections import deque
    rng=np.random.default_rng(33)
    e=events(np.arange(200)*20,latitudes=rng.normal(-30,.1,200))
    e['longitude']=rng.normal(-60,.1,200)
    f,a=group_flashes(e)
    anchors=[];expected=[]
    for row in e.sort_values(['time','event_id']).itertuples():
        candidate=None
        for j,g in enumerate(anchors):
            if row.time.value-g.time.value<=700_000_000 and haversine_km(row.latitude,row.longitude,g.latitude,g.longitude)<=10:
                candidate=j;break
        if candidate is None:candidate=len(anchors);anchors.append(row)
        expected.append(anchors[candidate].event_id)
    assert list(a.flash_id.str.split(':').str[-1])==expected

def test_missing_coverage_and_half_open_bins():
    f,_=group_flashes(events([0]))
    r=aggregate(f,'2018-12-13T23:45Z','2018-12-14T00:30Z',[-1,-1,1,1],[True,False,True])
    assert r['count'].iloc[0]==1 and np.isnan(r['count'].iloc[1]) and r['count'].iloc[2]==0
    assert np.isnan(r.cumulative.iloc[-1])

def prediction_fixture():
    times=pd.date_range('2018-12-13',periods=12,freq='15min',tz='UTC')
    rows=[{'time':t,'method':'IxG','configuration':c,'member':m,'predictor':float(i+1)*mult}
        for i,t in enumerate(times) for c,mult in [('TH_C0',1),('TH_CF',2)] for m in [0,1,2]]
    return pd.DataFrame(rows),pd.DataFrame({'time':times,'count':2*np.arange(1,13.)})

def test_fit_zero_is_unrecoverable():
    scale,status=fit_scale([0,0],[1,2]); assert np.isnan(scale) and status=='identically_zero_unrecoverable'

def test_zero_forecast_distinct_from_missing_and_unidentified_nonzero():
    p,o=prediction_fixture();p['predictor']=0.
    c,h,m,*_=evaluate(p,o,'TH_C0','2018-12-13T01:00Z','2018-12-13T01:30Z')
    assert set(c.status)=={'identically_zero_unrecoverable'}
    assert (h.prediction_A==0).all() and (m.n>0).all()
    p.loc[p.time>=pd.Timestamp('2018-12-13T01:30Z'),'predictor']=1.
    _,h,*_=evaluate(p,o,'TH_C0','2018-12-13T01:00Z','2018-12-13T01:30Z')
    assert h.prediction_A.isna().all()

def test_holdout_no_leakage_fixed_reference_and_members():
    p,o=prediction_fixture()
    args=('TH_C0','2018-12-13T01:00Z','2018-12-13T01:30Z')
    c,h,m,contrast,physical=evaluate(p,o,*args)
    changed=o.copy(); changed.loc[6:,'count']*=100
    c2,*_=evaluate(p,changed,*args)
    np.testing.assert_array_equal(c.scale,c2.scale)
    assert set(c.scale)=={1.,2.}
    assert set(h[h.configuration=='TH_CF'].reference_scale)=={2.}
    assert (h.time>=pd.Timestamp(args[-1])).all()
    assert h.groupby(['time','configuration']).scale.nunique().max()==1
    assert (m[m.comparison=='A'].rmse==0).all()
    assert (physical[physical.configuration=='TH_CF'].physical_change>0).all()
    p.loc[p.member==0,'predictor']*=100
    _,_,_,_,physical=evaluate(p,o,*args)
    ref=physical[physical.configuration=='TH_C0'].merge(o,on='time')
    np.testing.assert_allclose(ref.reference_prediction,ref['count'])

def test_joint_eda_covariance_paired_reproducibility():
    rng=np.random.default_rng(1); t=rng.normal(280,1,(10,4)); r=np.exp((t-280)*.1-5)
    a,w=joint_perturbations(t,r); b,_=joint_perturbations(t,r)
    np.testing.assert_array_equal(a,b)
    np.testing.assert_allclose(a[:,1],.1*a[:,0],atol=1e-14)
    assert a.shape==(20,2,4)

def test_registry_and_experiment_matrix():
    from leona.experiments import CONFIGURATIONS,namelist
    assert len(CONFIGURATIONS)==6
    assert 'use_mynn=.true.' in namelist('TH_C0','input.nc','out.nc')
    assert 'couple_mynn_condensate=.false.' in namelist('TH_C0','input.nc','out.nc')
    assert 'mynn_condensate_local_only=.true.' in namelist('WS_CL','input.nc','out.nc')
    registry=json.loads(Path('configs/parameterizations.json').read_text())
    assert len(registry['methods'])==11

def test_native_bin_integral_and_unavailable_gap():
    from leona.postprocess import integrate_bins
    t=pd.date_range('2018-12-14',periods=31,freq='1min',tz='UTC')
    f=pd.DataFrame({'time':t,'method':'IxG','value':2.})
    r=integrate_bins(f,t[0],t[-1])
    np.testing.assert_allclose(r.predictor,[1800,1800])
    f.loc[4,'value']=np.nan
    r=integrate_bins(f,t[0],t[-1])
    assert np.isnan(r.predictor.iloc[0]) and r.predictor.iloc[1]==1800

def test_native_diagnostics_exact_mass_basis_and_read_only(tmp_path):
    from leona.postprocess import native_diagnostics
    from leona.provenance import sha256
    p=np.array([[100000.,80000.]]*3);t=np.array([[300.,280.]]*3);rv=.01
    ds=xr.Dataset(coords={'time':[0.,900.,1800.],'lev':[250.,750.]})
    for name,value in {'theta':t/(p/1e5)**(2/7),'pressure':p,'qv':np.full((3,2),rv),
        'qi':np.full((3,2),1e-4),'qg':np.full((3,2),2e-4),'mf_w':np.full((3,2),100.)}.items():
        ds[name]=(('time','lev'),value)
    path=tmp_path/'native.nc';ds.to_netcdf(path);before=sha256(path)
    f=native_diagnostics(path,'2018-12-13T00:00Z',top=1000.)
    rho=p[0]/(287*t[0]*(1+rv/(287/461.6)))
    np.testing.assert_allclose(f[f.method=='IxG'].value,np.sum(rho**2*1e-4*2e-4*500))
    assert f[f.method=='PW10'].value.isna().all()
    assert sha256(path)==before

def test_incomplete_member_support_rejected():
    p,o=prediction_fixture()
    with pytest.raises(ValueError,match='Member support'):
        evaluate(p.iloc[1:],o,'TH_C0','2018-12-13T01:00Z','2018-12-13T01:30Z')
