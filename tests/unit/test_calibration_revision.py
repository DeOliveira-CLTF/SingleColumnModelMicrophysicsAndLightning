import copy
import numpy as np
import pandas as pd
import pytest
from scm_lightning.config import load_config
from scm_lightning.provenance import ROOT
from scm_lightning.calibration.revision import METHODS,fit_controls,apply,normalize,freeze,load_frozen
from scm_lightning.forcing.contracts import area_km2


def inputs():
    c=load_config(ROOT/'configs/cases/argentina_20191028.yaml')
    clocks=[pd.date_range(c['start'],c['end'],freq=f'{n}min',inclusive='left').astype(str) for n in [1,5]]
    p=pd.DataFrame([dict(time=t,method=m,configuration=s+'_C0',member=0,predictor=float((i%17)+1)*(1 if s=='TH' else 2))
        for s in ['TH','WS'] for m in METHODS for i,t in enumerate(clocks[0])])
    p5=pd.DataFrame([dict(time=t,method='McCaulF3',configuration=s+'_C0',member=0,predictor=float((i%13)+1)*(1 if s=='TH' else 2))
        for s in ['TH','WS'] for i,t in enumerate(clocks[1])])
    obs=[pd.DataFrame({'time':clock,'count':np.arange(len(clock))%23,'covered':True}) for clock in clocks]
    return c,p,obs[0],p5,obs[1]


def test_six_coefficients_use_all_valid_intervals_and_mccaul_is_OLS_not_peak():
    c,p,o,p5,o5=inputs();o5.loc[7,'count']=np.nan
    coef,_=fit_controls(c,p,o,p5,o5)
    assert len(coef)==6 and set(coef.fit_method)=={'nonnegative_proportional_OLS_origin_all_valid_intervals'}
    f=normalize(p5,c);f=f[f.configuration=='TH_C0']
    x=f.predictor.to_numpy()*area_km2(c['area']);y=o5['count'].to_numpy();valid=np.isfinite(y)
    actual=coef[(coef.method==METHODS[2])&(coef.microphysics=='TH')].iloc[0]
    assert actual.n_valid==575 and actual.fit_interval_minutes==5
    assert actual.scale==pytest.approx(np.dot(x[valid],y[valid])/np.dot(x[valid],x[valid]))
    assert actual.scale!=pytest.approx(y[valid].max()/x[valid].max())
    th=coef[coef.microphysics=='TH'].set_index('method').scale
    ws=coef[coef.microphysics=='WS'].set_index('method').scale
    np.testing.assert_allclose(ws,th/2)


@pytest.mark.parametrize('where',['observation','predictor','event','member','coupling','dates'])
def test_no_evaluation_event_member_or_coupling_can_calibrate(where):
    c,p,o,p5,o5=inputs()
    if where=='observation':o.loc[0,'time']='2018-12-13T00:00Z'
    if where=='predictor':p.loc[0,'time']='2018-12-13T00:00Z'
    if where=='event':c['event_id']='argentina_20181213'
    if where=='member':p.loc[0,'member']=1
    if where=='coupling':p.loc[0,'configuration']='TH_CL'
    if where=='dates':c['start']='2018-12-13T00:00Z'
    with pytest.raises(ValueError):fit_controls(c,p,o,p5,o5)


def test_each_microphysics_shares_coefficients_across_members_and_couplings(tmp_path):
    c,p,o,p5,o5=inputs();coef,baseline=fit_controls(c,p,o,p5,o5)
    frame=pd.DataFrame([dict(time='2018-12-13T00:00Z',method=m,configuration=s+'_'+mode,member=member,predictor=2.)
        for m in METHODS for s in ['TH','WS'] for mode in ['C0','CL','CF'] for member in range(41)])
    result=apply(frame,coef,123.)
    assert result.groupby(['method','microphysics']).scale.nunique().eq(1).all()
    assert result.groupby('method').scale.nunique().eq(2).all()
    common=apply(frame,coef,123.,'common_TH_same_MUCAPE_OLS')
    assert common.groupby('method').scale.nunique().eq(1).all()
    assert not result.prediction.equals(common.prediction)
    freeze(tmp_path,c,coef,baseline,[]);saved,meta=load_frozen(tmp_path)
    assert len(saved)==6 and not meta['evaluation_observations_used']
    with pytest.raises(FileExistsError):freeze(tmp_path,c,coef,baseline,[])


def test_unidentifiable_and_missing_inputs_are_not_fabricated():
    c,p,o,p5,o5=inputs();p.loc[p.method==METHODS[0],'predictor']=0
    p.loc[p.method==METHODS[1],'predictor']=np.nan
    coef,_=fit_controls(c,p,o,p5,o5)
    assert coef[coef.method==METHODS[0]].status.eq('identically_zero_unrecoverable').all()
    assert coef[coef.method==METHODS[1]].status.eq('no_valid_calibration_data').all()
    assert coef[coef.method!=METHODS[2]].scale.isna().all()
