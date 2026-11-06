import numpy as np
import pytest
from scm_lightning.ensemble.updraft_sensitivity import pulse
from scm_lightning.evaluation.updraft_sensitivity import raw_proxies


def test_historical_pulse_support_peak_and_integral():
    z=np.array([2000.,3000.,4000.]);t=np.array([99.,100.,600.,1100.,1101.])
    w=pulse(z,t,100.,1000.,3000.,1000.,2.)
    np.testing.assert_allclose(w[[0,1,3,4]],0,atol=1e-15)
    np.testing.assert_allclose(w[2],[2*np.exp(-.5),2.,2*np.exp(-.5)])
    timeline=np.linspace(100.,1100.,1001)
    integral=np.trapz(pulse([3000.],timeline,100.,1000.,3000.,1000.,2.)[:,0],timeline)
    assert integral==pytest.approx(1000.)
    np.testing.assert_array_equal(pulse(z,t,100.,1000.,3000.,1000.,0.),np.zeros((5,3)))
    with pytest.raises(ValueError):pulse(z,t,100.,0.,3000.,1000.,1.)


def state():
    shape=(2,4)
    return {'temperature':np.array([[280.,270.,255.,240.]]*2),'pressure':np.array([[100000.,70000.,44000.,20000.]]*2),
            'z':np.array([0.,2000.,5000.,8000.]),'dz':np.array([1000.,2000.,3000.,3000.]),
            'rho_d':np.ones(shape),'P':np.array([.001,.002]),
            **{q:np.full(shape,.001) for q in ['qc','qr','qi','qs','qg']}}


def test_velocity_changes_do_not_modify_state_or_ice_product():
    s=state();before=s['qg'].copy()
    zero={'w':np.zeros((2,4)),'z':s['z'],'active':None}
    fast={'w':np.full((2,4),11.),'z':s['z'],'active':None}
    low,high=raw_proxies(s,zero),raw_proxies(s,fast)
    np.testing.assert_array_equal(s['qg'],before)
    np.testing.assert_array_equal(low['IxG'],high['IxG'])
    np.testing.assert_array_equal(low['McCaul_F2'],high['McCaul_F2'])
    np.testing.assert_array_equal(low['PW10_1D'],[0.,0.])
    np.testing.assert_array_equal(high['PW10_1D'],s['P'])
    assert np.all(high['IFluxT_1D']>low['IFluxT_1D'])
    exact={'w':np.full((2,4),10.),'z':s['z'],'active':None}
    np.testing.assert_array_equal(raw_proxies(s,exact)['PW10_1D'],[0.,0.])


def test_zero_hydrometeors_are_valid_zero_predictors():
    s=state()
    for q in ['qc','qr','qi','qs','qg']:s[q][:]=0.
    values=raw_proxies(s,{'w':np.full((2,4),2.),'z':s['z'],'active':None})
    for name in ['IxG','VII','McCaul_F2','IFluxT_1D','McCaul_F1_1D','McCaul_F3_1D','LPI_1D']:
        np.testing.assert_array_equal(values[name],[0.,0.])
    assert np.all(values['PR92_wmax_1D']>0)
