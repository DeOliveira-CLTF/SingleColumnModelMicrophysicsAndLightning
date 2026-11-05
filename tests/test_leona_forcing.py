"""Integration with actual audited dephycf, using analytic fields only."""
import copy
import json
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from leona.forcing import build,load_dephy

def test_official_dephy_roundtrip_and_full_time_forcing(tmp_path):
    c=json.loads(Path('configs/leona.json').read_text())
    c.update(start='2018-12-13T00:00:00Z',end='2018-12-13T01:00:00Z',area=[-63,-32,-61,-30])
    c['area_selection']['status']='fixed_before_lightning_evaluation'
    coords={'time':pd.date_range('2018-12-12T23:00',periods=5,freq='1h'),
        'level':[950,850,700,500],'latitude':[-33.,-32.,-31.,-30.,-29.],'longitude':[-64.,-63.,-62.,-61.,-60.]}
    shape=(5,4,5,5)
    pl=xr.Dataset(coords=coords)
    for name,profile in {'t':[290,280,270,250],'q':[.01,.008,.004,.001],
        'z':np.array([500,1500,3000,5500])*9.80665,'u':[5,5,5,5],'v':[0,0,0,0],'w':[0,0,0,0]}.items():
        pl[name]=(('time','level','latitude','longitude'),np.broadcast_to(np.array(profile)[None,:,None,None],shape))
    sl=xr.Dataset(coords={k:v for k,v in coords.items() if k!='level'})
    for name,value in {'sp':100000.,'skt':300.,'t2m':295.,'d2m':285.,'u10':5.,'v10':0.,'z':0.,'avg_ishf':-10.,'avg_slhtf':-20.}.items():
        sl[name]=(('time','latitude','longitude'),np.ones((5,5,5))*value)
    pl.to_netcdf(tmp_path/'pl.nc');sl.to_netcdf(tmp_path/'sl.nc')
    output=build(c,tmp_path/'pl.nc',tmp_path/'sl.nc',tmp_path/'output',np.array([0.,500.,1000.,2000.]))
    with xr.open_dataset(output,decode_times=False) as scm:
        assert scm.time.values[0]==0 and scm.time.values[-1]==3600
        assert '2018-12-13' in scm.time.attrs['units']
        assert {'tntheta_adv','tnrv_adv','pa_forc','theta','rv'}<=set(scm.variables)
        for name in ['theta','rv','ta','pa','ua','va']:
            assert scm[name].dims==('t0','lev')
        for name in ['tntheta_adv','tnrv_adv','pa_forc']:
            assert scm[name].dims==('time','lev')
        np.testing.assert_allclose(scm.tntheta_adv,0,atol=1e-12)
        np.testing.assert_allclose(scm.tnrv_adv,0,atol=1e-12)
        assert np.isfinite(scm.pa_forc).all()
        np.testing.assert_allclose(scm.pa.values.squeeze(),[100000,95000,90000,80000],atol=.02)
    manifest=json.loads((tmp_path/'output/forcing_manifest.json').read_text())
    assert manifest['dephy']['commit']==c['dephy_expected_commit']
    from leona.experiments import perturb_scm
    member=tmp_path/'member.nc'
    perturb_scm(output,member,np.array([[.2,.2,.2,.2],[.01,.01,.01,.01]]))
    with xr.open_dataset(output,decode_times=False) as control, xr.open_dataset(member,decode_times=False) as changed:
        # Official DEPHY writer stores float32 temperature near 300 K.
        assert changed.ta.dims==control.ta.dims
        np.testing.assert_allclose(changed.ta.values-control.ta.values,.2,atol=3e-5,rtol=0)
        np.testing.assert_allclose(changed.rv.values/control.rv.values,np.exp(.01))
        np.testing.assert_array_equal(changed.pa,control.pa)
        np.testing.assert_array_equal(changed.pa_forc,control.pa_forc)
    from leona.ensemble import collocate
    edapl=pl.expand_dims(number=np.arange(10)).copy()
    edasl=sl.expand_dims(number=np.arange(10)).copy()
    edapl['t']=edapl.t+xr.DataArray(np.arange(10)*.1,dims='number',coords={'number':np.arange(10)})
    edapl.to_netcdf(tmp_path/'edapl.nc');edasl.to_netcdf(tmp_path/'edasl.nc')
    collocate(c,output,tmp_path/'edapl.nc',tmp_path/'edasl.nc',tmp_path/'eda_profiles.npz')
    eda=np.load(tmp_path/'eda_profiles.npz')
    assert eda['temperature'].shape==(10,4)
    np.testing.assert_allclose(eda['temperature'][1,1:]-eda['temperature'][0,1:],.1)
