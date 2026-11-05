"""Read-only native-timestep diagnostics and complete-bin integration."""
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from .proxies import diagnose

def native_diagnostics(path, origin, top=20000.):
    with xr.open_dataset(path,decode_times=False) as source:
        ds=source.load()
    z=ds.lev.values
    zw=np.zeros(len(z)+1)
    for k,center in enumerate(z): zw[k+1]=2*center-zw[k]
    if np.any(np.diff(zw)<=0) or not np.isclose(zw[-1],top):
        raise ValueError('Layer geometry inconsistent with driver grid')
    dz=np.diff(zw)
    # Constants from src/mpas_constants.F: cp=7*Rd/2, kappa=2/7.
    t=ds.theta.values*(ds.pressure.values/1.e5)**(2./7.)
    rho_d=ds.pressure.values/(287.*t*(1+ds.qv.values/(287./461.6)))
    times=pd.Timestamp(origin)+pd.to_timedelta(ds.time.values,unit='s')
    rows=[]
    for i,stamp in enumerate(times):
        fields={'rho_d':rho_d[i],'qi':ds.qi.values[i],'qg':ds.qg.values[i],'dz':dz}
        for row in diagnose(fields,{'mass_basis':'dry_air'}):
            rows.append({'time':stamp,**row})
    return pd.DataFrame(rows)

def integrate_bins(native, start, end, minutes=15):
    """Integral of native linear snapshots, refusing gaps and missing endpoints.

    Raw proxy-unit seconds, not labeled as flashes until scale calibration.
    NaNs are unavailable. No silent integration across missing values.
    """
    rows=[]
    edges=pd.date_range(start,end,freq=f'{minutes}min')
    for method,f in native.groupby('method'):
        f=f.sort_values('time')
        times=pd.DatetimeIndex(f.time).as_unit('ns').asi8
        values=f.value.to_numpy()
        if len(times)<2:raise ValueError('At least two native snapshots required')
        if np.any(np.diff(times)<=0):raise ValueError('Duplicate/nonmonotonic time')
        cadence=np.median(np.diff(times))
        for a,b in zip(edges[:-1],edges[1:]):
            left,right=np.searchsorted(times,[a.value,b.value])
            lo=max(0,left-1);hi=min(len(times),right+1)
            ts=times[lo:hi];vs=values[lo:hi]
            available=(len(ts)>=2 and ts[0]<=a.value and ts[-1]>=b.value and
                np.isfinite(vs).all() and np.max(np.diff(ts))<=1.5*cadence)
            value=np.nan
            if available:
                selected=np.r_[a.value,ts[(ts>a.value)&(ts<b.value)],b.value]
                seconds=(selected-a.value)/1e9
                interpolated=np.interp(seconds,(ts-a.value)/1e9,vs)
                value=float(np.trapz(interpolated,seconds))
            rows.append({'time':a,'method':method,'predictor':value,'status':'available' if available else 'unavailable'})
    return pd.DataFrame(rows)
