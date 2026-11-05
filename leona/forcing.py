"""ERA5 regional fields -> official dephycf.Case -> convert2SCM.

Never writes a substitute DEPHY schema. Rejects missing support instead of
extrapolating temperature/moisture through terrain or beyond available times.
"""
import inspect
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from .contracts import horizontal_advection, KAPPA
from .provenance import ROOT,git,sha256,write_json

def load_dephy():
    path = ROOT/'era5_forcing/DEPHY-SCM'
    sys.path.insert(0,str(path))
    from dephycf.Case import Case
    origin = Path(inspect.getfile(Case)).resolve()
    return Case,{'imported_path':str(origin),'sha256':sha256(origin),'commit':git('rev-parse','HEAD',cwd=path),'local_diff':git('diff','--stat',cwd=path)}

def normalize(ds):
    ren = {a:b for a,b in [('valid_time','time'),('pressure_level','level')] if a in ds.dims}
    ds = ds.rename(ren).sortby('time').sortby('latitude').sortby('longitude')
    if 'expver' in ds.dims:
        raise ValueError('Resolve ERA5 expver explicitly before forcing conversion')
    return ds

def common_grid(case,z):
    """Retain the single SCM height axis after official add_* calls."""
    from dephycf.Axis import Axis
    shared=Axis('lev',np.asarray(z),name='height',units='m')
    visited=set()
    def normalize_variable(variable,initial):
        if id(variable) in visited:return
        visited.add(id(variable))
        tid='t0' if initial else 'time'
        if variable.time is not None:variable.time.id=tid
        if variable.level is not None:
            if variable.data.shape[-1] != len(z):
                raise ValueError(f'Inconsistent SCM grid for {variable.id}')
            variable.set_level(lev=shared)
            variable.set_coordinates(tid,'zh' if initial else 'zh_forc','lat','lon')
        elif variable.time is not None:
            variable.set_coordinates(tid,'lat','lon')
        # This audited Variable caches axis IDs in axlist; renaming the Axis
        # alone leaves stale time_<variable> dimensions in nested coordinates.
        variable.axlist=([tid] if variable.time is not None else [])+(['lev'] if variable.level is not None else [])
        if variable.height is not None:
            variable.height.id='zh' if initial else 'zh_forc'
            normalize_variable(variable.height,initial)
        if variable.pressure is not None:
            variable.pressure.id='pa' if initial else 'pa_forc'
            normalize_variable(variable.pressure,initial)
    for name,variable in case.variables.items():
        normalize_variable(variable,name in case.var_init_list)

def build(config, pl_path, sl_path, output, zgrid):
    if config['area'] is None or config['area_selection']['status'] != 'fixed_before_lightning_evaluation':
        raise ValueError('Meteorological area selection must be recorded before forcing/calibration')
    Case,version = load_dephy()
    if version['commit'] != config['dephy_expected_commit']:
        raise ValueError('DEPHY commit differs from configured audited version')
    pl,sl = normalize(xr.load_dataset(pl_path)),normalize(xr.load_dataset(sl_path))
    start,end = pd.Timestamp(config['start']).tz_localize(None),pd.Timestamp(config['end']).tz_localize(None)
    times = pd.date_range(start,end,freq='10min') # both endpoints: no frozen final forcing
    w,s,e,n = config['area']
    # Gradients calculated on regional grid, before selecting fixed area.
    if pl.longitude.min() >= w or pl.longitude.max() <= e or pl.latitude.min() >= s or pl.latitude.max() <= n:
        raise ValueError('Regional pressure fields need gradient halo beyond the area')
    for ds in [pl,sl]:
        if pd.Timestamp(ds.time.values[0]) > start or pd.Timestamp(ds.time.values[-1]) < end:
            raise ValueError('Full UTC time coverage required')
    p = pl.level*100.
    # Preserve xarray labels for differentiating pressure-surface quantities.
    theta = pl.t*(1.e5/p)**KAPPA
    rv = pl.q/(1-pl.q)
    advth,advr = horizontal_advection(theta,pl.u,pl.v),horizontal_advection(rv,pl.u,pl.v)
    # Dry-air mass basis, moist equation of state.
    rho_d = p/(287.04*pl.t*(1+rv/(287.04/461.5)))
    rho = rho_d*(1+rv)
    fields = xr.Dataset({'theta':theta,'rv':rv,'u':pl.u,'v':pl.v,
        'w':-pl.w/(rho*9.80665),'advth':advth,'advr':advr,'height':pl.z/9.80665,'pressure':p.broadcast_like(pl.t)})
    fields = fields.sel(longitude=slice(w,e),latitude=slice(s,n)).interp(time=times)
    surf = sl.sel(longitude=slice(w,e),latitude=slice(s,n)).interp(time=times)
    ps = surf.sp
    above = (fields.pressure < ps) & (fields.height > surf.z/9.80665)
    # Same pressure levels must be above terrain in every selected grid cell/time.
    valid = above.all(('time','latitude','longitude'))
    fields = fields.sel(level=fields.level[valid])
    weights = np.cos(np.deg2rad(fields.latitude))
    f = fields.weighted(weights).mean(('latitude','longitude'))
    sf = surf.weighted(np.cos(np.deg2rad(surf.latitude))).mean(('latitude','longitude'))
    orog = float(sf.z.mean()/9.80665)
    z = np.asarray(zgrid,float)
    if np.any(np.diff(z)<=0) or z[0] != 0:
        raise ValueError('SCM output grid must increase from 0 m AGL')
    # Surface is constrained by actual ERA5 2-m/10-m fields. Horizontal surface
    # gradients use the same regional-grid contract, not vertically advected twice.
    for name in ['t2m','d2m','u10','v10','skt']:
        if name not in sf:
            raise ValueError(f'Required surface field: {name}')
    es = 611.2*np.exp(17.67*(sl.d2m-273.15)/(sl.d2m-29.65))
    r_surface = (287.04/461.5)*es/(sl.sp-es)
    th_surface = sl.t2m*(1.e5/sl.sp)**KAPPA
    def surface_mean(x):
        return x.sel(longitude=slice(w,e),latitude=slice(s,n)).interp(time=times).weighted(np.cos(np.deg2rad(surf.latitude))).mean(('latitude','longitude')).values
    base = {'theta':surface_mean(th_surface),'rv':surface_mean(r_surface),'u':sf.u10.values,'v':sf.v10.values,
        'pressure':sf.sp.values,'w':np.zeros(len(times)),
        'advth':surface_mean(horizontal_advection(th_surface,sl.u10,sl.v10)),
        'advr':surface_mean(horizontal_advection(r_surface,sl.u10,sl.v10))}
    profiles = {name:[] for name in base}
    for it in range(len(times)):
        hz = f.height.isel(time=it).values-orog
        ix = np.argsort(hz)
        hz = hz[ix]
        if hz[0] <= 0 or hz[-1] < z[-1]:
            raise ValueError('Terrain/top outside native profile; no extrapolation allowed')
        for name in profiles:
            values = f[name].isel(time=it).values[ix]
            profiles[name].append(np.interp(z,np.r_[0,hz],np.r_[base[name][it],values]))
    profiles = {k:np.asarray(v) for k,v in profiles.items()}
    if any(not np.isfinite(v).all() for v in profiles.values()):
        raise ValueError('Nonfinite forcing after interpolation')
    seconds = (times-times[0]).total_seconds().to_numpy()
    case = Case(config['case'],lat=(s+n)/2,lon=(w+e)/2,startDate=start.strftime('%Y%m%d%H%M%S'),endDate=end.strftime('%Y%m%d%H%M%S'),surfaceType='land',zorog=orog)
    kw = {'lev':z,'levtype':'altitude'}
    case.add_init_ps(float(sf.sp.values[0]))
    # Supply complete initial aliases on explicit p/z coordinates. In this
    # audited DEPHY version, incomplete height-only state triggers hydrostatic
    # pressure reconstruction, overriding independently supplied ERA5 pressure.
    case.add_init_pressure(profiles['pressure'][0],**kw)
    pvar=case.variables['pa']
    ikw=dict(kw,pressure=pvar)
    case.add_init_height(z,**ikw)
    case.add_init_theta(profiles['theta'][0],**ikw)
    case.add_init_rv(profiles['rv'][0],**ikw)
    case.add_init_wind(u=profiles['u'][0],v=profiles['v'][0],**ikw)
    case.add_init_variable('ta',case.compute_temp(pressure=pvar),**ikw)
    case.add_init_variable('thetal',profiles['theta'][0],**ikw)
    for name in ['rl','ri','ql','qi','tke']:
        case.add_init_variable(name,np.zeros_like(z),**ikw)
    for name in ['rt']:
        case.add_init_variable(name,profiles['rv'][0],**ikw)
    for name in ['qv','qt']:
        case.add_init_variable(name,profiles['rv'][0]/(1+profiles['rv'][0]),**ikw)
    case.add_init_variable('hur',case.compute_hur(),**ikw)
    fkw = dict(kw,time=seconds)
    case.add_pressure_forcing(profiles['pressure'],**fkw)
    case.add_surface_pressure_forcing(sf.sp.values,time=seconds)
    case.add_theta_advection(profiles['advth'],include_rad=False,**fkw)
    case.add_rv_advection(profiles['advr'],**fkw)
    case.add_vertical_velocity(w=profiles['w'],**fkw)
    case.add_wind_nudging(unudg=profiles['u'],vnudg=profiles['v'],timescale=3600.,**fkw)
    # Mean fluxes are previous-hour means: shift to interval midpoint. No extrapolation.
    flux_names = [('avg_ishf','avg_slhtf'),('msshf','mslhf')]
    pair = next((a for a in flux_names if set(a)<=set(sl)),None)
    if pair is None:
        raise ValueError('Require mean surface fluxes (W/m2), not accumulated J/m2')
    flux = sl[list(pair)].assign_coords(time=sl.time-np.timedelta64(30,'m')).interp(time=times)
    sens,lat = [surface_mean(-flux[name]) for name in pair]
    if not np.isfinite(sens).all() or not np.isfinite(lat).all():
        raise ValueError('Missing interval-centered flux support')
    case.add_surface_fluxes(sens,lat,time=seconds,forc_wind='z0',z0=.035)
    output = Path(output); output.mkdir(parents=True,exist_ok=True)
    case.write(str(output/'LEONA_DEF.nc'))
    scm = case.convert2SCM(time=seconds,lev=z,levtype='altitude')
    # convert2SCM's altitude interpolation drops pressure-coordinate links and
    # its completion step reconstructs pa/pa_forc hydrostatically. Restore the
    # supplied ERA5 pressure through the official Case API, then relink all
    # variables. No library edit or substitute NetCDF writer is involved.
    scm.add_init_pressure(profiles['pressure'][0],**kw)
    scm.add_pressure_forcing(profiles['pressure'],time=seconds,**kw)
    for name in scm.var_init_list:
        variable=scm.variables[name]
        if variable.level is not None:
            variable.pressure=scm.variables['pa']
    for name in scm.var_forcing_list:
        variable=scm.variables[name]
        if variable.level is not None:
            variable.pressure=scm.variables['pa_forc']
    # Temperature aliases and tendencies follow the restored pressure.
    initial=dict(kw,pressure=scm.variables['pa'])
    scm.add_init_temp(scm.compute_temp(pressure=scm.variables['pa']),**initial)
    if 'tnta_adv' in scm.variables:
        scm.add_temp_advection(scm.compute_tnta_adv(),include_rad=False,time=seconds,**kw)
    scm.add_surface_temp(sf.skt.values,time=seconds,timeid='time')
    common_grid(scm,z)
    scm.write(str(output/'LEONA_SCM.nc'))
    write_json(output/'forcing_manifest.json',{'dephy':version,'pl_sha256':sha256(pl_path),'sl_sha256':sha256(sl_path),
        'configuration':config,'origin_utc':str(start),'surface_policy':'2m thermodynamics and 10m wind mapped to 0m AGL; w=0 at surface',
        'pressure_policy':'ERA5 prescribed pressure restored via Case API after convert2SCM hydrostatic completion; no exact hydrostatic claim','top_m':float(z[-1]),
        'advection':'horizontal regional gradients only; driver adds vertical transport',
        'SCM_sha256':sha256(output/'LEONA_SCM.nc')})
    return output/'LEONA_SCM.nc'
