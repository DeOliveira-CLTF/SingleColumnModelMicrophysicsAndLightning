#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GOAMAZON/SINGLE - definicao original do caso a partir do ERA5.
Substitui sing_sounding_278_5.txt, sing_forcing.nc e sing_sfc_flux_278_5.txt.
Gera GOAMAZON_SINGLE_DEF_driver.nc (e era5_ts.npz com a Ts para o driver_SCM.py).
"""
import numpy as np
import xarray as xr
from dephycf.Case import Case

g, Rd, cp, Re = 9.80665, 287.04, 1004., 6.371e6
kappa = Rd / cp
lat0, lon0 = -8., -63.
t0 = np.datetime64('2014-10-05T12:00')
t1 = np.datetime64('2014-10-06T00:00')
lplot = True

# ---- opcoes -------------------------------------------------------------
DT_MIN = 10              # passo de tempo das forcantes (min); 10 -> 72 passos, 60 -> 12 passos
INTERP_FLUX = 'linear'   # interpolacao temporal dos fluxos: 'linear' ou 'cubic'
# --------------------------------------------------------------------------


def load(fn, tend):
    ds = xr.open_dataset(fn)
    ren = {k: v for k, v in {'valid_time': 'time', 'pressure_level': 'level'}.items()
           if k in ds.dims or k in ds.coords}
    ds = ds.rename(ren).sortby('latitude').sortby('longitude')
    return ds.sel(time=slice(t0, tend))


# sl vai ate 01 UTC: msshf/mslhf (media da hora anterior) sao centrados em t-30 min
pl = load('era5_pl.nc', t1)
sl = load('era5_sl.nc', t1 + np.timedelta64(1, 'h'))
box = ('latitude', 'longitude')

# ---------------------------------------------------------------- grade temporal
# 12 h a partir de t0, com passo DT_MIN (10 min -> 72 passos: 12:00 ... 23:50 UTC)
nsteps = int(12 * 60 / DT_MIN)
tnew = t0 + np.arange(nsteps) * np.timedelta64(DT_MIN * 60, 's')
tsec = ((tnew - t0) / np.timedelta64(1, 's')).astype(np.float64)

# fluxos: msshf/mslhf ja sao taxa media em W/m2 (ECMWF param table 235),
# associada ao ponto medio da hora anterior -- ao contrario de sshf/slhf
# (acumulados em J/m2, param table 128), nao precisam dividir por 3600.
# Sinal: ERA5 usa positivo para baixo; convencao DEPHY/CMOR (hfss/hfls) e
# positivo para cima -- daqui o sinal trocado.
#PK flux = xr.Dataset({'hfss': -sl.msshf, 'hfls': -sl.mslhf})
flux = xr.Dataset({'hfss': -sl.avg_ishf, 'hfls': -sl.avg_slhtf})
flux = flux.assign_coords(time=flux.time - np.timedelta64(30, 'm'))
flux = flux.interp(time=tnew, method=INTERP_FLUX,
                   kwargs={'fill_value': 'extrapolate'})

# demais campos (instantaneos): interpolacao linear no tempo
pl = pl.interp(time=tnew)
sl = sl[['sp', 'skt', 'z']].interp(time=tnew)

# ---------------------------------------------------------------- campos derivados
P = pl.level * 100.                                   # Pa
theta = pl.t * (1.e5 / P) ** kappa
coslat = np.cos(np.deg2rad(pl.latitude))
fac = 180. / np.pi


def hgrad(f):
    dfdx = f.differentiate('longitude') * fac / (Re * coslat)
    dfdy = f.differentiate('latitude') * fac / Re
    return dfdx, dfdy


dthx, dthy = hgrad(theta)
dqx, dqy = hgrad(pl.q)
adv_theta = -(pl.u * dthx + pl.v * dthy)              # K/s
adv_q = -(pl.u * dqx + pl.v * dqy)                    # kg/kg/s

rho = P / (Rd * pl.t * (1. + 0.608 * pl.q))
w = -pl.w / (rho * g)                                 # omega (Pa/s) -> w (m/s)

# altura acima do solo (m)
zs = (sl.z / g).isel(time=0).mean(box)
zagl = (pl.z / g).mean(box) - zs                      # (time, level)

# medias na caixa
m = lambda f: f.mean(box)
T_adv, q_adv, W = m(adv_theta), m(adv_q), m(w)
U, V = m(pl.u), m(pl.v)
TH, Q = m(theta), m(pl.q)
ps = float(m(sl.sp).isel(time=0))

# niveis validos (acima do solo em todos os tempos), ordenados por altitude crescente
hz = zagl.mean('time').values
valid = np.where((zagl.min('time').values > 0.))[0]
idx = valid[np.argsort(hz[valid])]
z = hz[idx].astype(np.float32)


def sel(f, t=None):
    a = f.transpose('time', 'level').values[:, idx] if t is None else f.isel(time=t).values[idx]
    return a.astype(np.float32)


################################################
# Case
################################################
case = Case('GOAMAZON/SINGLE', lat=lat0, lon=lon0,
            startDate="20141005120000", endDate="20141006000000",
            surfaceType='land', zorog=0.)
case.set_title("Forcing and initial conditions for the single pulse case - ERA5-based definition")
case.set_reference("ERA5 reanalysis (Hersbach et al. 2020), horizontal average over a ~2 deg box")
case.set_author("P. Kubota")
case.set_script("driver_DEF_era5.py")

# estado inicial (primeiro instante)
case.add_init_ps(ps)
case.add_init_wind(u=sel(U, 0), v=sel(V, 0), lev=z, levtype='altitude')
case.add_init_theta(sel(TH, 0), lev=z, levtype='altitude')
case.add_init_qv(sel(Q, 0), lev=z, levtype='altitude')

# forcantes (mesma altura para todos os tempos, como no original)
case.add_theta_advection(sel(T_adv), include_rad=False, lev=z, levtype='altitude', time=tsec)
case.add_rv_advection(sel(q_adv), lev=z, levtype='altitude', time=tsec)
case.add_vertical_velocity(w=sel(W), lev=z, levtype='altitude', time=tsec)
case.add_wind_nudging(unudg=sel(U), vnudg=sel(V), timescale=3600.,
                      time=tsec, lev=z, levtype='altitude')

# fluxos de superficie: ERA5 acumula em J/m2 na hora anterior, positivo para baixo
# -> W/m2 positivo para cima (convencao DEPHY)
hfss = flux.hfss.mean(box).values.astype(np.float32)
hfls = flux.hfls.mean(box).values.astype(np.float32)
case.add_surface_fluxes(hfss, hfls, time=tsec, forc_wind='z0', z0=0.035)

# Ts (skin temperature) para o driver_SCM.py
np.savez('era5_ts.npz', time=tsec, ts=m(sl.skt).values.astype(np.float32))

case.write('GOAMAZON_SINGLE_DEF_driver.nc')

if lplot:
    case.plot(rep_images='./images/driver_DEF/', timeunits='hours')
