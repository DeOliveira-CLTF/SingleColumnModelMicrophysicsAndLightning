"""Shared, inspectable definitions for the unmodified SCM experiments."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'resultados'
INPUT = ROOT / 'data_goamazon/GOAMAZON_SINGLE_SCM_driver.nc'
IDS = ['TH_C0', 'TH_CL', 'TH_CF', 'WS_C0', 'WS_CL', 'WS_CF']
SPECIES = ['qc', 'qr', 'qi', 'qs', 'qg']
RD, CP, RV, P0 = 287.0, 1004.5, 461.6, 100000.0
FIXED = dict(use_mynn=True, mynn_edmf=1, mynn_edmf_type=0,
             mynn_uedmf_mp=0, apply_rt_adv=False, nz=60, ztop=20000., dz0=50., dt_microp=30.)

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def parse_nml(path):
    text = Path(path).read_text(encoding='utf-8-sig')
    text = re.sub(r'!.*', '', text)
    out = {}
    for key, val in re.findall(r'(\w+)\s*=\s*("[^"]*"|\x27[^\x27]*\x27|[^,\n/]+)', text):
        val = val.strip()
        if val[0] in '\x27"':
            val = val[1:-1].strip()
        elif val.lower() in ('t', 'f', '.true.', '.false.'):
            val = val.lower() in ('t', '.true.')
        else:
            val = int(val) if re.fullmatch(r'[+-]?\d+',val) else float(re.sub('[dD]', 'e', val))
        out[key.lower()] = val
    return out

def nml_path(eid, directory=RESULTS):
    for p in (directory/'namelists'/f'{eid}.nml', directory/f'{eid}.nml'):
        if p.exists():
            return p
    raise FileNotFoundError(f'Namelist ausente: {eid}')

def expected(eid):
    scheme, coupling = eid.split('_')
    return dict(FIXED, microp_scheme='thompson' if scheme == 'TH' else 'wsm6',
                couple_mynn_condensate=coupling != 'C0', mynn_condensate_local_only=coupling == 'CL')

def config_errors(cfg, eid, reference=None):
    errors = [f'{k}: {cfg.get(k)!r} != {v!r}' for k,v in expected(eid).items() if cfg.get(k) != v]
    if reference:
        exempt = {'path_input', 'path_output', 'microp_scheme', 'couple_mynn_condensate', 'mynn_condensate_local_only'}
        errors += [f'Parâmetro de referência divergente: {k}' for k in set(reference)|set(cfg)
                   if k not in exempt and cfg.get(k) != reference.get(k)]
    if Path(str(cfg.get('path_input',''))).name != INPUT.name:
        errors.append('Arquivo de entrada incompatível')
    declared=str(cfg.get('path_input','')).replace('\\','/')
    allowed={INPUT.as_posix(),'/mnt/'+INPUT.drive[0].lower()+INPUT.as_posix()[2:] if INPUT.drive else INPUT.as_posix(),
             '../data_goamazon/'+INPUT.name,'data_goamazon/'+INPUT.name}
    # Historical namelists remain evidence, not active run configurations.
    # Accept their retired path only when the original input hash matches.
    if '/data/' in declared or declared.startswith('data/'):
        baseline=RESULTS/'proveniencia/execution_20261102T214208919408Z.json'
        if baseline.exists() and sha(INPUT)==json.loads(baseline.read_text(encoding='utf-8'))['input_sha256']:
            allowed.add(declared)
    if declared not in allowed: errors.append('Caminho de entrada não reconhecido; resolver e conferir hash antes de aceitar')
    return errors

def load(path):
    with xr.open_dataset(path, decode_times=False, engine='netcdf4') as ds:
        return ds.load()

def interfaces(z, top):
    # Exact inversion of z[k]=(zw[k]+zw[k+1])/2, zw[0]=0 in the driver.
    zw = np.zeros(len(z)+1)
    for k, center in enumerate(z):
        zw[k+1] = 2*center-zw[k]
    if np.any(np.diff(zw)<=0) or not np.isclose(zw[-1],top,rtol=0,atol=1e-7):
        raise ValueError('Centros incompatíveis com interfaces do driver/topo declarado')
    return zw

def temperature(ds):
    return ds.theta.values*(ds.pressure.values/P0)**(RD/CP)

def thermo(ds):
    t = temperature(ds)
    rho = ds.pressure.values/(RD*t*(1+0.61*ds.qv.values))
    return t, rho

def metrics(ds, cfg):
    t, rho = thermo(ds)
    time = ds.time.values.astype(float)
    dt = np.diff(time)
    if not np.isfinite(time).all() or np.any(dt<=0):
        raise ValueError('Tempo não finito ou não crescente')
    for v in ds.variables.values():
        if np.issubdtype(v.dtype,np.number) and not np.isfinite(v.values).all():
            raise ValueError('Dados não finitos: análise interrompida')
    z = ds.lev.values
    dz = np.diff(interfaces(z,cfg['ztop']))
    mixed = (t>=233.15)&(t<=273.15)
    rainrate = np.r_[np.nan, ds.rainncv.values[1:]*3600/dt]
    flash = ds.flash_rate.values
    active = flash>0
    out = dict(duration_h=(time[-1]-time[0])/3600, last_time_s=time[-1], samples=len(time),
               rain_final_mm=float(ds.rainnc.values[-1]),
               rain_35h_mm=float(np.interp(126000,time,ds.rainnc.values)) if time[0]<=126000<=time[-1] else None,
               rain_rate_max_mm_h=float(np.nanmax(rainrate)),
               rain_peak_h=float(time[np.nanargmax(rainrate)]/3600) if np.nanmax(rainrate)>0 else None,
               rain_onset_h=float(time[np.flatnonzero(ds.rainncv.values>0)[0]]/3600) if (ds.rainncv.values>0).any() else None,
               flash_active=bool(active.any()), flash_max_min=float(flash.max()),
               flash_integral_nominal=float(np.sum(flash[1:]*dt/60)),
               flash_onset_h=float(time[np.flatnonzero(active)[0]]/3600) if active.any() else None,
               flash_peak_h=float(time[np.argmax(flash)]/3600) if active.any() else None,
               flash_last_h=float(time[np.flatnonzero(active)[-1]]/3600) if active.any() else None,
               flash_active_duration_h=float(np.sum(dt[active[1:]])/3600),
               flash_end_censored=bool(active[-1]))
    series = {'time_h':time/3600,'rain_rate_mm_h':rainrate,'rain_mm':ds.rainnc.values,'flash_min':flash}
    for name in SPECIES:
        a = ds[name].values
        imax = np.unravel_index(np.argmax(a),a.shape)
        b = np.where(mixed,a,0)
        mass = (rho*b*dz).sum(axis=1)
        out.update({f'{name}_max_g_kg':float(a.max()*1000),
                    f'{name}_peak_h':float(time[imax[0]]/3600) if a.max()>0 else None,
                    f'{name}_peak_z_m':float(z[imax[1]]) if a.max()>0 else None,
                    f'{name}_mixed_max_g_kg':float(b.max()*1000),
                    f'{name}_mixed_mass_max_kg_m2':float(mass.max()),
                    f'{name}_mixed_mass_peak_h':float(time[np.argmax(mass)]/3600) if mass.max()>0 else None})
        series[f'{name}_mixed_kg_m2'] = mass
        series[f'{name}_mixed_max_g_kg'] = b.max(axis=1)*1000
    def correlation(a,b):
        return float(np.corrcoef(a,b)[0,1]) if np.ptp(a)>0 and np.ptp(b)>0 else None
    out['corr_flash_qg_mixed_mass'] = correlation(flash,series['qg_mixed_kg_m2'])
    out['corr_flash_qg_mixed_max'] = correlation(flash,series['qg_mixed_max_g_kg'])
    for a,b in [('flash_peak_h','qg_mixed_mass_peak_h'),('rain_peak_h','flash_peak_h'),('flash_peak_h','qg_peak_h')]:
        out[f'lag_{a}_minus_{b}'] = out[a]-out[b] if out[a] is not None and out[b] is not None else None
    coexist = mixed & (ds.qc.values>0)&(ds.qs.values>0)&(ds.qg.values>0)
    out['coexist_qc_qs_qg_same_cell_duration_h'] = float(np.sum(dt[np.any(coexist,axis=1)[1:]])/3600)
    out['qg_mixed_positive_duration_h'] = float(np.sum(dt[(series['qg_mixed_kg_m2']>0)[1:]])/3600)
    out['flash_with_same_cell_coexist_fraction'] = float(np.sum(dt[(active & coexist.any(axis=1))[1:]])/np.sum(dt[active[1:]])) if active[1:].any() else None
    return out,series

def save_json(path, obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False,default=str)+'\n',encoding='utf-8')
