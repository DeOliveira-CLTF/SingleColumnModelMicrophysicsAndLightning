"""Weather-only event selection and the existing official DEPHY forcing builder."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from scm_lightning.config import load_config
from scm_lightning.forcing import normalize, build
from scm_lightning.forcing.spatial import subset
from scm_lightning.provenance import ROOT, sha256, write_json


def weather_support(config_path):
    c = load_config(config_path)
    if c['area'] is not None:
        raise ValueError('Area already frozen; never reselect using lightning errors')
    cache = ROOT / c['paths']['era5']
    days = pd.date_range(c['start'], c['end'], freq='D', inclusive='left')
    paths = [cache / f'reanalysis_{kind}_{day:%Y%m%d}.nc' for kind in ['pl', 'sl'] for day in days]
    datasets = []
    for kind in ['pl', 'sl']:
        datasets.append(xr.concat([normalize(xr.load_dataset(p)) for p in paths if f'_{kind}_' in p.name], dim='time').sortby('time'))
    pl, sl = datasets
    expected = pd.date_range(c['start'], c['end'], freq='h', inclusive='left').tz_localize(None)
    for ds in datasets:
        np.testing.assert_array_equal(ds.time.values, expected.values)
    if 'tp' not in sl:
        raise ValueError('Independent precipitation field required before area selection')
    # Negative omega is ascent. Hourly regional fields, no ENTLN enters this function.
    ascent = (-pl.w.sel(level=[500, 700])).clip(min=0).mean('level')
    rain = sl.tp.clip(min=0)
    rain = rain.interp(latitude=ascent.latitude, longitude=ascent.longitude)
    mean_ascent, mean_rain = ascent.mean('time'), rain.mean('time')
    score = mean_ascent * mean_rain
    west, south, east, north = c['meteorological_selection']['search_region_wsen']
    half = c['meteorological_selection']['box_width_degrees'] / 2
    candidates = score.sel(longitude=slice(west+half, east-half), latitude=slice(south+half, north-half))
    if not np.isfinite(candidates).all() or float(candidates.max()) <= 0:
        raise ValueError('No identifiable joint ascent/precipitation support')
    # Area-mean ranking, rather than selecting a single anomalous pixel.
    ranked = []
    for lat in candidates.latitude.values:
        for lon in candidates.longitude.values:
            box = [float(lon-half), float(lat-half), float(lon+half), float(lat+half)]
            selected, weights = subset(score.to_dataset(name='score'), box)
            ranked.append((float(selected.score.weighted(weights).mean()), float(lat), float(lon)))
    value, lat, lon = sorted(ranked, key=lambda r: (-r[0], r[1], r[2]))[0]
    area = [lon-half, lat-half, lon+half, lat+half]
    out = ROOT / c['paths']['results'] / 'support'
    out.mkdir(parents=True, exist_ok=True)
    xr.Dataset({'positive_ascent': ascent, 'precipitation_m': rain, 'selection_score': score}).to_netcdf(out/'meteorology.nc')
    tracks = []
    for i, stamp in enumerate(expected):
        joint = ascent.isel(time=i)*rain.isel(time=i)
        region = joint.sel(longitude=slice(west, east), latitude=slice(south, north))
        if float(region.max()) > 0:
            iy, ix = np.unravel_index(np.argmax(region.values), region.shape)
            tracks.append({'time':str(stamp), 'latitude':float(region.latitude[iy]), 'longitude':float(region.longitude[ix]),
                           'max_joint_score':float(region.values[iy, ix]), 'regional_rain_max_mm':float(1000*rain.isel(time=i).max())})
    pd.DataFrame(tracks).to_csv(out/'hourly_meteorological_evolution.csv', index=False)
    write_json(out/'selection.json', {'area_wsen':area, 'criterion':'area mean of time-mean positive omega ascent times time-mean hourly ERA5 precipitation',
        'omega_levels_hPa':[500,700], 'score':value, 'period':[c['start'],c['end']], 'search_region_wsen':[west,south,east,north],
        'lightning_used':False, 'evaluation_2018_used':False, 'inputs':[{'path':str(p), 'sha256':sha256(p)} for p in paths],
        'track_is':'hourly joint-field maximum; not a tracked individual convective cell',
        'boundary_review':'hourly evolution retained; extend event explicitly if meteorological evidence requires it'})
    c['area'] = area
    c['area_selection'].update(status='fixed_before_lightning_evaluation', decision_sha256=sha256(out/'selection.json'))
    write_json(config_path, c)
    return area


def assemble_and_build(config_path):
    c = load_config(config_path)
    days = pd.date_range(c['acquisition']['start'], c['acquisition']['end'], freq='D')
    out = ROOT/'data/intermediate/era5'/c['event_id']
    out.mkdir(parents=True, exist_ok=True)
    for kind in ['pl', 'sl']:
        required=pd.date_range(c['start'],c['end'],freq='D')
        for day in required:
            if not (ROOT/c['paths']['era5']/f'reanalysis_{kind}_{day:%Y%m%d}.nc').exists():
                raise FileNotFoundError(f'Missing required {kind} field {day:%Y%m%d}')
        paths = [p for day in days if (p:=ROOT/c['paths']['era5']/f'reanalysis_{kind}_{day:%Y%m%d}.nc').exists()]
        target = out/f'reanalysis_{kind}.nc'
        if not target.exists():
            w,s,e,n=c['area']
            parts=[]
            for p in paths:
                with xr.open_dataset(p) as raw:
                    parts.append(normalize(raw).sel(latitude=slice(s-2,n+2),longitude=slice(w-2,e+2)).load())
            ds = xr.concat(parts, dim='time').sortby('time')
            if ds.time.to_index().duplicated().any():
                raise ValueError('Duplicate ERA5 hours')
            ds.to_netcdf(target)
            write_json(target.with_suffix('.json'), {'inputs':[{'path':str(p),'sha256':sha256(p)} for p in paths], 'sha256':sha256(target),
                'gradient_halo_degrees':2,'native_regional_grid_retained':True})
    return build(c, out/'reanalysis_pl.nc', out/'reanalysis_sl.nc', ROOT/c['paths']['forcing'], np.arange(0,20001,10))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['select','build']);ap.add_argument('--config',required=True)
    a=ap.parse_args()
    print(weather_support(a.config) if a.stage=='select' else assemble_and_build(a.config))


if __name__=='__main__':main()
