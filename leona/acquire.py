"""Separate cached acquisitions, checksums and explicit availability states."""
from pathlib import Path
from io import StringIO
import json
import time
import pandas as pd
import requests
from .provenance import verified_copy, write_json, sha256

def cached_get(url, params, target):
    target = Path(target)
    meta = target.with_suffix(target.suffix + '.json')
    if target.exists() and meta.exists():
        info = json.loads(meta.read_text())
        if sha256(target) != info['sha256']:
            raise IOError(f'Corrupt cache: {target}')
        return target.read_text(encoding='utf-8')
    error = None
    for attempt in range(3):
        try:
            r = requests.get(url, params=params, timeout=40)
            r.raise_for_status()
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(r.text, encoding='utf-8')
            write_json(meta, {'url':r.url,'sha256':sha256(target),'retrieved_utc':pd.Timestamp.now(tz='UTC').isoformat()})
            return r.text
        except requests.RequestException as exc:
            error = str(exc)
            if attempt < 2:
                time.sleep(attempt+1)
    raise ConnectionError(error)

def sounding_inventory(date, cache):
    text = cached_get('https://weather.uwyo.edu/cgi-bin/sounding_json',
        {'datetime':pd.Timestamp(date).strftime('%Y-%m-%d %H:%M:%S')}, Path(cache)/f'{pd.Timestamp(date):%Y%m%d%H}.json')
    data = json.loads(text)
    if 'stations' not in data or 'datetime' not in data:
        raise ValueError('Unexpected Wyoming inventory schema')
    if pd.Timestamp(data['datetime']) != pd.Timestamp(date).tz_localize(None):
        raise ValueError('Wyoming returned a different inventory datetime')
    return data

def wyoming_profile(date, station, source, cache):
    return cached_get('https://weather.uwyo.edu/wsgi/sounding',
        {'datetime':pd.Timestamp(date).strftime('%Y-%m-%d %H:%M:%S'), 'id':station,'src':source,'type':'TEXT:CSV'},
        Path(cache)/f'{pd.Timestamp(date):%Y%m%d%H}-{station}-{source}.csv')

def sounding_quality(path):
    f = pd.read_csv(path)
    required = ['pressure_hPa','geopotential height_m','temperature_C']
    if not set(required) <= set(f):
        raise ValueError('Unexpected sounding columns')
    trajectory = {'time','longitude','latitude'} <= set(f)
    if 'time' in f:
        f['time'] = pd.to_datetime(f.time, utc=True)
    for c in f.columns.difference(['time']):
        f[c] = pd.to_numeric(f[c], errors='coerce')
    # No derived atmospheric vertical wind from balloon ascent.
    top = float(f['geopotential height_m'].max())
    return f, {'file':str(path),'rows':len(f),'actual_start':str(f.time.min()) if trajectory else None,'actual_end':str(f.time.max()) if trajectory else None,
        'top_m':top,'min_pressure_hPa':float(f.pressure_hPa.min()),
        'status':'partial' if f.pressure_hPa.min() > 100 else 'upper_air_available',
        'missing_cells':int(f.isna().sum().sum()),'trajectory_preserved':trajectory,'sha256':sha256(path),
        'nominal_time_only':not trajectory}

def local_inputs(config, root):
    root = Path(root)
    records = []
    # Adjacent days protect midnight and spatial/time grouping boundaries.
    for day in ['2018-11-25','2018-12-12','2018-12-13','2018-12-14','2018-12-15']:
        src = Path(config['entln_source'])/(day+'.csv')
        try:
            records.append(verified_copy(src, root/'raw/entln'/src.name))
        except OSError as exc:
            records.append({'source':str(src),'status':'unavailable','reason':str(exc)})
    qualities = []
    for src in sorted(Path(config['sounding_source']).glob('*.csv')):
        records.append(verified_copy(src, root/'raw/soundings'/src.name))
        f,q = sounding_quality(root/'raw/soundings'/src.name)
        target = root/'processed/soundings'/src.name
        target.parent.mkdir(parents=True, exist_ok=True)
        f.to_csv(target,index=False)
        qualities.append(q)
    write_json(root/'manifest.json',records)
    write_json(root/'processed/soundings/quality.json',qualities)
    return records,qualities

def era5_request(config, product='reanalysis', single=False):
    # Domain includes a gradient halo; not a final observation mask.
    return {'product_type':[product], 'year':['2018'],'month':['12'],
        'day':['12','13','14','15'],'time':[f'{h:02}:00' for h in range(24)],
        'area':[-27,-68,-40,-50],'data_format':'netcdf','download_format':'unarchived',
        **({'variable':['surface_pressure','skin_temperature','geopotential','2m_temperature','2m_dewpoint_temperature','10m_u_component_of_wind','10m_v_component_of_wind',
            'mean_surface_sensible_heat_flux','mean_surface_latent_heat_flux']} if single else
           {'variable':['temperature','specific_humidity','u_component_of_wind','v_component_of_wind','vertical_velocity','geopotential'],
            'pressure_level':['1000','975','950','925','900','875','850','825','800','775','750','700','650','600','550','500','450','400','350','300','250','225','200','175','150','125','100','70','50']})}

def acquire_era5(config, cache, execute=False, selection=None):
    records = []
    inventory=Path(cache)/('inventory.json' if selection is None else f'inventory_{selection[0]}_{"sl" if selection[1] else "pl"}.json')
    for product in ['reanalysis','ensemble_members']:
        for single in [False,True]:
            if selection is not None and (product,single)!=selection:
                continue
            dataset = 'reanalysis-era5-' + ('single-levels' if single else 'pressure-levels')
            request = era5_request(config,product,single)
            if product == 'ensemble_members':
                request['time'] = [f'{h:02}:00' for h in range(0,24,3)]
            # Daily chunks stay below CDS cost limits; failure of one chunk does
            # not suppress independent levels/products/dates.
            for day in request['day']:
                target = Path(cache)/f'{product}_{"sl" if single else "pl"}_201812{day}.nc'
                daily = dict(request,day=[day])
                write_json(target.with_suffix('.request.json'),{'dataset':dataset,'request':daily})
                if execute:
                    import cdsapi
                    import xarray as xr
                    import zipfile
                    try:
                        tmp = target.with_suffix('.download')
                        sidecar = target.with_suffix('.nc.json')
                        if target.exists() and sidecar.exists():
                            previous=json.loads(sidecar.read_text())
                            if sha256(target)!=previous['sha256']:
                                raise IOError('Cached ERA5 checksum differs from original download')
                        if not target.exists():
                            cdsapi.Client().retrieve(dataset,daily,str(tmp))
                            if zipfile.is_zipfile(tmp):
                                # Read archive members as bytes, never extract arbitrary paths.
                                import tempfile
                                import shutil
                                with tempfile.TemporaryDirectory(dir=target.parent) as td, zipfile.ZipFile(tmp) as archive:
                                    parts=[]
                                    for i,name in enumerate(archive.namelist()):
                                        if name.endswith('.nc'):
                                            local=Path(td)/f'{i}.nc';local.write_bytes(archive.read(name))
                                            parts.append(xr.load_dataset(local))
                                    if not parts: raise ValueError('No NetCDF in CDS archive')
                                    xr.merge(parts,join='outer',compat='no_conflicts').to_netcdf(target)
                            else:
                                tmp.replace(target)
                        with xr.open_dataset(target) as ds:
                            if not ds.sizes: raise ValueError('Empty ERA5 file')
                        record={'path':str(target),'sha256':sha256(target),'status':'downloaded','dataset':dataset,'request':daily}
                        write_json(sidecar,record)
                        records.append(record)
                    except Exception as exc:
                        records.append({'path':str(target),'status':'acquisition_failed','reason':str(exc)})
                else:
                    records.append({'path':str(target),'status':'request_prepared_not_downloaded'})
                write_json(inventory,records)
    write_json(inventory,records)
    return records

def acquire_era5_parallel(config,cache):
    """At most three independent CDS streams, each with its own cache manifest."""
    from concurrent.futures import ThreadPoolExecutor
    selections=[(product,single) for product in ['reanalysis','ensemble_members'] for single in [False,True]]
    records=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(acquire_era5,config,cache,True,selection) for selection in selections]
        for future in futures:
            records.extend(future.result())
    write_json(Path(cache)/'inventory.json',records)
    return records
