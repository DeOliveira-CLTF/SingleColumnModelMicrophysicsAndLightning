"""External CHIRPS2 precipitation comparison, on its native daily support only.

CHIRPS2 daily is CFSv2-disaggregated pentad gauge/satellite precipitation. It is
independent of this SCM/ERA5 forcing, but not independent daily gauge truth.
It cannot validate 15-minute convective timing. Preserve that distinction.
"""
import argparse
import gzip
from io import BytesIO
from pathlib import Path
import numpy as np
import pandas as pd
import requests
import xarray as xr
from PIL import Image
from scm_lightning.config import load_config
from scm_lightning.forcing.spatial import subset
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.ensemble.experiments import control_directory

BASE='https://data.chc.ucsb.edu/products/CHIRPS/v2.0/'
CACHE=ROOT/'data/raw/precipitation/chirps_v2'


def acquire_era5_precipitation(config):
    """ERA5 hourly totals, in metres over the hour ending at valid_time."""
    import cdsapi
    import zipfile
    cache=ROOT/'data/raw/precipitation/era5'/config['event_id'];cache.mkdir(parents=True,exist_ok=True)
    target=cache/'total_precipitation.nc'
    if target.exists():
        meta=__import__('json').loads(target.with_suffix('.json').read_text())
        if sha256(target)!=meta['sha256']:raise ValueError('ERA5 precipitation cache checksum mismatch')
        return target
    days=pd.date_range(config['start'],config['end'],freq='D')
    w,s,e,n=config['area']
    request={'product_type':['reanalysis'],'variable':['total_precipitation'],
        'year':[str(days[0].year)],'month':[f'{days[0].month:02}'],'day':[f'{d.day:02}' for d in days],
        'time':[f'{h:02}:00' for h in range(24)],'area':[n+.5,w-.5,s-.5,e+.5],
        'data_format':'netcdf','download_format':'unarchived'}
    write_json(cache/'request.json',{'dataset':'reanalysis-era5-single-levels','request':request})
    pending=target.with_suffix('.download')
    cdsapi.Client().retrieve('reanalysis-era5-single-levels',request,str(pending))
    if zipfile.is_zipfile(pending):
        with zipfile.ZipFile(pending) as archive:
            names=[n for n in archive.namelist() if n.endswith('.nc')]
            if len(names)!=1:raise ValueError('Unexpected ERA5 precipitation archive')
            target.write_bytes(archive.read(names[0]))
    else:pending.replace(target)
    with xr.open_dataset(target) as d:
        if 'tp' not in d or d.tp.attrs.get('units')!='m':raise ValueError('Expected ERA5 tp in metres')
    write_json(target.with_suffix('.json'),{'source':'ERA5 CDS','sha256':sha256(target),'request':request,
        'temporal_support':'one hour ending at valid_time','independent_observation':False})
    return target


def era5_hourly(config):
    from scm_lightning.forcing import normalize
    source=ROOT/'data/raw/precipitation/era5'/config['event_id']/'total_precipitation.nc'
    with xr.open_dataset(source) as original:
        selected,weights=subset(normalize(original),config['area'])
        rain=1000*selected.tp.weighted(weights).mean(('latitude','longitude')).load()
    start,end=pd.Timestamp(config['start']),pd.Timestamp(config['end'])
    valid=pd.DatetimeIndex(rain.time.values).tz_localize('UTC')
    wanted=pd.date_range(start+pd.Timedelta(hours=1),end,freq='h')
    f=pd.DataFrame({'valid_time':valid,'ERA5_mm':rain.values}).set_index('valid_time').reindex(wanted)
    if f.ERA5_mm.isna().any() or (f.ERA5_mm<0).any():raise ValueError('Incomplete or negative ERA5 hourly precipitation')
    f['time']=wanted-pd.Timedelta(hours=1);f['end_time']=wanted
    rows=[]
    for configuration in ['TH_C0','TH_CL','TH_CF','WS_C0','WS_CL','WS_CF']:
        path=control_directory(configuration)/'output.nc'
        with xr.open_dataset(path,decode_times=False) as d:
            for row in f.itertuples():
                a=(row.time-start).total_seconds();b=a+3600
                amount=daily_amount(d.time.values,d.rainnc.values,a,b)
                rows.append({'time':str(row.time),'end_time':str(row.end_time),'configuration':configuration,
                    'simulation_mm':amount,'ERA5_mm':row.ERA5_mm,'difference_mm':amount-row.ERA5_mm})
    out=ROOT/'results'/config['event_id']/'meteorological_validation';out.mkdir(parents=True,exist_ok=True)
    result=pd.DataFrame(rows);result.to_csv(out/'precipitation_comparison_ERA5_hourly.csv',index=False)
    write_json(out/'ERA5_manifest.json',{'source_sha256':sha256(source),'area_wsen':config['area'],
        'units':'mm per hour; numerical rate mm/h','ERA5_time':'valid_time labels preceding hour; table time labels interval start',
        'support':'same fixed SCM area, spherical area weights, exact complete UTC hourly intervals',
        'simulation':'rainnc differences only; no snownc addition','independent_observation':False,
        'interpretation':'comparison with forcing reanalysis, not independent precipitation validation',
        'temporal_convention_source':'https://confluence.ecmwf.int/pages/viewpage.action?pageId=111150967'})
    return result


def acquire(days):
    CACHE.mkdir(parents=True,exist_ok=True)
    for day in days:
        stamp=pd.Timestamp(day)
        filename=f'chirps-v2.0.{stamp:%Y.%m.%d}.tif.gz'
        target=CACHE/filename;url=BASE+f'global_daily/tifs/p05/{stamp:%Y}/'+filename
        if not target.exists():
            r=requests.get(url,timeout=60);r.raise_for_status()
            pending=target.with_suffix('.pending');pending.write_bytes(r.content)
            with Image.open(BytesIO(gzip.decompress(pending.read_bytes()))) as img:img.verify()
            pending.replace(target)
        write_json(target.with_suffix('.json'),{'url':url,'sha256':sha256(target),'product':'CHIRPS v2 final daily',
            'units':'mm/day','day_label':'UTC calendar day; underlying station windows heterogeneous',
            'daily_disaggregation':'CFSv2; gauge/satellite constraint at pentad/monthly support',
            'documentation':BASE+'README-CHIRPS2-Daily.txt'})


def regional_field(path,area):
    with Image.open(BytesIO(gzip.decompress(Path(path).read_bytes()))) as img:
        scale=img.tag_v2[33550];tie=img.tag_v2[33922]
        dx,dy=float(scale[0]),float(scale[1]);x0,y0=float(tie[3]),float(tie[4])
        w,s,e,n=area
        left=max(0,int(np.floor((w-x0)/dx))-1);right=min(img.width,int(np.ceil((e-x0)/dx))+1)
        top=max(0,int(np.floor((y0-n)/dy))-1);bottom=min(img.height,int(np.ceil((y0-s)/dy))+1)
        a=np.asarray(img.crop((left,top,right,bottom)),float)
        lon=x0+(np.arange(left,right)+.5)*dx;lat=y0-(np.arange(top,bottom)+.5)*dy
    a[a<0]=np.nan
    return xr.DataArray(a,dims=('latitude','longitude'),coords={'latitude':lat,'longitude':lon},name='precipitation_mm').sortby('latitude')


def daily_amount(seconds,accumulated,start_second,end_second):
    seconds=np.asarray(seconds,float);a=np.asarray(accumulated,float)
    if np.any(np.diff(seconds)<=0) or np.any(np.diff(a)<-1e-8):raise ValueError('Invalid cumulative total precipitation')
    if seconds[0]>start_second or seconds[-1]<end_second:raise ValueError('Incomplete day; no extrapolation')
    return float(np.interp(end_second,seconds,a)-np.interp(start_second,seconds,a))


def build(config_path):
    c=load_config(config_path)
    if c['area'] is None:raise ValueError('Freeze weather-only area first')
    event=c.get('event_id','argentina_20181213')
    out=ROOT/'results'/event/'meteorological_validation';out.mkdir(parents=True,exist_ok=True)
    days=pd.date_range(c['start'],c['end'],freq='D',inclusive='left')
    observed=[];grids=[]
    for day in days:
        path=CACHE/f'chirps-v2.0.{day:%Y.%m.%d}.tif.gz'
        field=regional_field(path,c['area']);selected,weights=subset(field.to_dataset(),c['area'])
        if not np.isfinite(selected.precipitation_mm).all():raise ValueError('Incomplete external precipitation area')
        observed.append({'time':str(day),'precipitation_mm':float(selected.precipitation_mm.weighted(weights).mean()),
            'product':'CHIRPS2 daily, CFSv2-disaggregated','source_sha256':sha256(path)})
        grids.append(field.expand_dims(time=[day.tz_localize(None)]))
    xr.concat(grids,dim='time').to_dataset().to_netcdf(out/'observed_precipitation_daily.nc')
    obs=pd.DataFrame(observed);obs.to_csv(out/'observed_precipitation_daily.csv',index=False)
    records=[]
    configs=['TH_C0','TH_CL','TH_CF','WS_C0','WS_CL','WS_CF'] if event=='argentina_20181213' else ['TH_C0']
    for cfg in configs:
        source=control_directory(cfg)/'output.nc' if event=='argentina_20181213' else ROOT/'runs'/event/'thompson_C0/control/output.nc'
        if not source.exists():continue
        digest=sha256(source)
        with xr.open_dataset(source,decode_times=False) as d:
            for i,day in enumerate(days):
                amount=daily_amount(d.time.values,d.rainnc.values,i*86400,(i+1)*86400)
                records.append({'time':str(day),'configuration':cfg,'simulation_mm':amount,'external_estimate_mm':observed[i]['precipitation_mm'],
                    'difference_mm':amount-observed[i]['precipitation_mm'],'simulation_sha256':digest})
    pd.DataFrame(records).to_csv(out/'precipitation_comparison_daily.csv',index=False)
    write_json(out/'manifest.json',{'event':event,'area_wsen':c['area'],'support':'calendar-day area means, mm/day',
        'independent_of_SCM_and_ERA5_forcing':True,'independent_daily_gauge_truth':False,'timing_validation_15min':False,
        'daily_CHIRPS2_method':'pentad/monthly gauge/IR blend, disaggregated with CFSv2; not independent daily timing',
        'documentation':BASE+'README-CHIRPS2-Daily.txt','documentation_sha256':sha256(CACHE/'README-CHIRPS2-Daily.txt'),
        'precipitation_simulation':'rainnc only; includes frozen-water equivalent; no snownc addition',
        'native_simulations_modified':False})


def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['download','build','era5-download','era5-build']);ap.add_argument('--config',default=str(ROOT/'configs/cases/argentina_20181213.yaml'));a=ap.parse_args()
    c=load_config(a.config)
    if a.stage=='download':acquire(pd.date_range(c['start'],c['end'],freq='D',inclusive='left'))
    elif a.stage=='era5-download':acquire_era5_precipitation(c)
    elif a.stage=='era5-build':era5_hourly(c)
    else:build(a.config)


if __name__=='__main__':main()
