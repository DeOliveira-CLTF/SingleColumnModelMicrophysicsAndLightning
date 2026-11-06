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
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['download','build']);ap.add_argument('--config',default=str(ROOT/'configs/cases/argentina_20181213.yaml'));a=ap.parse_args()
    c=load_config(a.config)
    if a.stage=='download':acquire(pd.date_range(c['start'],c['end'],freq='D',inclusive='left'))
    else:build(a.config)


if __name__=='__main__':main()
