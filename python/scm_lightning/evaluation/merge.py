"""Hourly MERGE/ERA5/SCM on identical spherical fixed-box support."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds
import xarray as xr
from scm_lightning.config import load_config
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.forcing.spatial import subset
from scm_lightning.forcing import normalize
from scm_lightning.ensemble.experiments import control_directory
from scm_lightning.evaluation.meteorology import daily_amount
from scm_lightning.download.merge import CACHE


def regional_field(path,area,expected_time):
    with rasterio.open(path) as d:
        indexes=[i for i in d.indexes if d.tags(i).get('GRIB_ELEMENT')=='PREC']
        if len(indexes)!=1:raise ValueError('Unambiguous PREC band required; station-count band is not rain')
        band=indexes[0];tags=d.tags(band)
        if tags.get('GRIB_UNIT')!='[kg/(m^2)]':raise ValueError('Unexpected MERGE accumulation unit')
        valid=pd.to_datetime(int(tags['GRIB_VALID_TIME']),unit='s',utc=True)
        if valid!=pd.Timestamp(expected_time):raise ValueError('MERGE valid-time mismatch')
        west,south,east,north=area
        window=from_bounds(west-.2,south-.2,east+.2,north+.2,d.transform).round_offsets().round_lengths()
        a=d.read(band,window=window).astype(float)
        transform=d.window_transform(window)
        lon=transform.c+(np.arange(a.shape[1])+.5)*transform.a
        lat=transform.f+(np.arange(a.shape[0])+.5)*transform.e
        a[(a<0)|(~np.isfinite(a))]=np.nan
        field=xr.DataArray(a,dims=['latitude','longitude'],coords={'latitude':lat,'longitude':lon},name='MERGE_mm').sortby('latitude')
    selected,weight=subset(field.to_dataset(),area)
    valid_area=float(weight.where(selected.MERGE_mm.notnull(),0).sum()/weight.sum())
    amount=float(selected.MERGE_mm.weighted(weight).mean()) if np.isclose(valid_area,1,rtol=0,atol=1e-12) else np.nan
    return amount,selected.MERGE_mm,{'GRIB':tags,'covered_area_fraction':valid_area,'grid_spacing_degrees':[abs(transform.a),abs(transform.e)]}


def era5(config):
    event=config['event_id'];paths=[];parts=[]
    if event=='argentina_20181213':
        paths=[ROOT/'data/raw/precipitation/era5'/event/'total_precipitation.nc']
    else:
        paths=[ROOT/config['paths']['era5']/f'reanalysis_sl_{d:%Y%m%d}.nc' for d in pd.date_range(config['start'],config['end'],freq='D')]
    for path in paths:
        with xr.open_dataset(path) as original:
            normalized=normalize(original)
            if normalized.tp.attrs.get('units')!='m':raise ValueError('ERA5 tp must be in metres')
            s,w=subset(normalized[['tp']],config['area'])
            if not np.isfinite(s.tp.values).all():raise ValueError('Incomplete ERA5 fixed-box precipitation')
            parts.append((1000*s.tp.weighted(w).mean(('latitude','longitude'))).load())
    rain=xr.concat(parts,dim='time').sortby('time')
    times=pd.DatetimeIndex(rain.time.values).tz_localize('UTC')
    f=pd.Series(rain.values,index=times)
    if f.index.duplicated().any():raise ValueError('Duplicate ERA5 hourly validity')
    return f,[{'path':str(p),'sha256':sha256(p)} for p in paths]


def metrics(table):
    ref=table.MERGE_mm.to_numpy(float)
    shared=np.isfinite(table[['MERGE_mm','ERA5_mm','TH_C0_mm','WS_C0_mm']].to_numpy(float)).all(axis=1)
    f=table[shared];y=f.MERGE_mm.to_numpy(float);rows=[]
    for name in ['ERA5','TH_C0','WS_C0']:
        p=f[name+'_mm'].to_numpy(float);e=p-y
        peak=np.nan
        if len(y)>1 and np.ptp(y)>0 and np.ptp(p)>0:
            clock=pd.to_datetime(f.time,utc=True);peak=(clock.iloc[np.argmax(p)]-clock.iloc[np.argmax(y)]).total_seconds()/3600
        rows.append({'series':name,'reference':'MERGE_CPTEC_INPE','n_shared_hourly_intervals':len(f),
            'RMSE_mm_h':float(np.sqrt(np.mean(e**2))) if len(e) else np.nan,
            'bias_mm_h':float(np.mean(e)) if len(e) else np.nan,
            'accumulated_difference_mm':float(e.sum()) if len(e) else np.nan,
            'peak_time_error_hours':peak,'peak_timing_resolution_hours':1,
            'full_period_accumulation_valid':bool(shared.all())})
    return pd.DataFrame(rows),shared


def build(event):
    c=load_config(ROOT/'configs/cases'/(event+'.yaml'));start=pd.Timestamp(c['start']);end=pd.Timestamp(c['end'])
    downloaded=json.loads((CACHE/event/'download_manifest.json').read_text());times=pd.date_range(start,end,freq='h',inclusive='right')
    directory=ROOT/'results'/event/'meteorological_validation';directory.mkdir(parents=True,exist_ok=True)
    era,era_sources=era5(c);rows=[];grids=[];sources=[]
    for valid in times:
        path=CACHE/event/f'MERGE_CPTEC_{valid:%Y%m%d%H}.grib2'
        descriptor=path.with_suffix('.ctl')
        if 'PREC' not in descriptor.read_text() or '[mm/hr]' not in descriptor.read_text():raise ValueError('Missing official hourly descriptor')
        amount,field,meta=regional_field(path,c['area'],valid)
        rows.append({'time':str(valid-pd.Timedelta(hours=1)),'end_time':str(valid),'MERGE_mm':amount,
            'ERA5_mm':float(era.get(valid,np.nan)),'MERGE_covered_area_fraction':meta['covered_area_fraction']})
        grids.append(field.expand_dims(time=[(valid-pd.Timedelta(hours=1)).tz_localize(None)]))
        sources.append({'path':str(path),'sha256':sha256(path),'descriptor_sha256':sha256(descriptor),'metadata':meta})
    table=pd.DataFrame(rows)
    native=[]
    for cfg,name in [('TH_C0','thompson'),('WS_C0','wsm6')]:
        path=control_directory(cfg)/'output.nc' if event=='argentina_20181213' else ROOT/'runs'/event/(name+'_C0')/'control/output.nc'
        before=sha256(path)
        with xr.open_dataset(path,decode_times=False) as d:
            if not np.isin(np.arange(len(table)+1)*3600,d.time.values).all():raise ValueError('SCM hourly bounds must be actual native snapshots')
            table[cfg+'_mm']=[daily_amount(d.time.values,d.rainnc.values,i*3600,(i+1)*3600) for i in range(len(table))]
        if sha256(path)!=before:raise IOError('Native output changed')
        native.append({'configuration':cfg,'path':str(path),'sha256':before})
    scores,shared=metrics(table);table['common_valid']=shared
    table.to_csv(directory/'precipitation_MERGE_ERA5_SCM_hourly.csv',index=False)
    scores.to_csv(directory/'precipitation_metrics_vs_MERGE.csv',index=False)
    data=xr.concat(grids,dim='time').to_dataset();data.MERGE_mm.attrs.update(units='mm',interval_seconds=3600)
    _,weights=subset(data,c['area']);data['cell_overlap_area_km2']=weights
    data.cell_overlap_area_km2.attrs.update(units='km2',long_name='spherical grid-cell intersection with fixed SCM box')
    data['area_mean_precipitation']=('time',table.MERGE_mm.values)
    data.area_mean_precipitation.attrs.update(units='mm',interval_seconds=3600,rate_convention='one-hour amount numerically mm/h')
    data['time_bounds']=(('time','bounds'),np.column_stack([data.time.values,data.time.values+np.timedelta64(1,'h')]))
    data.time.attrs['bounds']='time_bounds'
    for name in ['time','time_bounds']:data[name].encoding.update(units='hours since '+start.strftime('%Y-%m-%d %H:%M:%S'),calendar='proleptic_gregorian')
    data.attrs.update(area_wsen=json.dumps(c['area']),product='MERGE CPTEC/INPE GPM HOURLY',temporal_support='hour ending at raw GRIB validity; exported time=interval start')
    data.to_netcdf(directory/'MERGE_fixed_box_hourly.nc')
    write_json(directory/'MERGE_validation_manifest.json',{'event':event,'area_wsen':c['area'],
        'spatial_aggregation':'exact spherical cell/box overlap area weights; missing-cell area is not renormalized',
        'time_support':[c['start'],c['end']],'interval_seconds':3600,'expected_intervals':len(times),
        'common_valid_intervals':int(shared.sum()),'MERGE_sources':sources,'ERA5_sources':era_sources,'SCM_sources':native,
        'MERGE_catalog_sha256':sha256(CACHE/'prec_merge_hourly-1.json'),'MERGE_download_manifest_sha256':sha256(CACHE/event/'download_manifest.json'),
        'MERGE_version':downloaded['archive_algorithm_version'],'STAC_collection_version':'1','resolution_degrees':.1,
        'independent_observational_reference':'MERGE blended satellite/gauge estimate; not asserted error-free gauge truth',
        'ERA5_independent_observation':False,'GRIB_amount_conversion':'1 kg/m2 water equivalent = 1 mm in each one-hour interval',
        'peak_timing_limit':'one-hour resolution, not minute timing','no_daily_to_hourly_interpolation':True,'native_outputs_modified':False,
        'versions':{'rasterio':rasterio.__version__,'GDAL':rasterio.__gdal_version__}})
    print(f'{event}: MERGE comparison {shared.sum()}/{len(times)} complete UTC hours.',flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--event',choices=['argentina_20181213','argentina_20191028','both'],default='both');a=p.parse_args()
    for event in ['argentina_20181213','argentina_20191028'] if a.event=='both' else [a.event]:build(event)


if __name__=='__main__':main()
