from scm_lightning.config import load_config
"""Freeze a weather-only fixed column support before any lightning evaluation.

The supplied visual activity sector bounds the candidate search, never the
camera field of view. Ranking uses only regional ERA5 synoptic ascent.
This selection does not certify that a forced SCM resolves storm convection.
"""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scm_lightning.forcing import normalize
from scm_lightning.provenance import ROOT,sha256,write_json

def main():
    path=ROOT/'configs/cases/argentina_20181213.yaml';c=load_config(path)
    if c['area_selection']['status']=='fixed_before_lightning_evaluation':
        raise FileExistsError('Meteorological support is already frozen; retain the original selection manifest')
    if c['area'] is not None:
        raise ValueError('Area already frozen; refusing reselection')
    files=[ROOT/'data/raw/era5'/f'reanalysis_pl_201812{day}.nc' for day in ['13','14']]
    data=xr.concat([normalize(xr.load_dataset(p)) for p in files],dim='time').sortby('time')
    start=pd.Timestamp(c['start']).tz_localize(None);end=pd.Timestamp(c['end']).tz_localize(None)
    data=data.sel(time=slice(start,end-pd.Timedelta('1ns')))
    if len(data.time)!=48:raise ValueError('Require complete 48 hourly ERA5 snapshots')
    half=c['meteorological_selection']['box_width_degrees']/2
    w,s,e,n=c['area_selection']['visual_reference_only']
    omega=data.w.sel(level=[500,700])
    score=(-omega).clip(min=0).mean(('time','level'))
    eligible=score.sel(longitude=slice(w+half,e-half),latitude=slice(s+half,n-half))
    # Deterministic ties: lowest latitude, then lowest longitude, by sorted grid.
    iy,ix=np.unravel_index(np.nanargmax(eligible.values),eligible.shape)
    lat,lon=float(eligible.latitude[iy]),float(eligible.longitude[ix])
    area=[lon-half,lat-half,lon+half,lat+half]
    decision={'method':c['meteorological_selection'],'selected_area_wsen':area,
        'score_Pa_per_s':float(eligible.values[iy,ix]),'search_region_wsen':[w,s,e,n],
        'inputs':[{'path':str(p),'sha256':sha256(p)} for p in files],
        'time_start':c['start'],'time_end_exclusive':c['end'],'snapshots':48,
        'created_utc':pd.Timestamp.now(tz='UTC').isoformat(),
        'lightning_error_used':False,'camera_fov_used':False,
        'representativeness':'Meteorological support only; no claim of resolved storm dynamics or exact hydrostatic balance'}
    out=ROOT/'results/argentina_20181213/support';out.mkdir(parents=True,exist_ok=True)
    score.to_dataset(name='mean_positive_ascent').to_netcdf(out/'meteorological_score.nc')
    write_json(out/'selection.json',decision)
    c['area']=area
    c['area_selection'].update(status='fixed_before_lightning_evaluation',criterion='ERA5 mean positive ascent at 500/700 hPa, full case; 1 degree box fixed before lightning evaluation',decision_sha256=sha256(out/'selection.json'))
    write_json(path,c)
    fig,axes=plt.subplots(1,3,figsize=(15,5))
    for ax,field,title in [(axes[0],score,'Mean positive ascent 500/700 hPa (Pa/s)'),
        (axes[1],data.t.sel(level=850).mean('time'),'Mean 850-hPa temperature (K)'),
        (axes[2],1000*data.q.sel(level=850).mean('time'),'Mean 850-hPa specific humidity (g/kg)')]:
        field.plot(ax=ax,cbar_kwargs={'label':title})
        ax.add_patch(Rectangle((area[0],area[1]),2*half,2*half,fill=False,edgecolor='red',linewidth=2))
        ax.set_xlim(c['plot_extent_requested'][:2]);ax.set_ylim(c['plot_extent_requested'][2:]);ax.set_title('ERA5: 13–14 Dec 2018 UTC')
    fig.suptitle('Fixed meteorological support (red); blank margins outside acquired ERA5 domain')
    fig.tight_layout();fig.savefig(out/'selection.png',dpi=160,bbox_inches='tight');plt.close(fig)
    print(json.dumps({'area_wsen':area,'selection_manifest':str(out/'selection.json')}))

if __name__=='__main__':main()
