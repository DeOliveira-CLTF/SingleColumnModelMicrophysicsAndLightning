from scm_lightning.config import load_config
"""Four main technical figures. No interpretation text is generated."""
from pathlib import Path
import argparse
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
from scm_lightning.ensemble.experiments import CONFIGURATIONS
from scm_lightning.provenance import ROOT,write_json,sha256

OUT=ROOT/'results/argentina_20181213/figures'

def save(fig,name):
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/name,dpi=180,bbox_inches='tight');plt.close(fig)


def boundaries(ax):
    """Draw cached country outlines and global first-level subdivisions."""
    root=ROOT/'data/raw/boundaries'
    def line(coordinates, **style):
        xy=np.asarray(coordinates)
        ax.plot(xy[:,0],xy[:,1],**style)
    for name in ['states','countries']:
        content=json.loads((root/(name+'.geojson')).read_text(encoding='utf-8'))
        style={'color':'0.18','linewidth':.65 if name=='countries' else .35,'alpha':.9,'zorder':3}
        for feature in content['features']:
            geometry=feature['geometry']
            if geometry is None:continue
            kind=geometry['type'];coordinates=geometry['coordinates']
            if kind=='LineString':line(coordinates,**style)
            elif kind in ['MultiLineString','Polygon']:
                for points in coordinates:line(points,**style)
            elif kind=='MultiPolygon':
                for polygon in coordinates:
                    for points in polygon:line(points,**style)


def case(c):
    ds=normalize(xr.load_dataset(ROOT/'data/raw/era5/reanalysis_pl_20181214.nc'))
    field=ds.w.sel(time='2018-12-14T06:00',level=500)
    if not np.isfinite(field).all():raise ValueError('Map contains missing ERA5 data')
    fig,axes=plt.subplots(1,3,figsize=(15,5.5),gridspec_kw={'width_ratios':[1.5,1,1]})
    limit=float(abs(field).quantile(.99))
    field.plot(ax=axes[0],cmap='RdBu_r',vmin=-limit,vmax=limit,cbar_kwargs={'label':'Omega (Pa/s); negative = ascent'})
    boundaries(axes[0])
    w,s,e,n=c['area'];axes[0].add_patch(Rectangle((w,s),e-w,n-s,fill=False,color='red',linewidth=2,label='Fixed SCM support'))
    station=c['station']
    axes[0].scatter([station['longitude']],[station['latitude']],marker='^',color='black',s=35,zorder=5,label='SACO / WMO 87344')
    axes[0].annotate('SACO\n31.32 S, 64.22 W',(station['longitude'],station['latitude']),
        xytext=(-2,20),textcoords='offset points',fontsize=7,ha='center',
        arrowprops={'arrowstyle':'-','color':'black'},zorder=5)
    extent=c['plot_extent_requested']
    # Crop to actual cell-center coverage. No blank frame or invented extrapolation.
    bounds=[max(extent[0],float(field.longitude.min())),min(extent[1],float(field.longitude.max())),
        max(extent[2],float(field.latitude.min())),min(extent[3],float(field.latitude.max()))]
    axes[0].set_xlim(bounds[:2]);axes[0].set_ylim(bounds[2:]);axes[0].set_aspect(1/np.cos(np.deg2rad(np.mean(bounds[2:]))))
    axes[0].set_title('ERA5 500 hPa: 14 Dec 2018, 06 UTC');axes[0].legend(fontsize=7,loc='lower left')
    axes[0].set(xlabel='Longitude',ylabel='Latitude')
    write_json(OUT/'01_map_contract.json',{'omega_sign':'negative is ascent; native omega plotted',
        'station':station,'actual_extent':bounds,'requested_extent':extent,
        'finite_field':True,'country_and_state_boundaries':'Natural Earth, cached sources in data/raw/boundaries/manifest.json'})
    for p in sorted((ROOT/'data/processed/profiles').glob('*.csv')):
        f=pd.read_csv(p)
        if 'time' in f:
            stamp=pd.Timestamp(f.time.iloc[0]);label=stamp.strftime('%d Dec %H:%M UTC')
            if '1403' in p.name:label+=' (partial)'
        else:label='14 Dec 12 UTC (nominal)'
        z=f['geopotential height_m']/1000
        axes[1].plot(f.temperature_C,z,label=label)
        axes[2].plot(f['mixing ratio_g/kg'],z,label=label)
    axes[1].set(xlabel='Temperature (C)',ylabel='Geopotential height (km)',title='Córdoba profiles')
    axes[2].set(xlabel='Vapor mixing ratio (g/kg)',ylabel='Geopotential height (km)',title='Córdoba profiles')
    for ax in axes[1:]:ax.grid(alpha=.2)
    handles,labels=axes[1].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
    fig.tight_layout(rect=[0,.08,1,1]);save(fig,'01_case_area_environment.png')

def microphysics(c):
    fig,axes=plt.subplots(6,2,figsize=(14,16),sharex=True)
    categories=['Deposição/sublimação','Congelamento','Riming (água de nuvem)','Fusão','Sedimentação']
    root=ROOT/'results/argentina_20181213/control_diagnostics'
    for row,configuration in enumerate(CONFIGURATIONS):
        p=pd.read_csv(root/(configuration+'_paths.csv'))
        time=pd.Timestamp(c['start'])+pd.to_timedelta(p.time_s,unit='s')
        for label,value in [('cloud liquid',p.qc),('rain',p.qr),('cloud ice + snow',p.qi+p.qs),('graupel',p.qg)]:
            axes[row,0].plot(time,value,label=label,linewidth=.9)
        r=pd.read_csv(root/(configuration+'_processes.csv.gz'))
        r['time']=pd.Timestamp(c['start'])+pd.to_timedelta(r.time_s,unit='s')
        r['ice_rate']=r.dqi+r.dqs+r.dqg
        pivot=r.pivot_table(index='time',columns='category',values='ice_rate',aggfunc='sum').reindex(time.iloc[1:],fill_value=0).fillna(0)
        for category in categories:
            values=pivot[category] if category in pivot else pd.Series(0.,index=pivot.index)
            values=values.resample('15min',closed='right',label='right').mean()
            axes[row,1].plot(values.index,values,label=category,linewidth=.9)
        axes[row,0].set_ylabel(configuration+'\nkg m-2');axes[row,1].set_ylabel('kg m-2 s-1')
        for ax in axes[row]:ax.grid(alpha=.2);ax.tick_params(axis='x',rotation=25,labelsize=7)
    axes[0,0].set_title('Instantaneous dry-air density weighted species paths')
    axes[0,1].set_title('Frozen-water process rates / initial dry-air cell mass')
    axes[-1,0].set_xlabel('UTC');axes[-1,1].set_xlabel('UTC')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',bbox_to_anchor=(.28,.005),ncol=2,fontsize=8)
    h,l=axes[0,1].get_legend_handles_labels();fig.legend(h,l,loc='lower center',bbox_to_anchor=(.76,.005),ncol=2,fontsize=8)
    fig.tight_layout(rect=[0,.06,1,1]);save(fig,'04_microphysics_processes.png')

def evaluation(c):
    root=ROOT/'results/argentina_20181213/evaluation'
    h=pd.read_csv(root/'heldout_predictions.csv');h['time']=pd.to_datetime(h.time,utc=True)
    fig,ax=plt.subplots(figsize=(13,5))
    obs=h.groupby('time')['count'].first();ax.plot(obs.index,obs,color='black',linewidth=1.8,label='Observed total flashes')
    colors={name:plt.get_cmap('tab10')(i) for i,name in enumerate(CONFIGURATIONS)}
    for (method,configuration),f in h.groupby(['method','configuration']):
        f=f[f.member!=0]
        g=f.groupby('time').prediction_A
        mean=g.mean()
        if not mean.notna().any():continue
        ax.plot(mean.index,mean,color=colors[configuration],linewidth=1.1,label=configuration+' / '+method)
        ax.fill_between(mean.index,g.quantile(.05),g.quantile(.95),color=colors[configuration],alpha=.10)
    ax.set(xlabel='UTC (held-out block)',ylabel='Flashes / 15 min',title='Held-out calibrated predictions; bands: initial-condition 5–95% spread')
    ax.grid(alpha=.2);ax.legend(ncol=3,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.15));fig.tight_layout();save(fig,'02_heldout_predictions.png')
    m=pd.read_csv(root/'metrics.csv');a=m[(m.comparison=='A')&(m.n>0)]
    contrast=pd.read_csv(root/'fixed_scale_skill_contrasts.csv');contrast=contrast[contrast.n>0]
    physical=pd.read_csv(root/'physical_trajectories.csv');physical=physical.groupby(['method','configuration']).physical_change.sum(min_count=1).dropna()
    fig,axes=plt.subplots(1,3,figsize=(15,6))
    axes[0].barh(a.method+'/'+a.configuration,a.rmse);axes[0].set_xlabel('A: count RMSE')
    axes[1].barh(contrast.method+'/'+contrast.configuration,contrast.rmse_change);axes[1].axvline(0,color='black');axes[1].set_xlabel('B: RMSE change from TH_C0')
    axes[2].barh(['/'.join(k) for k in physical.index],physical.values);axes[2].axvline(0,color='black');axes[2].set_xlabel('B: accumulated count change from TH_C0')
    from matplotlib.ticker import MaxNLocator
    for ax in axes:
        ax.grid(axis='x',alpha=.2);ax.tick_params(axis='y',labelsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(4));ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    fig.tight_layout();save(fig,'03_metrics_contrasts.png')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--controls-only',action='store_true');args=ap.parse_args()
    c=load_config(ROOT/'configs/cases/argentina_20181213.yaml')
    case(c);microphysics(c)
    if not args.controls_only:evaluation(c)
    write_json(OUT/'manifest.json',{'configuration':c,'figures':[{'path':p.name,'sha256':sha256(p)} for p in sorted(OUT.glob('*.png'))],
        'process_normalization':'initial dry-air cell mass; not exact instantaneous mass flux','prediction_band_scope':'initial-condition spread only'})

if __name__=='__main__':main()
