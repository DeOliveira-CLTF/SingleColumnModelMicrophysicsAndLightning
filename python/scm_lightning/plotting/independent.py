"""Principal figures bound to frozen independent-event CSVs; legends inside axes."""
import argparse
import json
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
from scm_lightning.config import load_config
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.evaluation.independent import OUT, EVALUATION_CONFIG, CALIBRATION_CONFIG
from scm_lightning.calibration.independent import METHODS
from scm_lightning.plotting.article import STYLE,BLUE,RED,save,utc_axis,environment,observations,microphysical_response
from scm_lightning.plotting.figures import boundaries

LABELS=dict(zip(METHODS,['CAPE × P','I × G','McCaul F3 — LS']))


def calibration_weather():
    c=load_config(CALIBRATION_CONFIG)
    if c['area'] is None:return
    from scm_lightning.forcing import normalize
    from scm_lightning.forcing.spatial import subset
    from scm_lightning.observations.profiles import selected
    from scm_lightning.plotting.profiles import draw
    cache=ROOT/c['paths']['era5'];start=pd.Timestamp(c['start']);end=pd.Timestamp(c['end'])
    omega=[];rain=[]
    for day in pd.date_range(start,end,freq='D'):
        with xr.open_dataset(cache/f'reanalysis_pl_{day:%Y%m%d}.nc') as source:
            omega.append(normalize(source).w.sel(level=[500,700]).mean('level').load())
        with xr.open_dataset(cache/f'reanalysis_sl_{day:%Y%m%d}.nc') as source:
            rain.append(normalize(source).tp.load())
    omega=xr.concat(omega,dim='time').sel(time=slice(start.tz_localize(None),end.tz_localize(None)-pd.Timedelta('1ns')))
    rain=xr.concat(rain,dim='time')
    # ERA5 tp is an hourly amount ending at valid_time. Align to the actual event.
    rain=rain.assign_coords(time=rain.time-np.timedelta64(1,'h')).sel(time=omega.time)*1000
    area=c['area'];west,south,east,north=area
    fig=plt.figure(figsize=(14,8.5));grid=fig.add_gridspec(2,6)
    axes=[fig.add_subplot(grid[0,i*2:(i+1)*2]) for i in range(3)]
    profiles=[fig.add_subplot(grid[1,:3]),fig.add_subplot(grid[1,3:])]
    for ax,field,unit in [(axes[0],omega.mean('time'),'Pa s$^{-1}$'),(axes[1],rain.sum('time'),'mm')]:
        if ax is axes[0]:
            limit=float(abs(field).max());field.plot(ax=ax,cmap='RdBu_r',vmin=-limit,vmax=limit,cbar_kwargs={'label':unit},add_labels=False)
        else:field.plot(ax=ax,cbar_kwargs={'label':unit},add_labels=False)
        boundaries(ax,color='white' if ax is axes[1] else '0.18');ax.add_patch(Rectangle((west,south),east-west,north-south,fill=False,ec='black',lw=1.5))
        ax.set(xlim=(float(field.longitude.min()),float(field.longitude.max())),ylim=(float(field.latitude.min()),float(field.latitude.max())),xlabel='Longitude',ylabel='Latitude')
    def mean_area(field):
        box,weights=subset(field.to_dataset(name='value'),area)
        return box.value.weighted(weights).mean(('latitude','longitude'))
    clock=pd.to_datetime(omega.time.values,utc=True)
    axes[2].plot(clock,mean_area(omega),c=BLUE,label=r'$\omega$')
    rainax=axes[2].twinx();rainax.plot(clock+pd.Timedelta(minutes=30),mean_area(rain),c=RED,label='Rain')
    axes[2].set_ylabel('Pa s$^{-1}$');rainax.set_ylabel('mm h$^{-1}$');axes[2].set_xlabel('Day / UTC')
    axes[2].xaxis.set_major_locator(mdates.HourLocator(byhour=[0,6,12,18]))
    axes[2].xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'))
    axes[2].tick_params(axis='x',labelrotation=90);axes[2].grid(alpha=.2)
    handles,labels=axes[2].get_legend_handles_labels();h,l=rainax.get_legend_handles_labels()
    axes[2].legend(handles+h,labels+l,loc='upper center',ncol=2,frameon=False)
    axes[0].set_title('(a) Mean $\omega$, 500/700 hPa\n'+start.strftime('%m/%d/%Y %H:%M')+' -\n'+clock[-1].strftime('%m/%d/%Y %H:%M UTC'),fontsize=8)
    axes[1].set_title('(b) Accumulated precipitation\n'+start.strftime('%m/%d/%Y %H:%M')+' -\n'+end.strftime('%m/%d/%Y %H:%M UTC'),fontsize=8)
    axes[2].set_title('(c)')
    f,metadata=draw(profiles,selected('argentina_20191028'),map_axes=axes[:2])
    profiles[0].set_title('(d) Temperature');profiles[1].set_title('(e) Water vapor')
    f.to_csv(OUT/'data/S2019_profiles.csv.gz',index=False)
    write_json(OUT/'data/S2019_metadata.json',{'profile_sources':metadata,'map_omega':'native signed omega; negative=ascent; mean 500/700 hPa',
        'omega_snapshot_range_UTC':[str(clock[0]),str(clock[-1])],'precipitation_interval_UTC':[c['start'],c['end']],
        'date_display':'MM/DD/YYYY HH:MM UTC','time_axis':'day and HH:MM; 90 degrees; no month',
        'sondes_used_as_initial_or_boundary_conditions':False})
    fig.tight_layout();save(fig,OUT/'supplementary/S2019_20191028_29_ERA5_fixed_support_ascent_precipitation')


def predictions(stats,observed,coupling='C0',supplement=False,symlog=False):
    fig,axes=plt.subplots(3,1,figsize=(9.6,7.2),sharex=True)
    rows=[]
    for i,(ax,method) in enumerate(zip(axes,METHODS)):
        ax.step(observed.time,observed['count'],where='post',c='black',lw=1.2,label='ENTLN')
        for cfg,color,label in [('TH_'+coupling,BLUE,'Thompson'),('WS_'+coupling,RED,'WSM6')]:
            f=stats[(stats.method==method)&(stats.configuration==cfg)].sort_values('time')
            ax.plot(f.time,f['mean'],c=color,lw=1.1,label=label)
            ax.fill_between(f.time,np.maximum(0,f['mean']-f['std']),f['mean']+f['std'],color=color,alpha=.18,lw=0)
            rows.append(f)
        ax.set(ylabel='Flashes min$^{-1}$',title=f'({chr(97+i)}) {LABELS[method]}');ax.set_ylim(bottom=0);utc_axis(ax)
        if symlog:ax.set_yscale('symlog',linthresh=1)
        ax.legend(loc='upper right',ncol=3,frameon=False)
    axes[-1].set_xlabel('UTC');fig.tight_layout()
    name=('supplementary/S03_'+coupling if supplement else 'figures/Figure_03')+'_20181213_14_independent_20191028_TH_C0_calibration_mean_std40'
    if symlog:name+='_symlog_linthresh1'
    save(fig,OUT/name)
    pd.concat(rows).to_csv(OUT/'data'/('Figure_03_'+coupling+'_plotted_statistics.csv'),index=False)


def performance():
    table=pd.read_csv(OUT/'tables/Table_04_independent_evaluation.csv')
    contrast=pd.read_csv(OUT/'tables/Table_05_microphysics_sensitivity.csv')
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.2))
    y=np.arange(3)
    for cfg,color,offset,label in [('TH_C0',BLUE,-.1,'Thompson'),('WS_C0',RED,.1,'WSM6')]:
        f=table[table.configuration==cfg].set_index('method').reindex(METHODS)
        axes[0,0].scatter(f.rmse,y+offset,c=color,label=label)
        for metric,marker,label2 in [('skill_zero','o','zero'),('skill_constant','^','constant')]:
            axes[0,1].scatter(f[metric],y+offset,c=color,marker=marker,label=label+' / '+label2)
    for mode,marker in [('C0','o'),('CL','^'),('CF','s')]:
        values=contrast[(contrast.metric=='rmse')&(contrast.configuration_a=='WS_'+mode)&(contrast.configuration_b=='TH_'+mode)].set_index('method').reindex(METHODS)
        axes[1,0].scatter(values.difference,y,marker=marker,label=mode)
    positions=[(m,mode) for m in METHODS for mode in ['CL','CF']]
    for scheme,color,offset,label in [('TH',BLUE,-.1,'Thompson'),('WS',RED,.1,'WSM6')]:
        values=[]
        for method,mode in positions:
            row=contrast[(contrast.method==method)&(contrast.metric=='rmse')&(contrast.configuration_a==scheme+'_'+mode)&(contrast.configuration_b==scheme+'_C0')]
            values.append(row.difference.iloc[0])
        axes[1,1].scatter(values,np.arange(6)+offset,c=color,label=label)
    for ax in [axes[0,0],axes[0,1],axes[1,0]]:ax.set_yticks(y,[LABELS[m] for m in METHODS])
    axes[1,1].set_yticks(np.arange(6),[LABELS[m]+' / '+mode+'−C0' for m,mode in positions])
    for ax,title,xlabel in zip(axes.flat,['(a)','(b)','(c) WSM6 − Thompson','(d)'],['RMSE (flashes/min)','MSE skill','ΔRMSE (flashes/min)','ΔRMSE (flashes/min)']):
        ax.set(title=title,xlabel=xlabel);ax.invert_yaxis();ax.grid(axis='x',alpha=.2);ax.legend(loc='best',frameon=False,fontsize=7)
        if ax is not axes[0,0]:ax.axvline(0,color='0.3',lw=.7)
    axes[0,0].set_xscale('log')
    reference=table[table.configuration.isin(['TH_C0','WS_C0'])]
    axes[0,0].set_xlim(reference.rmse.min()*.7,reference.rmse.max()*1.6)
    axes[0,1].set_xscale('symlog',linthresh=.01)
    skills=reference[['skill_zero','skill_constant']].to_numpy()
    axes[0,1].set_xlim(min(-.01,np.nanmin(skills)*1.6),max(.05,np.nanmax(skills)+.05))
    for ax in axes[1]:ax.set_xscale('symlog',linthresh=1);ax.margins(x=.08)
    fig.tight_layout();save(fig,OUT/'figures/Figure_04_20181213_14_fixed_2019_coefficients_RMSE_skills_physics_contrasts')
    table.to_csv(OUT/'data/Figure_04_plotted_metrics.csv',index=False)
    contrast.to_csv(OUT/'data/Figure_04_plotted_contrasts.csv',index=False)


def precipitation():
    path=ROOT/'results/argentina_20181213/meteorological_validation/precipitation_comparison_ERA5_hourly.csv'
    if not path.exists():return
    f=pd.read_csv(path);f=f[f.configuration.isin(['TH_C0','WS_C0'])]
    f.time=pd.to_datetime(f.time,utc=True);f.end_time=pd.to_datetime(f.end_time,utc=True)
    fig,axes=plt.subplots(2,1,figsize=(9.6,6.2),sharex=True)
    ref=f[f.configuration=='TH_C0'].sort_values('time')
    axes[0].step(ref.time,ref.ERA5_mm,where='post',c='black',label='ERA5')
    axes[1].plot(pd.concat([ref.time.iloc[:1],ref.end_time]),np.r_[0,ref.ERA5_mm.cumsum()],c='black',label='ERA5')
    for cfg,color,label in [('TH_C0',BLUE,'Thompson'),('WS_C0',RED,'WSM6')]:
        s=f[f.configuration==cfg].sort_values('time')
        axes[0].step(s.time,s.simulation_mm,where='post',color=color,label=label)
        axes[1].plot(pd.concat([s.time.iloc[:1],s.end_time]),np.r_[0,s.simulation_mm.cumsum()],c=color,label=label)
    for ax,title,unit in zip(axes,['(a)','(b)'],['mm h$^{-1}$','mm']):
        ax.set(title=title,ylabel=unit);ax.legend(loc='upper left',frameon=False,ncol=3);utc_axis(ax)
    axes[1].set_xlabel('UTC');fig.tight_layout()
    save(fig,OUT/'figures/Figure_06_20181213_14_precipitation_ERA5_hourly_C0')
    f.to_csv(OUT/'data/Figure_06_plotted_precipitation.csv',index=False)


def cape_figure():
    minutes=1 if (OUT/'minute_mode.json').exists() else 15
    path=OUT/('data/CAPE_precipitation_controls_1min.csv.gz' if minutes==1 else 'data/CAPE_precipitation_uncalibrated_15min.csv.gz')
    if not path.exists():return
    f=pd.read_csv(path);f.time=pd.to_datetime(f.time,utc=True);f=f[(f.member==0)&f.configuration.isin(['TH_C0','WS_C0'])]
    fig,axes=plt.subplots(3,1,figsize=(9.6,6.7),sharex=True)
    for cfg,color,label in [('TH_C0',BLUE,'Thompson'),('WS_C0',RED,'WSM6')]:
        s=f[f.configuration==cfg].sort_values('time')
        for ax,column in zip(axes,['cape_J_kg','precipitation_amount_kg_m2','CAPExP_J_m2']):ax.plot(s.time,s[column],c=color,label=label,lw=1.1)
    for ax,title,unit in zip(axes,['(a) CAPE','(b) Precipitation','(c) CAPE × P'],['J kg$^{-1}$',f'kg m$^{{-2}}$ / {minutes} min',f'J m$^{{-2}}$ / {minutes} min']):
        ax.set(title=title,ylabel=unit);utc_axis(ax);ax.legend(loc='upper right',frameon=False)
    fig.tight_layout();save(fig,OUT/f'supplementary/S_CAPE_20181213_14_{minutes}min_uncalibrated_most_cape_pseudoadiabatic')


def main():
    active=OUT/'independent_manifest.json'
    if active.exists() and json.loads(active.read_text()).get('version')==3:
        raise ValueError('Historical plotting driver cannot relabel v3; use scripts/plot.py article')
    ap=argparse.ArgumentParser();ap.add_argument('--available-only',action='store_true');ap.add_argument('--skip-supplementary',action='store_true');a=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'data').mkdir(parents=True,exist_ok=True)
    c=load_config(EVALUATION_CONFIG)
    if (OUT/'minute_mode.json').exists():c['interval_minutes']=1
    minutes=c['interval_minutes']
    o=pd.read_csv(ROOT/f'results/argentina_20181213/observations/counts_{minutes}min.csv');o.time=pd.to_datetime(o.time,utc=True)
    with plt.rc_context(STYLE):
        environment(c);observations(c,o);microphysical_response(c);precipitation();calibration_weather();cape_figure()
        if (OUT/'tables/Table_04_independent_evaluation.csv').exists():
            stats=pd.read_csv(OUT/'data/ensemble_statistics.csv');stats.time=pd.to_datetime(stats.time,utc=True)
            predictions(stats,o);performance()
            if not a.skip_supplementary:
                for mode in ['CL','CF']:predictions(stats,o,coupling=mode,supplement=True)
                predictions(stats,o,supplement=True,symlog=True)
        elif not a.available_only:
            raise FileNotFoundError('Independent 2019 calibration/evaluation missing; Figures 3 and 4 are not generated from old coefficients')
    write_json(OUT/'figure_manifest.json',{'calibration_event':'argentina_20191028','evaluation_event':'argentina_20181213',
        'legends':'inside axes','band':'mean +/- sample standard deviation, members 1..40; lower edge clipped for display only',
        'principal_methods':METHODS,'files':[{'path':str(p.relative_to(OUT)),'sha256':sha256(p)} for p in (OUT/'figures').glob('*') if p.is_file()],
        'Figure_04_axes':'RMSE log; skills symlog linthresh .01; signed RMSE contrasts symlog linthresh 1; no values transformed in CSV',
        'interval_minutes':minutes,'rate_units':'flashes/min','Figure_06_reference':'ERA5 hourly precipitation, same fixed area; not independent observations',
        'missing_figures_3_4':not (OUT/'tables/Table_04_independent_evaluation.csv').exists()})


if __name__=='__main__':main()
