"""Principal figures bound to frozen independent-event CSVs; legends inside axes."""
import argparse
import json
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
    root=ROOT/'results/argentina_20191028/support'
    with xr.open_dataset(root/'meteorology.nc') as source:d=source.load()
    area=c['area'];west,south,east,north=area
    fig,axes=plt.subplots(1,3,figsize=(12.5,4.2))
    for ax,field,unit in [(axes[0],d.positive_ascent.mean('time'),'Pa s$^{-1}$'),(axes[1],1000*d.precipitation_m.sum('time'),'mm')]:
        field.plot(ax=ax,cbar_kwargs={'label':unit},add_labels=False)
        boundaries(ax);ax.add_patch(Rectangle((west,south),east-west,north-south,fill=False,ec='black',lw=1.5))
        ax.set(xlim=(float(field.longitude.min()),float(field.longitude.max())),ylim=(float(field.latitude.min()),float(field.latitude.max())),xlabel='Longitude',ylabel='Latitude')
    tracks=pd.read_csv(root/'hourly_meteorological_evolution.csv')
    axes[1].plot(tracks.longitude,tracks.latitude,'k.',ms=2)
    selected=d.sel(longitude=slice(west,east),latitude=slice(south,north))
    clock=pd.to_datetime(d.time.values,utc=True)
    axes[2].plot(clock,selected.positive_ascent.mean(('latitude','longitude')),c=BLUE,label='Ascent')
    rainax=axes[2].twinx();rainax.plot(clock,1000*selected.precipitation_m.mean(('latitude','longitude')),c=RED,label='Rain')
    axes[2].set_ylabel('Pa s$^{-1}$');rainax.set_ylabel('mm h$^{-1}$');utc_axis(axes[2])
    handles,labels=axes[2].get_legend_handles_labels();h,l=rainax.get_legend_handles_labels()
    axes[2].legend(handles+h,labels+l,loc='upper left',frameon=False)
    for i,ax in enumerate(axes):ax.set_title('('+chr(97+i)+')')
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
        ax.set(ylabel='Flashes / 15 min',title=f'({chr(97+i)}) {LABELS[method]}');ax.set_ylim(bottom=0);utc_axis(ax)
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
    for ax,title,xlabel in zip(axes.flat,['(a)','(b)','(c) WSM6 − Thompson','(d)'],['RMSE (flashes / 15 min)','MSE skill','ΔRMSE (flashes / 15 min)','ΔRMSE (flashes / 15 min)']):
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
    path=ROOT/'results/argentina_20181213/meteorological_validation/precipitation_comparison_daily.csv'
    if not path.exists():return
    f=pd.read_csv(path);f=f[f.configuration.isin(['TH_C0','WS_C0'])]
    days=sorted(f.time.unique());x=np.arange(len(days));fig,ax=plt.subplots(figsize=(7.6,4.1))
    observed=f.groupby('time').external_estimate_mm.first().reindex(days)
    ax.bar(x-.24,observed,width=.24,color='black',label='CHIRPS2')
    for cfg,color,offset,label in [('TH_C0',BLUE,0,'Thompson'),('WS_C0',RED,.24,'WSM6')]:
        ax.bar(x+offset,f[f.configuration==cfg].set_index('time').simulation_mm.reindex(days),width=.24,color=color,label=label)
    ax.set(xticks=x,xticklabels=[pd.Timestamp(day).strftime('%d %b') for day in days],ylabel='mm / day')
    ax.legend(loc='upper right',frameon=False);ax.grid(axis='y',alpha=.2);fig.tight_layout()
    save(fig,OUT/'figures/Figure_06_20181213_14_daily_precipitation_CHIRPS2_CFSv2_disaggregated_C0')
    f.to_csv(OUT/'data/Figure_06_plotted_precipitation.csv',index=False)


def cape_figure():
    path=OUT/'data/CAPE_precipitation_uncalibrated_15min.csv.gz'
    if not path.exists():return
    f=pd.read_csv(path);f.time=pd.to_datetime(f.time,utc=True);f=f[(f.member==0)&f.configuration.isin(['TH_C0','WS_C0'])]
    fig,axes=plt.subplots(3,1,figsize=(9.6,6.7),sharex=True)
    for cfg,color,label in [('TH_C0',BLUE,'Thompson'),('WS_C0',RED,'WSM6')]:
        s=f[f.configuration==cfg].sort_values('time')
        for ax,column in zip(axes,['cape_J_kg','precipitation_amount_kg_m2','CAPExP_J_m2']):ax.plot(s.time,s[column],c=color,label=label,lw=1.1)
    for ax,title,unit in zip(axes,['(a) CAPE','(b) Precipitation','(c) CAPE × P'],['J kg$^{-1}$','kg m$^{-2}$ / 15 min','J m$^{-2}$ / 15 min']):
        ax.set(title=title,ylabel=unit);utc_axis(ax);ax.legend(loc='upper right',frameon=False)
    fig.tight_layout();save(fig,OUT/'supplementary/S_CAPE_20181213_14_uncalibrated_most_cape_pseudoadiabatic')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--available-only',action='store_true');ap.add_argument('--skip-supplementary',action='store_true');a=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'data').mkdir(parents=True,exist_ok=True)
    c=load_config(EVALUATION_CONFIG)
    o=pd.read_csv(ROOT/'results/argentina_20181213/observations/counts_15min.csv');o.time=pd.to_datetime(o.time,utc=True)
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
        'missing_figures_3_4':not (OUT/'tables/Table_04_independent_evaluation.csv').exists()})


if __name__=='__main__':main()
