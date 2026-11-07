"""Versioned scientific figures; own/common coefficients never pooled."""
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.calibration.revision import METHODS,MU_METHOD
from scm_lightning.evaluation.revision import OUT,PREVIOUS,diagnostics_root
from scm_lightning.plotting.article import STYLE,BLUE,RED,save,utc_axis

LABELS={MU_METHOD:'MUCAPE × P','IxG':'I × G','McCaulF3_LS_1D_adapted':'McCaul F3 — LS'}


def predictions(base,comparison,coupling='C0',symlog=False):
    stats=pd.read_csv(base/'data/ensemble_statistics.csv');stats.time=pd.to_datetime(stats.time,utc=True)
    if set(stats.method)!=set(METHODS):raise ValueError('v3 figures require validated MUCAPE/OLS statistics')
    o=pd.read_csv(base/'data/observations_1min.csv');o.time=pd.to_datetime(o.time,utc=True)
    fig,axes=plt.subplots(3,1,figsize=(9.6,7.3),sharex=True);rows=[]
    for ax,method,letter in zip(axes,METHODS,'abc'):
        ax.step(o.time,o['count'],where='post',color='black',lw=1.1,label='ENTLN')
        for scheme,color,label in [('TH',BLUE,'Thompson'),('WS',RED,'WSM6')]:
            f=stats[(stats.method==method)&(stats.configuration==scheme+'_'+coupling)].sort_values('time')
            ax.plot(f.time,f['mean'],c=color,label=label,lw=1.1)
            ax.fill_between(f.time,np.maximum(0,f['mean']-f['std']),f['mean']+f['std'],color=color,alpha=.18,lw=0)
            rows.append(f)
        ax.set(ylabel='Flashes min$^{-1}$',title=f'({letter}) {LABELS[method]}');ax.set_ylim(bottom=0);utc_axis(ax)
        if symlog:ax.set_yscale('symlog',linthresh=1)
        ax.legend(loc='upper right',ncol=3,frameon=False)
    fig.suptitle(comparison,fontsize=10);axes[-1].set_xlabel('UTC');fig.tight_layout()
    prefix='figures/Figure_03' if coupling=='C0' and not symlog else 'supplementary/S03_'+coupling
    name=prefix+'_20181213_14_MUCAPE_OLS_mean_std40'+('_symlog' if symlog else '')
    save(fig,base/name);pd.concat(rows).to_csv(base/f'data/Figure_03_{coupling}_plotted_statistics.csv',index=False)


def contrast_performance(base,comparison):
    table=pd.read_csv(base/'tables/Table_04_independent_evaluation.csv')
    from scm_lightning.evaluation.revision import contrasts
    contrast=contrasts(table,comparison)
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.2));y=np.arange(3)
    for cfg,color,offset,label in [('TH_C0',BLUE,-.1,'Thompson'),('WS_C0',RED,.1,'WSM6')]:
        f=table[table.configuration==cfg].set_index('method').reindex(METHODS)
        axes[0,0].scatter(f.rmse,y+offset,c=color,label=label)
        for key,marker,label2 in [('skill_zero','o','zero'),('skill_constant','^','constant')]:
            axes[0,1].scatter(f[key],y+offset,c=color,marker=marker,label=label+' / '+label2)
    for mode,marker in [('C0','o'),('CL','^'),('CF','s')]:
        f=contrast[(contrast.metric=='rmse')&(contrast.configuration_a=='WS_'+mode)&(contrast.configuration_b=='TH_'+mode)].set_index('method').reindex(METHODS)
        axes[1,0].scatter(f.difference,y,marker=marker,label=mode)
    positions=[(m,mode) for m in METHODS for mode in ['CL','CF']]
    for scheme,color,offset,label in [('TH',BLUE,-.1,'Thompson'),('WS',RED,.1,'WSM6')]:
        values=[contrast[(contrast.method==m)&(contrast.metric=='rmse')&(contrast.configuration_a==scheme+'_'+mode)&(contrast.configuration_b==scheme+'_C0')].difference.iloc[0] for m,mode in positions]
        axes[1,1].scatter(values,np.arange(6)+offset,c=color,label=label)
    for ax in [axes[0,0],axes[0,1],axes[1,0]]:ax.set_yticks(y,[LABELS[m] for m in METHODS])
    axes[1,1].set_yticks(np.arange(6),[LABELS[m]+' / '+mode+'−C0' for m,mode in positions])
    for ax,title,label in zip(axes.flat,['(a)','(b)','(c) WSM6 − Thompson','(d) Coupling − C0'],['RMSE (flashes/min)','MSE skill','ΔRMSE (flashes/min)','ΔRMSE (flashes/min)']):
        ax.set(title=title,xlabel=label);ax.invert_yaxis();ax.grid(axis='x',alpha=.2);ax.legend(loc='best',frameon=False,fontsize=7)
        if ax is not axes[0,0]:ax.axvline(0,c='0.3',lw=.7)
    axes[0,0].set_xscale('log');axes[0,1].set_xscale('symlog',linthresh=.01)
    for ax in axes[1]:ax.set_xscale('symlog',linthresh=1);ax.margins(x=.08)
    fig.suptitle(comparison,fontsize=10);fig.tight_layout()
    save(fig,base/'supplementary/S04_RMSE_skills_microphysics_coupling_contrasts')
    table.to_csv(base/'data/Figure_04_plotted_metrics.csv',index=False);contrast.to_csv(base/'data/Figure_04_plotted_contrasts.csv',index=False)


def performance(base,comparison):
    table=pd.read_csv(base/'tables/Table_04_independent_evaluation.csv')
    panels=[('rmse','RMSE (flashes/min)','log'),('mae','MAE (flashes/min)','log'),
        ('correlation','Correlation','linear'),('bias','Bias (flashes/min)','symlog'),
        ('skills','MSE skill','symlog'),('relative_peak_error','Relative peak amplitude error','symlog'),
        ('peak_time_error_minutes','Peak time error (min)','linear'),('accumulated_error','Accumulated error (flashes)','symlog')]
    fig,axes=plt.subplots(4,2,figsize=(11.8,12.3));y=np.arange(3)
    for ax,(metric,label,scale),letter in zip(axes.flat,panels,'abcdefgh'):
        for cfg,color,offset,name in [('TH_C0',BLUE,-.1,'Thompson'),('WS_C0',RED,.1,'WSM6')]:
            f=table[table.configuration==cfg].set_index('method').reindex(METHODS)
            keys=[('pearson','o','Pearson'),('spearman','^','Spearman')] if metric=='correlation' else (
                [('skill_zero','o','zero'),('skill_constant','^','constant')] if metric=='skills' else [(metric,'o','')])
            for key,marker,suffix in keys:ax.scatter(f[key],y+offset,c=color,marker=marker,label=name+(' / '+suffix if suffix else ''))
        ax.set_yticks(y,[LABELS[m] for m in METHODS]);ax.invert_yaxis();ax.set(title=f'({letter})',xlabel=label)
        ax.grid(axis='x',alpha=.2);ax.legend(loc='upper left' if metric in ['correlation','accumulated_error'] else 'best',fontsize=6.5,frameon=False)
        if scale=='symlog':ax.set_xscale(scale,linthresh=.01 if metric=='skills' else (.1 if metric=='relative_peak_error' else 1))
        else:ax.set_xscale(scale)
        if metric not in ['rmse','mae']:ax.axvline(0,c='0.4',lw=.7)
        ax.margins(x=.15,y=.15)
        if metric=='correlation':ax.set_xlim(-1.05,1.05)
    fig.suptitle(comparison,fontsize=10);fig.tight_layout()
    save(fig,base/'figures/Figure_04_20181213_14_MUCAPE_OLS_performance')
    contrast_performance(base,comparison)


def precipitation(event,destination):
    folder=ROOT/'results'/event/'meteorological_validation';f=pd.read_csv(folder/'precipitation_MERGE_ERA5_SCM_hourly.csv')
    f.time=pd.to_datetime(f.time,utc=True);f.end_time=pd.to_datetime(f.end_time,utc=True)
    fig,axes=plt.subplots(2,1,figsize=(9.6,6.4),sharex=True)
    for name,color,style,label in [('MERGE','black','-','MERGE CPTEC/INPE'),('ERA5','#756bb1',':','ERA5'),('TH_C0',BLUE,'-','Thompson C0'),('WS_C0',RED,'--','WSM6 C0')]:
        values=f[name+'_mm'];axes[0].step(pd.concat([f.time,f.end_time.iloc[-1:]]),np.r_[values,values.iloc[-1]],where='post',c=color,ls=style,label=label,lw=1.4)
        accumulated=np.r_[0,values.cumsum(skipna=False)]
        axes[1].plot(pd.concat([f.time.iloc[:1],f.end_time]),accumulated,c=color,ls=style,label=label,lw=1.4)
    for ax,title,unit in zip(axes,['(a) Hourly precipitation','(b) Accumulated precipitation'],['mm h$^{-1}$','mm']):
        ax.set(title=title,ylabel=unit);utc_axis(ax);ax.set_xlim(f.time.iloc[0],f.end_time.iloc[-1])
    fig.legend(*axes[0].get_legend_handles_labels(),ncol=4,loc='upper center',frameon=False,bbox_to_anchor=(.5,1.015))
    axes[1].set_xlabel('UTC');fig.tight_layout(rect=(0,0,1,.96));save(fig,destination)
    f.to_csv(destination.parent/(destination.name+'_plotted.csv'),index=False)


def cape_figures():
    for cfg in ['TH_C0','WS_C0']:
        path=ROOT/'results/argentina_20181213/thermodynamic_revision'/f'CAPE_profile_transition_diagnostics_{cfg}.csv.gz'
        f=pd.read_csv(path);f.time=pd.to_datetime(f.time,utc=True)
        if 'SBCAPE_J_kg' not in f:raise ValueError('Mandatory control definition comparison missing')
        fig,axes=plt.subplots(6,1,figsize=(10,11),sharex=True)
        for label,color in [('SBCAPE','#636363'),('previous_most_cape','#e6550d'),('MUCAPE',BLUE)]:
            cape='raw_cape_J_kg' if label=='MUCAPE' else label+'_J_kg'
            cin='cin_J_kg' if label=='MUCAPE' else ('SBCIN_J_kg' if label=='SBCAPE' else 'previous_most_cin_J_kg')
            axes[0].plot(f.time,f[cape],c=color,label=label+' (raw)' if label=='MUCAPE' else label,lw=1)
            axes[1].plot(f.time,f[cin],c=color,label=label,lw=1)
            product='CAPExP_J_m2' if label=='MUCAPE' else label+'_P_J_m2'
            axes[5].plot(f.time,f[product],c=color,label=label,lw=1)
            axes[3].plot(f.time,f[label+'_LFC_Pa']/100,c=color,label=label+' LFC',lw=.8)
            axes[3].plot(f.time,f[label+'_EL_Pa']/100,c=color,ls='--',label=label+' EL',lw=.8)
        invalid=f.raw_cape_J_kg.lt(0)
        axes[0].scatter(f.loc[invalid,'time'],f.loc[invalid,'raw_cape_J_kg'],s=9,c='black',marker='x',label='Invalid for prediction',zorder=5)
        axes[2].plot(f.time,f.parcel_pressure_Pa/100,c=BLUE,label='MU parcel pressure')
        height=axes[2].twinx();height.plot(f.time,f.parcel_height_AGL_m,c='0.5',ls=':',lw=.8);height.set_ylabel('Parcel height (m AGL)')
        axes[4].plot(f.time,f.precipitation_amount_kg_m2,c='black',lw=1)
        for ax,unit,title in zip(axes,['J kg$^{-1}$','J kg$^{-1}$','hPa','hPa','mm / min','J m$^{-2}$ / min'],
            ['(a) CAPE','(b) CIN','(c) Selected MUCAPE parcel','(d) LFC / EL','(e) Precipitation','(f) CAPE × P']):
            ax.set(ylabel=unit,title=title);utc_axis(ax);ax.axvline(pd.Timestamp('2018-12-13T21:00Z'),color='0.4',ls=':',lw=.8)
        for ax in [axes[0],axes[1],axes[5]]:ax.legend(loc='upper right',ncol=3,fontsize=6,frameon=False)
        axes[3].legend(loc='upper center',bbox_to_anchor=(.5,-.12),ncol=3,fontsize=6,frameon=False)
        fig.suptitle(cfg+' — native profiles; definitions selected without lightning observations',fontsize=10)
        axes[-1].set_xlabel('UTC');fig.tight_layout();save(fig,OUT/f'supplementary/S_CAPE_definition_comparison_{cfg}')
        zoom=f[(f.time>=pd.Timestamp('2018-12-13T20:00Z'))&(f.time<=pd.Timestamp('2018-12-13T22:00Z'))]
        fig,ax=plt.subplots(figsize=(9.6,3.8))
        for name,color,label in [('SBCAPE_J_kg','0.4','SBCAPE'),('previous_most_cape_J_kg','#e6550d','Previous most_cape'),('raw_cape_J_kg',BLUE,'MUCAPE (raw)')]:
            ax.plot(zoom.time,zoom[name],c=color,label=label,lw=1.2)
        ax.set(ylabel='J kg$^{-1}$',title=cfg+' — 12/13/2018 20:00–22:00 UTC');ax.legend(frameon=False);ax.grid(alpha=.2)
        import matplotlib.dates as mdates
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'));fig.tight_layout();save(fig,OUT/f'supplementary/S_CAPE_21UTC_zoom_{cfg}')
        f.to_csv(OUT/f'data/CAPE_definition_comparison_{cfg}.csv.gz',index=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('--meteorology-only',action='store_true');a=p.parse_args()
    with plt.rc_context(STYLE):
        if a.meteorology_only:
            for event in ['argentina_20181213','argentina_20191028']:
                precipitation(event,ROOT/'results'/event/'meteorological_validation/Figure_precipitation_MERGE_ERA5_SCM')
            return
        manifest=OUT/'evaluation_manifest.json'
        if not manifest.exists() or json.loads(manifest.read_text()).get('version')!=3:
            raise ValueError('Complete v3 evaluation before generating principal lightning figures')
        predictions(OUT,'Microphysics-specific 2019 calibration — six frozen coefficients')
        performance(OUT,'Microphysics-specific 2019 calibration — six frozen coefficients')
        for mode in ['CL','CF']:predictions(OUT,'Microphysics-specific 2019 calibration',coupling=mode)
        predictions(OUT,'Microphysics-specific 2019 calibration',symlog=True)
        base=OUT/'supplementary/common_TH_same_MUCAPE_OLS'
        predictions(base,'Common Thompson C0 2019 coefficients — MUCAPE / OLS definitions')
        performance(base,'Common Thompson C0 2019 coefficients — MUCAPE / OLS definitions')
        precipitation('argentina_20181213',OUT/'figures/Figure_06_20181213_14_MERGE_ERA5_precipitation_C0')
        precipitation('argentina_20191028',OUT/'supplementary/S2019_MERGE_ERA5_precipitation_C0')
        cape_figures()
    # Superseded principal graphics live in the hash-verified v2 archive.
    for pattern in ['Figure_03*independent_20191028*','Figure_04*fixed_2019*','Figure_06*ERA5_hourly*']:
        for path in (OUT/'figures').glob(pattern):
            if not (PREVIOUS/'figures'/path.name).exists():raise IOError('Old principal figure not preserved')
            path.unlink()
    write_json(OUT/'figure_manifest.json',{'version':3,'interval_minutes':1,'rate_units':'flashes/min',
        'principal_methods':METHODS,'main_coefficients':'microphysics-specific six coefficients fitted only with 2019',
        'supplementary_comparisons':['common_TH_same_MUCAPE_OLS','legacy_original_TH_v2'],
        'Figure_06_reference':'MERGE hourly observational estimate; ERA5 remains in both panels as additional comparison',
        'Figure_05_native_results_preserved':True,'band':'mean +/- sample std of 40 members; member0 excluded; lower band clipped only for display',
        'files':[{'path':str(f.relative_to(OUT)),'sha256':sha256(f)} for folder in ['figures','supplementary'] for f in (OUT/folder).rglob('*') if f.suffix in ['.png','.pdf']]})


if __name__=='__main__':main()
