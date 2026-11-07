"""Publication exports backed by audited CSVs; no physical/model changes."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr
from scm_lightning.config import load_config
from scm_lightning.provenance import ROOT,sha256,write_json,git
from scm_lightning.forcing import normalize
from scm_lightning.forcing.contracts import area_km2
from scm_lightning.ensemble.experiments import CONFIGURATIONS,control_directory
from scm_lightning.evaluation.article import OUT,CASE,METHODS
from scm_lightning.plotting.figures import boundaries,microphysics

BLUE='#2166ac';RED='#b2182b'
LABELS={'CAPExP':'CAPExP','IxG':'IxG','McCaulF3':'McCaulF3 (1D, prescribed w)'}
STYLE={'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':10,
       'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,
       'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','timezone':'UTC'}


def save(fig,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    candidates=[]
    for extension in ['.pdf','.png']:
        temporary=path.with_suffix('.pending'+extension)
        fig.savefig(temporary,bbox_inches='tight',**({'dpi':300} if extension=='.png' else {}))
        candidates.append((temporary,path.with_suffix(extension)))
    for temporary,destination in candidates:
        if destination.exists() and sha256(temporary)==sha256(destination):
            temporary.unlink();continue
        for attempt in range(6):
            try:temporary.replace(destination);break
            except OSError:
                if attempt==5:raise
                time.sleep(.2)
    plt.close(fig)


def utc_axis(ax):
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=[0,6,12,18]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b\n%H:%M'))
    ax.grid(axis='y',alpha=.2)


def load_data():
    audit=json.loads((OUT/'audit/current_statistics.json').read_text())
    if audit['status']!='passed':raise ValueError('Statistical audit must pass before figure generation')
    stats=pd.read_csv(OUT/'data/ensemble_statistics_A.csv');stats.time=pd.to_datetime(stats.time,utc=True)
    obs=pd.read_csv(OUT/'data/observations_15min.csv');obs.time=pd.to_datetime(obs.time,utc=True)
    return load_config(CASE),stats,obs


def environment(c):
    d=normalize(xr.load_dataset(ROOT/'data/raw/era5/reanalysis_pl_20181214.nc'))
    w=d.w.sel(time='2018-12-14T06:00',level=500)
    assert np.isfinite(w).all()
    w.to_dataframe(name='omega_Pa_s').reset_index().to_csv(OUT/'data/Figure_01_omega.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(12.8,5),gridspec_kw={'width_ratios':[1.45,1,1]})
    limit=float(abs(w).max())
    mesh=axes[0].pcolormesh(w.longitude,w.latitude,w,cmap='RdBu_r',vmin=-limit,vmax=limit,shading='auto')
    boundaries(axes[0])
    west,south,east,north=c['area']
    axes[0].add_patch(Rectangle((west,south),east-west,north-south,fill=False,ec='black',lw=1.7,label='Fixed SCM area',zorder=6))
    station=c['station'];lon,lat=station['longitude'],station['latitude']
    axes[0].scatter([lon],[lat],marker='^',c='black',s=30,zorder=7,label='SACO / WMO 87344')
    axes[0].annotate('SACO',(lon,lat),xytext=(-1,10),textcoords='offset points',ha='center',fontsize=8,zorder=7)
    axes[0].set(xlim=(float(w.longitude.min()),float(w.longitude.max())),ylim=(float(w.latitude.min()),float(w.latitude.max())),
                xlabel='Longitude (°E)',ylabel='Latitude (°N)',title='(a) ERA5, 14 Dec 2018 06 UTC')
    axes[0].set_aspect(1/np.cos(np.deg2rad((south+north)/2)))
    cb=fig.colorbar(mesh,ax=axes[0],shrink=.87,pad=.025);cb.set_label(r'$\omega$ (Pa s$^{-1}$); negative = ascent')
    if cb.solids is not None:cb.solids.set_rasterized(False)
    axes[0].legend(loc='lower left',fontsize=7)
    labels=[]
    from scm_lightning.observations.profiles import selected
    from scm_lightning.plotting.profiles import draw
    records=selected('argentina_20181213')
    frames,profile_metadata=draw(axes[1:],records,map_axes=[axes[0]])
    frames.to_csv(OUT/'data/Figure_01_campaign_profiles.csv.gz',index=False)
    axes[1].set(title='(b) Temperature');axes[2].set(title='(c) Water vapor')
    for ax in axes[1:]:ax.grid(alpha=.2)
    from scm_lightning.observations import haversine_km
    distance=float(haversine_km(lat,lon,(south+north)/2,(west+east)/2))
    handles,labels=axes[1].get_legend_handles_labels()
    axes[1].legend(handles,labels,loc='upper right',frameon=False,fontsize=6.5,ncol=2)
    fig.tight_layout();save(fig,OUT/'figures/Figure_01_case_environment')
    write_json(OUT/'data/Figure_01_metadata.json',{'SACO_to_box_center_km':distance,'station':station,'profile_labels':labels,
        'ERA5_source_sha256':sha256(ROOT/'data/raw/era5/reanalysis_pl_20181214.nc'),
        'omega_sign':'native Pa/s; negative=ascent','full_data_extent':[float(w.longitude.min()),float(w.longitude.max()),float(w.latitude.min()),float(w.latitude.max())],
        'sounding_scope':'spatially distinct reference; not observations collocated with the model box',
        'profile_sources':profile_metadata,'preferred_source':'CSU_RELAMPAGO_in_projetoCHUVA; Wyoming retained in normalized inventory',
        'used_for_initial_or_boundary_conditions':False})


def observations(c,o):
    minutes=c.get('interval_minutes',15)
    bytype=pd.read_csv(ROOT/f'results/argentina_20181213/observations/counts_by_type_{minutes}min.csv')
    bytype.time=pd.to_datetime(bytype.time,utc=True)
    assert bytype.covered.all() and c['flag_mapping_confirmed']
    sub=bytype[bytype.type.isin(['IC','CG'])].pivot(index='time',columns='type',values='count').reindex(o.time)
    np.testing.assert_allclose(sub.sum(axis=1),o['count'])
    fig,axes=plt.subplots(2,1,figsize=(9.5,6),sharex=True)
    axes[0].step(o.time,o['count'],where='post',c='black',lw=1.3,label='Total')
    for name,color in [('IC','#666666'),('CG','#b4b4b4')]:
        axes[0].step(o.time,sub[name].values,where='post',c=color,lw=.8,label=name)
        axes[1].plot(o.time+pd.Timedelta(minutes=minutes),sub[name].cumsum(),c=color,lw=.9,label=name)
    peak=int(np.argmax(o['count'].values));stamp=o.time.iloc[peak]
    axes[0].scatter([stamp],[o['count'].iloc[peak]],c='black',s=22,zorder=5)
    axes[0].annotate(stamp.strftime('Peak interval: %d Dec %H:%M UTC'),(stamp,o['count'].iloc[peak]),xytext=(10,-20),textcoords='offset points',fontsize=8)
    axes[1].plot(o.time+pd.Timedelta(minutes=minutes),o['count'].cumsum(),c='black',lw=1.4,label='Total')
    axes[0].set(title='(a) Reconstructed ENTLN flashes',ylabel='Flashes min$^{-1}$' if minutes==1 else 'Flashes / 15 min')
    axes[1].set(title='(b) Accumulated activity',ylabel='Cumulative flashes',xlabel='UTC')
    for ax in axes:utc_axis(ax)
    axes[0].legend(ncol=3,frameon=False,loc='upper left')
    fig.tight_layout();save(fig,OUT/'figures/Figure_02_observed_lightning')
    data=o.copy()
    for name in ['IC','CG']:data[name]=sub[name].values;data['cumulative_'+name]=sub[name].cumsum().values
    data['cumulative_total']=data['count'].cumsum()
    data.to_csv(OUT/'data/Figure_02_observations.csv',index=False)


def predictions(c,stats,o,coupling='C0',full=False,percentiles=False):
    start=pd.Timestamp(c['start'] if full else c['split']['evaluation_start']);end=pd.Timestamp(c['end'])
    source=stats[(stats.configuration.isin(['TH_'+coupling,'WS_'+coupling]))&(stats.time>=start)&(stats.time<end)]
    observed=o[(o.time>=start)&(o.time<end)]
    upper=np.nanmax(np.r_[source['mean'].values+source['std'].values,observed['count'].values])
    if percentiles:upper=max(upper,float(source.p95.max()))
    fig,axes=plt.subplots(3,1,figsize=(10.2,7.6),sharex=True,sharey=True)
    for ax,method,panel in zip(axes,METHODS,'abc'):
        ax.step(observed.time,observed['count'],where='post',color='black',lw=1.1,label='ENTLN total')
        for scheme,color,ls,label in [('TH',BLUE,'-','Thompson'),('WS',RED,'--','WSM6')]:
            f=source[(source.method==method)&(source.configuration==scheme+'_'+coupling)].sort_values('time')
            assert len(f)==len(observed) and (f.n==40).all()
            low=f.p05.values if percentiles else np.maximum(0,f['mean'].values-f['std'].values)
            high=f.p95.values if percentiles else f['mean'].values+f['std'].values
            ax.fill_between(f.time,low,high,color=color,alpha=.17,linewidth=0)
            ax.plot(f.time,f['mean'],color=color,ls=ls,lw=1.25,label=label)
        ax.set(title=f'({panel}) {LABELS[method]}',ylabel='Flashes / 15 min',ylim=(0,upper*1.06),xlim=(start,end))
        if full:
            ax.axvspan(pd.Timestamp(c['start']),pd.Timestamp(c['split']['development_end']),color='.7',alpha=.16)
            ax.axvspan(pd.Timestamp(c['split']['development_end']),pd.Timestamp(c['split']['evaluation_start']),color='.3',alpha=.25)
        utc_axis(ax)
    axes[-1].set_xlabel('UTC')
    handles,labels=axes[0].get_legend_handles_labels()
    band='Bands: member 5–95 percentiles.' if percentiles else 'Bands: ±1 sample SD across 40 perturbed members; lower edge clipped at zero for display only.'
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,frameon=False)
    period='Shaded: development; dark strip: excluded gap.' if full else 'Evaluation: 14 Dec 01:00 to 15 Dec 00:00 UTC (92 intervals).'
    fig.text(.5,.037,f'{coupling}. {period}',ha='center',fontsize=8)
    fig.text(.5,.009,band,ha='center',fontsize=7.8)
    fig.tight_layout(rect=[0,.13,1,1])
    if coupling=='C0' and not full and not percentiles:path=OUT/'figures/Figure_03_ensemble_predictions'
    else:path=OUT/'supplementary'/('S03_ensemble_'+coupling+('_full' if full else '_evaluation')+('_percentiles' if percentiles else '_SD'))
    save(fig,path)
    source.to_csv(OUT/'data'/(path.name+'_statistics.csv'),index=False)


def performance_figure(protocol='A'):
    metric=pd.read_csv(OUT/('tables/Table_03_observational_performance.csv' if protocol=='A' else 'data/fixed_scale_B_performance.csv'))
    contrasts=pd.read_csv(OUT/'data/rmse_contrasts.csv');contrasts=contrasts[contrasts.calibration==protocol]
    fig,axes=plt.subplots(1,3,figsize=(13.8,5.5),gridspec_kw={'width_ratios':[1.1,1,1.1]})
    positions={(method,mode):i for i,(method,mode) in enumerate((m,k) for m in METHODS for k in ['C0','CL','CF'])}
    for scheme,color,offset,label in [('Thompson',BLUE,-.12,'Thompson'),('WSM6',RED,.12,'WSM6')]:
        f=metric[metric.microphysics==scheme]
        axes[0].scatter(f.rmse,[positions[(r.method,r.coupling)]+offset for r in f.itertuples()],c=color,s=25,label=label)
    labels=[method+' / '+mode for method,mode in positions]
    axes[0].set_yticks(range(9),labels);axes[0].set_xlim(left=0);axes[0].set_title('(a) Held-out RMSE')
    f=contrasts[contrasts.kind=='microphysics']
    axes[1].scatter(f.delta_rmse,[positions[(r.method,r.contrast)] for r in f.itertuples()],c='0.2',s=25)
    axes[1].set_yticks(range(9),labels);axes[1].set_title('(b) WSM6 − Thompson')
    position={(method,mode):i for i,(method,mode) in enumerate((m,k) for m in METHODS for k in ['CL','CF'])}
    for scheme,color,offset in [('TH',BLUE,-.12),('WS',RED,.12)]:
        f=contrasts[(contrasts.kind=='coupling')&contrasts.configuration_a.str.startswith(scheme)]
        axes[2].scatter(f.delta_rmse,[position[(r.method,r.configuration_a.split('_')[1])]+offset for r in f.itertuples()],c=color,s=25)
    axes[2].set_yticks(range(6),[m+' / '+k+'−C0' for m,k in position]);axes[2].set_title('(c) MYNN coupling contrasts')
    for i,ax in enumerate(axes):
        ax.invert_yaxis();ax.grid(axis='x',alpha=.2);ax.set_xlabel(('RMSE' if i==0 else 'ΔRMSE')+' (flashes / 15 min)')
        if i:ax.axvline(0,color='black',lw=.7)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.035),ncol=2,frameon=False)
    fig.text(.5,.005,'Calibration '+protocol+': '+('one development-mean coefficient per method/configuration.' if protocol=='A' else 'TH_C0 development-mean coefficients fixed across configurations.')+' All contrasts use the same 92 evaluation intervals.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.11,1,1]);save(fig,OUT/('figures/Figure_04_performance_sensitivity' if protocol=='A' else 'supplementary/S04_fixed_reference_B_performance'))


def microphysical_response(c):
    fig,axes=plt.subplots(3,1,figsize=(9.6,7.1),sharex=True)
    frames=[]
    for configuration,color,ls,label in [('TH_C0',BLUE,'-','Thompson C0'),('WS_C0',RED,'--','WSM6 C0')]:
        f=pd.read_csv(ROOT/'results/argentina_20181213/control_diagnostics'/(configuration+'_paths.csv'))
        with xr.open_dataset(control_directory(configuration)/'output.nc',decode_times=False) as d:
            np.testing.assert_array_equal(f.time_s,d.time.values);rain=d.rainnc.values.copy()
        clock=pd.Timestamp(c['start'])+pd.to_timedelta(f.time_s,unit='s')
        for ax,values in zip(axes,[f.qg,f.qi+f.qs,rain]):ax.plot(clock,values,color=color,ls=ls,lw=1.2,label=label)
        frames.append(pd.DataFrame({'time':clock,'configuration':configuration,'graupel_kg_m2':f.qg,'cloud_ice_snow_kg_m2':f.qi+f.qs,'precipitation_mm':rain}))
    for ax,title,unit in zip(axes,['(a) Graupel path','(b) Cloud ice + snow path','(c) Accumulated precipitation'],[r'kg m$^{-2}$',r'kg m$^{-2}$','mm']):
        ax.set(title=title,ylabel=unit);utc_axis(ax)
    axes[-1].set_xlabel('UTC');axes[0].legend(loc='upper left',ncol=2,frameon=False)
    fig.tight_layout();save(fig,OUT/'figures/Figure_05_microphysical_response')
    pd.concat(frames).to_csv(OUT/'data/Figure_05_control_paths.csv',index=False)


def table_exports(c):
    with xr.open_dataset(control_directory('TH_C0')/'output.nc',decode_times=False) as d:
        z=d.lev.values;dt=float(np.diff(d.time.values)[0])
    edges=np.zeros(len(z)+1)
    for k,value in enumerate(z):edges[k+1]=2*value-edges[k]
    ensemble=ROOT/c['ensemble']['active_directory'];plan=json.loads((ensemble/'manifest.json').read_text())
    n=len({r['member'] for r in plan['runs'] if r['member']!=0});controls=len([r for r in plan['runs'] if r['member']==0])
    rows=[('Case','Argentina, 13-14 December 2018','UTC'),('Simulation period',c['start']+' to '+c['end'],'end exclusive for 15-min bins'),
        ('SCM box (W,S,E,N)',', '.join(str(x) for x in c['area']),'degrees; fixed independently of lightning'),
        ('SCM support area',area_km2(c['area']),'km2'),('Vertical grid',len(z),'layers'),('Model top',edges[-1],'m AGL'),
        ('Layer thickness range',str(np.diff(edges).min())+' to '+str(np.diff(edges).max()),'m; stretched'),
        ('Integration step',dt,'s'),('Microphysics','Thompson; WSM6','two physical schemes'),
        ('MYNN condensate treatment','C0; CL; CF','condensate tendencies off; local only; full'),
        ('MYNN thermodynamic mixing','Enabled in all six configurations','C0 does not disable MYNN'),
        ('Initial-condition members',n,'paired perturbed members per physical configuration'),
        ('Deterministic controls',controls,'one per configuration; excluded from member statistics'),
        ('Primary lightning methods','CAPExP; IxG; McCaulF3','separately calibrated; no method pooling'),
        ('Observed flash aggregation',c['interval_minutes'],'min; total IC+CG'),
        ('Development block',c['start']+' to '+c['split']['development_end'],'96 intervals; end exclusive'),
        ('Excluded gap',c['split']['development_end']+' to '+c['split']['evaluation_start'],'4 intervals; end exclusive'),
        ('Evaluation block',c['split']['evaluation_start']+' to '+c['end'],'92 intervals; end exclusive'),
        ('Primary calibration','A: each method/configuration development ensemble mean','one nonnegative multiplicative scale; same scale on all members'),
        ('Fixed reference (supplement)','TH_C0, predefined','protocol B; no choice by held-out performance')]
    t1=pd.DataFrame(rows,columns=['parameter','value','units_or_convention'])
    t1.to_csv(OUT/'tables/Table_01_experimental_design.csv',index=False)
    t2=pd.DataFrame([
        dict(method='CAPExP',reference='Romps (2019), eq. 1',equation='CAPE * P',variables='Simulated pressure, theta, vapor; total rainnc increment',original_units='W m-2',
             SCM_adaptation='MetPy liquid pseudoadiabatic parcel in lowest simulated cell at 15-min bin start; exact following precipitation amount. Not original isentropic CAPE.',integrated_predictor_units='J m-2',coefficient_units='flash m2 J-1',final_units='flashes per 15 min'),
        dict(method='IxG',reference='Romps (2019), eq. 4',equation='integral (rho_d * qi) * (rho_d * qg) dz',variables='Cloud ice qi; graupel qg; dry-air density rho_d; layer thickness',original_units='kg2 m-5',
             SCM_adaptation='Full-column product of mass concentrations; mixing ratios per kg dry air. Native snapshots trapezoid over each 15-min bin.',integrated_predictor_units='kg2 s m-5',coefficient_units='flash m5 kg-2 s-1',final_units='flashes per 15 min'),
        dict(method='McCaulF3',reference='McCaul et al. (2009), eqs. 3-5',equation='F3 = 0.95 F1 + 0.05 F2; F1 = 0.042 [max(w_LS,0) qg(g/kg)] at 258.15 K; F2 = 0.20 integral rho_d(qi+qs+qg) dz',
             variables='Prescribed ERA5 w_LS; graupel at -15 C; vertically integrated cloud ice+snow+graupel',original_units='Published flash-origin score per grid cell per 5 min',
             SCM_adaptation='1D prescribed large-scale w, not original resolved-convective formulation. Factors interpolated separately to lowest descending isotherm; no rho in F1; F2 uses dry-air density. Score integrated over 900 s once.',
             integrated_predictor_units='published score s',coefficient_units='flash (published score s)-1',final_units='flashes per 15 min')])
    t2.to_csv(OUT/'tables/Table_02_lightning_parameterizations.csv',index=False)
    tablepaths=list((OUT/'tables').glob('Table_*.csv'))
    for path in tablepaths:
        table=pd.read_csv(path)
        latex_table(table,path.with_suffix('.tex'))
    write_json(OUT/'tables/metadata.json',{'Table_03':'Metrics of 40-member ensemble MEAN, evaluation only; deterministic control excluded; coefficient fitted on development ensemble mean.',
        'metric_units':{'rmse':'flashes/15min','mae':'flashes/15min','bias':'flashes/15min; prediction-observation','pearson':'dimensionless; undefined retained as NaN','spearman':'dimensionless; undefined retained as NaN',
                        'peak_time_error_minutes':'minutes; predicted minus observed; earliest peak on ties','accumulated_error':'flashes; predicted minus observed over identical 92 bins'},
        'CSV_precision':'full float64 serialization; LaTeX is formatted for display','source_case_sha256':sha256(CASE),
        'native_control_sha256':sha256(control_directory('TH_C0')/'output.nc'),'ensemble_manifest_sha256':sha256(ensemble/'manifest.json')})


def latex_table(table,path):
    """Self-contained table fragment; avoids optional pandas/Jinja rendering dependencies."""
    def escape(value):
        if pd.isna(value):return '--'
        if isinstance(value,(float,np.floating)):return f'{value:.6g}'
        mapping={'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_',
                 '{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}
        return ''.join(mapping.get(char,char) for char in str(value))
    equations=None
    if path.name.startswith('Table_01'):
        headers=['Parameter','Value','Units / convention'];columns=list(table.columns)
        fmt=r'p{.23\linewidth}p{.42\linewidth}p{.25\linewidth}'
    elif path.name.startswith('Table_02'):
        columns=['method','reference','equation','variables','original_units','SCM_adaptation','final_units']
        headers=['Method','Reference','Equation','Variables','Original units','SCM adaptation','Calibrated units']
        fmt=r'p{.07\linewidth}p{.10\linewidth}p{.25\linewidth}p{.11\linewidth}p{.10\linewidth}p{.20\linewidth}p{.10\linewidth}'
        equations=[r'$\mathrm{CAPE}\,P$',r'$\int (\rho_d r_i)(\rho_d r_g)\,dz$',
            r'$\begin{gathered}F_3=0.95F_1+0.05F_2\\F_1=0.042(w^+ r_g^{\rm g/kg})_{258.15\,K}\\F_2=0.20\int\rho_d(r_i+r_s+r_g)\,dz\end{gathered}$']
    else:
        columns=['method','microphysics','coupling','coefficient','rmse','mae','bias','pearson','spearman','peak_time_error_minutes','accumulated_error']
        headers=['Method','MP','MYNN',r'$a$','RMSE','MAE','Bias',r'$r$',r'$r_s$',r'$\Delta t_{\rm peak}$ (min)',r'$\Delta N$']
        fmt='lllrrrrrrrr'
    header=' & '.join(headers)+r' \\'
    lines=[r'% Requires booktabs, longtable, array, amsmath.',r'\begingroup\small\setlength{\tabcolsep}{2pt}',
           r'\begin{longtable}{'+fmt+'}',r'\toprule',header,r'\midrule\endfirsthead',r'\toprule',header,r'\midrule\endhead']
    for i,(_,row) in enumerate(table.iterrows()):
        values=[equations[i] if equations is not None and col=='equation' else escape(row[col]) for col in columns]
        lines.append(' & '.join(values)+r' \\')
    lines.extend([r'\bottomrule',r'\end{longtable}',r'\endgroup'])
    if path.name.startswith('Table_03'):
        lines.append(r'% RMSE, MAE, bias: flashes per 15 min; Delta N: flashes; correlations undefined -> --. Evaluation only.')
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')


def supplementary(c):
    microphysics(c,output_dir=OUT/'supplementary',filename='S05_six_configuration_processes')
    # Same scenario-only tables and exact closing decomposition already verified;
    # regenerate vector graphics in a new destination without touching old exports.
    from scm_lightning.evaluation.updraft_sensitivity import make_figures
    root=ROOT/'results/argentina_20181213/updraft_sensitivity'
    manifest=json.loads((root/'manifest.json').read_text());assert manifest['run_count']==18
    decomposition=pd.read_csv(root/'tables/response_decomposition_15min.csv.gz')
    valid=decomposition[['X00','Xr0','X0r','Xrr']].notna().all(axis=1)
    np.testing.assert_allclose((decomposition.state_response+decomposition.velocity_on_control_state+decomposition.interaction)[valid],decomposition.total_response[valid],atol=1e-10,rtol=1e-10)
    make_figures(pd.read_csv(root/'tables/predictors_15min.csv.gz'),decomposition,OUT/'supplementary/vertical_sensitivity',load_config(ROOT/'configs/experiments/updraft_sensitivity.yaml'),dpi=300,vector_pdf=True)
    write_json(OUT/'supplementary/source_metadata.json',{'vertical_scenarios':'18 deterministic cases only; never initial-condition ensemble members',
        'vertical_manifest_sha256':sha256(root/'manifest.json'),'decomposition_closure_test':'passed',
        'processes':'Actual frozen-water species tendencies, not exclusively graupel; native rate convention per initial dry-air cell mass retained',
        'process_source_manifest_paths':[str(ROOT/'results/argentina_20181213/control_diagnostics'/(key+'_manifest.json')) for key in CONFIGURATIONS]})


def main():
    if load_config(CASE).get('role')=='evaluation_only':
        from scm_lightning.plotting.independent import main as independent_main
        return independent_main()
    ap=argparse.ArgumentParser();ap.add_argument('--skip-supplementary',action='store_true');ap.add_argument('--supplementary-only',action='store_true');args=ap.parse_args()
    c,stats,o=load_data()
    with plt.rc_context(STYLE):
        if not args.supplementary_only:
            environment(c);observations(c,o);predictions(c,stats,o);performance_figure();microphysical_response(c)
            for coupling in ['C0','CL','CF']:
                if coupling!='C0':predictions(c,stats,o,coupling=coupling)
                predictions(c,stats,o,coupling=coupling,full=True)
            predictions(c,stats,o,full=True,percentiles=True);performance_figure('B');table_exports(c)
        if not args.skip_supplementary:supplementary(c)
    write_manifest()
    print('Saved five principal figures, three tables and separate supplements (PDF vector + PNG 300 dpi).')


def write_manifest():
    write_json(OUT/'manifest.json',{'methods':METHODS,'primary_reference_coupling':'C0','perturbed_members':40,'control_in_statistics':False,
        'technical_method_ids':{'CAPExP':'CAPExP_SBCAPE15_native_snapshot_1D_adaptado',
            'IxG':'Romps2019_eq4_full_column_mass_concentrations',
            'McCaulF3':'McCaul_F3_LS_native_snapshot_1D_adaptado'},
        'limitations':'CAPExP uses MetPy liquid pseudoadiabatic, not original isentropic CAPE. McCaulF3 uses prescribed large-scale ERA5 w, not resolved convective plume w; published score is empirically calibrated to interval counts.',
        'references':{'Romps2019':'https://doi.org/10.1029/2019GL085748','McCaul2009':'https://doi.org/10.1175/2008WAF2222152.1'},
        'model_reruns':0,'new_ensembles':0,'original_NetCDF_modified':False,'original_graphics_overwritten':False,
        'git_commit':git('rev-parse','HEAD'),
        'source_code':[{'path':str(p.relative_to(ROOT)),'sha256':sha256(p)} for p in [
            ROOT/'python/scm_lightning/evaluation/article.py',ROOT/'python/scm_lightning/plotting/article.py']],
        'MC_spread':'initial-condition spread only; no method/physics/pulse pooling',
        'supplementary_state_response_scope':'The state component includes thermodynamics, hydrometeors, density, precipitation and isotherm-height changes; it is not a purely microphysical attribution.',
        'statistics_audit_sha256':sha256(OUT/'audit/current_statistics.json'),'legacy_audit_sha256':sha256(OUT/'audit/legacy_audit.json'),
        'files':[{'path':str(p.relative_to(ROOT)),'sha256':sha256(p),'bytes':p.stat().st_size} for folder in ['figures','tables','supplementary','audit','data'] for p in (OUT/folder).rglob('*') if p.is_file()]})


if __name__=='__main__':main()
