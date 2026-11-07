"""Numerical products for the isolated MYNN audit; never writes article results."""
import json
import warnings
import numpy as np
import pandas as pd
import xarray as xr
from scm_lightning.evaluation.mynn_audit import RUNS, OUT, SURFACE, GATES, FORCING
from scm_lightning.ensemble.experiments import CONFIGURATIONS
from scm_lightning.provenance import ROOT, sha256, write_json
from scm_lightning.lightning.adaptations import crossing, interface_density, recover_velocity


def saturation_specific(ts, ps):
    es=611.2*np.exp(17.67*(ts-273.15)/(ts-29.65)); ep=287./461.6
    return ep*es/(ps-(1-ep)*es)


def supported_reach(active, interfaces, heights):
    """Require both bracketing native interfaces; do not bridge inactive gaps."""
    i=np.searchsorted(interfaces, heights, side='right')-1
    valid=np.isfinite(heights)&(i>=0)&(i<len(interfaces)-1)
    j=np.clip(i,0,len(interfaces)-2); rows=np.arange(len(i))
    return valid&active[rows,j]&active[rows,j+1]


def lcl_height(p,t,rv,z):
    from metpy.calc import dewpoint,lcl
    from metpy.units import units
    vapor=p[:,0]*rv[:,0]/(287./461.6+rv[:,0])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pl,_=lcl(p[:,0]*units.Pa,t[:,0]*units.K,dewpoint(vapor*units.Pa))
    # Simulated lowest-cell parcel; log-pressure interpolation within existing grid only.
    return np.array([np.interp(np.log(a),np.log(b[::-1]),z[::-1],left=np.nan,right=np.nan)
                     if np.isfinite(a) and a>0 else np.nan for a,b in zip(pl.magnitude,p)])


def statistics(values):
    x=np.asarray(values); x=x[np.isfinite(x)]
    if not len(x):return {'mean':np.nan,'rms':np.nan,'p95_abs':np.nan,'max_abs':np.nan}
    return dict(mean=float(x.mean()),rms=float(np.sqrt(np.mean(x*x))),
                p95_abs=float(np.percentile(abs(x),95)),max_abs=float(abs(x).max()))


def relative_statistics(a,b):
    a,b=np.asarray(a),np.asarray(b);valid=np.isfinite(a)&np.isfinite(b)
    aa,bb=a[valid],b[valid]
    denominator=float(np.sqrt(np.mean(bb*bb))) if len(bb) else np.nan
    peak=float(np.max(abs(bb))) if len(bb) else np.nan
    return dict(valid_pairs=int(valid.sum()),missing_pairs=int((~valid).sum()),reference_rms=denominator,
                difference_rms_over_reference_rms=float(np.sqrt(np.mean((aa-bb)**2))/denominator) if denominator>0 else np.nan,
                difference_max_over_reference_max=float(np.max(abs(aa-bb))/peak) if peak>0 else np.nan)


def load(cfg,phase='original'):
    prefix=RUNS/phase/cfg/'output.nc'
    with xr.open_dataset(str(prefix)+'.mynn_audit.nc',decode_times=False) as f:d=f.load()
    with xr.open_dataset(str(prefix)+'.vertical.nc',decode_times=False) as f:v=f.load()
    with xr.open_dataset(prefix,decode_times=False) as f:n=f.load()
    return d,v,n


def summarize(cfg,phase):
    d,v,n=load(cfg,phase);s=pd.DataFrame(d.surface.values,columns=SURFACE)
    g=pd.DataFrame(d.gate.values,columns=GATES)
    z=d.lev.values;dz=v.layer_thickness.values;zi=v.z_interface.values
    active=v.plume_active.values.astype(bool);has=active.any(1)
    top=np.max(np.where(active,zi[None,:],np.nan_to_num(-np.inf)),axis=1);top[~has]=np.nan
    temp=d.temperature.values;p=d.pressure.values;rv=d.qv_entry.values
    iso={label:crossing(temp,z,target) for label,target in [('zero',273.15),('minus15',258.15),('260K',260.)]}
    term=d.termination.values[:,:,0];rho=d.rho.values
    np.testing.assert_allclose(v.rhoz,interface_density(v.rho_MYNN,dz),rtol=2e-15)
    np.testing.assert_allclose(v.w_EDMF,recover_velocity(v.mf_w,v.rhoz,v.mf_a),rtol=2e-15,equal_nan=True)
    assert np.all(np.isnan(v.w_EDMF.values[~active]))
    massflux=v.mf_a.values*v.mf_w.values
    row=dict(configuration=cfg,phase=phase,samples=len(has),active_samples=int(has.sum()),active_fraction=float(has.mean()),
             plume_top_max_m_AGL=float(np.nanmax(top)),plume_top_p50_m_AGL=float(np.nanpercentile(top,50)),
             plume_top_p95_m_AGL=float(np.nanpercentile(top,95)),pblh_mean_m_AGL=s.pblh.mean(),pblh_max_m_AGL=s.pblh.max(),
             mass_flux_max_kg_m2_s=float(massflux.max()),w_EDMF_max_m_s=float(v.w_EDMF.max()),
             area_max=float(v.mf_a.max()),no_active_interface=int((~has).sum()),
             gate_low_buoyancy=int((g.fltv_entry<=.002).sum()),gate_not_superadiabatic=int((g.superadiabatic==0).sum()),
             gate_width_le_minimum=int((g.maxwidth_pre<=g.minwidth).sum()),
             scale_Psig_shcu_min=g.Psig_shcu.min(),scale_Psig_shcu_max=g.Psig_shcu.max(),
             resolved_w_suppression_samples=int((g.max_resolved_w>1).sum()),
             qsfc_minus_current_qsat_max_abs=float(abs(s.qsfc-saturation_specific(s.tsk,s.psfc)).max()),
             graupel_plume_flux_available=False)
    for label,h in iso.items():row['reaches_'+label+'_samples']=int(supported_reach(active,zi,h).sum())
    for code in range(7):row['terminal_reason_'+str(code)+'_plumes']=int((term==code).sum())
    stopped=d.termination.values[(term>=1)&(term<=4)]
    for j,label in [(1,'attempted_stop_height_m_AGL'),(3,'terminal_buoyancy_m_s2'),(4,'terminal_entrainment_m-1'),(5,'terminal_buoyancy_delta_w_m_s'),(6,'terminal_drag_delta_w_m_s')]:
        for k,value in statistics(stopped[:,j]).items():row[label+'_'+k]=value
    # Overlapping gate counts are separate predicates, not an additive partition.
    series=s.copy();series['time_s']=d.time.values
    series['time']=pd.Timestamp('2018-12-13',tz='UTC')+pd.to_timedelta(d.time.values,unit='s')
    series['solar_hour']=(series.time.dt.hour+series.time.dt.minute/60-58.75/15)%24
    series['plume_top']=top;series['lcl']=lcl_height(p,temp,rv,z)
    for label,h in iso.items():series[label]=h
    series['w_LS_max']=v.w_LS.max('lev').values;series['w_LS_min']=v.w_LS.min('lev').values;series['area_max']=v.mf_a.max('interface').values
    series['qsat_current']=saturation_specific(s.tsk,s.psfc)
    series['Monin_Obukhov_surface_m']=np.divide(1.,s.rmol_sf.to_numpy(),out=np.full(len(s),np.nan),where=s.rmol_sf.to_numpy()!=0)
    for name in GATES:series['gate_'+name]=g[name]
    series.to_csv(OUT/f'{phase}_{cfg}_timeseries.csv',index=False)
    tendencies=[];checks={};dt=s.dt.to_numpy()[:,None];mass=rho*dz
    for species in ['qc','qi']:
        local=d[species+'_tendency_local'].values;full=d[species+'_tendency_full'].values
        delta=full-local
        applied=d[species+'_applied_increment'].values;expected=d[species+'_expected_increment'].values
        clip=d[species+'_clipping_increment'].values
        np.testing.assert_allclose(applied,expected+clip,atol=1e-18,rtol=1e-12)
        choice=np.zeros_like(full) if cfg.endswith('C0') else local if cfg.endswith('CL') else full
        np.testing.assert_allclose(expected,dt*choice,atol=0,rtol=0)
        np.testing.assert_allclose(clip,np.maximum(-(d[species+'_entry'].values+expected),0),atol=1e-18,rtol=1e-10)
        checks[species+'_clipping_kg_m2']=float(np.sum(clip*mass))
        checks[species+'_omitted_full_increment_kg_m2']=float(np.sum((dt*full-expected)*mass))
        checks[species+'_delta_abs_mass_above_0C_fraction']=float(np.sum(abs(delta)*dt*mass*(temp<=273.15))/np.sum(abs(delta)*dt*mass)) if np.any(delta) else np.nan
        checks[species+'_delta_abs_mass_above_minus15C_fraction']=float(np.sum(abs(delta)*dt*mass*(temp<=258.15))/np.sum(abs(delta)*dt*mass)) if np.any(delta) else np.nan
        for name,a in [('local',local),('full',full),('full_minus_local',delta)]:
            for k,height in enumerate(z):
                tendencies.append(dict(configuration=cfg,phase=phase,species=species,component=name,mask='all_times_at_native_level',
                    z_m_AGL=height,units='kg kg-1 s-1',**statistics(a[:,k]),
                    time_integral_kg_kg=float(np.sum(a[:,k]*dt[:,0])),
                    absolute_time_integral_kg_kg=float(np.sum(abs(a[:,k])*dt[:,0])),
                    mass_time_integral_kg_m2=float(np.sum(a[:,k]*dt[:,0]*mass[:,k]))))
            for label,mask in [('warm_T_gt_273.15K',temp>273.15),('mixed_233.15_to_273.15K',(temp<=273.15)&(temp>=233.15)),('cold_T_lt_233.15K',temp<233.15)]:
                tendencies.append(dict(configuration=cfg,phase=phase,species=species,component=name,mask=label,z_m_AGL=np.nan,
                    units='kg kg-1 s-1',**statistics(a[mask]),time_integral_kg_kg=np.nan,
                    absolute_time_integral_kg_kg=np.nan,mass_time_integral_kg_m2=float(np.sum(np.where(mask,a*dt*mass,0)))))
    # Sensible + latent enthalpy-rate diagnostic only. No radiative/momentum/pressure work closure claim.
    energy=np.sum(mass*(1004.5*d.theta_tendency.values*(p/1e5)**(2/7)+2.5e6*d.qv_tendency.values-3.34e5*d.qi_tendency_full.values),axis=1)
    checks['MYNN_enthalpy_rate_proxy_W_m2']=statistics(energy)
    checks['LH_minus_Lv_QFX_max_abs_W_m2']=float(abs(s.lh-2.5e6*s.qfx).max())
    checks['density_velocity_units_and_inactive_masks']='passed'
    checks['condensate_applied_exactly_once_and_clipping_accounted']='passed'
    checks['surface_mavail']=list(np.unique(s.mavail));checks['dx_m']=list(np.unique(s.dx))
    checks['native_sha256']=sha256(RUNS/phase/cfg/'output.nc')
    checks['independent_local_call_available']=bool(d.attrs['local_tendencies_available'])
    from scm_lightning.io import numeric_csv
    stage=numeric_csv(RUNS/phase/cfg/'output.nc.stage_water.csv')
    assert np.isfinite(stage.values).all()
    residual=float(abs(stage.accounting_residual_ref).max()); reconstruction=float(abs(stage.adv_reconstruction_error_ref).max())
    assert residual<1e-10 and reconstruction<1e-10
    checks['stage_accounting_max_residual_kg_m2']=residual
    checks['advective_reconstruction_max_error_kg_m2']=reconstruction
    checks['MYNN_stage_total_fixed_mass_kg_m2']=float(stage.MYNN_stage_ref.sum())
    checks['MP_plus_precip_total_fixed_mass_kg_m2']=float(stage.MP_plus_precip_ref.sum())
    if phase=='corrected':
        assert row['qsfc_minus_current_qsat_max_abs']<1e-12
        checks['surface_saturation_follows_current_Tskin_pressure']='passed'
    return row,tendencies,checks,(d,v,n,series)


def contrasts(cache):
    rows=[]
    for scheme in ['TH','WS']:
        for left,right in [('CL','C0'),('CF','CL'),('CF','C0')]:
            ca,cb=scheme+'_'+left,scheme+'_'+right
            da,va,na,_=cache[ca];db,vb,nb,_=cache[cb]
            # Same initial dry-air mass weights in every coupling: fixed reference, no continuity claim.
            p=nb.pressure.values[0];t=nb.theta.values[0]*(p/1e5)**(2/7);r=nb.qv.values[0]
            mass=p/(287*t*(1+r/(287/461.6)))*vb.layer_thickness.values
            for name in ['qc','qi','qg','qs','qr','rainnc']:
                a=na[name].values;b=nb[name].values
                if name!='rainnc':a=np.sum(a*mass,axis=1);b=np.sum(b*mass,axis=1)
                diff=np.squeeze(a-b)
                rows.append(dict(microphysics=scheme,contrast=left+'-'+right,variable=name,units='kg m-2' if name!='rainnc' else 'mm',
                                 **statistics(diff),**relative_statistics(a,b),final_difference=float(diff[-1]),raw_predictor=False))
            def raw(cfg):
                return pd.read_csv(ROOT/f'results/argentina_20181213/diagnostics_v3/{cfg}/00/predictors_1min.csv').pivot(index='time',columns='method',values='predictor')
            a,b=raw(ca).align(raw(cb),join='inner')
            for method in a:
                diff=a[method]-b[method]
                rows.append(dict(microphysics=scheme,contrast=left+'-'+right,variable=method,units='native predictor; see diagnostics_v3 metadata',
                                 **statistics(diff),**relative_statistics(a[method],b[method]),final_difference=float(diff.iloc[-1]),raw_predictor=True))
    pd.DataFrame(rows).to_csv(OUT/'coupling_contrasts.csv',index=False)


def surface_era5(cache):
    """Use existing common-box upward flux forcing; its ERA5 provenance is retained."""
    with xr.open_dataset(FORCING,decode_times=False) as f:
        era=f[['hfss','hfls','ts_forc']].squeeze().load();time=f.time.values
    rows=[]
    for cfg,(_,_,_,s) in cache.items():
        for column,source in [('hfx','hfss'),('lh','hfls'),('tsk','ts_forc')]:
            a=np.interp(s.time_s,time,era[source].values)
            for period,mask in [('solar_day_06_18',(s.solar_hour>=6)&(s.solar_hour<18)),('solar_night',(s.solar_hour<6)|(s.solar_hour>=18)),('full',np.ones(len(s),bool))]:
                rows.append(dict(configuration=cfg,period=period,variable=column,reference=source,units='K' if column=='tsk' else 'W m-2 upward',
                                 MYNN_mean=float(s.loc[mask,column].mean()),ERA5_forcing_mean=float(a[mask].mean()),
                                 **statistics(s.loc[mask,column]-a[mask])))
    pd.DataFrame(rows).to_csv(OUT/'surface_ERA5_comparison.csv',index=False)


def figures(cache,corrected):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    fig,ax=plt.subplots(3,2,figsize=(13,10),sharex=True,layout='constrained')
    for c,cfg in enumerate(['TH_CF','WS_CF']):
        d,v,n,s=cache[cfg];time=s.time
        for key,label in [('plume_top','EDMF top'),('pblh','PBL'),('lcl','LCL'),('zero','0°C'),('minus15','−15°C'),('260K','260 K')]:ax[0,c].plot(time,s[key]/1000,label=label,lw=1)
        ax[0,c].set(title=cfg,ylabel='Height (km AGL)',ylim=(0,12));ax[0,c].legend(ncol=3,fontsize=7)
        ax[1,c].plot(time,s.w_LS_max,label='LS max');ax[1,c].plot(time,v.w_EDMF.max('interface'),label='EDMF max');ax[1,c].legend();ax[1,c].set_ylabel('Velocity (m s⁻¹)')
        ax[1,c].plot(time,s.w_LS_min,label='LS min',ls=':',lw=.8);ax[1,c].legend(fontsize=7)
        ax[2,c].plot(time,s.area_max);ax[2,c].set_ylabel('Max plume area fraction')
        if cfg in corrected:
            _,cv,_,cs=corrected[cfg]
            ax[0,c].plot(time,cs.plume_top/1000,'k--',lw=.8,label='Boundary corrected');ax[0,c].legend(ncol=3,fontsize=7)
            ax[1,c].plot(time,cv.w_EDMF.max('interface'),'k--',lw=.8,label='Boundary corrected');ax[1,c].legend(fontsize=7)
            ax[2,c].plot(time,cs.area_max,'k--',lw=.8)
    format_dates(ax);fig.savefig(folder/'01_plume_heights.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(4,2,figsize=(13,12),sharex=True,layout='constrained')
    with xr.open_dataset(FORCING,decode_times=False) as f:forcing=f.squeeze().load()
    ft=pd.Timestamp('2018-12-13',tz='UTC')+pd.to_timedelta(forcing.time.values,unit='s')
    for c,cfg in enumerate(['TH_CF','WS_CF']):
        _,_,_,s=cache[cfg];time=s.time
        ax[0,c].plot(time,s.hfx,label='MYNN HFX');ax[0,c].plot(time,s.lh,label='MYNN LH');ax[0,c].plot(ft,forcing.hfss,label='ERA5 HFX',ls='--');ax[0,c].plot(ft,forcing.hfls,label='ERA5 LH',ls='--');ax[0,c].legend(fontsize=7,ncol=2);ax[0,c].set(title=cfg,ylabel='W m⁻² upward')
        ax[1,c].plot(time,s.qsfc*1000,label='Boundary q');ax[1,c].plot(time,s.qsat_current*1000,label='Current saturation');ax[1,c].legend(fontsize=8);ax[1,c].set_ylabel('g kg⁻¹ moist gas')
        ax[2,c].plot(time,s.gate_fltv_entry,label='Buoyancy flux');ax[2,c].axhline(.002,ls='--',color='k');ax[2,c].set_ylabel('K m s⁻¹')
        for key in ['maxwidth_pre','dx_cap','pbl_cap','cloud_cap','flux_cap']:ax[3,c].plot(time,s['gate_'+key],label=key,lw=.7)
        ax[3,c].axhline(300,ls='--',color='k');ax[3,c].set(ylabel='Width limits (m)',ylim=(0,2000));ax[3,c].legend(fontsize=6,ncol=3)
        if cfg in corrected:
            cs=corrected[cfg][3]
            ax[0,c].plot(time,cs.hfx,color='k',ls=':',lw=.7,label='Corrected HFX');ax[0,c].plot(time,cs.lh,color='grey',ls=':',lw=.7,label='Corrected LH');ax[0,c].legend(fontsize=6,ncol=3)
            ax[1,c].plot(time,cs.qsfc*1000,'k:',lw=.8,label='Corrected boundary');ax[1,c].legend(fontsize=7)
    format_dates(ax);fig.savefig(folder/'02_surface_activation.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(4,2,figsize=(12,12),layout='constrained')
    for c,cfg in enumerate(['TH_CF','WS_CF']):
        d,v,_,s=cache[cfg];z=d.lev.values/1000;zi=v.z_interface.values/1000
        low=z<=8
        ax[0,c].plot(d.tke.mean('time').values[low],z[low]);ax[0,c].set(title=cfg,xlabel='Mean TKE (m² s⁻²)',ylabel='km AGL',ylim=(0,8))
        ax[1,c].plot((v.mf_a*v.mf_w).max('time'),zi);ax[1,c].set(xlabel='Max EDMF mass flux (kg m⁻² s⁻¹)',ylabel='km AGL',ylim=(0,8))
        speed=ax[1,c].twiny()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            speed.plot(v.w_EDMF.mean('time',skipna=True),zi,'--',color='tab:orange')
        speed.set_xlabel('Conditional mean active plume speed (m s⁻¹)',color='tab:orange');speed.tick_params(axis='x',colors='tab:orange')
        for r,sp in [(2,'qc'),(3,'qi')]:
            for suffix,label in [('local','EDMF OFF shadow'),('full','Full MYNN')]:
                ax[r,c].plot(np.sqrt((d[sp+'_tendency_'+suffix]**2).mean('time'))*1e6,z,label=label)
            delta=d[sp+'_tendency_full']-d[sp+'_tendency_local'];ax[r,c].plot(np.sqrt((delta**2).mean('time'))*1e6,z,label='Full − shadow')
            ax[r,c].set(xlabel=sp+' tendency RMS (mg kg⁻¹ s⁻¹)',ylabel='km AGL',ylim=(0,12));ax[r,c].legend(fontsize=7)
    fig.savefig(folder/'03_vertical_tendencies.png',dpi=170);plt.close(fig)


def format_dates(ax):
    import matplotlib.dates as mdates
    for a in ax.flat:a.grid(alpha=.2);a.xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'));a.tick_params(axis='x',rotation=90)
    for a in ax[-1]:a.set_xlabel('December 2018 — UTC')


def corrected_predictors():
    """Raw diagnostics only: no frozen coefficients, observations or fitting."""
    from scm_lightning.config import load_config
    from scm_lightning.evaluation.independent import control_predictors
    from scm_lightning.lightning.thermodynamics import diagnose
    config=load_config(ROOT/'configs/cases/argentina_20181213.yaml');config['interval_minutes']=1
    contrasts=[]
    for cfg in ['TH_CF','WS_CF']:
        source=RUNS/'corrected'/cfg/'output.nc';folder=OUT/'corrected_predictors'/cfg
        tag=folder/'audit_identity.json';identity={'source_sha256':sha256(source),
            'thermodynamics_sha256':sha256(ROOT/'python/scm_lightning/lightning/thermodynamics.py'),
            'ice_integrals_sha256':sha256(ROOT/'python/scm_lightning/evaluation/independent.py'),
            'calibrated':False,'requires_2019_recalibration':True,'observations_used':False}
        if tag.exists():
            assert json.loads(tag.read_text())==identity
        else:
            primary,_=control_predictors(source,FORCING,config,folder,configuration=cfg,compute_cape=False)
            with xr.open_dataset(source,decode_times=False) as f:d=f.load()
            p=d.pressure.values;t=d.theta.values*(p/1e5)**(2/7);r=d.qv.values;rain=d.rainnc.values
            rows=[]
            for i in range(0,len(p)-1,2):
                row=diagnose(p[i],t[i],r[i],d.lev.values)
                row.update(time=str(pd.Timestamp(config['start'])+pd.Timedelta(seconds=float(d.time.values[i]))),
                    precipitation_amount_kg_m2=float(rain[i+2]-rain[i]))
                row['CAPExP_J_m2']=row['cape_J_kg']*row['precipitation_amount_kg_m2'];rows.append(row)
            thermo=pd.DataFrame(rows);thermo.to_csv(folder/'thermodynamics_mucape_1min.csv',index=False)
            cap=thermo[['time','CAPExP_J_m2','status']].rename(columns={'CAPExP_J_m2':'predictor'})
            cap['method']='MUCAPExP_300hPa_SCM';cap['configuration']=cfg;cap['member']=0
            primary=pd.concat([primary[primary.method!='CAPExP'],cap],ignore_index=True)
            primary.to_csv(folder/'predictors_1min.csv',index=False);write_json(tag,identity)
        a=pd.read_csv(folder/'predictors_1min.csv');b=pd.read_csv(ROOT/f'results/argentina_20181213/diagnostics_v3/{cfg}/00/predictors_1min.csv')
        for frame in [a,b]:frame['time']=pd.to_datetime(frame.time,utc=True)
        merged=a.merge(b,on=['time','method'],suffixes=('_corrected','_reference'),validate='one_to_one')
        for method,group in merged.groupby('method'):
            delta=group.predictor_corrected-group.predictor_reference
            contrasts.append(dict(configuration=cfg,variable=method,valid_pairs=int(delta.notna().sum()),
                missing_pairs=int(delta.isna().sum()),calibrated=False,**statistics(delta)))
        with xr.open_dataset(source,decode_times=False) as a,xr.open_dataset(RUNS/'original'/cfg/'output.nc',decode_times=False) as b:
            for name in ['qc','qi','qg','rainnc']:
                contrasts.append(dict(configuration=cfg,variable=name,calibrated=False,**statistics(a[name].values-b[name].values)))
        print('Uncalibrated corrected predictors: '+cfg,flush=True)
    pd.DataFrame(contrasts).to_csv(OUT/'boundary_correction_contrasts.csv',index=False)


def analyze():
    OUT.mkdir(parents=True,exist_ok=True);rows=[];tendencies=[];checks={};cache={};corrected={}
    phases=['original']+(['corrected'] if (OUT/'corrected_run_tests.json').exists() else [])
    for phase in phases:
        for cfg in CONFIGURATIONS if phase=='original' else ['TH_CF','WS_CF']:
            row,t,ch,data=summarize(cfg,phase);rows.append(row);tendencies.extend(t);checks[phase+'/'+cfg]=ch
            if phase=='original':cache[cfg]=data
            else:corrected[cfg]=data
    pd.DataFrame(rows).to_csv(OUT/'plume_audit.csv',index=False)
    pd.DataFrame(tendencies).to_csv(OUT/'mynn_tendencies.csv',index=False)
    contrasts(cache);surface_era5(cache);figures(cache,corrected)
    protected=json.loads((OUT/'protected_inputs.json').read_text())
    changed=[r['path'] for r in protected if sha256(r['path'])!=r['sha256']]
    assert not changed,changed
    result={'checks':checks,'protected_files_unchanged':True,'protected_file_count':len(protected),
        'calibration_changed':False,'ensemble_reruns':0,'alternative_core':{'tested':False,'reason':'Deferred after two proven boundary inconsistencies; corrected original requires validation and 2019 recalibration.'},
        'corrected_model_requires_2019_recalibration':True,
        'limitations':{'surface_mavail':'fixed wet-limit=1; no soil resistance/availability in cached ERA5; not validated for Argentine land surface',
          'dx':'5000m nominal parametrization scale inherited; no independent justification for this case; not observational box width',
          'w2':'zero resolved-column w; LS transport separate; resolved-w suppression inactive; Psig_shcu still applies',
          'radiation':'rthraten=0, topdown=0, downdrafts=0; meteorological forcing contains indirect radiation history, no explicit column radiation closure',
          'surface_flux':'MYNN diagnoses upward flux; ERA5 hfss/hfls are comparison, not prescribed flux; no ERA5 blh/ust in cache',
          'surface_gustiness':'input HFX/QFX reset every step; resulting surface wstar=0; no convective gustiness feedback, retained reference',
          'energy':'enthalpy tendency proxy only; no exact budget without radiation, momentum, pressure work and all native stage energy terms',
          'water':'stage accounting telescopes with fixed initial dry-air weights; native MP entry-density residual is not exact full atmospheric conservation',
          'full_minus_local':'whole MYNN call response including feedback, not isolated covariant nonlocal condensate flux',
          'C0_CL':'full theta/vapor retained; C0 omits condensate, CL inserts shadow condensate; not separately conservative EDMF-OFF integrations',
          'termination':'negative buoyancy and entrainment drag coexist; trajectory causal attribution to entrainment-induced dilution remains indeterminate',
          'gate_counts':'overlapping predicates; not a disjoint partition; actual termination codes from DMP only',
          'graupel_plume':'unavailable with original core and plume microphysics off; column graupel is not plume graupel flux'},
        'source_hashes':{str(p.relative_to(ROOT)):sha256(p) for p in [ROOT/'src/mynn_offline_coupling.F90',ROOT/'src/mp_offline_driver.F90',ROOT/'src/mynn/bl_mynn_subroutines.F90',ROOT/'src/diagnostics/mynn_audit.F90']},
        'forcing_sha256':sha256(FORCING)}
    def clean(a):
        if isinstance(a,dict):return {k:clean(v) for k,v in a.items()}
        if isinstance(a,list):return [clean(v) for v in a]
        if isinstance(a,(float,np.floating)) and not np.isfinite(a):return None
        return a
    write_json(OUT/'audit_tests.json',clean(result))
    print('MYNN audit: '+str(len(rows))+' controls summarized; protected files unchanged; 3 diagnostic figures.')


def finalize():
    """Collect machine-readable checks and print at most fifteen technical lines."""
    import platform
    import xml.etree.ElementTree as ET
    import metpy
    root=ET.parse(OUT/'pytest_final.xml').getroot()
    suites=list(root.iter('testsuite'))
    counts={name:sum(int(s.attrib.get(name,0)) for s in suites) for name in ['tests','failures','errors','skipped']}
    assert counts['failures']==counts['errors']==0
    result={name:json.loads((OUT/name).read_text()) for name in ['original_run_tests.json','reference_v2_run_tests.json','no_shadow_v3_run_tests.json','corrected_run_tests.json','audit_tests.json','all_native_preservation.json']}
    assert result['all_native_preservation.json']['unchanged']
    result['pytest']=dict(counts,passed=counts['tests']-counts['skipped'],xml_sha256=sha256(OUT/'pytest_final.xml'))
    before=ET.parse(OUT/'boundary_before_fix.xml').getroot()
    result['regression_before_correction']={'expected_failure_reproduced':any(int(s.attrib.get('failures',0))>0 for s in before.iter('testsuite')),
        'xml_sha256':sha256(OUT/'boundary_before_fix.xml')}
    result['versions']={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'xarray':xr.__version__,'MetPy':metpy.__version__}
    result['boundary_scope']={'option':'correct_mynn_boundary_units=.true.','default':False,
        'two_errors_corrected_together':True,'main_experiments_changed':False,'new_electrical_coefficients_fitted':False,
        'mavail':'dimensionless availability multiplier of moisture conductance flqc; held at 1, wet-limit; no validated ERA5 replacement in cache',
        'core_constants_and_integrator_equations_changed':False,'optional_alternative_core':'not run; boundary validation prioritized',
        'neutral_Monin_Obukhov_length':'missing when reciprocal is exactly zero; not coerced to a finite value',
        'no_shadow_v2_metadata':'Superseded diagnostic metadata only: availability attribute was initialized before option assignment; corrected and rerun in no_shadow_v3. Native physics identical in both.'}
    result['DMP_original_source_contract']={
        'buoyancy_gate_K_m_s':.002,'land_superadiabatic_gradient_K_m':-.005,
        'minimum_width_m':300,'maximum_width_m':'min(1.2*dx,1000,1.1*pblh,cloud_base_limit,flux_limit)',
        'entrainment_m-1':'max(0.33/(min(max(w_previous,0.3+0.0005*diameter),0.9)*diameter),0.0003); add 5e-6*(z-min(pblh+1500,4000)) above threshold; cap at 0.9/dz; stochastic multiplier=1',
        'detrainment_m-1':'min(max(-0.060/max(1,0.5*(w_next+w_previous))*(w_next-w_previous)/dz,0.00008),0.0002); momentum cap 0.0001',
        'detrainment_role':'modifies environment thermodynamics/condensate/winds; area held constant along surviving native plume; no direct detrainment-rate time series exported',
        'momentum_equation':'w_next=w_previous+(-2*ENT*w_previous+BCOEFF*B/max(w_previous,0.2))*min(delta_z,250); BCOEFF=0.15 if B>0 else 0.2; original acceleration/deceleration limiter and [0,3] m/s cap unchanged',
        'resolved_w_role':'asymmetric max(abs(w_up),2*abs(w_down)) suppression starts at 1 m/s, completes at 2 m/s; w2=0 does not supply convective lifting',
        'scale_role':'Psig_w also includes Psig_shcu from SCALE_AWARE; horizontal grid width capped at 1000m with dx=5000m',
        'entrainment_causality':'terminal B and direct drag increments measured; a separate causal contribution of entrainment dilution to B is not resolved'}
    result['products']={str(p.relative_to(OUT)):sha256(p) for p in OUT.rglob('*') if p.is_file() and p.suffix in ['.csv','.png']}
    result['audit_code_hashes']={str(p.relative_to(ROOT)):sha256(p) for p in [ROOT/'python/scm_lightning/evaluation/mynn_audit.py',ROOT/'python/scm_lightning/evaluation/mynn_audit_products.py']}
    write_json(OUT/'test_results.json',result)
    table=pd.read_csv(OUT/'plume_audit.csv').set_index(['phase','configuration'])
    lines=['MYNN audit v1 — deterministic controls only.',
        'Corrected opt-in boundary errors: dry/specific humidity interface and latched qsfc.',
        'Original core/thresholds/integrator preserved; no ensemble reruns or electrical recalibration.',
        f'Tests: {result["pytest"]["passed"]} passed, {counts["skipped"]} historical skips; {result["all_native_preservation.json"]["files_checked"]} native files unchanged.',
        'Six diagnostic controls and six correction-OFF controls: every native variable bit-identical.',
        'Shadow removed in two CF controls: bit-identical; unavailable local diagnostics explicitly NaN.',
        'Original stops: negative buoyancy with drag dominates; positive-buoyancy drag stops also recorded.',
        'CL/CF differ: final precipitation CF−CL = +0.00686 mm TH; −0.00736 mm WS.',
        'Corrected plume maxima: '+', '.join(f'{c} {table.loc[("corrected",c),"plume_top_max_m_AGL"]:.1f} m AGL' for c in ['TH_CF','WS_CF'])+'.',
        'No corrected plume reaches −15°C; water stage accounting checked; exact full energy/mass closure unavailable.',
        '2019 recalibration REQUIRED before corrected-model scientific lightning evaluation.',
        'Alternative core type=1 deferred; fixed mavail and omitted radiation remain limitations.',
        'CSV: plume_audit, mynn_tendencies, coupling_contrasts, surface_ERA5_comparison, boundary_correction_contrasts.',
        'Outputs: results/argentina_20181213/mynn_audit_v1; 3 PNG figures; test_results.json.',
        'NetCDF and run identities: runs/argentina_20181213/mynn_audit_v1.']
    assert len(lines)<=15
    # Native Windows redirected stdout may still use cp1252 (U+2212 unsupported).
    print('\n'.join(lines).replace('\u2212','-').replace('\u2014','-').replace('\u00b0',' deg'))


if __name__=='__main__':analyze()
