"""NetCDF, tables and figures for the frozen deterministic pulse experiment.

No ENTLN access, fitting, ranking or interpretation report is performed here.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from scm_lightning.config import load_config
from scm_lightning.ensemble.updraft_sensitivity import RUNS, OUT, CONFIG, FORCING
from scm_lightning.ensemble.experiments import CONFIGURATIONS
from scm_lightning.lightning.adaptations import (
    validate_units, interface_density, recover_velocity, crossing, at_height, lpi_1d, surface_cape,
)
from scm_lightning.io import numeric_csv
from scm_lightning.provenance import ROOT, sha256, write_json

UNITS={'IxG':'kg2 m-5','VII':'kg m-2','McCaul_F2':'published_grid_score',
       'cold_precip_ice_path_1D':'kg m-2','cold_ice_path_no_radar_1D':'kg m-2',
       'PR92_wmax_1D':'published_nominal_flash_min-1_score','PW10_1D':'kg m-2 s-1',
       'IFluxT_1D':'kg m-2 s-1','IFluxP440_1D':'kg m-2 s-1',
       'McCaul_F1_1D':'published_grid_score','McCaul_F3_1D':'published_grid_score',
       'graupel_mass_flux258_1D':'kg m-2 s-1','LPI_1D':'J kg-1','updraft_depth5_1D':'m',
       'CAPExP_SBCAPE15_1D':'W m-2'}


def load_run(regime, configuration):
    folder=RUNS/regime/configuration
    with xr.open_dataset(folder/'output.nc.vertical.nc',decode_times=False) as source:d=source.load()
    with xr.open_dataset(folder/'output.nc',decode_times=False) as source:n=source.load()
    validate_units(d)
    np.testing.assert_allclose(d.rhoz,interface_density(d.rho_MYNN.values,d.layer_thickness.values),rtol=2e-14,atol=0)
    np.testing.assert_array_equal(d.w_EDMF.values,recover_velocity(d.mf_w.values,d.rhoz.values,d.mf_a.values))
    np.testing.assert_array_equal(d.time.values,n.time.values[:-1])
    state={key:d[key].values for key in ['temperature','pressure','qv','qc','qr','qi','qs','qg']}
    state.update(z=d.lev.values,dz=d.layer_thickness.values,time=d.time.values,
                 rho_d=d.pressure.values/(287.*d.temperature.values*(1+d.qv.values/(287./461.6))),
                 P=np.maximum(0,np.diff(n.rainnc.values)/np.diff(n.time.values)))
    return folder,d,n,state


def velocity(diag, representation):
    if representation=='prescribed':
        return {'w':diag.w_effective.values,'z':diag.lev.values,'active':None}
    return {'w':np.nan_to_num(diag.w_EDMF.values,nan=0.),'z':diag.z_interface.values,'active':diag.mf_a.values>0}


def isotherm_velocity(state, speed, target):
    height=crossing(state['temperature'],state['z'],target)
    w=at_height(speed['w'],speed['z'],height)
    if speed['active'] is not None:
        k=np.clip(np.searchsorted(speed['z'],height),1,len(speed['z'])-1)
        rows=np.arange(len(height));active=speed['active'][rows,k-1]&speed['active'][rows,k]
        w=np.where(active,w,np.where(np.isfinite(height),0.,np.nan))
    return height,w


def raw_proxies(state, speed):
    """All evaluable mathematical column predictors; never a storm-volume substitute."""
    t,z,dz,rho=state['temperature'],state['z'],state['dz'],state['rho_d']
    qc,qr,qi,qs,qg=(state[q] for q in ['qc','qr','qi','qs','qg'])
    wmax=np.maximum(0,np.max(speed['w'],axis=-1))
    w=np.asarray([np.interp(z,speed['z'],row,left=0.,right=0.) for row in speed['w']])
    vii=np.sum(rho*(qi+qs+qg)*dz,axis=-1)
    raw={'IxG':np.sum(rho*rho*qi*qg*dz,axis=-1),'VII':vii,'McCaul_F2':.20*vii,
         'cold_precip_ice_path_1D':np.sum(np.where(t<268.15,rho*(qs+qg)*dz,0.),axis=-1),
         'cold_ice_path_no_radar_1D':np.sum(np.where(t<263.15,rho*(qi+qs+qg)*dz,0.),axis=-1),
         'PR92_wmax_1D':5.e-6*wmax**4.5,'PW10_1D':state['P']*(wmax>10.),
         'LPI_1D':lpi_1d(t,dz,qc,qr,qi,qs,qg,w),
         'updraft_depth5_1D':np.sum(((t<268.15)&(w>5.))*dz,axis=-1)}
    h,wi=isotherm_velocity(state,speed,260.)
    raw['IFluxT_1D']=at_height(rho*qi,z,h)*wi*(wi>1.)
    h,wi=isotherm_velocity(state,speed,258.15)
    raw['McCaul_F1_1D']=.042*at_height(1000*qg,z,h)*np.maximum(0,wi)
    raw['McCaul_F3_1D']=.95*raw['McCaul_F1_1D']+.05*raw['McCaul_F2']
    raw['graupel_mass_flux258_1D']=at_height(rho*qg,z,h)*np.maximum(0,wi)
    height=np.array([np.interp(44000.,p[::-1],z[::-1],left=np.nan,right=np.nan) for p in state['pressure']])
    wi=at_height(w,z,height)
    raw['IFluxP440_1D']=at_height(rho*qi,z,height)*wi*(wi>1.)
    return raw


def bin_integrals(raw):
    result={}
    for key,values in raw.items():
        blocks=np.asarray(values).reshape(192,30)
        result[key]=np.where(np.isfinite(blocks).all(axis=1),np.sum(blocks*30.,axis=1),np.nan)
    return result


def cape_product(folder, state, source_sha):
    cache=folder/'cape_15min.csv';tag=folder/'cape_identity.json'
    if cache.exists() and tag.exists() and json.loads(tag.read_text())['source_sha256']==source_sha:
        return pd.read_csv(cache)
    rows=[]
    for i in range(0,5760,30):
        cape=surface_cape(state['pressure'][i],state['temperature'][i],state['qv'][i])
        rows.append({'time_s':float(state['time'][i]),'SBCAPE_J_kg':cape,
                     'precip_bin_kg_m2':float(np.sum(state['P'][i:i+30]*30)),
                     'CAPExP_integral':cape*float(np.sum(state['P'][i:i+30]*30))})
    result=pd.DataFrame(rows);result.to_csv(cache,index=False)
    write_json(tag,{'source_sha256':source_sha,'parcel':'lowest simulated cell; MetPy liquid pseudoadiabatic; 15min bin-start sample'})
    return result


def paths_and_budget(folder,d,n,state):
    rho=n.pressure.values/(287.*n.theta.values*(n.pressure.values/1.e5)**(2/7)*(1+n.qv.values/(287./461.6)))
    paths={'time_s':n.time.values}
    for q in ['qv','qc','qr','qi','qs','qg']:paths[q+'_kg_m2']=np.sum(rho*n[q].values*state['dz'],axis=-1)
    paths['rainnc_kg_m2']=n.rainnc.values
    paths['EDMF_sample_time_s']=np.r_[np.nan,d.time.values]
    paths['EDMF_active_depth_m']=np.r_[np.nan,np.sum((d.mf_a.values>0)*state['dz'][:-1],axis=1)]
    paths['EDMF_max_speed_m_s']=np.r_[np.nan,np.max(np.nan_to_num(d.w_EDMF.values,nan=0.),axis=1)]
    paths['EDMF_max_area']=np.r_[np.nan,np.max(d.mf_a.values,axis=1)]
    for name in ['w_LS','w_perturbation','w_effective']:
        paths[name+'_max_m_s']=np.r_[np.nan,np.max(d[name].values,axis=1)]
    frame=pd.DataFrame(paths)
    stage=numeric_csv(folder/'output.nc.stage_water.csv')
    assert len(stage)==5760
    assert np.isfinite(stage.values).all()
    assert stage.accounting_residual_ref.abs().max()<2.e-11
    assert stage.adv_reconstruction_error_ref.abs().max()<2.e-11
    np.testing.assert_allclose(stage.precip_kg_m2,np.diff(n.rainnc),rtol=1.e-8,atol=1.e-10)
    # Unbalanced prescribed pulse implies an additional mass-divergence demand.
    # Diagnose this demand, do not invent a compensating horizontal tendency.
    mass_flux=state['rho_d']*d.w_perturbation.values
    required_hdiv=-np.gradient(mass_flux,state['z'],axis=1)/state['rho_d']
    stage['pulse_continuity_demand_max_abs_s_1']=np.max(abs(required_hdiv),axis=1)
    wtotal=np.sum(frame[[q+'_kg_m2' for q in ['qv','qc','qr','qi','qs','qg']]],axis=1)
    stage['instant_density_inventory_change_kg_m2']=np.diff(wtotal)
    # Label the diagnostic density/scheme discrepancy separately, not a closure proof.
    stage['density_weighting_difference_kg_m2']=np.diff(wtotal)-stage.water_change_ref.values
    return frame,stage,required_hdiv


def main():
    cfg=load_config(CONFIG);checks=json.loads((OUT/'run_checks.json').read_text())
    assert checks['run_count']==18 and all(r['status']=='passed_basic_numerical_checks' for r in checks['runs'])
    tables=OUT/'tables';figures=OUT/'figures';netcdf=OUT/'netcdf'
    for folder in [tables,figures,netcdf]:folder.mkdir(parents=True,exist_ok=True)
    all_bins=[];decomposition=[];summaries=[];availability=[];reach_tables=[]
    for configuration in CONFIGURATIONS:
        basefolder,bd,bn,bs=load_run('ERA5',configuration)
        basecape=cape_product(basefolder,bs,sha256(basefolder/'output.nc.vertical.nc'))
        basebins={rep:bin_integrals(raw_proxies(bs,velocity(bd,rep))) for rep in ['prescribed','EDMF']}
        for rep in basebins:basebins[rep]['CAPExP_SBCAPE15_1D']=basecape.CAPExP_integral.values
        for regime in cfg['regimes']:
            folder,d,n,state=load_run(regime,configuration)
            np.testing.assert_array_equal(d.w_LS.values,bd.w_LS.values)
            np.testing.assert_array_equal(n.pressure.values,bn.pressure.values)
            before=n.time.values<=cfg['pulse']['start_s']
            for name in ['theta','qv','qc','qr','qi','qs','qg','rainnc']:
                np.testing.assert_array_equal(n[name].values[before],bn[name].values[before])
            cape=cape_product(folder,state,sha256(folder/'output.nc.vertical.nc'))
            paths,stage,continuity=paths_and_budget(folder,d,n,state)
            paths['configuration']=configuration;paths['regime']=regime
            paths.to_csv(tables/(regime+'_'+configuration+'_paths.csv.gz'),index=False)
            stage.to_csv(tables/(regime+'_'+configuration+'_water_budget.csv.gz'),index=False)
            ts=pd.Timestamp('2018-12-13T00:00:00Z')+pd.to_timedelta(np.arange(192)*900,unit='s')
            for rep in ['prescribed','EDMF']:
                current=bin_integrals(raw_proxies(state,velocity(d,rep)))
                fixed=bin_integrals(raw_proxies(state,velocity(bd,rep)))
                velocity_only=bin_integrals(raw_proxies(bs,velocity(d,rep)))
                current['CAPExP_SBCAPE15_1D']=cape.CAPExP_integral.values
                fixed['CAPExP_SBCAPE15_1D']=cape.CAPExP_integral.values
                velocity_only['CAPExP_SBCAPE15_1D']=basecape.CAPExP_integral.values
                for method,xrr in current.items():
                    x00=basebins[rep][method];xr0=fixed[method];x0r=velocity_only[method]
                    frame=pd.DataFrame({'time':ts,'configuration':configuration,'regime':regime,'velocity':rep,
                                        'method':method,'predictor_integral':xrr,
                                        'status':np.where(np.isfinite(xrr),'available','unavailable'),
                                        'raw_units':UNITS[method],'time_support':'native 30s step-begin integral / 15min'})
                    all_bins.append(frame)
                    diff=pd.DataFrame({'time':ts,'configuration':configuration,'regime':regime,'velocity':rep,'method':method,
                        'X00':x00,'Xr0':xr0,'X0r':x0r,'Xrr':xrr,
                        'state_response':xr0-x00,'velocity_on_control_state':x0r-x00,
                        'interaction':xrr-xr0-x0r+x00,'velocity_on_perturbed_state':xrr-xr0,'total_response':xrr-x00})
                    valid=diff[['X00','Xr0','X0r','Xrr']].notna().all(axis=1)
                    np.testing.assert_allclose((diff.state_response+diff.velocity_on_control_state+diff.interaction)[valid],diff.total_response[valid],rtol=1e-10,atol=1e-10)
                    decomposition.append(diff)
                    availability.append({'regime':regime,'configuration':configuration,'velocity':rep,'method':method,
                        'valid_bins':int(np.isfinite(xrr).sum()),'zero_bins':int((xrr==0).sum()),'total_bins':192,
                        'scope':'SCM exploratory' if method not in ['IxG','VII'] else 'physical column integral'})
            reach={'time_s':state['time'],'regime':regime,'configuration':configuration}
            for target,label in [(258.15,'minus15C'),(260.,'260K')]:
                h,w=isotherm_velocity(state,velocity(d,'EDMF'),target)
                reach[label+'_height_m_AGL']=h;reach[label+'_EDMF_speed_m_s']=w
                reach[label+'_EDMF_ascent_present']=np.isfinite(w)&(w>0)
            reach_tables.append(pd.DataFrame(reach))
            # Derived NetCDF contains the separated velocities and needed state/budget context.
            derived=d.copy()
            derived['pulse_required_horizontal_divergence']=(('time','lev'),continuity,{'units':'s-1','description':'-d(rho_d*w_perturbation)/dz/rho_d, dry-air continuity demand; diagnosed, NOT applied'})
            derived['rho_d']=(('time','lev'),state['rho_d'],{'units':'kg m-3','description':'Dry-air density for mass concentrations and dry-air continuity demand'})
            derived['precipitation_rate']=(('time',),state['P'],{'units':'kg m-2 s-1','cell_methods':'time: mean','description':'Exact rainnc increment over the following 30-second step; no snownc addition'})
            derived['precipitation_amount']=(('time',),state['P']*30.,{'units':'kg m-2','cell_methods':'time: sum'})
            derived['rainnc_step_start']=(('time',),n.rainnc.values[:-1],{'units':'kg m-2'})
            derived['time_bounds']=(('time','bounds'),np.column_stack([d.time.values,d.time.values+30.]),
                                    {'units':'seconds since 2018-12-13 00:00:00 +00:00','calendar':'gregorian'})
            derived.attrs.update(regime=regime,configuration=configuration,origin_utc='2018-12-13T00:00:00Z',
                forcing_sha256=cfg['forcing_sha256'],design_sha256=sha256(CONFIG),
                velocity_contract='w_effective=w_LS+w_perturbation; w_EDMF remains separate; MYNN receives w2=0',
                hydrometeor_contract='dry-air grid-column mixing ratios; EDMF products assume plume representativeness, not species covariance')
            derived.time.attrs.update(units='seconds since 2018-12-13 00:00:00 +00:00',calendar='gregorian',bounds='time_bounds')
            derived.to_netcdf(netcdf/(regime+'_'+configuration+'.nc'))
            summaries.append({'configuration':configuration,'regime':regime,'precipitation_kg_m2':float(n.rainnc.values[-1]),
                'pre_pulse_state_bit_identical':True,'prescribed_pressure_bit_identical':True,
                'peak_graupel_path_kg_m2':float(paths.qg_kg_m2.max()),'peak_cloud_ice_path_kg_m2':float(paths.qi_kg_m2.max()),
                'peak_frozen_path_kg_m2':float((paths.qi_kg_m2+paths.qs_kg_m2+paths.qg_kg_m2).max()),
                'EDMF_active_samples':int((d.mf_a.values>0).any(axis=1).sum()),
                'EDMF_peak_speed_m_s':float(paths.EDMF_max_speed_m_s.max()),
                'horizontal_water_ref_kg_m2':float(stage.horizontal_ref.sum()),
                'LS_vertical_water_ref_kg_m2':float(stage.LS_vertical_ref.sum()),
                'pulse_vertical_water_ref_kg_m2':float(stage.pulse_vertical_ref.sum()),
                'MYNN_water_stage_ref_kg_m2':float(stage.MYNN_stage_ref.sum()),
                'advection_clipping_water_ref_kg_m2':float(stage.adv_clip_ref.sum()),
                'MP_plus_precip_reference_residual_kg_m2':float(stage.MP_plus_precip_ref.sum()),
                'accounting_max_residual_ref_kg_m2':float(stage.accounting_residual_ref.abs().max()),
                'pulse_continuity_demand_max_abs_s_1':float(np.max(abs(continuity)))})
            print('Postprocessed',regime,configuration,flush=True)
    predictions=pd.concat(all_bins,ignore_index=True);differences=pd.concat(decomposition,ignore_index=True)
    predictions.to_csv(tables/'predictors_15min.csv.gz',index=False)
    differences.to_csv(tables/'response_decomposition_15min.csv.gz',index=False)
    pd.DataFrame(availability).to_csv(tables/'method_availability.csv',index=False)
    pd.DataFrame(summaries).to_csv(tables/'run_summary.csv',index=False)
    pd.concat(reach_tables).to_csv(tables/'isotherm_reach.csv.gz',index=False)
    unavailable=[{'method':name,'status':'unavailable_original','reason':reason} for name,reason in {
        'Barthe_PM':'Storm-volume precipitating-ice mass not determined by a column path',
        'Petersen_IWP18':'Validated reflectivity >18 dBZ mask unavailable; separate no-radar column integral exported',
        'Deierling_flux_product':'Terminal fall speeds and horizontal divergence of nonprecipitating ice unavailable',
        'updraft_volume':'Resolved 3D storm volume unavailable; separate column depth exported',
        'PR92_echo_top':'20 dBZ echo top not diagnosed',
        'original_CRM_velocity_schemes':'Original resolved-convective support not supplied by LS or area-mean EDMF',
        'uedmf_graupel_flux':'type=0 / explicit plume MP=0; module zeros are not available species fluxes'}.items()]
    pd.DataFrame(unavailable).to_csv(tables/'original_method_requirements.csv',index=False)
    make_figures(predictions,differences,figures,cfg)
    write_json(OUT/'manifest.json',{'case':'argentina_20181213','regimes':cfg['regimes'],'run_count':18,
        'design_sha256':sha256(CONFIG),'run_checks_sha256':sha256(OUT/'run_checks.json'),'ENTLN_accessed':False,'calibration_performed':False,
        'ensemble_expanded':False,'integrator_changed':False,'physics_changed':False,
        'transport_scope':'Unmodified principal integrator applies -w*d/dz to theta and vapor only. No resolved hydrometeor advection was added; electrical transport products are exploratory diagnostics.',
        'budget':'Exact fixed-initial-dry-mass stage accounting; MP-plus-precipitation residual and variable-density discrepancy retained, not certified full atmospheric closure.',
        'forcing_consistency':'Identical full forcing SHA for all runs; pressure/horizontal/surface forcing unchanged. Added w is a deliberately unbalanced advective sensitivity; required compensating divergence is diagnosed but not applied.',
        'decomposition':'Four corners X00, Xr0, X0r, Xrr: state + velocity_on_control + interaction = total. Velocity_on_perturbed is Xrr-Xr0. EDMF changes are diagnosed responses, not direct forcing activation.',
        'CAPE':'MetPy liquid pseudoadiabatic parcel from lowest simulated cell at bin start; 15min CAPE*P adaptation',
        'units':'Raw proxy units integrated over seconds. Published coefficients are numeric scores, not absolute SCM flash rates.',
        'source_files':[{'path':str(p.relative_to(ROOT)),'sha256':sha256(p),'bytes':p.stat().st_size} for p in list(tables.iterdir())+list(figures.iterdir())+list(netcdf.iterdir())]})


def make_figures(predictors,decomposition,out,cfg,dpi=150,vector_pdf=False):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def emit(fig,path):
        if vector_pdf:
            from matplotlib.ticker import FixedLocator,FuncFormatter
            for ax in fig.axes:
                if ax.get_xlim()[1]>24 and ax.get_xlim()[0]<1:
                    ax.set_xlim(0,48)
                    ax.xaxis.set_major_locator(FixedLocator([0,12,24,36,48]))
                    ax.xaxis.set_major_formatter(FuncFormatter(lambda hour,pos:(pd.Timestamp('2018-12-13T00:00Z')+pd.Timedelta(hours=hour)).strftime('%d %b\n%H:%M')))
                    if ax.get_xlabel():ax.set_xlabel('UTC')
            fig.subplots_adjust(bottom=max(.075,fig.subplotpars.bottom))
            if not any(ax.get_xlabel()=='UTC' for ax in fig.axes):fig.supxlabel('UTC',y=.005,fontsize=9)
        fig.savefig(path,dpi=dpi)
        if vector_pdf:fig.savefig(path.with_suffix('.pdf'))
    colors={'ERA5':'black','weak':'#167aaf','strong':'#ce4936'}
    basis=pd.read_csv(OUT/'design_basis.csv');hours=(pd.to_datetime(basis.time)-pd.Timestamp('2018-12-13')).dt.total_seconds()/3600
    from scm_lightning.ensemble.updraft_sensitivity import pulse
    z=np.linspace(0,20000,401);time=np.arange(0,172801,30.);p=cfg['pulse']
    fig,axes=plt.subplots(1,3,figsize=(15,5))
    axes[0].plot(hours,basis.ascent_score_m_s,color='black',label='ERA5 ascent score')
    for regime,amplitude in cfg['regimes'].items():
        value=pulse(z,time,p['start_s'],p['duration_s'],p['center_m_AGL'],p['sigma_m'],amplitude)
        axes[0].plot(time/3600,value.max(axis=1),color=colors[regime],label=regime+' increment')
        axes[1].plot(amplitude*np.exp(-.5*((z-p['center_m_AGL'])/p['sigma_m'])**2),z/1000,color=colors[regime])
    axes[0].set(xlabel='Hours since case start',ylabel='m s-1');axes[0].legend(fontsize=7)
    axes[1].set(xlabel='Peak imposed increment (m s-1)',ylabel='Height AGL (km)')
    axes[2].plot(hours,basis.z_minus15C_m_AGL/1000,label='ERA5 -15 C')
    axes[2].plot(hours,basis.z_minus5C_m_AGL/1000,label='ERA5 -5 C')
    axes[2].plot(hours,basis.z_minus20C_m_AGL/1000,label='ERA5 -20 C')
    axes[2].axhline(p['center_m_AGL']/1000,color='black',ls='--',label='Pulse center')
    axes[2].set(xlabel='Hours since case start',ylabel='Height AGL (km)');axes[2].legend(fontsize=7)
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();emit(fig,out/'00_frozen_meteorological_design.png');plt.close(fig)
    fig,axes=plt.subplots(6,3,figsize=(14,16),sharex=True)
    for i,c in enumerate(CONFIGURATIONS):
        for regime in cfg['regimes']:
            f=pd.read_csv(OUT/'tables'/(regime+'_'+c+'_paths.csv.gz'))
            hours=f.time_s/3600
            axes[i,0].plot(hours,f.qg_kg_m2,label=regime,color=colors[regime])
            axes[i,1].plot(hours,f.qi_kg_m2+f.qs_kg_m2,color=colors[regime])
            axes[i,2].plot(hours,f.rainnc_kg_m2,color=colors[regime])
        axes[i,0].set_ylabel(c+'\nkg m-2')
        for ax in axes[i]:ax.grid(alpha=.2);ax.axvspan(cfg['pulse']['start_s']/3600,(cfg['pulse']['start_s']+cfg['pulse']['duration_s'])/3600,color='.8',alpha=.2)
    for ax,title in zip(axes[0],['Graupel path','Cloud ice + snow path','Accumulated precipitation']):ax.set_title(title)
    for ax in axes[-1]:ax.set_xlabel('Hours since 13 Dec 2018 00 UTC')
    axes[0,0].legend(fontsize=8);fig.tight_layout();emit(fig,out/'01_hydrometeor_response.png');plt.close(fig)
    for rep in ['prescribed','EDMF']:
        methods=['IxG','CAPExP_SBCAPE15_1D','PR92_wmax_1D','IFluxT_1D','McCaul_F3_1D','LPI_1D']
        fig,axes=plt.subplots(6,6,figsize=(19,15),sharex=True)
        for i,c in enumerate(CONFIGURATIONS):
            for j,method in enumerate(methods):
                for regime in cfg['regimes']:
                    f=predictors[(predictors.configuration==c)&(predictors.regime==regime)&(predictors.velocity==rep)&(predictors.method==method)]
                    axes[i,j].plot(np.arange(192)/4,f.predictor_integral,color=colors[regime],label=regime,linewidth=.8)
                axes[i,j].grid(alpha=.2);axes[i,j].ticklabel_format(axis='y',style='sci',scilimits=(0,0));axes[i,j].tick_params(labelsize=6)
                if i==0:axes[i,j].set_title(method,fontsize=8)
                if j==0:axes[i,j].set_ylabel(c,fontsize=8)
        axes[0,0].legend(fontsize=6);fig.suptitle(rep+' velocity; raw proxy-unit seconds / 15 min; zero retained')
        fig.tight_layout(rect=[0,0,1,.97]);emit(fig,out/('02_predictors_'+rep+'.png'));plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for ax,c in zip(axes.flat,CONFIGURATIONS):
        for field,style in [('state_response','-'),('velocity_on_control_state','--'),('interaction',':')]:
            f=decomposition[(decomposition.configuration==c)&(decomposition.regime=='strong')&(decomposition.velocity=='prescribed')&(decomposition.method=='McCaul_F3_1D')]
            ax.plot(np.arange(192)/4,f[field],style,label=field)
        ax.set_title(c);ax.grid(alpha=.2);ax.set_xlabel('Hours since case start');ax.set_ylabel('F3 score s / 15 min');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    axes[0,0].legend(fontsize=7);fig.suptitle('Four-corner diagnostic decomposition: McCaul F3 score (strong pulse)')
    fig.tight_layout(rect=[0,0,1,.95]);emit(fig,out/'03_direct_velocity_and_state_components.png');plt.close(fig)
    fig,axes=plt.subplots(6,3,figsize=(14,16),sharex=True)
    for i,c in enumerate(CONFIGURATIONS):
        for regime in cfg['regimes']:
            f=pd.read_csv(OUT/'tables'/(regime+'_'+c+'_paths.csv.gz'));time=f.EDMF_sample_time_s/3600
            axes[i,0].plot(time,f.w_effective_max_m_s,color=colors[regime],label=regime)
            axes[i,1].plot(time,f.EDMF_max_speed_m_s,color=colors[regime])
            axes[i,2].plot(time,f.EDMF_active_depth_m/1000,color=colors[regime])
        axes[i,0].set_ylabel(c+'\nm s-1')
        for ax in axes[i]:ax.grid(alpha=.2)
    for ax,title in zip(axes[0],['Prescribed maximum w_LS + increment','Diagnosed EDMF maximum (m/s)','Active EDMF depth (km)']):ax.set_title(title,fontsize=9)
    axes[0,0].legend(fontsize=7)
    for ax in axes[-1]:ax.set_xlabel('Hours since case start')
    fig.tight_layout();emit(fig,out/'04_prescribed_velocity_and_EDMF.png');plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for ax,c in zip(axes.flat,CONFIGURATIONS):
        for regime in cfg['regimes']:
            f=pd.read_csv(OUT/'tables'/(regime+'_'+c+'_water_budget.csv.gz'))
            ax.plot(f.time_s/3600,f.pulse_vertical_ref.cumsum(),color=colors[regime],label=regime+' pulse transport')
            ax.plot(f.time_s/3600,f.adv_clip_ref.cumsum(),color=colors[regime],ls=':',label=regime+' clipping')
        ax.set(title=c,xlabel='Hours since case start',ylabel='kg m-2 / initial dry-mass weights');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=6);fig.tight_layout();emit(fig,out/'05_water_stage_accounting.png');plt.close(fig)


def plot_main():
    """Regenerate figures from saved tables without model runs or recalibration."""
    root=OUT/'tables';figures=OUT/'figures'
    make_figures(pd.read_csv(root/'predictors_15min.csv.gz'),pd.read_csv(root/'response_decomposition_15min.csv.gz'),figures,load_config(CONFIG))
    path=OUT/'manifest.json'
    manifest=json.loads(path.read_text(encoding='utf-8'))
    manifest['transport_scope']='Unmodified main integrator applies vertical transport to theta and vapor only; no resolved hydrometeor advection added.'
    manifest['pulse_time_definition']='Duration is full cosine support; temporal FWHM is half of this duration.'
    prefix=figures.relative_to(ROOT).as_posix()+'/'
    records=[record for record in manifest['source_files'] if not record['path'].replace('\\','/').startswith(prefix)]
    records.extend({'path':str(p.relative_to(ROOT)),'sha256':sha256(p),'bytes':p.stat().st_size} for p in figures.iterdir() if p.is_file())
    manifest['source_files']=records
    write_json(path,manifest)
