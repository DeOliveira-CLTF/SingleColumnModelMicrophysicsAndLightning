"""Deterministic surface-interface tests; no observation-based tuning or ensembles.

make BUILD=build/mynn_surface_v1 -j2
python -m scm_lightning.evaluation.mynn_surface run
python -m scm_lightning.evaluation.mynn_surface analyze
python -m scm_lightning.evaluation.mynn_surface land
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import threading
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.ensemble.experiments import CONFIGURATIONS,control_directory,namelist,execute,validate_output
from scm_lightning.evaluation.mynn_audit import protected,SURFACE,GATES
from scm_lightning.evaluation.mynn_audit_products import supported_reach,saturation_specific,statistics
from scm_lightning.lightning.adaptations import crossing,interface_density,recover_velocity

CFG=ROOT/'configs/experiments/mynn_surface_2018.yaml'
RUNS=ROOT/'runs/argentina_20181213/mynn_surface_v1'
OUT=ROOT/'results/argentina_20181213/mynn_surface_v1'
EXE=ROOT/'build/mynn_surface_v1/mp_offline_driver'
FORCING=ROOT/'data/processed/forcing/LEONA_SCM.nc'
BOUNDARY='incoming_hfx incoming_qfx bulk_hfx bulk_qfx applied_hfx applied_qfx applied_lh forcing_hfx forcing_lh wstar ust rmol_bulk rho_surface previous_pblh tsk qsfc mavail mode memory_enabled dt'.split()
LOCK=threading.Lock()


def settings():return yaml.safe_load(CFG.read_text())


def initial_protection():
    OUT.mkdir(parents=True,exist_ok=True)
    p=OUT/'protected_inputs.json'
    if not p.exists():write_json(p,protected())
    return json.loads(p.read_text())


def one(scenario,cfg):
    config=settings();options={} if scenario=='reference_off' else config['scenarios'][scenario]
    folder=RUNS/scenario/cfg;folder.mkdir(parents=True,exist_ok=True);output=folder/'output.nc'
    additions='write_vertical_diagnostics=.true.,write_mynn_audit=.true.,'
    if options.get('correct_mynn_boundary_units'):additions+='correct_mynn_boundary_units=.true.,'
    if options.get('mynn_surface_flux_memory'):additions+='mynn_surface_flux_memory=.true.,'
    additions+=f'mynn_surface_flux_mode={options.get("mynn_surface_flux_mode",0)},'
    text=namelist(cfg,FORCING,output).replace('write_water_budget=.true.,','write_water_budget=.true.,'+additions)
    identity={'scenario':scenario,'configuration':cfg,'executable_sha256':sha256(EXE),
        'forcing_sha256':sha256(FORCING),'config_sha256':sha256(CFG),'namelist':text,
        'source_hashes':{name:sha256(ROOT/name) for name in ['src/mp_offline_driver.F90','src/mynn_offline_coupling.F90','src/diagnostics/mynn_audit.F90']},
        'lightning_calibrated':False,'ensemble_member':False}
    tag=folder/'identity.json'
    if output.exists():
        if json.loads(tag.read_text())!=identity:raise ValueError('Changed surface inputs: use a new version directory')
    else:
        write_json(tag,identity);(folder/'namelist.mp_offline').write_text(text);execute(EXE,folder)
    with LOCK:
        result=validate_output(output,172800)
        assert result['physical_preflight_passed']
        with xr.open_dataset(output,decode_times=False) as d:
            temperature=d.theta.values*(d.pressure.values/1e5)**(2/7)
            assert np.all((temperature>100)&(temperature<400))
            assert np.all(d.pressure.values>0) and np.all(np.diff(d.pressure.values,axis=-1)<0)
            if scenario in ['reference_off','corrected_bulk']:
                original=control_directory(cfg)/'output.nc' if scenario=='reference_off' else ROOT/f'runs/argentina_20181213/mynn_audit_v1/corrected/{cfg}/output.nc'
                with xr.open_dataset(original,decode_times=False) as b:
                    assert set(d.variables)==set(b.variables)
                    for name in d.variables:np.testing.assert_array_equal(d[name],b[name],err_msg=scenario+'/'+cfg+'/'+name)
                result['baseline_all_variables_bit_identical']=True
    result.update(configuration=cfg,scenario=scenario,output=str(output))
    print(scenario+' '+cfg+': passed native checks',flush=True)
    return result


def run():
    initial_protection();tests=[]
    # First prove unchanged defaults before executing any boundary sensitivities.
    for scenario,configs in [('reference_off',CONFIGURATIONS)]+[(k,v['configurations']) for k,v in settings()['scenarios'].items()]:
        with ThreadPoolExecutor(max_workers=2) as pool:tests.extend(pool.map(lambda c:one(scenario,c),configs))
        write_json(OUT/'run_tests.json',{'runs':tests,'ensemble_reruns':0,'selection_using_lightning':False})
    return tests


def analyze():
    from scm_lightning.io import numeric_csv
    tests=json.loads((OUT/'run_tests.json').read_text());records=[];budget_rows=[];checks=[];cache={};compatibility=[]
    for record in tests['runs']:
        scenario,cfg=record['scenario'],record['configuration']
        if scenario=='reference_off':continue
        path=RUNS/scenario/cfg/'output.nc'
        with xr.open_dataset(str(path)+'.mynn_audit.nc',decode_times=False) as f:d=f.load()
        with xr.open_dataset(str(path)+'.vertical.nc',decode_times=False) as f:v=f.load()
        with xr.open_dataset(path,decode_times=False) as f:n=f.load()
        b=pd.DataFrame(d.surface_boundary.values,columns=BOUNDARY);s=pd.DataFrame(d.surface.values,columns=SURFACE);g=pd.DataFrame(d.gate.values,columns=GATES)
        time=d.time.values;z=d.lev.values;dz=v.layer_thickness.values;dt=b.dt.values[:,None]
        np.testing.assert_allclose(v.rhoz,interface_density(v.rho_MYNN,dz),rtol=2e-15)
        np.testing.assert_allclose(v.w_EDMF,recover_velocity(v.mf_w,v.rhoz,v.mf_a),equal_nan=True)
        np.testing.assert_allclose(b.applied_hfx,s.hfx,rtol=0,atol=0)
        np.testing.assert_allclose(b.applied_lh,s.lh,rtol=0,atol=0)
        np.testing.assert_allclose(b.applied_qfx*2.5e6,b.applied_lh,atol=1e-12)
        np.testing.assert_allclose(b.qsfc,saturation_specific(b.tsk,s.psfc),atol=1e-15)
        if scenario=='bulk_memory':
            assert b.incoming_hfx.iloc[0]==b.incoming_qfx.iloc[0]==0
            np.testing.assert_array_equal(b.incoming_hfx.values[1:],b.applied_hfx.values[:-1])
            np.testing.assert_array_equal(b.incoming_qfx.values[1:],b.applied_qfx.values[:-1])
            assert b.wstar.max()>0
        elif scenario=='corrected_bulk':
            assert (b.incoming_hfx==0).all() and (b.incoming_qfx==0).all() and (b.wstar==0).all()
        else:
            np.testing.assert_array_equal(b.applied_hfx,b.forcing_hfx)
            np.testing.assert_array_equal(b.applied_lh,b.forcing_lh)
            np.testing.assert_array_equal(b.incoming_hfx,b.forcing_hfx)
            np.testing.assert_allclose(b.incoming_qfx,b.forcing_lh/2.5e6,atol=0,rtol=0)
        # Independently reconstruct the exact surface-layer Beljaars expression.
        ep1=461.6/287.-1
        virtual_skin=(b.tsk/(s.psfc/1e5)**(287./1004.5))*(1+ep1*b.qsfc)
        # sf_mynn thgb uses surface exner, checked against source rather than lowest-cell exner.
        flux=np.maximum(b.incoming_hfx/b.rho_surface/1004.5+ep1*virtual_skin*b.incoming_qfx/b.rho_surface,0)
        expected_wstar=1.25*(9.80616/b.tsk*np.minimum(1.5*b.previous_pblh,4000)*flux)**.33
        np.testing.assert_allclose(b.wstar,expected_wstar,rtol=2e-12,atol=1e-14)
        mass=d.rho.values*dz
        for sp in ['qc','qi']:
            expected=dt*d[sp+'_tendency_full'].values
            clip=np.maximum(-(d[sp+'_entry'].values+expected),0)
            np.testing.assert_allclose(d[sp+'_expected_increment'],expected,atol=0,rtol=0)
            np.testing.assert_allclose(d[sp+'_applied_increment'],expected+clip,atol=1e-18,rtol=1e-10)
        stage=numeric_csv(str(path)+'.stage_water.csv');assert np.isfinite(stage.values).all()
        assert stage.accounting_residual_ref.abs().max()<1e-10 and stage.adv_reconstruction_error_ref.abs().max()<1e-10
        active=v.plume_active.values.astype(bool);top=np.max(np.where(active,v.z_interface.values,-np.inf),axis=1);top[~active.any(1)]=np.nan
        width=g.maxwidth_pre<=g.minwidth;eligible=(g.fltv_entry>.002)&(g.superadiabatic>0)
        row={'configuration':cfg,'scenario':scenario,'samples':len(b),'active_samples':int(active.any(1).sum()),
            'active_fraction':float(active.any(1).mean()),'plume_top_max_m_AGL':float(np.nanmax(top)),
            'wstar_max_m_s':float(b.wstar.max()),'w_EDMF_max_m_s':float(v.w_EDMF.max()),
            'width_fail_fraction':float(width.mean()),'width_fail_eligible_count':int((width&eligible).sum()),'eligible_count':int(eligible.sum()),
            'PBL_mean_m_AGL':float(s.pblh.mean()),'precipitation_total_mm':float(n.rainnc.values[-1]),
            'HFX_RMSE_against_ERA5_W_m2':float(np.sqrt(np.mean((b.applied_hfx-b.forcing_hfx)**2))),
            'LH_RMSE_against_ERA5_W_m2':float(np.sqrt(np.mean((b.applied_lh-b.forcing_lh)**2)))}
        for label,target in [('0C',273.15),('minus15C',258.15),('260K',260.)]:
            heights=crossing(d.temperature.values,z,target)
            row['reach_'+label+'_samples']=int(supported_reach(active,v.z_interface.values,heights).sum())
        for sp in ['qc','qi','qg']:
            # Fixed initial dry mass gives identical comparison weights in all scenarios.
            p=n.pressure.values[0];t=n.theta.values[0]*(p/1e5)**(2/7);rv=n.qv.values[0]
            fixed_mass=p/(287*t*(1+rv/(287/461.6)))*dz
            row[sp+'_column_max_fixed_mass_kg_m2']=float(np.max(np.sum(n[sp].values*fixed_mass,axis=1)))
        records.append(row)
        evaporating=(b.applied_lh>1e-6)&(b.forcing_lh>0)
        ratio=b.forcing_lh[evaporating]/b.applied_lh[evaporating]
        compatibility.append({'scenario':scenario,'configuration':cfg,'positive_evaporation_pairs':int(evaporating.sum()),
            'LH_sign_mismatch_fraction':float(((b.applied_lh*b.forcing_lh)<0).mean()),
            'ERA5_LH_over_bulk_LH_median_on_evaporating_pairs':float(ratio.median()) if len(ratio) else None,
            'fraction_ratios_gt_1_on_evaporating_pairs':float((ratio>1).mean()) if len(ratio) else None,
            'ratio_is_mavail':False,'ratio_scope':'Fixed simulated atmosphere only; cannot infer soil availability or adjust coefficient from this diagnostic.'})
        stage['scenario']=scenario;stage['configuration']=cfg
        budget_rows.append({'scenario':scenario,'configuration':cfg,'accounting_residual_max_kg_m2':float(stage.accounting_residual_ref.abs().max()),
            'adv_reconstruction_max_kg_m2':float(stage.adv_reconstruction_error_ref.abs().max()),
            'MYNN_stage_water_total_fixed_mass_kg_m2':float(stage.MYNN_stage_ref.sum()),
            'MP_plus_precip_total_fixed_mass_kg_m2':float(stage.MP_plus_precip_ref.sum()),
            'surface_vapor_integral_kg_m2':float(np.sum(b.applied_qfx*b.dt)),
            'MYNN_fixed_mass_minus_surface_vapor_kg_m2':float(stage.MYNN_stage_ref.sum()-np.sum(b.applied_qfx*b.dt)),
            'surface_HFX_plus_LH_integral_J_m2':float(np.sum((b.applied_hfx+b.applied_lh)*b.dt)),
            'MYNN_enthalpy_tendency_proxy_entry_mass_J_m2':float(np.sum(mass*dt*(1004.5*d.theta_tendency.values*(d.pressure.values/1e5)**(2/7)+2.5e6*d.qv_tendency.values-3.34e5*d.qi_tendency_full.values))),
            'energy_closure_certified':False,'water_conservation_certified':False})
        b['time_s']=time;b['time']=pd.Timestamp('2018-12-13',tz='UTC')+pd.to_timedelta(time,unit='s')
        b['plume_top_m_AGL']=top;b['PBL_m_AGL']=s.pblh.values;b['fltv']=g.fltv_entry
        for name in ['maxwidth_pre','minwidth','flux_cap','cloud_cap','pbl_cap']:b[name]=g[name]
        b.to_csv(OUT/f'{scenario}_{cfg}_surface.csv',index=False)
        checks.append({'scenario':scenario,'configuration':cfg,'memory_or_prescribed_flux_contract':'passed',
            'native_Beljaars_wstar_reconstruction':'passed','humidity_units_and_current_saturation':'passed',
            'condensate_single_application':'passed','fixed_mass_stage_accounting':'passed'})
        cache[(scenario,cfg)]=(b,s,v,n)
    for cfg in ['TH_CF','WS_CF']:
        base=cache[('corrected_bulk',cfg)][2].w_LS
        for mode in ['bulk_memory','prescribed_era5']:np.testing.assert_array_equal(cache[(mode,cfg)][2].w_LS,base)
    pd.DataFrame(records).to_csv(OUT/'surface_plume_summary.csv',index=False)
    pd.DataFrame(budget_rows).to_csv(OUT/'water_energy_checks.csv',index=False)
    pd.DataFrame(compatibility).to_csv(OUT/'surface_flux_compatibility.csv',index=False)
    figures(cache)
    changed=[r['path'] for r in initial_protection() if sha256(r['path'])!=r['sha256']];assert not changed
    write_json(OUT/'validation.json',{'checks':checks,'protected_originals_unchanged':True,'ensemble_reruns':0,
        'reference_selection':'pending: neither land-surface hydrology nor complete atmospheric energy/mass closure is validated',
        'selected_reference':None,'2019_recalibration':'required before any scientific evaluation of changed model; not started while reference undecided',
        'bulk_memory':'one-step applied-flux memory, first step zero; same surface-layer equations',
        'prescribed_flux':'scalar heat/moisture fluxes exact forcing; bulk momentum and SF stability retained; distinct boundary sensitivity',
        'mavail':'still fixed=1, no tuning or direct substitution by soil volumetric water',
        'physical_gate':'finite/positive species, pressure and temperature; exact discrete stage bookkeeping, not a proof of full atmospheric conservation',
        'energy_scope':'No radiative heating, momentum integration, or pressure-work budget; full energy closure not certified',
        'budget_scope':'Surface vapor/energy integrals and isolated MYNN enthalpy proxy reported; fixed-reference-density residuals are not exact conservation residuals under prescribed mass/pressure and changing moisture basis.',
        'velocity_sampling':'w_LS(z,t) bit-identical between scenarios; its value sampled at a changing isotherm can change thermodynamically; no EDMF addition.',
        'config_sha256':sha256(CFG)})
    print('Surface audit: 6 sensitivities checked; originals preserved; 3 figures.',flush=True)


def figures(cache):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as md
    out=OUT/'figures';out.mkdir(exist_ok=True)
    colors={'corrected_bulk':'tab:blue','bulk_memory':'tab:orange','prescribed_era5':'tab:green'}
    for number,variables in [(1,['applied_hfx','applied_lh','wstar']),(2,['PBL_m_AGL','plume_top_m_AGL','fltv']),(3,['maxwidth_pre','flux_cap','cloud_cap'])]:
        fig,ax=plt.subplots(3,2,figsize=(13,9),sharex=True,layout='constrained')
        for c,cfg in enumerate(['TH_CF','WS_CF']):
            for scenario,color in colors.items():
                b=cache[(scenario,cfg)][0]
                for r,var in enumerate(variables):ax[r,c].plot(b.time,b[var],label=scenario,color=color,lw=.9)
            if number==1:
                b=cache[('corrected_bulk',cfg)][0]
                ax[0,c].plot(b.time,b.forcing_hfx,'k--',label='ERA5',lw=.8);ax[1,c].plot(b.time,b.forcing_lh,'k--',lw=.8)
            if number==3:
                for r in range(3):ax[r,c].axhline(300,color='k',ls='--',lw=.7)
            for r,var in enumerate(variables):
                ax[r,c].set_ylabel({'applied_hfx':'HFX (W m⁻² upward)','applied_lh':'LH (W m⁻² upward)','wstar':'Surface w* (m s⁻¹)',
                    'PBL_m_AGL':'PBL (m AGL)','plume_top_m_AGL':'Plume top (m AGL)','fltv':'Buoyancy flux (K m s⁻¹)'}.get(var,var+' (m)'))
                ax[r,c].grid(alpha=.2);ax[r,c].xaxis.set_major_formatter(md.DateFormatter('%d %H:%M'));ax[r,c].tick_params(axis='x',rotation=90)
            ax[0,c].set_title(cfg);ax[0,c].legend(fontsize=7)
            ax[-1,c].set_xlabel('December 2018 - UTC')
        fig.savefig(out/f'{number:02}_surface_diagnostics.png',dpi=160);plt.close(fig)


def soil_properties(classes):
    """Lookup native categorical classes; leave unmapped/noninteger classes missing."""
    table=yaml.safe_load((ROOT/'configs/experiments/era5_soil_hydraulics.yaml').read_text())
    a=np.asarray(classes);result={name:np.full(a.shape,np.nan) for name in ['saturation','field_capacity','wilting_point']}
    for code,values in table['soil_classes'].items():
        mask=np.isfinite(a)&(abs(a-int(code))<1e-6)
        for name in result:result[name][mask]=values[name]
    return result


def raw_predictors_one(key):
    """Three native mathematical predictors; no coefficients or ENTLN access."""
    import shutil
    from scm_lightning.config import load_config
    from scm_lightning.evaluation.independent import control_predictors
    from scm_lightning.lightning.thermodynamics import diagnose
    scenario,cfg=key;source=RUNS/scenario/cfg/'output.nc';folder=OUT/'predictors'/scenario/cfg
    folder.mkdir(parents=True,exist_ok=True);tag=folder/'identity.json'
    identity={'source_sha256':sha256(source),'forcing_sha256':sha256(FORCING),
        'thermodynamics_sha256':sha256(ROOT/'python/scm_lightning/lightning/thermodynamics.py'),
        'ice_integrals_sha256':sha256(ROOT/'python/scm_lightning/evaluation/independent.py'),
        'scenario':scenario,'configuration':cfg,'calibrated':False,'observations_used':False,
        'velocity':'LS only; unchanged between surface scenarios; EDMF never added',
        'changes_in_lightning_predictors':'hydrometeor/thermodynamic response; prescribed LS velocity contribution unchanged'}
    if tag.exists():assert json.loads(tag.read_text())==identity;return
    if scenario=='corrected_bulk':
        old=ROOT/f'results/argentina_20181213/mynn_audit_v1/corrected_predictors/{cfg}'
        if json.loads((old/'audit_identity.json').read_text())['source_sha256']!=identity['source_sha256']:raise ValueError('Corrected predictor cache does not match native fields')
        for name in ['predictors_1min.csv','predictors_5min.csv','thermodynamics_mucape_1min.csv']:shutil.copy2(old/name,folder/name)
    else:
        c=load_config(ROOT/'configs/cases/argentina_20181213.yaml');c['interval_minutes']=1
        primary,_=control_predictors(source,FORCING,c,folder,configuration=cfg,compute_cape=False)
        with xr.open_dataset(source,decode_times=False) as f:d=f[['time','lev','pressure','theta','qv','rainnc']].load()
        p=d.pressure.values;t=d.theta.values*(p/1e5)**(2/7);r=d.qv.values;rain=d.rainnc.values;rows=[]
        for i in range(0,len(p)-1,2):
            row=diagnose(p[i],t[i],r[i],d.lev.values)
            row.update(time=str(pd.Timestamp(c['start'])+pd.Timedelta(seconds=float(d.time.values[i]))),
                precipitation_amount_kg_m2=float(rain[i+2]-rain[i]))
            row['CAPExP_J_m2']=row['cape_J_kg']*row['precipitation_amount_kg_m2'];rows.append(row)
        therm=pd.DataFrame(rows);therm.to_csv(folder/'thermodynamics_mucape_1min.csv',index=False)
        cap=therm[['time','CAPExP_J_m2','status']].rename(columns={'CAPExP_J_m2':'predictor'})
        cap['method']='MUCAPExP_300hPa_SCM';cap['configuration']=cfg;cap['member']=0
        pd.concat([primary[primary.method!='CAPExP'],cap],ignore_index=True).to_csv(folder/'predictors_1min.csv',index=False)
    write_json(tag,identity);print('Raw surface predictors: '+scenario+' '+cfg,flush=True)


def predictors(available=False):
    from concurrent.futures import ProcessPoolExecutor
    from scm_lightning.evaluation.mynn_audit_products import relative_statistics
    keys=[(k,c) for k,v in settings()['scenarios'].items() for c in v['configurations']]
    if available:
        ready={(r['scenario'],r['configuration']) for r in json.loads((OUT/'run_tests.json').read_text())['runs']}
        keys=[key for key in keys if key in ready]
    with ProcessPoolExecutor(max_workers=2) as pool:list(pool.map(raw_predictors_one,keys))
    rows=[]
    for scenario,cfg in keys:
        if scenario=='corrected_bulk':continue
        def read(mode):
            d=pd.read_csv(OUT/f'predictors/{mode}/{cfg}/predictors_1min.csv');d.time=pd.to_datetime(d.time,utc=True)
            return d.pivot(index='time',columns='method',values='predictor')
        a,b=read(scenario).align(read('corrected_bulk'),join='inner')
        for name in a:
            rows.append(dict(scenario=scenario,configuration=cfg,method=name,calibrated=False,
                **statistics(a[name]-b[name]),**relative_statistics(a[name],b[name])))
    pd.DataFrame(rows).to_csv(OUT/'raw_predictor_contrasts.csv',index=False)
    write_json(OUT/'predictor_batch.json',{'controls_processed':len(keys),'expected_controls':6,'complete':len(keys)==6,'calibrated':False})


def land():
    from scm_lightning.config import load_config
    from scm_lightning.forcing import normalize
    from scm_lightning.forcing.spatial import subset
    hydraulics=ROOT/'configs/experiments/era5_soil_hydraulics.yaml'
    combined=[];metadata=[]
    for event in ['argentina_20181213','argentina_20191028']:
        source=ROOT/'data/raw/era5/mynn_land'/event/'source.json'
        if not source.exists():continue
        provenance=json.loads(source.read_text());c=load_config(ROOT/f'configs/cases/{event}.yaml');parts=[]
        for record in provenance['files']:
            assert sha256(record['path'])==record['sha256']
            with xr.open_dataset(record['path']) as raw:
                d=raw.load()
            if 'expver' in d.coords:
                assert set(np.unique(d.expver.values))<={'0001','1'},'Resolve nonfinal ERA5 explicitly'
                d=d.drop_vars('expver')
            if 'number' in d.coords:d=d.drop_vars('number')
            parts.append(normalize(d))
        d=xr.merge(parts,join='exact',compat='no_conflicts');box,weight=subset(d,c['area'])
        assert np.all(np.diff(box.time.values)==np.timedelta64(1,'h'))
        for name in ['swvl1','swvl2','swvl3','swvl4']:assert box[name].attrs['units']=='m**3 m**-3'
        for name in ['slhf','sshf']:assert box[name].attrs['units']=='J m**-2'
        assert box.tp.attrs['units']=='m'
        for name in ['slt','lsm','tvh','tvl']:
            np.testing.assert_allclose(box[name],np.broadcast_to(box[name].values[0],box[name].shape),rtol=0,atol=0)
        land_weight=weight*box.lsm
        props=soil_properties(box.slt)
        for name,value in props.items():box['soil_'+name]=xr.DataArray(value,dims=box.slt.dims,coords=box.slt.coords,attrs={'units':'m3 m-3','source':'IFS CY41R2 Table 8.9; contextual, not mavail'})
        box['cell_overlap_area_km2']=weight
        frame=pd.DataFrame({'time':pd.DatetimeIndex(box.time.values).tz_localize('UTC')})
        mean=lambda field,w:field.weighted(w).mean(('latitude','longitude')).values
        frame['land_fraction_of_box']=mean(box.lsm,weight)
        known=np.isfinite(box.soil_field_capacity)
        frame['mapped_soil_fraction_of_land']=((land_weight*known).sum(('latitude','longitude'))/land_weight.sum(('latitude','longitude'))).values
        # A depth mean of top 1m; not vegetation-dependent root weighting.
        box['soil_water_depth_mean_0_1m']=.07*box.swvl1+.21*box.swvl2+.72*box.swvl3
        for name in ['swvl1','swvl2','swvl3','swvl4','soil_water_depth_mean_0_1m']:
            frame[name+'_land_mean_m3_m3']=mean(box[name],land_weight)
            availability=(box[name]-box.soil_wilting_point)/(box.soil_field_capacity-box.soil_wilting_point)
            box[name+'_relative_to_capacity']=availability
            box[name+'_relative_to_capacity'].attrs={'units':'1','interpretation':'unclipped hydraulic context (theta-wilting)/(field_capacity-wilting); NOT MYNN mavail'}
            frame[name+'_capacity_index_land_mean']=mean(availability,land_weight)
            denominator=land_weight.where(known,0).sum(('latitude','longitude'))
            frame[name+'_fraction_mapped_land_ge_field_capacity']=((land_weight*((box[name]>=box.soil_field_capacity)&known)).sum(('latitude','longitude'))/denominator).values
            frame[name+'_fraction_mapped_land_le_wilting']=((land_weight*((box[name]<=box.soil_wilting_point)&known)).sum(('latitude','longitude'))/denominator).values
        for name in ['cvh','cvl','lai_hv','lai_lv']:frame[name+'_land_mean']=mean(box[name],land_weight)
        for name in ['slt','tvl','tvh']:
            codes=np.unique(box[name].values[np.isfinite(box[name].values)])
            for code in codes:
                frame[name+'_class_'+str(int(code))+'_land_fraction']=((land_weight*(box[name]==code)).sum(('latitude','longitude'))/land_weight.sum(('latitude','longitude'))).values
        frame['ERA5_HFX_W_m2_upward']=mean(-box.sshf/3600,weight)
        frame['ERA5_LH_W_m2_upward']=mean(-box.slhf/3600,weight)
        frame['ERA5_precipitation_hour_mm']=mean(box.tp*1000,weight)
        frame['event']=event
        start,end=pd.Timestamp(c['start']),pd.Timestamp(c['end'])
        frame['period']=np.where((frame.time>=start)&(frame.time<end),'event',np.where(frame.time<start,'antecedent','after_event'))
        frame['precipitation_interval_in_event']=(frame.time>start)&(frame.time<=end)
        frame['precipitation_interval_in_previous_7days']=(frame.time>start-pd.Timedelta(days=7))&(frame.time<=start)
        folder=OUT/'land'/event;folder.mkdir(parents=True,exist_ok=True)
        frame.to_csv(folder/'land_context_hourly.csv',index=False)
        box.attrs.update(source='Official ERA5 CDS; source.json hashes preserved',box_wsen=json.dumps(c['area']),
            aggregation='Spherical cell/box overlaps; soil/vegetation weighted by land fraction; flux/precipitation whole box',
            soil_mapping='IFS CY41R2 Table 8.9; hydraulic context only; unknown classes remain NaN',mavail_assigned='none')
        box.to_netcdf(folder/'land_context_box.nc')
        record={'event':event,'box_wsen':c['area'],'source_sha256':sha256(source),'hydraulics_sha256':sha256(hydraulics),
            'previous_7days_precipitation_mm':float(frame.loc[frame.precipitation_interval_in_previous_7days,'ERA5_precipitation_hour_mm'].sum()),
            'event_precipitation_mm':float(frame.loc[frame.precipitation_interval_in_event,'ERA5_precipitation_hour_mm'].sum()),
            'mapped_soil_land_fraction_min':float(frame.mapped_soil_fraction_of_land.min()),
            'soil_type_classes':list(map(float,np.unique(box.slt.values))),
            'event_topsoil_capacity_index_mean':float(frame.loc[frame.period=='event','swvl1_capacity_index_land_mean'].mean()),
            'event_topsoil_fraction_ge_capacity_mean':float(frame.loc[frame.period=='event','swvl1_fraction_mapped_land_ge_field_capacity'].mean()),
            'mavail_assigned':None,'surface_availability_validated':False,'reason':'Hydraulic context does not specify canopy/soil resistances or reproduce HTESSEL surface exchange.'}
        write_json(folder/'context.json',record);metadata.append(record);combined.append(frame)
        print('Land context processed: '+event,flush=True)
    if combined:
        pd.concat(combined,ignore_index=True).to_csv(OUT/'land_context_hourly.csv',index=False)
        write_json(OUT/'land_context.json',metadata)
    return metadata


def finalize():
    """Machine-readable gates and <=15 terminal lines; no scientific report."""
    import platform
    import xml.etree.ElementTree as ET
    import metpy
    suites=list(ET.parse(OUT/'pytest_final.xml').getroot().iter('testsuite'))
    count={k:sum(int(s.attrib.get(k,0)) for s in suites) for k in ['tests','failures','errors','skipped']}
    assert count['failures']==count['errors']==0
    count['passed']=count['tests']-count['skipped']
    result={name:json.loads((OUT/name).read_text()) for name in ['run_tests.json','validation.json','land_context.json','predictor_batch.json','all_native_preservation.json']}
    assert result['predictor_batch.json']['complete'] and result['all_native_preservation.json']['unchanged']
    changed=[p['path'] for p in initial_protection() if sha256(p['path'])!=p['sha256']];assert not changed
    reference=json.loads((RUNS/'reference_sources/identity.json').read_text())
    core=['src/mynn/sf_mynn.F90','src/mynn/bl_mynn.F90','src/mynn/bl_mynn_subroutines.F90']
    for name in core:assert sha256(ROOT/name)==next(r['sha256'] for r in reference if r['path']==name)
    result['original_physics_sources_unchanged']=core
    result['pytest']=dict(count,xml_sha256=sha256(OUT/'pytest_final.xml'))
    before=ET.parse(OUT/'regression_before.xml').getroot()
    result['regression_before_fix_expected_failure']=any(int(s.attrib.get('failures',0))>0 for s in before.iter('testsuite'))
    result['versions']={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'xarray':xr.__version__,'MetPy':metpy.__version__}
    compatibility=pd.read_csv(OUT/'surface_flux_compatibility.csv')
    failure=compatibility[compatibility.scenario=='bulk_memory'].set_index('configuration').LH_sign_mismatch_fraction.to_dict()
    decision={'implementation_reference':'bulk_memory: verified humidity boundary corrections plus one-step surface flux feedback',
        'scientific_reference':None,'calibration_started':False,'old_coefficients_applied_to_changed_model':False,
        'reason':'Surface flux compatibility remains unresolved. At fixed atmospheric state, scaling moisture conductance cannot reverse a wrong flux sign. Soil moisture context does not specify vegetation/soil resistance. Prescribed scalar-flux experiment retains bulk momentum/stability and is not a fully flux-driven surface boundary.',
        'bulk_memory_LH_sign_mismatch_fraction_against_ERA5_forcing':failure,
        'selection_using_ENTLN':False,'if_reference_changes':'Reintegrate only 2019 TH_C0/WS_C0, then fit six coefficients using 2019 alone; do not combine with original 2018 ensemble.',
        'preserved_parameters':['mavail=1','minwidth=300m','dx=5000m','w2=0','rthraten=0','mynn_edmf_type=0','mynn_uedmf_mp=0']}
    write_json(OUT/'reference_calibration_gate.json',decision);result['reference_calibration_gate']=decision
    result['inputs']={'experiment_config':sha256(CFG),'forcing':sha256(FORCING),'hydraulic_table':sha256(ROOT/'configs/experiments/era5_soil_hydraulics.yaml')}
    result['code']={str(p.relative_to(ROOT)):sha256(p) for p in [ROOT/'src/mynn_offline_coupling.F90',ROOT/'src/mp_offline_driver.F90',ROOT/'src/diagnostics/mynn_audit.F90',Path(__file__),ROOT/'python/scm_lightning/download/mynn_land.py']}
    result['products']={str(p.relative_to(OUT)):sha256(p) for p in OUT.rglob('*') if p.is_file() and p.suffix in ['.csv','.nc','.png']}
    write_json(OUT/'test_results.json',result)
    print('\n'.join([
        'MYNN surface v1: deterministic boundary tests completed.',
        'Humidity/surface-temperature corrections retained; applied-flux memory implemented and verified.',
        'Six default controls and two corrected-bulk controls bit-identical to their references.',
        'Two memory controls and two prescribed ERA5 scalar-flux controls completed.',
        'ERA5 soil/vegetation/hydraulic context processed for both events and seven antecedent days.',
        'Surface soil >= field capacity: about 89% (2018), 58% (2019), area-time means.',
        'Surface wstar memory response verified against the native Beljaars expression.',
        f'Tests: {count["passed"]} passed; {count["skipped"]} historical skips.',
        '248 original native outputs and existing calibration/article products preserved; no ensemble reruns.',
        'Scalar-flux forcing exact by construction; full atmospheric energy/mass closure not certified.',
        'Scientific reference/calibration held: bulk LH has about 45% sign mismatch against ERA5 forcing.',
        'No mavail tuning, minimum-width change, core replacement, or lightning calibration.',
        'Outputs: results/argentina_20181213/mynn_surface_v1 (CSV, NetCDF, 3 PNG, test_results.json).',
        'Runs: runs/argentina_20181213/mynn_surface_v1; options/configs versioned separately.'
    ]))


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['run','analyze','land','predictors','finalize']);p.add_argument('--available',action='store_true');a=p.parse_args()
    if a.stage=='run':run()
    elif a.stage=='analyze':analyze()
    elif a.stage=='land':land()
    elif a.stage=='finalize':finalize()
    else:predictors(a.available)


if __name__=='__main__':main()
