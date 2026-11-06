"""Water perturbation audit, vertical EP diagnostics and crossed resolution study.

No new figures, no parameterization/input edits. Uses the verified ensemble binary.
"""
from datetime import datetime, timezone
import subprocess
import pandas as pd
import xarray as xr
from scm_artigo import *
from executar_experimentos import linux_command, linux_path, write_nml
from executar_ensemble_goamazon import check

BASE=RESULTS/'ensemble_goamazon'
OUT=BASE/'diagnostico_cientifico'
CASES=['wsm6_media_w3','thompson_media_w3','thompson_alta_w3','wsm6_alta_w3']


def initial_water(ds,cfg):
    src=load(INPUT); z=ds.lev.values; dz=np.diff(interfaces(z,cfg['ztop']))
    old={k:np.interp(z,src.lev.values,src[v].values.squeeze()) for k,v in
         [('qv','rv'),('qc','rl'),('qi','ri'),('theta','theta'),('p','pa')]}
    t=old['theta']*(old['p']/P0)**(RD/CP)
    # Exact ideal-gas dry-air density for mixing ratios per kg dry air.
    dry=old['p']/(t*(RD+RV*old['qv']))
    model=old['p']/(RD*t*(1+.61*old['qv']))
    delta={k:ds[k].values[0]-old[k] for k in ['qv','qc','qi']}
    dq=sum(delta.values()); weight=dry*dz
    added=float(np.dot(weight,dq))
    rho_after=ds.pressure.values[0]/(t*(RD+RV*ds.qv.values[0]))
    rt0=old['qv']+old['qc']+old['qi']; rt1=sum(ds[k].values[0] for k in ['qv','qc','qi'])
    rain=float(ds.rainnc.values[-1]); mask=(z>=cfg['sat_z_lo'])&(z<=cfg['sat_z_hi'])
    return {'water_added_fixed_dry_air_mm':added,
            'water_added_model_density_mm':float(np.dot(model*dz,dq)),
            'water_inventory_delta_eulerian_mm':float(np.dot(rho_after*rt1-dry*rt0,dz)),
            'water_initial_original_fixed_dry_mm':float(np.dot(weight,rt0)),
            'water_initial_perturbed_fixed_dry_mm':float(np.dot(weight,rt1)),
            'vapor_change_fixed_dry_mm':float(np.dot(weight,delta['qv'])),
            'cloud_liquid_change_fixed_dry_mm':float(np.dot(weight,delta['qc'])),
            'cloud_ice_change_fixed_dry_mm':float(np.dot(weight,delta['qi'])),
            'latent_energy_proxy_delta_MJ_m2':float(np.dot(weight,2.5e6*delta['qv']-3.5e5*delta['qi'])/1e6),
            'rain_final_mm':rain,'added_over_rain_percent':100*added/rain if rain>0 else None,
            'rain_over_added':rain/added if added>0 else None,
            'saturated_cells':int(mask.sum()),'saturated_center_min_m':float(z[mask].min()),
            'saturated_center_max_m':float(z[mask].max()),
            'saturated_thickness_m':float(dz[mask].sum())},weight


def reconstruct_ep(ds,cfg):
    """Literal postprocessing of existing EP, including its endpoint-OR mask.
    w at output time t_n is the w used at t_(n-1), before advancing the step.
    """
    z=ds.lev.values; time=ds.time.values; source=load(INPUT)
    background=np.interp(z,source.lev.values,source.wa.values[0])
    used_time=np.maximum(time-cfg['dt_microp'],0)
    shape=np.where(used_time<=cfg['pulse_duration_s'],.5*(1-np.cos(2*3.14159265358979*used_time/cfg['pulse_duration_s'])),0)
    w=background[None,:]+cfg['pulse_w_max']*shape[:,None]*np.exp(-.5*((z-(cfg['sat_z_lo']+cfg['sat_z_hi'])/2)/cfg['pulse_sigma_z'])**2)
    T,rho=thermo(ds); tc=T-273.15; qg=ds.qg.values; qs=ds.qs.values; qc=ds.qc.values
    pi=3.14159265358979
    n0=np.clip(200/np.maximum(qg,1e-12),1e4,5e6)
    lam=(6*pi*400*n0/(rho*np.maximum(qg,1e-30)))**.25
    ng=np.where(qg>1e-12,n0/lam,0); dg=np.where(qg>1e-12,1/lam,.002)
    ni=rho*qs/((pi/6)*900*(2e-4)**3)
    vg=3*np.maximum(dg*1000,.001)**.15+1.5
    trev=np.clip(-15+6*qc*1000*rho,-20,-3)
    collision=np.exp(-.5*((tc+15)/8)**2)
    density=ng*ni*np.abs(vg-.5)*pi*(.5*dg+.0001)**2*collision*np.tanh(np.abs(tc-trev)/5)*(1+np.maximum(w,0)/5)
    mixed=(tc>=-40)&(tc<=0)
    edges=mixed[:,1:]|mixed[:,:-1]
    intervals=.5*(density[:,1:]+density[:,:-1])*np.diff(z)*edges
    ep=intervals.sum(axis=1); flash=.001*np.sqrt(ep)
    np.testing.assert_allclose(flash[1:],ds.flash_rate.values[1:],rtol=5e-11,atol=1e-12)
    outside=.5*(density[:,1:]*(~mixed[:,1:])+density[:,:-1]*(~mixed[:,:-1]))*np.diff(z)*edges
    return T,w,density,intervals,ep,outside.sum(axis=1)


def detailed(eid,ds,cfg):
    folder=OUT/'perfis'; folder.mkdir(exist_ok=True)
    T,w,ep_density,segments,ep,outside=reconstruct_ep(ds,cfg)
    z=ds.lev.values; dz=np.diff(interfaces(z,cfg['ztop'])); t,rho=thermo(ds)
    water,weight=initial_water(ds,cfg); m,series=metrics(ds,cfg)
    mixed=(T>=233.15)&(T<=273.15); step=np.diff(ds.time.values)
    series['ep_reconstructed']=ep; series['ep_outside_center_temperature']=outside
    series['w_max_used_m_s']=w.max(axis=1)
    rt=sum(ds[q].values for q in ['qv']+SPECIES)
    series['water_column_fixed_initial_dry_mm']=(rt*weight).sum(axis=1)
    series['apparent_water_residual_fixed_dry_mm']=series['water_column_fixed_initial_dry_mm']-series['water_column_fixed_initial_dry_mm'][0]+ds.rainnc.values
    # Reapply the exact, linear upwind update to the PREVIOUS saved vapor.
    # This diagnoses imposed transport under a fixed initial dry-air weight;
    # it is not a tracer or a closed Eulerian mass budget.
    advection=np.zeros(len(ds.time))
    for i in range(1,len(ds.time)):
        old=ds.qv.values[i-1].copy(); advected=old.copy(); speed=w[i]
        up=np.arange(len(z))-1; up[speed<0]=np.arange(len(z))[speed<0]+1
        valid=(up>=0)&(up<len(z)); indices=np.flatnonzero(valid)
        courant=np.abs(speed[valid])*cfg['dt_microp']/np.abs(z[valid]-z[up[valid]])
        nsub=max(1,int(np.ceil(courant.max()/.8))); c=courant/nsub
        for _ in range(nsub):
            prev=advected.copy(); advected[indices]=(1-c)*prev[indices]+c*prev[up[valid]]
        advection[i]=np.dot(weight,advected-old)
    series['water_advection_cumulative_fixed_dry_mm']=np.cumsum(advection)
    series['water_remainder_after_advection_mm']=series['apparent_water_residual_fixed_dry_mm']-np.cumsum(advection)
    for q in ['qc','qi','qs','qg']:
        a=ds[q].values; series[q+'_full_column_kg_m2']=(a*rho*dz).sum(axis=1)
        series[q+'_global_max_g_kg']=a.max(axis=1)*1000
    pd.DataFrame(series).to_csv(folder/f'{eid}_series.csv',index=False,float_format='%.17g')
    derived=ds[['qc','qi','qs','qg','qr','qv','theta','pressure','flash_rate','rainnc']].copy()
    derived['temperature']=(('time','lev'),T,{'units':'K'})
    derived['w_used']=(('time','lev'),w,{'units':'m s-1','description':'reconstructed from initial wa plus pulse at previous timestep'})
    derived['ep_local_density']=(('time','lev'),ep_density,{'units':'nominal','description':'literal EP integrand, before mixed-layer interval mask'})
    derived.attrs.update(source_profile='Original archived ensemble member '+eid,
                         warning='Idealized experiment; w reconstructed, not a saved prognostic output; no process budgets')
    derived.to_netcdf(folder/f'{eid}_diagnostico.nc',engine='netcdf4')
    snapshots=[]; eventrows=[]
    times={'inicial':0.,'pico_pulso':900.,'fim_pulso':1800.,'final':7200.}
    for name,key in [('inicio_EP','flash_onset_h'),('pico_EP','flash_peak_h'),('pico_qg','qg_peak_h'),('pico_chuva','rain_peak_h')]:
        if m[key] is not None: times[name]=m[key]*3600
    for event,sec in times.items():
        i=int(np.argmin(np.abs(ds.time.values-sec)))
        snapshots.append(pd.DataFrame({'event':event,'time_s':float(ds.time.values[i]),'z_m':z,
                                      'temperature_C':T[i]-273.15,'w_used_m_s':w[i],
                                      **{q+'_g_kg':ds[q].values[i]*1000 for q in ['qc','qi','qs','qg']},
                                      'ep_local_density':ep_density[i],'flash_min':float(ds.flash_rate.values[i])}))
        j=int(np.argmax(ep_density[i])); active=ds.qg.values[i]>1e-12
        eventrows.append({'experiment':eid,'event':event,'time_s':float(ds.time.values[i]),'flash_min':float(ds.flash_rate.values[i]),
                          'ep_dominant_z_m':float(z[j]) if ep[i]>0 else None,'ep_dominant_temperature_C':float(T[i,j]-273.15) if ep[i]>0 else None,
                          'ep_dominant_w_m_s':float(w[i,j]) if ep[i]>0 else None,
                          **{q+'_at_ep_dominant_g_kg':float(ds[q].values[i,j]*1000) if ep[i]>0 else None for q in ['qc','qi','qs','qg']},
                          'ep_outside_centers_percent':float(100*outside[i]/ep[i]) if ep[i]>0 else None,
                          'qg_extent_min_m':float(z[active].min()) if active.any() else None,
                          'qg_extent_max_m':float(z[active].max()) if active.any() else None})
    pd.concat(snapshots).to_csv(folder/f'{eid}_eventos_verticais.csv',index=False,float_format='%.17g')
    coexist=mixed&(ds.qg.values>1e-12)&(ds.qs.values>1e-12)
    m.update(water)
    m['qs_qg_same_mixed_cell_duration_min']=float(np.sum(step[np.any(coexist,axis=1)[1:]])/60)
    m['apparent_water_residual_final_mm']=float(series['apparent_water_residual_fixed_dry_mm'][-1])
    m['water_advection_cumulative_final_mm']=float(np.sum(advection))
    m['water_remainder_after_advection_final_mm']=float(series['water_remainder_after_advection_mm'][-1])
    return dict(experiment=eid,**m),eventrows


def main():
    OUT.mkdir(exist_ok=True); stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    original=json.loads((BASE/'manifest.json').read_text(encoding='utf-8'))
    binary=ROOT/original['binary']
    if sha(binary)!=original['binary_sha256'] or sha(INPUT)!=original['identity']['input_sha256']: raise ValueError('Baseline hashes differ')
    for p,h in original['identity']['sources'].items():
        if sha(ROOT/p)!=h: raise ValueError('Source changed: '+p)
    waterrows=[]; detailrows=[]; events=[]
    for eid,record in original['members'].items():
        path=ROOT/record['directory']/'profile_1.nc'
        if sha(path)!=record['output_sha256']: raise ValueError('Output hash changed: '+eid)
        ds=load(path); cfg=record['config']; wr,_=initial_water(ds,cfg)
        waterrows.append(dict(experiment=eid,**wr))
        if eid in CASES:
            row,ev=detailed(eid,ds,cfg); detailrows.append(row); events+=ev
    pd.DataFrame(waterrows).to_csv(OUT/'agua_introduzida_12_membros.csv',index=False,float_format='%.17g')
    pd.DataFrame(detailrows).to_csv(OUT/'diagnosticos_quatro_membros.csv',index=False,float_format='%.17g')
    pd.DataFrame(events).to_csv(OUT/'eventos_EP.csv',index=False,float_format='%.17g')
    previous=json.loads((OUT/'manifest.json').read_text(encoding='utf-8')) if (OUT/'manifest.json').exists() else {}
    manifest={'baseline_manifest_sha256':sha(BASE/'manifest.json'),'binary_sha256':sha(binary),'input_sha256':sha(INPUT),
              'sources':original['identity']['sources'],'created_utc':stamp,'runs':{},
              'design':'4 cases x nz60/120/240 x dt20/10/5; dz0=3000/nz, fixed top 20000; no new physics'}
    rows=[]
    for eid in CASES:
        orig=original['members'][eid]
        for nz in [60,120,240]:
            for dt in [20.,10.,5.]:
                label=f'{eid}_nz{nz}_dt{int(dt)}'
                cfg=dict(orig['config'],nz=nz,dz0=3000./nz,dt_microp=dt,path_goamazon=linux_path(INPUT))
                if nz==60 and dt==20:
                    run=ROOT/orig['directory']; action='baseline_reused'
                else:
                    prev=previous.get('runs',{}).get(label,{})
                    compatible=(previous.get('binary_sha256')==sha(binary) and previous.get('input_sha256')==sha(INPUT) and prev.get('config')==cfg)
                    if compatible and (ROOT/prev['directory']/'profile_1.nc').exists() and sha(ROOT/prev['directory']/'profile_1.nc')==prev['output_sha256']:
                        run=ROOT/prev['directory']; action='reused'
                    else:
                        run=OUT/f'{label}_{stamp}'; run.mkdir()
                        write_nml(run/'ensemble.nml',cfg,'ensemble')
                        try:
                            result=linux_command([linux_path(binary)],run)
                            (run/'execution.log').write_text(result.stdout+result.stderr,encoding='utf-8')
                        except subprocess.CalledProcessError as e:
                            (run/'execution.log').write_text((e.stdout or '')+(e.stderr or ''),encoding='utf-8'); raise
                        action='executed'
                path=run/'profile_1.nc'
                try:
                    ds,m,series,qa=check(path,cfg)
                    wr,_=initial_water(ds,cfg); m.update(wr)
                    z=ds.lev.values; dz=np.diff(interfaces(z,cfg['ztop']))
                    layer=(z>=cfg['sat_z_lo'])&(z<=cfg['sat_z_hi'])
                    m['dz_in_saturation_min_m']=float(dz[layer].min()); m['dz_in_saturation_max_m']=float(dz[layer].max())
                    row=dict(case=eid,nz=nz,dt_s=dt,status='valid',**m)
                except (ValueError,AssertionError) as error:
                    qa={'valid':False,'error':str(error)}; row=dict(case=eid,nz=nz,dt_s=dt,status='invalid',error=str(error))
                save_json(run/('validation_convergence.json'),qa)
                record={'directory':str(run.relative_to(ROOT)),'action':action,'config':cfg,'output_sha256':sha(path),
                        'namelist_sha256':sha(run/'ensemble.nml'),'quality':qa}
                manifest['runs'][label]=record; rows.append(row)
                save_json(OUT/'manifest.json',manifest); save_json(OUT/f'manifest_{stamp}.json',manifest)
                pd.DataFrame(rows).to_csv(OUT/'sensibilidade_resolucao.csv',index=False,float_format='%.17g')
                print(label,action,row['status'],flush=True)
    print(pd.DataFrame(rows)[['case','nz','dt_s','rain_final_mm','qg_max_g_kg','flash_integral_nominal']].to_string(index=False))


if __name__=='__main__': main()
