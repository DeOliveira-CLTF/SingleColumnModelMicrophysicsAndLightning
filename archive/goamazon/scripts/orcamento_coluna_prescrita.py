"""Stage-observed water budget, retaining the prescribed advective experiment.

No closure term is silently assigned to a physical process. Native microphysics
residual includes sedimentation/limiters; process attribution remains separate.
"""
from datetime import datetime, timezone
import pandas as pd
from scm_artigo import *
from executar_experimentos import linux_command, linux_path, write_nml
from controlar_perturbacao_inicial import project
from diagnosticar_ensemble_goamazon import CASES, reconstruct_ep

OUT = RESULTS/'ensemble_goamazon/orcamento_prescrito'


def transport(q, w, z, edges, mass, dt):
    """Independent discrete product-rule audit, accumulated every CFL substep.

    Interface mass flux is linearly interpolated rho_ref*w, with endpoint
    extension. Boundary water donor is zero-gradient, matching the driver.
    The stencil correction is explicitly numerical, not horizontal forcing.
    """
    old=q.copy(); n=len(z); up=np.arange(n)-1
    up[w<0]=np.arange(n)[w<0]+1
    ix=np.flatnonzero((up>=0)&(up<n))
    rate=np.abs(w[ix])/np.abs(z[ix]-z[up[ix]])
    ns=max(1,int(np.ceil(np.max(rate)*dt/.8))); ds=dt/ns
    air=np.interp(edges,z,mass/np.diff(edges)*w)
    bottom=top=divergence=stencil=0.
    for _ in range(ns):
        prev=q.copy(); tendency=np.zeros(n)
        tendency[ix]=rate*(prev[up[ix]]-prev[ix])
        donor=np.r_[prev[0],np.where(air[1:-1]>=0,prev[:-1],prev[1:]),prev[-1]]
        flux=air*donor
        bottom+=ds*flux[0]; top+=ds*flux[-1]
        divergence+=ds*np.sum(prev*np.diff(air))
        stencil+=ds*np.sum(mass*tendency+np.diff(flux)-prev*np.diff(air))
        c=rate*ds; q[ix]=(1-c)*prev[ix]+c*prev[up[ix]]
    actual=float(np.dot(mass,q-old))
    return dict(bottom_in_mm=bottom,top_out_mm=top,
                prescribed_divergence_mm=divergence,stencil_correction_mm=stencil,
                reconstructed_transport_mm=actual,
                transport_closure_mm=actual-(bottom-top+divergence+stencil))


def audit(folder,cfg):
    ds=load(folder/'profile_1.nc'); frame=pd.read_csv(folder/'profile_water_budget_1.csv')
    projection=pd.read_csv(folder/'projection.csv'); mass=projection.reference_dry_mass_kg_m2.values
    edges=np.r_[projection.bottom_m.values,projection.top_m.values[-1]]
    _,w,*_=reconstruct_ep(ds,cfg)
    flux=pd.DataFrame([transport(ds.qv.values[i-1].copy(),w[i],ds.lev.values,edges,mass,cfg['dt_microp'])
                       for i in range(1,len(ds.time))])
    np.testing.assert_allclose(flux.reconstructed_transport_mm,frame.transport_ref_mm,atol=2e-13,rtol=1e-10)
    frame=pd.concat([frame,flux],axis=1)
    frame['micro_native_residual_mm']=frame.micro_native_mm+frame.precip_mm
    frame['bookkeeping_residual_mm']=(frame.storage_after_mm-frame.storage_before_mm+frame.precip_mm
        -frame.transport_ref_mm-frame.micro_native_residual_mm-frame.weight_conversion_mm-frame.density_remap_mm)
    assert np.max(np.abs(frame.bookkeeping_residual_mm))<1e-11
    assert np.max(np.abs(frame.transport_closure_mm))<1e-11
    assert np.isfinite(frame.values).all()
    np.testing.assert_allclose(frame.precip_mm,np.diff(ds.rainnc.values),atol=1e-13)
    frame.to_csv(folder/'budget_audit.csv',index=False,float_format='%.17g')
    summary={key:float(frame[key].sum()) for key in frame if key.endswith('_mm') and not key.startswith('storage_')}
    summary['storage_change_mm']=float(frame.storage_after_mm.iloc[-1]-frame.storage_before_mm.iloc[0])
    summary['max_step_bookkeeping_residual_mm']=float(frame.bookkeeping_residual_mm.abs().max())
    summary['max_step_native_residual_mm']=float(frame.micro_native_residual_mm.abs().max())
    summary['native_residual_over_precip_percent']=100*summary['micro_native_residual_mm']/summary['precip_mm']
    return summary


def main():
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    root=OUT/stamp; root.mkdir(parents=True)
    build=root/'build'
    run=linux_command(['make','-j4','BUILD='+build.relative_to(ROOT).as_posix()],ROOT)
    (root/'build.log').write_text(run.stdout+run.stderr,encoding='utf-8')
    binary=build/'ensemble_driver'; source=load(INPUT); input_hash=sha(INPUT)
    baseline=json.loads((RESULTS/'ensemble_goamazon/manifest.json').read_text(encoding='utf-8'))
    rows=[]; records={}
    for eid in CASES:
        cfg=dict(baseline['members'][eid]['config'])
        folder=root/eid
        perturbation,projection=project(source,60,cfg['sat_z_lo'],cfg['sat_z_hi'],folder)
        cfg.update(path_goamazon=linux_path(INPUT),path_initial_perturbation=linux_path(perturbation),write_water_budget=True)
        write_nml(folder/'ensemble.nml',cfg,'ensemble')
        run=linux_command([linux_path(binary)],folder)
        (folder/'execution.log').write_text(run.stdout+run.stderr,encoding='utf-8')
        ds=load(folder/'profile_1.nc')
        for var in ds.variables: assert np.isfinite(ds[var].values).all(),var
        for var in ['qv','qc','qr','qi','qs','qg']: assert float(ds[var].min())>=0,var
        row=dict(experiment=eid,**audit(folder,cfg),**projection)
        # All prognostics must be bitwise equal with diagnostics disabled.
        off=folder/'diagnostics_off'; off.mkdir()
        write_nml(off/'ensemble.nml',dict(cfg,write_water_budget=False),'ensemble')
        run=linux_command([linux_path(binary)],off)
        (off/'execution.log').write_text(run.stdout+run.stderr,encoding='utf-8')
        other=load(off/'profile_1.nc')
        for var in ds.variables: np.testing.assert_array_equal(ds[var].values,other[var].values)
        row['diagnostic_toggle_bit_identical']=True
        rows.append(row)
        records[eid]=dict(directory=str(folder.relative_to(ROOT)),config=cfg,output_sha256=sha(folder/'profile_1.nc'))
        print(eid,'native residual mm:',row['micro_native_residual_mm'],flush=True)
    assert sha(INPUT)==input_hash
    pd.DataFrame(rows).to_csv(root/'summary.csv',index=False,float_format='%.17g')
    manifest=dict(directory=str(root.relative_to(ROOT)),created_utc=stamp,input_sha256=input_hash,
                  binary_sha256=sha(binary),sources={p.relative_to(ROOT).as_posix():sha(p) for p in (ROOT/'src').rglob('*') if p.is_file()},
                  runs=records,design='Prescribed column, controlled dose; diagnostic-only stage snapshots and independent transport reconstruction',
                  caution='Bookkeeping closure is not physical conservation; native microphysics residual remains explicit.')
    save_json(root/'manifest.json',manifest); save_json(OUT/'latest.json',manifest)
    print(root,flush=True)


if __name__=='__main__': main()
