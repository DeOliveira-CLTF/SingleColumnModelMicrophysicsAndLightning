"""Project a single native-grid water-dose profile onto arbitrary model cells.

Only initialization smoke tests (zero duration). Subsequent transport design
must be settled before scientific integrations or a new resolution study.
"""
from datetime import datetime,timezone
import subprocess
import pandas as pd
from scm_artigo import *
from executar_experimentos import linux_command,linux_path,write_nml
from executar_ensemble_goamazon import saturation,LAYERS

OUT=RESULTS/'ensemble_goamazon/controle_inicial'


def grid(nz,dz0,top):
    lo,hi=1.,1.2
    for _ in range(60):
        mid=(lo+hi)/2
        s=dz0*nz if abs(mid-1)<1e-12 else dz0*(mid**nz-1)/(mid-1)
        if s<top: lo=mid
        else: hi=mid
    edges=np.r_[0,np.cumsum(dz0*((lo+hi)/2)**np.arange(nz))]
    return (edges[1:]+edges[:-1])/2,edges


def integrate_linear(x,y,a,b):
    if b<=a: return 0.
    nodes=np.r_[a,x[(x>a)&(x<b)],b]
    vals=np.interp(nodes,x,y)
    return float(np.sum((vals[1:]+vals[:-1])*.5*np.diff(nodes)))


def project(source,nz,lo,hi,directory):
    z,edge=grid(nz,3000/nz,20000.)
    x=source.lev.values; p=source.pa.values.squeeze(); th=source.theta.values.squeeze(); rv=source.rv.values.squeeze()
    T=th*(p/100000)**(RD/CP); dry=p/(T*(RD+RV*rv))
    # Define once on native 10m grid, then integrate the same piecewise linear
    # volumetric source with exact physical layer boundaries on every grid.
    native=dry*np.maximum(saturation(p,T)-rv,0.)
    target=integrate_linear(x,native,lo,hi)
    pg=np.interp(z,x,p); tg=np.interp(z,x,th)*(pg/100000)**(RD/CP); rvg=np.interp(z,x,rv)
    drymass=pg/(tg*(RD+RV*rvg))*np.diff(edge)
    dose=np.array([integrate_linear(x,native,max(a,lo),min(b,hi)) for a,b in zip(edge[:-1],edge[1:])])
    overlap=np.maximum(0,np.minimum(edge[1:],hi)-np.maximum(edge[:-1],lo))
    frac=overlap/np.diff(edge); delta=dose/drymass
    np.testing.assert_allclose(dose.sum(),target,rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(overlap.sum(),hi-lo,rtol=0,atol=1e-8)
    assert not np.any(delta[frac==0])
    directory.mkdir(parents=True,exist_ok=True)
    file=directory/'perturbation.txt'
    lines=[f'{nz} {target:.17g} {2.5e6*target:.17g}']
    lines += [' '.join(f'{v:.17g}' for v in row) for row in zip(z,delta,drymass,frac)]
    file.write_text('\n'.join(lines)+'\n',encoding='ascii')
    pd.DataFrame(dict(z_m=z,bottom_m=edge[:-1],top_m=edge[1:],overlap_m=overlap,fraction=frac,
                      reference_dry_mass_kg_m2=drymass,water_added_kg_m2=dose,delta_rv_kg_kg=delta)).to_csv(directory/'projection.csv',index=False,float_format='%.17g')
    return file,dict(nz=nz,layer_bottom_m=lo,layer_top_m=hi,target_water_mm=target,
                    projected_water_mm=float(dose.sum()),physical_support_m=float(overlap.sum()),
                    latent_energy_external_MJ_m2=2.5*target,max_delta_rv=float(delta.max()),
                    partially_covered_cells=int(((frac>0)&(frac<1)).sum()))


def main():
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    root=OUT/stamp; root.mkdir(parents=True)
    source=load(INPUT); source_before=sha(INPUT)
    build=root/'build'
    try:
        result=linux_command(['make','-j4','BUILD='+str(build.relative_to(ROOT)).replace('\\','/')],ROOT)
        (root/'build.log').write_text(result.stdout+result.stderr,encoding='utf-8')
    except subprocess.CalledProcessError as e:
        (root/'build.log').write_text((e.stdout or '')+(e.stderr or ''),encoding='utf-8'); raise
    binary=build/'ensemble_driver'
    original=json.loads((RESULTS/'ensemble_goamazon/manifest.json').read_text(encoding='utf-8'))
    base=original['members']['thompson_media_w3']['config']
    rows=[]
    for layer,lo,hi in LAYERS:
        for nz in [60,120,240]:
            folder=root/f'{layer}_nz{nz}'
            file,record=project(source,nz,lo,hi,folder)
            cfg=dict(base,path_goamazon=linux_path(INPUT),path_initial_perturbation=linux_path(file),
                     nz=nz,dz0=3000/nz,sat_z_lo=lo,sat_z_hi=hi,satconfig_label=layer,run_duration_s=0.)
            write_nml(folder/'ensemble.nml',cfg,'ensemble')
            run=linux_command([linux_path(binary)],folder)
            (folder/'execution.log').write_text(run.stdout+run.stderr,encoding='utf-8')
            ds=load(folder/'profile_1.nc'); data=pd.read_csv(folder/'projection.csv')
            z=ds.lev.values; rv0=np.interp(z,source.lev,source.rv.values.squeeze())
            dose=float(np.sum((ds.qv.values[0]-rv0)*data.reference_dry_mass_kg_m2.values))
            np.testing.assert_allclose(dose,record['target_water_mm'],rtol=1e-12,atol=1e-12)
            np.testing.assert_allclose(ds.theta.values[0],np.interp(z,source.lev,source.theta.values.squeeze()),rtol=1e-13)
            for q,src in [('qc','rl'),('qi','ri')]:
                np.testing.assert_allclose(ds[q].values[0],np.interp(z,source.lev,source[src].values.squeeze()),rtol=1e-13,atol=1e-14)
            record.update(layer=layer,applied_water_mm=dose,run_directory=str(folder.relative_to(ROOT)),output_sha256=sha(folder/'profile_1.nc'))
            rows.append(record); print(layer,nz,dose,'verified',flush=True)
    # Regression: empty optional path must preserve the historical integration.
    reg=root/'legacy_regression'; reg.mkdir()
    write_nml(reg/'ensemble.nml',dict(base,path_goamazon=linux_path(INPUT)),'ensemble')
    run=linux_command([linux_path(binary)],reg)
    (reg/'execution.log').write_text(run.stdout+run.stderr,encoding='utf-8')
    a=load(reg/'profile_1.nc'); b=load(ROOT/original['members']['thompson_media_w3']['directory']/'profile_1.nc')
    for key in b.variables: np.testing.assert_array_equal(a[key].values,b[key].values)
    assert sha(INPUT)==source_before
    frame=pd.DataFrame(rows); frame.to_csv(root/'validacao_projecao.csv',index=False,float_format='%.17g')
    for _,group in frame.groupby('layer'): assert np.ptp(group.applied_water_mm)<1e-11
    record={'directory':str(root.relative_to(ROOT)),'input_sha256':source_before,'binary_sha256':sha(binary),
            'source_sha256':sha(ROOT/'src/ensemble_driver.F90'),'legacy_regression_bit_identical':True,
            'initializations_verified':9,'scientific_runs_controlled':0,'energy_treatment':'Fixed temperature, externally supplied vapor; latent source Lv*dM declared, not energy-conserving adjustment',
            'method':'Native 10m piecewise-linear dry-density weighted positive saturation deficit, integrated over exact physical overlap; model-cell mean increments',
            'pending':'Transport/mass-continuity methodological choice before repeating resolution or process attribution'}
    save_json(root/'manifest.json',record); save_json(OUT/'latest.json',record)
    print(frame[['layer','nz','target_water_mm','applied_water_mm','physical_support_m','latent_energy_external_MJ_m2']].to_string(index=False))


if __name__=='__main__': main()
