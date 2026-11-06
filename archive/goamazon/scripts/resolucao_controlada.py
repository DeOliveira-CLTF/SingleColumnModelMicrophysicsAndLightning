"""Repeat the crossed grid/time study with identical projected initial water."""
from datetime import datetime,timezone
import pandas as pd
from scm_artigo import *
from executar_experimentos import linux_command,linux_path,write_nml
from controlar_perturbacao_inicial import project
from orcamento_coluna_prescrita import audit,OUT as BUDGET
from diagnosticar_ensemble_goamazon import reconstruct_ep

OUT=RESULTS/'ensemble_goamazon/resolucao_controlada'


def main():
    baseline=json.loads((BUDGET/'latest.json').read_text(encoding='utf-8'))
    binary=ROOT/baseline['directory']/'build/ensemble_driver'
    assert sha(binary)==baseline['binary_sha256']
    assert sha(INPUT)==baseline['input_sha256']
    for name,digest in baseline['sources'].items(): assert sha(ROOT/name)==digest,name
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    root=OUT/stamp; root.mkdir(parents=True); source=load(INPUT)
    rows=[]; records={}
    for eid,base in baseline['runs'].items():
        for nz in [60,120,240]:
            for dt in [20.,10.,5.]:
                label=f'{eid}_nz{nz}_dt{int(dt)}'
                if nz==60 and dt==20:
                    folder=ROOT/base['directory']; cfg=base['config']
                    assert sha(folder/'profile_1.nc')==base['output_sha256']
                else:
                    folder=root/label; cfg=dict(base['config'],nz=nz,dz0=3000/nz,dt_microp=dt)
                    file,_=project(source,nz,cfg['sat_z_lo'],cfg['sat_z_hi'],folder)
                    cfg['path_initial_perturbation']=linux_path(file)
                    write_nml(folder/'ensemble.nml',cfg,'ensemble')
                    run=linux_command([linux_path(binary)],folder)
                    (folder/'execution.log').write_text(run.stdout+run.stderr,encoding='utf-8')
                ds=load(folder/'profile_1.nc')
                for key in ds.variables: assert np.isfinite(ds[key].values).all(),key
                for key in ['qv']+SPECIES: assert float(ds[key].min())>=0,key
                np.testing.assert_allclose(ds.time.values,np.arange(0,7200+dt/2,dt),atol=1e-10)
                projection=pd.read_csv(folder/'projection.csv')
                dose=float(projection.water_added_kg_m2.sum())
                rv0=np.interp(ds.lev.values,source.lev,source.rv.values.squeeze())
                applied=float(np.dot(ds.qv.values[0]-rv0,projection.reference_dry_mass_kg_m2))
                np.testing.assert_allclose(applied,dose,atol=1e-11,rtol=0)
                budget=audit(folder,cfg); metric,_=metrics(ds,cfg)
                rows.append(dict(experiment=eid,nz=nz,dt_s=dt,initial_added_mm=applied,**budget,**metric))
                records[label]=dict(directory=str(folder.relative_to(ROOT)),config=cfg,output_sha256=sha(folder/'profile_1.nc'))
                pd.DataFrame(rows).to_csv(root/'summary.csv',index=False,float_format='%.17g')
                save_json(root/'runs.json',records)
                print(label,'OK',flush=True)
    table=pd.DataFrame(rows)
    for _,g in table.groupby('experiment'): assert np.ptp(g.initial_added_mm)<1e-11
    manifest=dict(directory=str(root.relative_to(ROOT)),runs=records,input_sha256=sha(INPUT),binary_sha256=sha(binary),
                  baseline_manifest_sha256=sha(BUDGET/'latest.json'),design='4 cases x 3 grids x 3 time steps; fixed native physical vapor source')
    save_json(root/'manifest.json',manifest); save_json(OUT/'latest.json',manifest)
    print(root,flush=True)


if __name__=='__main__': main()
