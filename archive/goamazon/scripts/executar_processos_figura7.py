"""Reproduce Figure 6 members with observational process hooks and archive rates."""
from datetime import datetime,timezone
import subprocess
import argparse
import pandas as pd
import xarray as xr
from scm_artigo import *
from executar_experimentos import linux_command,linux_path,write_nml
from diagnosticar_ensemble_goamazon import initial_water
from instrumentar_processos import PROV,strip_hooks

OUT=RESULTS/'ensemble_goamazon/processos'
GROUPS={
 'Condensação/evap. de nuvem':['th_prw_vcd','ws_pcond'],
 'Deposição/sublimação':['th_pri_ide','th_prs_ide','th_prs_sde','th_prg_gde','ws_pidep','ws_psdep','ws_pgdep'],
 'Nucleação de gelo':['th_pri_inu','th_pri_iha','ws_pigen'],
 'Congelamento':['th_pri_wfz','th_pri_rfz','th_prg_rfz','th_homogeneous_freezing','ws_pihmf','ws_pfrzdtc','ws_pfrzdtr'],
 'Riming (água de nuvem)':['th_prs_scw','th_prg_scw','th_prg_gcw','ws_paacw_snow','ws_paacw_graupel'],
 'Coleta chuva–gelo':['th_collection_rci','th_collection_rcs','th_collection_rcg','ws_piacr','ws_praci','ws_psacr','ws_pracs','ws_pgacr'],
 'Coleta gelo–gelo':['th_prs_sci','ws_psaci','ws_pgaci','ws_pgacs'],
 'Autoconversão/agregação':['th_prr_wau','th_prs_iau','ws_praut','ws_psaut','ws_pgaut'],
 'Coleta quente':['th_prr_rcw','ws_pracw','ws_paacw_warm'],
 'Fusão':['th_prr_sml','th_prr_gml','th_instant_melting','ws_pimlt','ws_psmlt','ws_pgmlt','ws_pseml','ws_pgeml'],
 'Evaporação precipitação':['th_prv_rev','ws_prevp','ws_psevp','ws_pgevp'],
 'Gelo secundário':['th_prs_ihm','th_prg_ihm'],
 'Sedimentação':['sedimentation'],
 'Ajustes/limiares':['tendency_adjustments','positivity_clipping','final_clipping'],
 'Resíduo não atribuído':['unattributed_residual'],
}
CATEGORY={p:g for g,ps in GROUPS.items() for p in ps}
QS=['qv','qc','qr','qi','qs','qg']


def analyze(folder,cfg,eid):
    ds=load(folder/'profile_1.nc'); _,mass=initial_water(ds,cfg)
    raw=pd.read_csv(folder/'profile_processes_1.csv')
    assert set(raw.process).issubset(CATEGORY),set(raw.process)-set(CATEGORY)
    assert np.isfinite(raw[['time_s','k']+['d'+q for q in QS]].values).all()
    # Full time/height/process archive of applied increments, not unlimited potentials.
    names=sorted(CATEGORY); nt=len(ds.time)-1; nz=len(ds.lev)
    values=np.zeros((nt,nz,len(names),6))
    if len(raw):
        ti=np.rint(raw.time_s.values/cfg['dt_microp']).astype(int)-1
        ki=raw.k.values.astype(int)-1
        pi=raw.process.map({p:i for i,p in enumerate(names)}).values
        assert not raw.duplicated(['time_s','k','process']).any()
        assert np.all((ti>=0)&(ti<nt)) and np.all((ki>=0)&(ki<nz))
        values[ti,ki,pi,:]=raw[['d'+q for q in QS]].values
    residual=values[:,:,names.index('unattributed_residual'),:]
    maxres=float(np.max(np.abs(residual)))
    weighted_abs=float(np.sum(np.abs(residual)*mass[None,:,None]))
    if maxres>1e-12 or weighted_abs>1e-7:
        raise ValueError(f'Unexplained process residual {eid}: {maxres}, {weighted_abs} mm')
    # Verify species closure independently against stored post-step hydrometeors.
    for j,q in enumerate(QS[1:],1):
        np.testing.assert_allclose(values[:,:,:,j].sum(axis=2),np.diff(ds[q].values,axis=0),atol=1e-12,rtol=1e-10)
    rates=xr.Dataset({'delta_q':(('time','lev','process','species'),values,{'units':'kg kg-1','description':'Applied increments accumulated over internal substeps; divide by dt for tendency'}),
                      'reference_dry_mass':(('lev',),mass,{'units':'kg m-2'}),
                      'process_name':(('process',),np.asarray(names,dtype='S48')),
                      'species_name':(('species',),np.asarray(QS,dtype='S2'))},
                     coords={'time':ds.time.values[1:],'lev':ds.lev.values,'process':np.arange(len(names)),'species':np.arange(6)})
    rates.attrs.update(experiment=eid,dt_s=cfg['dt_microp'],input_sha256=sha(INPUT),
                       normalization='Fixed pre-perturbation dry-air cell mass; not instantaneous native density',
                       process_categories=json.dumps(CATEGORY,ensure_ascii=False))
    rates.to_netcdf(folder/'process_tendencies.nc',engine='netcdf4',encoding={'delta_q':{'zlib':True,'complevel':4}})
    weighted=values*mass[None,:,None,None]
    rows=[]; series=[]
    for j,name in enumerate(names):
        for k,q in enumerate(QS):
            data=weighted[:,:,j,k]
            rows.append(dict(experiment=eid,process=name,category=CATEGORY[name],species=q,
                             net_mm=float(data.sum()),positive_mm=float(np.maximum(data,0).sum()),negative_mm=float(np.minimum(data,0).sum())))
        for it,t in enumerate(ds.time.values[1:]):
            series.append(dict(time_s=t,process=name,category=CATEGORY[name],
                               **{q+'_mm':float(weighted[it,:,j,k].sum()) for k,q in enumerate(QS)}))
    pd.DataFrame(rows).to_csv(folder/'process_integrals.csv',index=False,float_format='%.17g')
    pd.DataFrame(series).to_csv(folder/'process_timeseries.csv',index=False,float_format='%.17g')
    return rows,dict(max_unattributed_kg_kg=maxres,absolute_unattributed_mm=weighted_abs,
                    original_bit_identical=True,profile_sha256=sha(folder/'profile_1.nc'),diagnostics_sha256=sha(folder/'process_tendencies.nc'))


def main():
    baseline=json.loads((RESULTS/'ensemble_goamazon/manifest.json').read_text(encoding='utf-8'))
    assert sha(INPUT)==baseline['identity']['input_sha256']
    for relative in ['src/module_mp_thompson.F','src/wsm6/mp_wsm6.F90']:
        assert strip_hooks((ROOT/relative).read_text(encoding='utf-8'))==(PROV/Path(relative).name).read_text(encoding='utf-8')
        source_hashes={p.replace('\\','/'):v for p,v in baseline['identity']['sources'].items()}
        assert sha(PROV/Path(relative).name)==source_hashes[relative]
    parser=argparse.ArgumentParser(); parser.add_argument('--resume'); args=parser.parse_args()
    stamp=args.resume or datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    root=OUT/stamp
    assert root.resolve().parent==OUT.resolve()
    root.mkdir(parents=True,exist_ok=bool(args.resume))
    build=root/'build'
    try: result=linux_command(['make','-j4','BUILD='+build.relative_to(ROOT).as_posix()],ROOT)
    except subprocess.CalledProcessError as e:
        (root/'build.log').write_text((e.stdout or '')+(e.stderr or ''),encoding='utf-8'); raise
    (root/'build.log').write_text(result.stdout+result.stderr,encoding='utf-8')
    binary=build/'ensemble_driver'; rows=[]
    manifest=dict(directory=str(root.relative_to(ROOT)),created_utc=stamp,binary_sha256=sha(binary),input_sha256=sha(INPUT),
                  baseline_manifest_sha256=sha(RESULTS/'ensemble_goamazon/manifest.json'),runs={},
                  sources={p.relative_to(ROOT).as_posix():sha(p) for p in (ROOT/'src').rglob('*') if p.is_file()},
                  design='Exact Figure 6 initializations and forcing; diagnostics only; no retuning')
    for eid,base in baseline['members'].items():
        folder=root/eid; folder.mkdir(exist_ok=bool(args.resume))
        cfg=dict(base['config'],path_goamazon=linux_path(INPUT),write_process_diagnostics=True)
        write_nml(folder/'ensemble.nml',cfg,'ensemble')
        result=linux_command([linux_path(binary)],folder)
        (folder/'execution.log').write_text(result.stdout+result.stderr,encoding='utf-8')
        current=load(folder/'profile_1.nc'); prior=ROOT/base['directory']/'profile_1.nc'
        assert sha(prior)==base['output_sha256']
        old=load(prior)
        for key in old.variables: np.testing.assert_array_equal(current[key].values,old[key].values,err_msg=eid+' '+key)
        memberrows,validation=analyze(folder,cfg,eid); rows.extend(memberrows)
        manifest['runs'][eid]=dict(directory=str(folder.relative_to(ROOT)),config=cfg,validation=validation)
        save_json(root/'manifest.json',manifest)
        pd.DataFrame(rows).to_csv(root/'process_integrals.csv',index=False,float_format='%.17g')
        print(eid,'bitwise identical; closure',validation['absolute_unattributed_mm'],'mm',flush=True)
    save_json(OUT/'latest.json',manifest)
    print(root,flush=True)


if __name__=='__main__': main()
