"""Prepare 12 Amazon-only members; run a diagnostic smoke test before acceptance.

Never changes Fortran or input. Rejects the extension if ERA5 specific humidity
is passed unchanged to microphysics expecting dry-air mixing ratio.
"""
from datetime import datetime, timezone
import argparse
import subprocess
from scm_artigo import *
from executar_experimentos import linux_command,linux_path,write_nml
from validar_resultados import validate_file

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--rerun-smoke',action='store_true'); args=ap.parse_args()
    dest=RESULTS/'ensemble_amazonia'; dest.mkdir(exist_ok=True)
    if (dest/'status.json').exists() and not args.rerun_smoke:
        print((dest/'status.json').read_text(encoding='utf-8')); return
    if (dest/'status.json').exists():
        import shutil
        shutil.copy2(dest/'status.json',dest/('status_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')+'.json'))
    pl=ROOT/'era5_forcing/era5_pl.nc'; sl=ROOT/'era5_forcing/era5_sl.nc'
    if not pl.exists() or not sl.exists():
        save_json(dest/'status.json',{'status':'pendente','reason':'Dados ERA5 locais ausentes'}); return
    p=load(pl); s=load(sl); lat,lon=-8.,-63.
    if not (float(p.latitude.min())<=lat<=float(p.latitude.max()) and float(p.longitude.min())<=lon<=float(p.longitude.max())):
        raise ValueError('Ponto fora da cobertura: não usar vizinho distante')
    if not np.array_equal(p.valid_time.values,s.valid_time.values): raise ValueError('Eixos de tempo PL/SL distintos')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    base=dict(npoints=1,nsatconfigs=1,point_label='amazonia_8S_63W',path_pl=linux_path(pl),path_sl=linux_path(sl),
              target_lat=lat,target_lon=lon,itime_pick=16,nz=60,ztop=20000.,dz0=50.,dt_microp=20.,
              run_duration_s=7200.,pulse_sigma_z=1000.,pulse_duration_s=1800.,sat_rh_target=1.,
              write_profiles=True,path_summary='summary.nc',path_profile_prefix='profile')
    for scheme in ['thompson','wsm6']:
        for label,lo,hi in [('baixa',500.,1500.),('media',3500.,5500.),('alta',7500.,9500.)]:
            for pulse in [0.,3.]:
                cfg=dict(base,microp_scheme=scheme,satconfig_label=label,sat_z_lo=lo,sat_z_hi=hi,pulse_w_max=pulse)
                file=dest/f'{scheme}_{label}_w{int(pulse)}.nml'
                if not file.exists(): write_nml(file,cfg,'ensemble')
    # A single smoke member establishes whether the unmodified path is acceptable.
    smoke=dest/f'smoke_{stamp}'; smoke.mkdir()
    cfg=dict(base,microp_scheme='thompson',satconfig_label='baixa',sat_z_lo=500.,sat_z_hi=1500.,pulse_w_max=0.)
    write_nml(smoke/'ensemble.nml',cfg,'ensemble')
    try:
        r=linux_command([linux_path(ROOT/'build/ensemble_driver')],smoke)
        (smoke/'execution.log').write_text(r.stdout+r.stderr,encoding='utf-8')
    except subprocess.CalledProcessError as e:
        (smoke/'execution.log').write_text((e.stdout or '')+(e.stderr or ''),encoding='utf-8')
        save_json(dest/'status.json',{'status':'pendente','reason':'Falha do driver sem alteração Fortran','run':str(smoke),'exit_code':e.returncode}); return
    path=smoke/'profile_1.nc'
    qc=validate_file(path,cfg,check_initial=False)
    ds=load(path)
    col=p.sel(latitude=lat,longitude=lon).isel(valid_time=15)
    # At the lowest level, outside the intentionally saturated 500-1500 m layer.
    z=col.z.values/9.80616; q=col.q.values
    order=np.argsort(z)
    q_interp=float(np.interp(float(ds.lev.values[0]),z[order],q[order]))
    r_interp=float(np.interp(float(ds.lev.values[0]),z[order],(q/(1-q))[order]))
    output=float(ds.qv.values[0,0])
    semantic_problem=bool(np.isclose(output,q_interp,rtol=1e-10,atol=1e-12) and not np.isclose(output,r_interp,rtol=1e-6,atol=1e-10))
    status={'status':'pendente_validacao_semantica' if semantic_problem or not qc['valid'] else 'smoke_aprovado',
            'smoke_path':str(smoke.relative_to(ROOT)),'quality':qc,'specific_humidity_input_kg_kg':q_interp,
            'mixing_ratio_required_kg_kg':r_interp,'qv_output_initial_kg_kg':output,
            'relative_deficit_percent':100*(output/r_interp-1),'humidity_semantics_failed':semantic_problem,
            'point':[lat,lon],'time_index_1based':16,'time_epoch_seconds':float(col.valid_time),
            'pl_sha256':sha(pl),'sl_sha256':sha(sl),'driver_sha256':sha(ROOT/'build/ensemble_driver'),
            'source_evidence':['src/era5_bridge.F90: col%qv=max(buf,1e-12) after reading q',
                               'src/ensemble_driver.F90: regrid qvcol=col%qv without conversion'],
            'decision':'Não publicar Tabela4/Figura6: umidade específica é entregue como razão de mistura; correção exige revisão do caminho Fortran, fora do escopo. 12 namelists preparados, não executados como resultados científicos.'}
    save_json(dest/'status.json',status)
    print(status)

if __name__=='__main__': main()
