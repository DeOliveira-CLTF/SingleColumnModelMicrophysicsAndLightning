"""Twelve idealized members using GOAMAZON initial profiles and original perturbations.

No input or parameterization is edited. Each run has an isolated directory.
Existing outputs are reusable only when configuration, source and input hashes match.
"""
import argparse
from datetime import datetime, timezone
import subprocess
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scm_artigo import *
from executar_experimentos import linux_command, linux_path, write_nml
from validar_resultados import validate_file

DEST=RESULTS/'ensemble_goamazon'
LAYERS=[('baixa',500.,1500.),('media',3500.,5500.),('alta',7500.,9500.)]


def saturation(p,t):
    """Independent evaluation of the unchanged RSLF/RSIF polynomials."""
    liquid=[611.583699,44.4606896,1.43177157,.0264224321,.000299291081,.00000203154182,7.02620698e-9,3.79534310e-12,-3.21582393e-14]
    ice=[609.868993,49.9320233,1.84672631,.0402737184,.000565392987,.00000521693933,3.07839583e-8,1.05785160e-10,1.61444444e-13]
    x=np.maximum(-80.,t-273.16)
    es=np.where(t>253.15,np.polynomial.polynomial.polyval(x,liquid),np.polynomial.polynomial.polyval(x,ice))
    es=np.minimum(es,p*.15)
    return .622*es/(p-es)


def check(path,cfg):
    qa=validate_file(path,cfg,check_initial=False)
    if not qa['valid']: raise ValueError(qa)
    ds=load(path); source=load(INPUT); z=ds.lev.values
    np.testing.assert_allclose(ds.time.values,np.arange(0,cfg['run_duration_s']+cfg['dt_microp']/2,cfg['dt_microp']),rtol=0,atol=1e-9)
    original={k:np.interp(z,source.lev.values,source[v].values.squeeze()) for k,v in
              [('theta','theta'),('pressure','pa'),('qv','rv'),('qc','rl'),('qi','ri')]}
    for k in ['theta','pressure']:
        np.testing.assert_allclose(ds[k][0],original[k],rtol=1e-12,atol=1e-12)
    mask=(z>=cfg['sat_z_lo'])&(z<=cfg['sat_z_hi'])
    t=original['theta']*(original['pressure']/100000)**(287/1004.5)
    target=saturation(original['pressure'],t)
    excess=np.maximum(original['qv']-target,0.)
    expected_qv=original['qv'].copy(); expected_qv[mask]=target[mask]
    expected_qc=original['qc']+np.where(mask&(t>253.15),excess,0)
    expected_qi=original['qi']+np.where(mask&(t<=253.15),excess,0)
    for k,a in [('qv',expected_qv),('qc',expected_qc),('qi',expected_qi)]:
        np.testing.assert_allclose(ds[k][0],a,rtol=2e-10,atol=1e-12)
    for k in ['qr','qs','qg']: np.testing.assert_array_equal(ds[k][0],np.zeros(len(z)))
    # Initial rain/flash are explicitly zero; pressure is the fixed initial profile.
    np.testing.assert_allclose(ds.pressure.values,np.broadcast_to(original['pressure'],ds.pressure.shape),rtol=1e-12)
    m,series=metrics(ds,cfg)
    summary=load(path.parent/'summary.nc')
    np.testing.assert_allclose(float(summary.flash_count_total.values[0]),m['flash_integral_nominal'],rtol=1e-11,atol=1e-12)
    np.testing.assert_allclose(float(summary.flash_rate_peak.values[0]),m['flash_max_min'],rtol=1e-12,atol=1e-12)
    for q in SPECIES:
        # Driver summary records post-step maxima, excluding the perturbed initial record.
        np.testing.assert_allclose(float(summary[q+'_peak'].values[0]),float(ds[q][1:].max()),rtol=1e-12,atol=1e-12)
    qa['warnings'].append('Camada superior artificial do GOAMAZON; ensemble idealizado, sem validação observacional')
    qa['initial_perturbation_checked']=True
    qa['input_sha256']=sha(INPUT)
    qa['actual_saturated_centers_m']=z[mask].tolist()
    qa['water_added_initial_kg_m2']=float(np.sum((expected_qv+expected_qc+expected_qi-original['qv']-original['qc']-original['qi'])*
                                     original['pressure']/(287*t*(1+.61*original['qv']))*np.diff(interfaces(z,cfg['ztop']))))
    return ds,m,series,qa


def products(rows, write_report=True):
    df=pd.DataFrame(rows)
    plt.rcParams.update({'font.size':20,'pdf.fonttype':42,'axes.titlesize':22,'axes.labelsize':20,'xtick.labelsize':18,'ytick.labelsize':18,'legend.fontsize':18})
    fig,axes=plt.subplots(2,2,figsize=(20,14),layout='constrained')
    keys=['rain_final_mm','qg_mixed_mass_max_kg_m2','flash_max_min','flash_integral_nominal']
    labels=['Accumulated precipitation (mm)','Maximum mixed-phase graupel mass (kg/m²)','Maximum electrical rate (flashes/min)','Nominal electrical integral (flashes)']
    cats=[(name,w) for name,_,_ in LAYERS for w in [0,3]]
    for ax,key,label in zip(axes.flat,keys,labels):
        for j,scheme in enumerate(['thompson','wsm6']):
            subset=df[df.microp_scheme==scheme].set_index(['layer','pulse_w_max'])
            vals=[float(subset.loc[c,key]) for c in cats]; x=np.arange(6)+(j-.5)*.36
            ax.bar(x,vals,width=.34,label=scheme.title() if j==0 else 'WSM6',color=['#0072B2','#D55E00'][j])
            for xx,v in zip(x,vals): ax.annotate(f'{v:.3g}',(xx,v),xytext=(0,8),textcoords='offset points',ha='center',fontsize=16)
        ax.set_xticks(range(6),[f'{c}\nw={w}' for c,w in cats]); ax.set_ylabel(label,fontsize=20)
        ax.set_ylim(0,max(float(df[key].max()),1e-10)*1.3); ax.legend(loc='upper left',ncol=2,fontsize=18); ax.grid(axis='y',alpha=.2)
    fig.suptitle('Idealized GOAMAZON ensemble · 2 h · without MYNN\nInitial layer saturation; additional pulse of 0 or 3 m/s',fontsize=28)
    out=ROOT/'figuras_artigo'; out.mkdir(exist_ok=True)
    for ext in ['png','pdf']:
        fig.savefig(out/f'Figura6_ensemble_goamazon.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)
    if not write_report:
        return
    (RESULTS/'tabelas').mkdir(exist_ok=True)
    df.to_csv(RESULTS/'tabelas/Tabela4_ensemble_goamazon.csv',index=False,float_format='%.17g')
    compact=df[['experiment','rain_final_mm','qg_max_g_kg','flash_max_min','flash_integral_nominal']]
    lines=['# Ensemble GOAMAZON — resultados','',
           'Doze membros aprovados nos controles implementados. Dados originais em `data_goamazon/`; nenhuma entrada ERA5 é lida. As rotinas de saturação, pulso, microfísica e EP permanecem inalteradas.',
           '', '## Desenho e interpretação','',
           'Perfil inicial (pa, theta, rv, rl, ri) e wa do primeiro registro do GOAMAZON. Pressão e wa de fundo fixas durante 2 h; theta e vapor evoluem pela advecção vertical e microfísica do ensemble original. Sem MYNN, tendências horizontais, atualização temporal das forçantes ou fluxos de superfície. Portanto, este é um experimento idealizado com o perfil GOAMAZON, distinto da matriz C0/CL/CF e da reprodução de 35,26 h.',
           '', 'Grade de 60 camadas, topo 20 km, primeira espessura 50 m, passo 20 s. Saturação aplicada uma vez nos centros entre 500–1500, 3500–5500 e 7500–9500 m, com RH alvo 1. A rotina original usa saturação sobre água acima de −20 °C e sobre gelo abaixo. Acrescenta vapor onde necessário, sem compensação energética: uma perturbação idealizada, não uma evolução fechada conservativa. Excesso preexistente vira condensado conforme a rotina original. As massas de água acrescentada e alturas efetivas estão nas validações individuais.',
           '', 'Pulso original: 0,5[1−cos(2πt/1800)] × exp[−0,5((z−zc)/1000)²] × amplitude, durante 30 min, centrado no meio de cada camada. Amplitude 0 ou 3 m/s. Pulso zero mantém wa inicial. A discretização vertical pode não amostrar exatamente o máximo de 3 m/s.',
           '', 'O transporte vertical no caminho GOAMAZON foi corrigido: upwind de primeira ordem com distâncias reais entre centros e subpassos CFL ≤0,8, mantendo a equação −w ∂φ/∂z. Contorno de entrada com gradiente nulo. O método centrado explícito anterior gerou temperaturas de 885 K antes de NaNs no teste de pulso baixo, em 1.260 s; essa tentativa foi preservada e excluída das tabelas. A mudança é numérica e introduz difusão de primeira ordem; não se afirma equivalência dos resultados com o transporte antigo nem convergência de resolução. As rotinas físicas e o pulso não foram alterados.',
           '', '## Síntese quantitativa','', '| Experimento | Chuva (mm) | qg máx. (g/kg) | Taxa elétrica máx. (min⁻¹) | Integral nominal |', '| --- | --- | --- | --- | --- |']
    for row in compact.itertuples(index=False): lines.append('| '+' | '.join([row[0]]+[f'{x:.6g}' for x in row[1:]])+' |')
    active=df[df.flash_max_min>0]
    lines+=['',f'{len(active)} de 12 membros apresentam EP positivo. Taxas e integrais são nominais, sem calibração. Horários de séries nulas permanecem ausentes no CSV. Indicador de censura identifica atividade no último passo.',
            '', 'Máximos globais da Tabela 4 incluem o estado inicial perturbado; o resumo nativo do driver considera somente os passos posteriores. A validação compara esse resumo com o mesmo subconjunto, evitando uma falsa discrepância.',
            '', '## Validação e proveniência','',
            'Verificados pressão estritamente decrescente, finitude, não negatividade, fechamento da chuva incremental, 361 registros em 2 h, estado inicial regradeado, saturação contra os polinômios RSLF/RSIF, identidade inicial entre esquemas/pulsos pareados e integral elétrica contra o resumo Fortran. Todos os membros usam exatamente o mesmo arquivo fonte. A camada superior artificial do GOAMAZON e a incerteza geográfica/temporal original permanecem limitações. Não há validação observacional.',
            '', 'Execuções e logs em subpastas individuais. `manifest.json` registra fontes, entrada, binário, compilação e hashes das saídas. O teste ERA5 antigo permanece apenas como histórico, fora destes resultados.',
            '', '**Figura 6.** Chuva final, máximo da massa de graupel entre −40 e 0 °C, taxa elétrica máxima e integral nominal nos 12 membros. Cores distinguem microfísicas; cada camada aparece com pulso adicional 0/3 m/s. Escalas lineares e zeros explícitos; sem inferência de habilidade observacional.',
            '', 'Reprodução: `python scripts/executar_ensemble_goamazon.py`. Reutiliza somente execuções com hashes compatíveis; `--rerun` cria novas pastas e preserva as anteriores.']
    if (DEST/'diagnostico_cientifico/RELATORIO_DIAGNOSTICO_CIENTIFICO.md').exists():
        lines[2:2]=['**Atualização científica:** o [diagnóstico de água, perfis e resolução](diagnostico_cientifico/RELATORIO_DIAGNOSTICO_CIENTIFICO.md) mostra aporte artificial relevante e convergência quantitativa ainda não demonstrada. A aprovação abaixo significa integridade numérica, não conservação ou convergência. A Figura 6 é exploratória.','']
    (DEST/'RELATORIO_ENSEMBLE_GOAMAZON.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(compact.to_string(index=False),flush=True)


def figure_only():
    manifest_path=DEST/'manifest.json'
    if not manifest_path.exists():
        raise FileNotFoundError(f'Manifesto existente não encontrado: {manifest_path}')
    manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    rows=[]
    for eid,record in manifest.get('members',{}).items():
        run=ROOT/record['directory']
        profile=run/'profile_1.nc'
        if not profile.exists():
            raise FileNotFoundError(f'Perfil existente não encontrado: {profile}')
        cfg=record['config']
        ds=load(profile)
        metrics_result,_=metrics(ds,cfg)
        scheme,layer,pulse=eid.split('_')
        rows.append(dict(experiment=eid,microp_scheme=scheme,layer=layer,
                         pulse_w_max=float(pulse.removeprefix('w')),
                         **metrics_result))
    if len(rows) != 12:
        raise ValueError(f'Esperados 12 membros existentes; encontrados {len(rows)}')
    products(rows,write_report=False)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--rerun',action='store_true')
    ap.add_argument('--figure-only',action='store_true',
                    help='gera somente a Figura 6 a partir dos resultados existentes')
    args=ap.parse_args()
    if args.figure_only:
        figure_only()
        return
    DEST.mkdir(exist_ok=True); stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    source_hashes={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'src').rglob('*') if p.suffix.lower() in ['.f','.f90']}
    identity={'input_sha256':sha(INPUT),'sources':source_hashes,'makefile_sha256':sha(ROOT/'Makefile')}
    previous=json.loads((DEST/'manifest.json').read_text(encoding='utf-8')) if (DEST/'manifest.json').exists() else {}
    reuse=previous.get('identity')==identity and not args.rerun
    if reuse:
        binary=ROOT/previous['binary']
        if not binary.exists() or sha(binary)!=previous['binary_sha256']: reuse=False
    if not reuse:
        build=DEST/f'build_{stamp}'; build.mkdir()
        log=DEST/f'build_{stamp}.log'
        try:
            r=linux_command(['make','-j4','BUILD='+str(build.relative_to(ROOT)).replace('\\','/')],ROOT)
            log.write_text(r.stdout+r.stderr,encoding='utf-8')
        except subprocess.CalledProcessError as e:
            log.write_text((e.stdout or '')+(e.stderr or ''),encoding='utf-8'); raise
        binary=build/'ensemble_driver'
    manifest={'identity':identity,'binary':str(binary.relative_to(ROOT)),'binary_sha256':sha(binary),
              'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'compiler':linux_command(['gfortran','--version']).stdout,'created_utc':stamp,'members':{}}
    rows=[]; initial_by_layer={}
    for scheme in ['thompson','wsm6']:
        for label,lo,hi in LAYERS:
            for pulse in [0.,3.]:
                eid=f'{scheme}_{label}_w{int(pulse)}'
                cfg=dict(input_source='goamazon',path_goamazon=linux_path(INPUT),npoints=1,nsatconfigs=1,
                         point_label='GOAMAZON_original',microp_scheme=scheme,satconfig_label=label,sat_z_lo=lo,sat_z_hi=hi,
                         sat_rh_target=1.,nz=60,ztop=20000.,dz0=50.,dt_microp=20.,run_duration_s=7200.,
                         pulse_sigma_z=1000.,pulse_duration_s=1800.,pulse_w_max=pulse,write_profiles=True,
                         path_summary='summary.nc',path_profile_prefix='profile')
                record=previous.get('members',{}).get(eid,{})
                cached=reuse and record.get('config')==cfg
                if cached:
                    run=ROOT/record['directory']; profile=run/'profile_1.nc'
                    cached=(profile.exists() and sha(profile)==record['output_sha256'] and
                            sha(run/'ensemble.nml')==record['namelist_sha256'] and sha(run/'summary.nc')==record['summary_sha256'])
                if not cached:
                    run=DEST/f'{eid}_{stamp}'; run.mkdir()
                    write_nml(run/'ensemble.nml',cfg,'ensemble')
                    try:
                        r=linux_command([linux_path(binary)],run)
                        (run/'execution.log').write_text(r.stdout+r.stderr,encoding='utf-8')
                    except subprocess.CalledProcessError as e:
                        (run/'execution.log').write_text((e.stdout or '')+(e.stderr or ''),encoding='utf-8'); raise
                    profile=run/'profile_1.nc'
                ds,m,series,qa=check(profile,cfg)
                state={k:ds[k].values[0] for k in ['theta','pressure','qv']+SPECIES}
                if label in initial_by_layer:
                    for k,v in state.items(): np.testing.assert_array_equal(v,initial_by_layer[label][k])
                else: initial_by_layer[label]=state
                save_json(run/'validation.json',qa)
                pd.DataFrame(series).to_csv(run/'series.csv',index=False,float_format='%.17g')
                record={'directory':str(run.relative_to(ROOT)),'config':cfg,'output_sha256':sha(profile),
                        'namelist_sha256':sha(run/'ensemble.nml'),'summary_sha256':sha(run/'summary.nc'),'exit_code':0,'quality':qa}
                manifest['members'][eid]=record
                save_json(DEST/f'manifest_{stamp}.json',manifest)
                rows.append(dict(experiment=eid,microp_scheme=scheme,layer=label,pulse_w_max=pulse,**m,
                                 status='aprovado_com_avisos',input_sha256=sha(INPUT),output_sha256=sha(profile)))
                print(eid,'reutilizado' if cached else 'executado','e aprovado',flush=True)
    save_json(DEST/'manifest.json',manifest)
    products(rows)


if __name__=='__main__': main()
