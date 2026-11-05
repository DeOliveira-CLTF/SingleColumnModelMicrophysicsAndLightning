"""Regenerate scientific tables from available, validated outputs."""
import argparse
import pandas as pd
from scm_artigo import *
from validar_resultados import validate_file

CONTRASTS=[('microfisica',f'WS_{c}',f'TH_{c}') for c in ['C0','CL','CF']]+[(kind,f'{s}_{a}',f'{s}_{b}') for kind,a,b in [('local','CL','C0'),('plumas_adicional','CF','CL'),('total','CF','C0')] for s in ['TH','WS']]
COMPARABLE=['rain_final_mm','rain_35h_mm','rain_rate_max_mm_h','flash_max_min','flash_integral_nominal','flash_active_duration_h']+[f'{q}_max_g_kg' for q in SPECIES]+[f'{q}_mixed_max_g_kg' for q in SPECIES]+[f'{q}_mixed_mass_max_kg_m2' for q in SPECIES]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--results',type=Path,default=RESULTS); args=ap.parse_args(); directory=args.results
    tab=directory/'tabelas'; tab.mkdir(parents=True,exist_ok=True)
    (directory/'series').mkdir(exist_ok=True)
    ref=parse_nml(nml_path('TH_CF',directory)); configs=[]; rows=[]; values={}
    for eid in IDS:
        path=directory/f'{eid}.nc'; cfg=dict(ref,**expected(eid))
        if path.exists(): cfg=parse_nml(nml_path(eid,directory))
        configs.append(dict(experiment=eid,available=path.exists(),**cfg))
        if not path.exists():
            rows.append(dict(experiment=eid,status='ausente')); continue
        quality=validate_file(path,cfg,eid,ref)
        if not quality['valid']: raise ValueError(f'{eid}: {quality["errors"]}')
        ds=load(path); m,series=metrics(ds,cfg); values[eid]=m
        rows.append(dict(experiment=eid,microp_scheme=cfg['microp_scheme'],coupling=eid.split('_')[1],
                         status='aprovado_com_avisos',**m,quality_warnings='; '.join(quality['warnings']),
                         output_sha256=sha(path)))
        pd.DataFrame(series).to_csv(directory/'series'/f'{eid}.csv',index=False,float_format='%.17g')
    pd.DataFrame(configs).to_csv(tab/'Tabela1_configuracoes.csv',index=False,float_format='%.17g')
    pd.DataFrame(rows).to_csv(tab/'Tabela2_resultados.csv',index=False,float_format='%.17g')
    contrasts=[]
    for kind,a,b in CONTRASTS:
        for metric in COMPARABLE:
            av=values.get(a,{}).get(metric); bv=values.get(b,{}).get(metric)
            ok=av is not None and bv is not None
            contrasts.append(dict(type=kind,contrast=f'{a} - {b}',metric=metric,a=av,b=bv,
                                  difference=av-bv if ok else None,
                                  relative_percent=100*(av-bv)/bv if ok and bv!=0 else None,
                                  status='disponivel' if ok else 'ausente',
                                  relative_status='definida' if ok and bv!=0 else 'indefinida'))
    for c in ['CL','CF']:
        for metric in COMPARABLE:
            ids=[f'WS_{c}','WS_C0',f'TH_{c}','TH_C0']; v=[values.get(i,{}).get(metric) for i in ids]
            ok=all(x is not None for x in v)
            contrasts.append(dict(type='interacao',contrast=f'(WS_{c}-WS_C0)-(TH_{c}-TH_C0)',metric=metric,
                                  difference=v[0]-v[1]-v[2]+v[3] if ok else None,status='disponivel' if ok else 'ausente',relative_status='nao_aplicavel'))
    pd.DataFrame(contrasts).to_csv(tab/'Tabela3_contrastes.csv',index=False,float_format='%.17g')
    save_json(tab/'metadados.json',{
        'time':'hours since t0, not reconciled civil time; earliest occurrence for tied maxima',
        'rain_rate':'rainncv[1:] * 3600 / diff(time); first rate absent',
        'electrical_integral':'sum(flash_rate[1:] * diff(time)/60), nominal flashes, uncalibrated',
        'activity':'strictly >0; duration sums positive-step intervals; censored if last rate >0',
        'mixed_phase':'233.15 <= T <= 273.15 K; entire layer selected by center temperature',
        'temperature':'theta*(pressure/100000)**(287/1004.5)',
        'density':'pressure/(287*T*(1+0.61*qv)), same virtual-temperature approximation as driver; recomputed at saved state',
        'interfaces':'zw[0]=0; zw[k+1]=2*lev[k]-zw[k]; exact inversion of driver centers; top checked to 1e-7 m',
        'mass':'sum(rho*q*dz) in kg/m2; center-cell quadrature at phase boundaries, not an ice flux',
        'correlation':'Pearson over full common window only for nonconstant series; no p-values; internal EP dependence, not independent validation',
        'lag':'difference of maxima times, not cross-correlation optimization',
        'missing':'empty CSV fields, never zero; absent experiments explicitly status=ausente',
        'constants_source':'src/mpas_constants.F, src/mp_offline_driver.F90',
        'source_attribute':'not used for scheme identification; legacy output attribute always mentions Thompson'})
    print(pd.DataFrame(rows)[['experiment','status','rain_final_mm','flash_max_min','flash_integral_nominal']].to_string(index=False))

if __name__=='__main__': main()
