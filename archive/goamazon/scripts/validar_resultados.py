"""Validate without silently correcting suspect fields. Missing != zero."""
import argparse
import json
import numpy as np
from scm_artigo import *

UNITS = {'time':'seconds since t0','lev':'m','theta':'K','pressure':'Pa',
         **{k:'kg kg-1' for k in ['qv']+SPECIES},'rainnc':'mm','rainncv':'mm',
         **{k:'flashes min-1' for k in ['flash_rate','flash_rate_ic','flash_rate_cg']},
         'snownc':'mm (liquid-water equivalent)','snowncv':'mm (liquid-water equivalent)',
         'pblh':'m','tke':'m2 s-2','mixing_length':'m','exch_h':'m2 s-1','exch_m':'m2 s-1',
         'cldfra_bl':'0-1','qc_bl':'kg kg-1','qi_bl':'kg kg-1','mf_a':'0-1','mf_w':'m s-1'}
PROFILES=['theta','pressure','qv']+SPECIES+['tke','mixing_length','exch_h','exch_m','cldfra_bl','qc_bl','qi_bl','mf_a','mf_w']

def validate(ds,cfg,eid=None,reference=None,check_initial=True):
    errors = config_errors(cfg,eid,reference) if eid else []
    warnings=[]; details={}
    for name,unit in UNITS.items():
        if name not in ds:
            errors.append(f'Variável ausente: {name}'); continue
        if ds[name].attrs.get('units')!=unit:
            errors.append(f'Unidade incompatível: {name}')
        dims = ('lev',) if name=='lev' else ('time','lev') if name in PROFILES else ('time',)
        if ds[name].dims!=dims: errors.append(f'Dimensões incompatíveis: {name}')
    for name,v in ds.variables.items():
        if np.issubdtype(v.dtype,np.number) and not np.isfinite(v.values).all():
            errors.append(f'Não finitos/ausentes: {name}')
    if errors: return {'valid':False,'errors':errors,'warnings':warnings,'details':details}
    time=ds.time.values; z=ds.lev.values; dt=np.diff(time)
    if len(time)<2 or np.any(dt<=0): errors.append('Tempo não crescente ou insuficiente')
    elif not np.allclose(dt,cfg['dt_microp'],rtol=0,atol=1e-8):
        warnings.append('Tempo irregular/diferente de DT_MICROP; métricas usam diff(time)')
    if len(z)!=int(cfg['nz']) or np.any(np.diff(z)<=0): errors.append('Grade vertical incompatível')
    try:
        dz=np.diff(interfaces(z,cfg['ztop']))
        if not np.isclose(dz[0],cfg['dz0'],atol=1e-8): errors.append('DZ0 incompatível')
        details['dz_min_max_m']=[float(dz.min()),float(dz.max())]
    except ValueError as e: errors.append(str(e))
    for name in ['qv']+SPECIES:
        if (ds[name].values<0).any(): errors.append(f'Razão de mistura negativa: {name}')
    if np.any(ds.pressure.values<=0) or np.any(np.diff(ds.pressure.values,axis=1)>=0): errors.append('Pressão inválida/não decrescente')
    if np.any(np.diff(ds.rainnc.values)<-1e-10): errors.append('Precipitação não monotônica')
    if np.any(ds.rainncv.values<0): errors.append('Incremento de precipitação negativo')
    if not np.allclose(np.diff(ds.rainnc.values),ds.rainncv.values[1:],rtol=1e-9,atol=1e-10): errors.append('rainnc/rainncv inconsistentes')
    if not np.allclose(ds.flash_rate,ds.flash_rate_ic+ds.flash_rate_cg,rtol=1e-12,atol=1e-14): errors.append('Partição elétrica inconsistente')
    if np.any(ds.flash_rate.values<0): errors.append('Taxa elétrica negativa')
    T=temperature(ds); tc=T-273.15
    # Bolton saturation over liquid, diagnostic approximation; exact driver Rd,cp,Rv.
    es=611.2*np.exp(17.67*tc/(tc+243.5))
    e=ds.qv.values*ds.pressure.values/(RD/RV+ds.qv.values)
    rh=100*e/es
    for label,sl in [('initial',slice(0,1)),('evolution',slice(1,None))]:
        mask=T[sl]>273.15; vals=rh[sl][mask]
        details[f'rh_warm_{label}_max_percent']=float(vals.max()) if vals.size else None
        details[f'rh_warm_{label}_above101_count']=int((vals>101).sum())
        if vals.size and vals.max()>101: warnings.append(f'RH quente {label} >101% (diagnóstico Bolton)')
    details['temperature_min_max_K']=[float(T.min()),float(T.max())]
    details['temperature_above16500m_min_K']=float(T[:,z>=16500].min()) if (z>=16500).any() else None
    if check_initial: warnings.append('Camada superior artificial; extremos frios não equivalem a validação atmosférica')
    details['samples']=len(time); details['last_time_s']=float(time[-1])
    if check_initial:
        forcing=load(INPUT)
        mapping={'theta':'theta','pressure':'pa','qv':'rv','qc':'rl','qi':'ri'}
        for name,src in mapping.items():
            ini=np.interp(z,forcing.lev.values,forcing[src].values.squeeze())
            if not np.allclose(ds[name].values[0],ini,rtol=1e-12,atol=1e-12): errors.append(f'Estado inicial incompatível com entrada: {name}')
        for name in ['qr','qs','qg']:
            if np.any(ds[name].values[0]!=0): errors.append(f'Inicial não nulo: {name}')
        end=float(forcing.time.values[-1]); expected_end=float(forcing.time.values[0])+np.ceil((end-float(forcing.time.values[0]))/cfg['dt_microp'])*cfg['dt_microp']
        if not np.isclose(time[-1],expected_end,atol=1e-8): errors.append('Integração incompleta/final incompatível')
    return {'valid':not errors,'errors':errors,'warnings':warnings,'details':details}

def validate_file(path,cfg,eid=None,reference=None,check_initial=True):
    try: return validate(load(path),cfg,eid,reference,check_initial)
    except Exception as e: return {'valid':False,'errors':[f'{type(e).__name__}: {e}'],'warnings':[],'details':{}}

def main():
    p=argparse.ArgumentParser(); p.add_argument('--results',type=Path,default=RESULTS); args=p.parse_args()
    reference=parse_nml(nml_path('TH_CF',args.results)); report={}; initial=None; grid=None
    for eid in IDS:
        path=args.results/f'{eid}.nc'
        if not path.exists(): report[eid]={'status':'ausente','valid':False}; continue
        cfg=parse_nml(nml_path(eid,args.results)); r=validate_file(path,cfg,eid,reference)
        if r['valid']:
            ds=load(path); state={k:ds[k].values[0] for k in ['theta','pressure','qv']+SPECIES}
            if grid is None: grid=(ds.time.values,ds.lev.values)
            elif not (np.array_equal(ds.time.values,grid[0]) and np.array_equal(ds.lev.values,grid[1])):
                r['valid']=False; r['errors'].append('Eixos temporal/vertical diferentes entre experimentos')
            if initial is None: initial=state
            elif any(not np.array_equal(state[k],initial[k]) for k in state):
                r['valid']=False; r['errors'].append('Estados iniciais não idênticos entre experimentos')
        r['status']='aprovado_com_avisos' if r['valid'] and r['warnings'] else 'aprovado' if r['valid'] else 'reprovado'
        r['sha256']=sha(path); report[eid]=r
    save_json(args.results/'validacao/qualidade.json',report)
    lines=['# Qualidade dos resultados','', 'Aprovação numérica não implica validação observacional. RH: fórmula de Bolton sobre água, com Rd/Rv do driver.','']
    for eid,r in report.items():
        lines += [f'## {eid}: {r["status"]}','',json.dumps(r,ensure_ascii=False,indent=2),'']
    (args.results/'validacao/RELATORIO_QUALIDADE.md').write_text('\n'.join(lines),encoding='utf-8')
    print({k:v['status'] for k,v in report.items()})
    if any(not r['valid'] and r['status']!='ausente' for r in report.values()): raise SystemExit(1)

if __name__=='__main__': main()
