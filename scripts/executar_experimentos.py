"""Run only missing valid cases. Windows Python orchestrates the existing WSL binary."""
import argparse
from datetime import datetime, timezone
import os
import platform
import shutil
import subprocess
import sys
from scm_artigo import *
from validar_resultados import validate_file

def linux_path(path):
    path=Path(path).resolve()
    if os.name!='nt': return str(path)
    return '/mnt/'+path.drive[0].lower()+'/'+path.as_posix()[3:]

def linux_command(argv,cwd=None,capture=True):
    if os.name=='nt':
        cmd=['wsl','-d','Ubuntu-24.04']
        if cwd: cmd+=['--cd',linux_path(cwd)]
        cmd+=['--']+argv
    else: cmd=argv
    return subprocess.run(cmd,cwd=cwd if os.name!='nt' else None,check=True,
                          capture_output=capture,text=True,encoding='utf-8',errors='replace')

def write_nml(path,cfg,group='mp_offline'):
    def value(v):
        if isinstance(v,bool): return '.true.' if v else '.false.'
        if isinstance(v,str): return "'"+v.replace("'","''")+"'"
        return str(v)
    path.write_text('&'+group+'\n'+''.join(f' {k} = {value(v)},\n' for k,v in cfg.items())+'/\n',encoding='ascii')

def provenance():
    def git(*args): return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()
    return {'recorded_utc':datetime.now(timezone.utc).isoformat(),'git_commit':git('rev-parse','HEAD'),
            'git_status':git('status','--short'),'git_diff_makefile':git('diff','--','Makefile'),
            'input':str(INPUT),'input_sha256':sha(INPUT),
            'physics_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'src').rglob('*')) if p.suffix.lower() in ['.f','.f90']},
            'executable_sha256':sha(ROOT/'build/mp_offline_driver'),
            'python':sys.version,'platform':platform.platform(),
            'compiler':linux_command(['gfortran','--version']).stdout,
            'netcdf':linux_command(['nf-config','--all']).stdout,
            'build_log':'resultados/build_TH_CF.log (clean successful build from earlier session)',
            'legacy_provenance_limit':'Original CF runs have namelists and data, but no contemporaneous hash manifest; hashes recorded now are not retrospective proof of their execution environment.'}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--experiments',nargs='+',choices=IDS,default=['TH_CF','WS_CF','TH_C0','TH_CL','WS_C0','WS_CL'])
    ap.add_argument('--rerun',action='store_true'); ap.add_argument('--rebuild',action='store_true'); args=ap.parse_args()
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    for folder in ['logs','namelists','validacao','proveniencia','runs','versions']: (RESULTS/folder).mkdir(parents=True,exist_ok=True)
    reference=parse_nml(nml_path('TH_CF'))
    pending=args.rerun or any(not (RESULTS/f'{e}.nc').exists() for e in args.experiments)
    binary=ROOT/'build/mp_offline_driver'
    # Timestamps cannot establish a clean build. Trust only a matching recorded
    # executable and complete source hash set; otherwise rebuild before new runs.
    trusted=False
    if binary.exists():
        current={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'src').rglob('*')) if p.suffix.lower() in ['.f','.f90']}
        make_diff=subprocess.check_output(['git','diff','--','Makefile'],cwd=ROOT,text=True).strip()
        for old in (RESULTS/'proveniencia').glob('execution_*.json'):
            record=json.loads(old.read_text(encoding='utf-8'))
            established=record.get('build_verified',False) or old.name=='execution_20261102T214208919408Z.json'
            if established and record.get('physics_sha256')==current and record.get('executable_sha256')==sha(binary) and record.get('git_diff_makefile')==make_diff:
                trusted=True
    new_build_log=None
    if args.rebuild or not binary.exists() or (pending and not trusted):
        # Preserve ALL existing build contents before the Makefile's rm -rf build.
        backup=RESULTS/'versions'/('build_'+stamp)
        if (ROOT/'build').exists(): shutil.copytree(ROOT/'build',backup)
        for command in [['make','clean'],['make','-j4']]:
            new_build_log=RESULTS/'logs'/f'build_{stamp}.log'
            try: r=linux_command(command,ROOT)
            except subprocess.CalledProcessError as e:
                with new_build_log.open('a',encoding='utf-8') as f: f.write((e.stdout or '')+(e.stderr or ''))
                raise
            with new_build_log.open('a',encoding='utf-8') as f: f.write(r.stdout+r.stderr)
        trusted=True
    manifest=provenance(); manifest['experiments']={}
    manifest['build_verified']=trusted
    if new_build_log: manifest['build_log']=str(new_build_log.relative_to(ROOT))
    manifest_path=RESULTS/'proveniencia'/f'execution_{stamp}.json'
    save_json(manifest_path,manifest)
    for eid in args.experiments:
        target=RESULTS/f'{eid}.nc'
        if target.exists() and not args.rerun:
            nml=nml_path(eid); cfg=parse_nml(nml); r=validate_file(target,cfg,eid,reference)
            if not r['valid']: raise RuntimeError(f'{eid}: existing result invalid; preserved. {r["errors"]}')
            dest=RESULTS/'namelists'/f'{eid}.nml'
            if dest!=nml and not dest.exists(): shutil.copy2(nml,dest)
            print(eid,'reutilizado',flush=True)
            manifest['experiments'][eid]={'action':'reused','output_sha256':sha(target),'namelist':cfg,'namelist_sha256':sha(nml),'quality':r}
            save_json(manifest_path,manifest)
            continue
        run=RESULTS/'runs'/f'{eid}_{stamp}'; run.mkdir()
        cfg=dict(reference,**expected(eid)); cfg['path_input']=linux_path(INPUT); cfg['path_output']=f'{eid}.nc'
        nml=run/'namelist.mp_offline'; write_nml(nml,cfg)
        log=RESULTS/'logs'/f'{eid}_{stamp}.log'
        try:
            result=linux_command([linux_path(ROOT/'build/mp_offline_driver')],run)
            log.write_text(result.stdout+result.stderr,encoding='utf-8')
        except subprocess.CalledProcessError as e:
            log.write_text((e.stdout or '')+(e.stderr or ''),encoding='utf-8'); raise
        produced=run/f'{eid}.nc'; quality=validate_file(produced,cfg,eid,reference)
        save_json(run/'validation.json',quality)
        if not quality['valid']: raise RuntimeError(f'{eid}: {quality["errors"]}; raw output retained in {run}')
        if target.exists():
            old=RESULTS/'versions'/f'{eid}_{stamp}'; old.mkdir()
            shutil.copy2(target,old/target.name)
            shutil.copy2(nml_path(eid),old/f'{eid}.nml')
        shutil.copy2(produced,target); shutil.copy2(nml,RESULTS/'namelists'/f'{eid}.nml')
        record={'action':'executed','exit_code':result.returncode,'output_sha256':sha(target),'namelist':cfg,'namelist_sha256':sha(nml),'log':str(log.relative_to(ROOT)),'quality':quality}
        manifest['experiments'][eid]=record
        save_json(RESULTS/'proveniencia'/f'{eid}_{stamp}.json',dict(manifest,experiments={eid:record}))
        print(eid,'executado e aprovado',flush=True)
        save_json(manifest_path,manifest)
    save_json(manifest_path,manifest)

if __name__=='__main__': main()
