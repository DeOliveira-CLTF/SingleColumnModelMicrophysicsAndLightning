"""Local-only provenance; manifests never contain credentials."""
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def git(*args, cwd=ROOT):
    return subprocess.check_output(['git', *args], cwd=cwd, text=True, encoding='utf-8').strip()

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')

def verified_copy(source, destination):
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    before = sha256(source)
    if destination.exists() and sha256(destination) != before:
        raise FileExistsError(f'Refusing to replace different cached file: {destination}')
    if not destination.exists():
        shutil.copy2(source, destination)
    after = sha256(destination)
    if before != after or before != sha256(source):
        raise IOError(f'Checksum mismatch: {source}')
    return {'source': str(source), 'cache': str(destination), 'sha256': after,
            'bytes': destination.stat().st_size, 'redistribution': 'restricted'}

def audit(destination):
    base = '5c1655c6d23514d9a993b61b208dc9454f0e36ef'
    files = []
    for name in git('ls-files').splitlines():
        p = ROOT / name
        if p.is_file():
            files.append({'path': name, 'bytes': p.stat().st_size, 'sha256': sha256(p)})
    write_json(destination, {'utc': datetime.now(timezone.utc).isoformat(), 'baseline': base,
        'head': git('rev-parse', 'HEAD'), 'status': git('status', '--short'),
        'baseline_diff': git('diff', '--stat', base), 'files': files,
        'python': platform.python_version(),
        'tools': {x: shutil.which(x) for x in ['gfortran', 'nf-config', 'make', 'wsl']},
        'dephy_head': git('rev-parse', 'HEAD', cwd=ROOT/'external/DEPHY-SCM'),
        'dephy_status': git('status', '--short', cwd=ROOT/'external/DEPHY-SCM')})


def resolve_path(value):
    """Resolve immutable pre-reorganization manifest paths into the current layout."""
    value = str(value).replace('\\', '/')
    root = ROOT.as_posix()
    if value.startswith(root + '/'):
        value = value[len(root)+1:]
    elif value.startswith('/mnt/c/git_repos/SingleColumnModelMicrophysicsAndLightning/'):
        value = value.split('SingleColumnModelMicrophysicsAndLightning/', 1)[1]
    for short, long in [('TH', 'thompson'), ('WS', 'wsm6')]:
        for mode in ['C0', 'CL', 'CF']:
            cfg = short+'_'+mode
            dest = 'runs/argentina_20181213/'+long+'_'+mode
            for old in ['work/leona/ensemble/', 'work/leona/ensemble40/']:
                prefix = old+cfg+'/'
                if value.startswith(prefix):
                    return ROOT/dest/'members'/value[len(prefix):]
            prefix = 'work/leona/controls/'+cfg
            if value == prefix or value.startswith(prefix+'/'):
                return ROOT/(dest+'/control'+value[len(prefix):])
    mapping = {'products/leona': 'results/argentina_20181213', 'data_private/leona/era5': 'data/raw/era5', 'data_private/leona/raw/entln': 'data/raw/entln', 'data_private/leona/raw/soundings': 'data/raw/radiosondes', 'data_private/leona/processed/soundings': 'data/processed/profiles', 'data_private/leona/wyoming': 'data/raw/radiosondes/wyoming', 'data_private/leona/forcing_source': 'data/intermediate/era5', 'data_private/leona': 'data/intermediate/input_cache', 'work/leona/case_forcing_v2': 'data/processed/forcing', 'work/leona/build': 'build/leona', 'work/leona': 'runs/argentina_20181213/_execution', 'era5_forcing/DEPHY-SCM': 'external/DEPHY-SCM', 'configs/leona.json': 'configs/cases/argentina_20181213.yaml', 'configs/parameterizations.json': 'configs/lightning/schemes.yaml', 'configs/ensemble_stability.json': 'configs/ensemble/stability.yaml', 'configs/budget_contract.json': 'configs/evaluation/budget_contract.yaml', 'configs/sounding_sources.json': 'configs/cases/sounding_sources.yaml'}
    for old, new in sorted(mapping.items(), key=lambda item: -len(item[0])):
        if value == old or value.startswith(old+'/'):
            return ROOT/(new+value[len(old):])
    return Path(value) if Path(value).is_absolute() else ROOT/value
