"""Local-only provenance; manifests never contain credentials."""
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]

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
        'dephy_head': git('rev-parse', 'HEAD', cwd=ROOT/'era5_forcing/DEPHY-SCM'),
        'dephy_status': git('status', '--short', cwd=ROOT/'era5_forcing/DEPHY-SCM')})
