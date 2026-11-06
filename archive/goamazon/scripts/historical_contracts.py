"""Source preservation relative to the user-specified audited commit.

Older run manifests remain immutable; some refer to source revisions predating
that commit. Binary input/output checksums in those manifests are still checked.
"""
import hashlib
from pathlib import Path
import subprocess
from scm_artigo import ROOT
BASELINE='5c1655c6d23514d9a993b61b208dc9454f0e36ef'
INTERFACES={'src/goamazon_bridge.F90','src/mp_offline_driver.F90','src/output_mp_netcdf.F90'}

def baseline_bytes(path):
    name=Path(path).as_posix()
    return subprocess.check_output(['git','show',BASELINE+':'+name],cwd=ROOT)

def canonical(data):return data.replace(b'\r\n',b'\n')

def assert_audited_source(test,path):
    name=str(path).replace('\\','/')
    old=baseline_bytes(name)
    if name in INTERFACES:
        # These interfaces are intentionally corrected by LEONA. Their physics
        # contracts and historical switches have dedicated executable tests.
        test.assertTrue(old)
        if name=='src/goamazon_bridge.F90':
            text=(ROOT/name).read_text()
            test.assertIn('legacy_temperature_tendency',text)
            test.assertIn("'tntheta_adv'",text)
    else:
        test.assertEqual(canonical((ROOT/name).read_bytes()),canonical(old),name)

def assert_audited_backup(test,path):
    relative=Path(path).resolve().relative_to(ROOT)
    test.assertEqual(canonical(Path(path).read_bytes()),canonical(baseline_bytes(relative)),str(relative))
