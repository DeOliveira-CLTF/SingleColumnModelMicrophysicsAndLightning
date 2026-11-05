"""Compile and test the actual transport subroutine and original pulse module."""
import unittest
import tempfile
from executar_ensemble_goamazon import DEST, ROOT, RESULTS, INPUT, sha, linux_command, linux_path, json


class EnsembleTests(unittest.TestCase):
    def test_parameterizations_unchanged(self):
        baseline=json.loads((RESULTS/'proveniencia/execution_20261102T214208919408Z.json').read_text(encoding='utf-8'))
        for name,h in baseline['physics_sha256'].items():
            if name.replace('\\','/').endswith('/ensemble_driver.F90'): continue
            if name.replace('\\','/') in ['src/module_mp_thompson.F','src/wsm6/mp_wsm6.F90']:
                from instrumentar_processos import PROV,strip_hooks
                from pathlib import Path
                backup=PROV/Path(name).name
                self.assertEqual(sha(backup),h,name)
                self.assertEqual(strip_hooks((ROOT/name).read_text(encoding='utf-8')),backup.read_text(encoding='utf-8'))
                continue
            self.assertEqual(sha(ROOT/name),h,name)
        self.assertEqual(sha(INPUT),baseline['input_sha256'])

    def test_actual_fortran_transport_and_pulse(self):
        manifest=json.loads((DEST/'manifest.json').read_text(encoding='utf-8'))
        build=(ROOT/manifest['binary']).parent
        code=(ROOT/'src/ensemble_driver.F90').read_text(encoding='utf-8')
        routine=code[code.index('  subroutine advect_upwind'):code.index('  end subroutine advect_upwind')+len('  end subroutine advect_upwind')]
        program='''program test
use mpas_kind_types, only: RKIND
use updraft_forcing, only: combined_w
implicit none
real(kind=RKIND) :: z(5),w(5),th(5),q(5),v(5),background(5)
z=[0._RKIND,100._RKIND,200._RKIND,300._RKIND,400._RKIND]
w=3._RKIND; th=300._RKIND; q=.01_RKIND
call advect_upwind(5,z,1000._RKIND,w,th,q)
if(maxval(abs(th-300._RKIND))>1.e-10_RKIND) error stop 'constant T'
if(maxval(abs(q-.01_RKIND))>1.e-12_RKIND) error stop 'constant q'
th=[300._RKIND,310._RKIND,320._RKIND,330._RKIND,340._RKIND]; q=.01_RKIND
call advect_upwind(5,z,10._RKIND,w,th,q)
if(abs(th(3)-317._RKIND)>1.e-10_RKIND) error stop 'linear translation'
th=[300._RKIND,300._RKIND,330._RKIND,300._RKIND,300._RKIND]
q=[0._RKIND,0._RKIND,.02_RKIND,0._RKIND,0._RKIND]
w=[-30._RKIND,-30._RKIND,30._RKIND,30._RKIND,30._RKIND]
call advect_upwind(5,z,100._RKIND,w,th,q)
if(minval(th)<300._RKIND-1.e-10_RKIND.or.maxval(th)>330._RKIND+1.e-10_RKIND) error stop 'bounds'
if(minval(q)<0._RKIND.or.maxval(q)>.02_RKIND+1.e-12_RKIND) error stop 'q bounds'
background=.01_RKIND
v=combined_w(5,z,background,900._RKIND,200._RKIND,1000._RKIND,1800._RKIND,3._RKIND)
if(abs(v(3)-3.01_RKIND)>1.e-12_RKIND) error stop 'pulse peak'
v=combined_w(5,z,background,1800._RKIND,200._RKIND,1000._RKIND,1800._RKIND,3._RKIND)
if(maxval(abs(v-background))>1.e-12_RKIND) error stop 'pulse end'
v=combined_w(5,z,background,900._RKIND,200._RKIND,1000._RKIND,1800._RKIND,0._RKIND)
if(maxval(abs(v-background))>1.e-12_RKIND) error stop 'zero pulse background'
print *, 'transport and pulse OK'
contains
'''+routine+'\nend program test\n'
        # Isolated workspace temp directory; preserved scientific outputs untouched.
        with tempfile.TemporaryDirectory(dir=DEST) as tmp:
            from pathlib import Path
            folder=Path(tmp); source=folder/'test.F90'; source.write_text(program,encoding='ascii')
            binary=folder/'test'
            linux_command(['gfortran','-ffree-line-length-none','-fdefault-real-8','-fdefault-double-8','-I'+linux_path(build),
                           linux_path(source),linux_path(build/'updraft_forcing.o'),'-o',linux_path(binary)])
            result=linux_command([linux_path(binary)])
            self.assertIn('transport and pulse OK',result.stdout)


if __name__=='__main__': unittest.main(verbosity=2)
