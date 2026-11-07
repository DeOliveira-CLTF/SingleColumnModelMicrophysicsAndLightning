"""Regression for the dry-mixing-ratio / MYNN specific-humidity boundary."""
from pathlib import Path
import os
import re
import shutil
import subprocess
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_boundary_adapter_is_explicit_and_opt_in():
    source = (ROOT/'src/mynn_offline_coupling.F90').read_text()
    assert 'sqv2(1,:) = qv/(1._RKIND+qv)' in source
    assert 'boundary_fix = .false.' in source
    assert 'if(boundary_fix) qsfc1=0._RKIND' in source


def test_exact_finite_step_moisture_conversion():
    rv=np.array([0.,.01,.025]); rc=np.array([0.,.001,.002]); dt=30.
    qv=rv/(1+rv); qc=rc/(1+rv)
    # Independent reference: after vapor-to-cloud conversion, total specific water stays fixed.
    transfer=np.array([0.,1e-5,2e-5]); qv1=qv-transfer; qc1=qc+transfer
    drc=(qc1/(1-qv1)-rc)/dt
    drv=(qv1/(1-qv1)-rv)/dt
    np.testing.assert_allclose((rc+dt*drc)/(1+rv+dt*drv),qc1,atol=1e-18)
    np.testing.assert_allclose((rv+dt*drv)/(1+rv+dt*drv),qv1,atol=1e-18)
    assert np.max(abs(drc-transfer/dt))>1e-10


def test_compiled_actual_fortran_boundary_blocks(tmp_path):
    """Execute the actual packing/inverse statements against an independent reference."""
    compiler=shutil.which('gfortran')
    if not compiler and not (os.name=='nt' and shutil.which('wsl')):pytest.skip('Fortran compiler unavailable')
    source=(ROOT/'src/mynn_offline_coupling.F90').read_text()
    blocks=re.findall(r'if\(boundary_fix\) then\n(.*?)\n    endif',source,re.S)
    assert len(blocks)==2
    code='''program check_adapter
    implicit none
    integer,parameter::RKIND=kind(1.d0)
    real(RKIND)::qv(3),qc(3),qi(3),qs(3),dt
    real(RKIND)::sqv2(1,3),sqc2(1,3),sqi2(1,3),sqs2(1,3)
    real(RKIND)::rqvblten2(1,3),rqcblten2(1,3),rqiblten2(1,3)
    qv=[0._RKIND,.01_RKIND,.025_RKIND];qc=[0._RKIND,.001_RKIND,.002_RKIND]
    qi=qc*.1_RKIND;qs=0.;dt=30.
    '''+blocks[0]+'''
    rqvblten2(1,:)=[0._RKIND,-1.e-5_RKIND,-2.e-5_RKIND]/dt
    rqcblten2=-rqvblten2;rqiblten2=0.
    '''+blocks[1]+'''
    print *,qv+dt*rqvblten2(1,:)
    print *,qc+dt*rqcblten2(1,:)
    print *,qi+dt*rqiblten2(1,:)
    end program
    '''
    path=tmp_path/'adapter.f90';path.write_text(code)
    if compiler:
        exe=tmp_path/'adapter';subprocess.run([compiler,'-ffree-line-length-none',str(path),'-o',str(exe)],check=True,capture_output=True)
        output=subprocess.check_output([str(exe)],text=True)
    else:
        pathlinux='/mnt/'+path.drive[0].lower()+path.as_posix()[2:]
        exe=pathlinux+'.exe'
        subprocess.run(['wsl','-d','Ubuntu-24.04','--','gfortran','-ffree-line-length-none',pathlinux,'-o',exe],check=True,capture_output=True)
        output=subprocess.check_output(['wsl','-d','Ubuntu-24.04','--',exe],text=True)
    actual=np.fromstring(output,sep=' ').reshape(3,3)
    rv=np.array([0.,.01,.025]);rc=np.array([0.,.001,.002]);rv1=rv/(1+rv)-np.array([0.,1e-5,2e-5])
    expected=np.array([rv1/(1-rv1),(rc/(1+rv)+np.array([0.,1e-5,2e-5]))/(1-rv1),.1*rc/(1+rv)/(1-rv1)])
    np.testing.assert_allclose(actual,expected,rtol=1e-13,atol=1e-18)
