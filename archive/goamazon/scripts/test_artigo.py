"""Boundary cases requested for the scientific workflow, using copies only."""
import unittest
import tempfile
import shutil
import subprocess
import sys
import pandas as pd
import netCDF4
import numpy as np
from scm_artigo import *
from validar_resultados import validate, validate_file

class ScientificTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.th=load(RESULTS/'TH_CF.nc'); cls.ws=load(RESULTS/'WS_CF.nc')
        cls.cfg=parse_nml(nml_path('TH_CF')); cls.wcfg=parse_nml(nml_path('WS_CF'))

    def test_reference_reproduction(self):
        # References rounded by user: half-unit rounding, with small FP margin.
        for ds,cfg,rain,graupel,flash,integral in [(self.th,self.cfg,.922617,0.,0.,0.),(self.ws,self.wcfg,2.706833,2.559396e-5,.02131269,2.939459)]:
            m,_=metrics(ds,cfg)
            self.assertAlmostEqual(m['rain_final_mm'],rain,delta=5.1e-7)
            self.assertAlmostEqual(m['qg_max_g_kg']/1000,graupel,delta=5.1e-12)
            self.assertAlmostEqual(m['flash_max_min'],flash,delta=5.1e-9)
            self.assertAlmostEqual(m['flash_integral_nominal'],integral,delta=5.1e-7)

    def test_zero_series(self):
        m,_=metrics(self.th,self.cfg)
        self.assertEqual(m['flash_integral_nominal'],0)
        for k in ['flash_onset_h','flash_peak_h','qg_peak_h','corr_flash_qg_mixed_mass']:
            self.assertIsNone(m[k])

    def test_censored(self):
        ds=self.th.isel(time=slice(0,4)).copy(deep=True)
        ds.flash_rate.values[:]=[0,1,0,2]
        m,_=metrics(ds,self.cfg)
        self.assertTrue(m['flash_end_censored'])
        self.assertAlmostEqual(m['flash_active_duration_h'],60/3600)

    def test_nonfinite(self):
        ds=self.th.copy(deep=True); ds.qg.values[1,0]=np.nan
        self.assertFalse(validate(ds,self.cfg,'TH_CF')['valid'])
        with self.assertRaises(ValueError): metrics(ds,self.cfg)

    def test_irregular_time(self):
        ds=self.th.isel(time=slice(0,4)).copy(deep=True).assign_coords(time=[0.,10.,40.,100.])
        ds.time.attrs['units']='seconds since t0'
        ds.flash_rate.values[:]=[0,1,2,3]
        ds.rainncv.values[:]=[0,1,2,3]; ds.rainnc.values[:]=[0,1,3,6]
        m,s=metrics(ds,self.cfg)
        self.assertAlmostEqual(m['flash_integral_nominal'],250/60)
        np.testing.assert_allclose(s['rain_rate_mm_h'][1:],[360,240,180])

    def test_missing_file(self):
        self.assertFalse(validate_file(RESULTS/'__absent__.nc',self.cfg)['valid'])

    def test_config_mismatch(self):
        self.assertFalse(validate(self.th,self.wcfg,'TH_CF')['valid'])
        cfg=dict(self.cfg,nz=59)
        self.assertFalse(validate(self.th,cfg,'TH_CF')['valid'])

    def test_initial_identity(self):
        for q in ['theta','qv','pressure']+SPECIES:
            np.testing.assert_array_equal(self.th[q][0],self.ws[q][0])

    def test_scalar_integer_roundtrip(self):
        from executar_experimentos import write_nml
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'case.nml'; write_nml(path,self.cfg)
            parsed=parse_nml(path)
            self.assertIsInstance(parsed['mynn_edmf'],int)
            self.assertEqual(parsed,self.cfg)

    def test_partial_pipeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            shutil.copy2(nml_path('TH_CF'),directory/'TH_CF.nml')
            shutil.copy2(RESULTS/'TH_CF.nc',directory/'TH_CF.nc')
            subprocess.run([sys.executable,str(ROOT/'scripts/analisar_experimentos.py'),'--results',tmp],check=True,capture_output=True)
            rows=pd.read_csv(directory/'tabelas/Tabela2_resultados.csv').set_index('experiment')
            self.assertEqual(rows.loc['WS_CF','status'],'ausente')
            self.assertTrue(pd.isna(rows.loc['WS_CF','rain_final_mm']))
            self.assertEqual(rows.loc['TH_CF','flash_max_min'],0.)
            contrasts=pd.read_csv(directory/'tabelas/Tabela3_contrastes.csv')
            missing=contrasts[contrasts.contrast=='WS_CF - TH_CF']
            self.assertTrue(missing.difference.isna().all())

    def test_independent_integrals_all_cases(self):
        rows=pd.read_csv(RESULTS/'tabelas/Tabela2_resultados.csv').set_index('experiment')
        for eid in IDS:
            path=RESULTS/f'{eid}.nc'
            if not path.exists(): continue
            with netCDF4.Dataset(path) as ds:
                t=ds['time'][:]; flash=ds['flash_rate'][:]
                total=sum(float(flash[i])*(float(t[i])-float(t[i-1]))/60 for i in range(1,len(t)))
                self.assertAlmostEqual(total,rows.loc[eid,'flash_integral_nominal'],delta=1e-12)
                self.assertAlmostEqual(float(sum(ds['rainncv'][1:])),float(ds['rainnc'][-1]-ds['rainnc'][0]),delta=1e-12)
                # Independent geometric-grid reconstruction, from grid endpoints.
                lo,hi=1.,2.
                for _ in range(100):
                    r=(lo+hi)/2
                    if sum(50*r**k for k in range(60))>20000: hi=r
                    else: lo=r
                dz=50*((lo+hi)/2)**np.arange(60)
                centers=np.cumsum(dz)-dz/2
                np.testing.assert_allclose(ds['lev'][:],centers,rtol=0,atol=1e-7)
                temp=ds['theta'][:]*(ds['pressure'][:]/100000)**(2/7)
                density=ds['pressure'][:]/(287*temp*(1+.61*ds['qv'][:]))
                mass=np.sum(density*ds['qg'][:]*dz*((temp>=233.15)&(temp<=273.15)),axis=1)
                self.assertAlmostEqual(float(mass.max()),rows.loc[eid,'qg_mixed_mass_max_kg_m2'],delta=1e-12)

if __name__=='__main__': unittest.main(verbosity=2)
