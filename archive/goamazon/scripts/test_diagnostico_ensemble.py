"""Independent arithmetic and provenance checks for scientific diagnostics."""
import unittest
import pandas as pd
from diagnosticar_ensemble_goamazon import *


class DiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.baseline=json.loads((BASE/'manifest.json').read_text(encoding='utf-8'))
        cls.results=pd.read_csv(OUT/'agua_introduzida_12_membros.csv').set_index('experiment')

    def test_water_mass_independent_scalar_sum(self):
        src=load(INPUT)
        for eid in CASES:
            r=self.baseline['members'][eid]; ds=load(ROOT/r['directory']/'profile_1.nc')
            z=ds.lev.values; edges=interfaces(z,20000); total=0.
            for k,zz in enumerate(z):
                p=float(np.interp(zz,src.lev,src.pa.values.squeeze()))
                th=float(np.interp(zz,src.lev,src.theta.values.squeeze()))
                rv=float(np.interp(zz,src.lev,src.rv.values.squeeze()))
                old=sum(float(np.interp(zz,src.lev,src[q].values.squeeze())) for q in ['rv','rl','ri'])
                new=sum(float(ds[q].values[0,k]) for q in ['qv','qc','qi'])
                T=th*(p/100000)**(2/7)
                # p/(Rd*T) is total pressure expressed as a dry-air density;
                # divide by 1+r*Rv/Rd to obtain actual dry-air density.
                density=p/(287*T)/(1+rv*461.6/287)
                total+=density*(new-old)*(edges[k+1]-edges[k])
            self.assertAlmostEqual(total,self.results.loc[eid,'water_added_fixed_dry_air_mm'],delta=1e-11)

    def test_ep_reconstruction_and_time_offset(self):
        for eid in CASES:
            r=self.baseline['members'][eid]; ds=load(ROOT/r['directory']/'profile_1.nc')
            T,w,density,intervals,ep,outside=reconstruct_ep(ds,r['config'])
            np.testing.assert_allclose(.001*np.sqrt(ep[1:]),ds.flash_rate.values[1:],rtol=5e-11,atol=1e-12)
            # Pulse is maximal at forcing time 900 s, saved at 920 s (dt20).
            self.assertEqual(int(np.argmax(w.max(axis=1))),46)
            np.testing.assert_allclose(w[0],w[-1],rtol=0,atol=1e-14)
            self.assertTrue(np.all(outside<=ep+1e-12))

    def test_matrix_completeness_and_initial_water_independence_of_dt(self):
        df=pd.read_csv(OUT/'sensibilidade_resolucao.csv')
        self.assertEqual(len(df),36)
        self.assertTrue((df.status=='valid').all())
        self.assertEqual(len(df[['case','nz','dt_s']].drop_duplicates()),36)
        for _,part in df.groupby(['case','nz']):
            self.assertEqual(part.water_added_fixed_dry_air_mm.nunique(),1)
        for case,part in df.groupby('case'):
            expected=case!='thompson_media_w3'
            self.assertTrue(((part.flash_integral_nominal>0)==expected).all())

    def test_input_sources_and_original_outputs_preserved(self):
        self.assertEqual(sha(INPUT),self.baseline['identity']['input_sha256'])
        from historical_contracts import assert_audited_source
        for path,digest in self.baseline['identity']['sources'].items(): assert_audited_source(self,path)
        for r in self.baseline['members'].values(): self.assertEqual(sha(ROOT/r['directory']/'profile_1.nc'),r['output_sha256'])


if __name__=='__main__': unittest.main(verbosity=2)
