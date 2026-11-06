import unittest
import numpy as np
from orcamento_coluna_prescrita import transport,OUT
from scm_artigo import ROOT,load,json
import pandas as pd


class BudgetTests(unittest.TestCase):
    def test_constant_tracer_with_divergent_prescribed_air(self):
        z=np.array([50.,175.,400.]); edges=np.array([0.,100.,250.,550.])
        mass=np.array([120.,130.,200.]); q=np.full(3,.01)
        result=transport(q,np.array([2.,-3.,1.]),z,edges,mass,1000.)
        self.assertAlmostEqual(result['reconstructed_transport_mm'],0.,places=13)
        self.assertAlmostEqual(result['transport_closure_mm'],0.,places=12)

    def test_uniform_grid_linear_profile_one_step(self):
        # Lower cell fixed by zero-gradient inflow; two interior cells lose .001.
        result=transport(np.array([.01,.02,.03]),np.ones(3),np.array([50.,150.,250.]),
                         np.array([0.,100.,200.,300.]),np.full(3,100.),10.)
        self.assertAlmostEqual(result['reconstructed_transport_mm'],-.2,places=13)
        self.assertAlmostEqual(result['transport_closure_mm'],0.,places=13)

    def test_saved_stage_budget_and_toggle(self):
        manifest=json.loads((OUT/'latest.json').read_text(encoding='utf-8'))
        for item in manifest['runs'].values():
            folder=ROOT/item['directory']; ds=load(folder/'profile_1.nc')
            other=load(folder/'diagnostics_off/profile_1.nc')
            for key in ds.variables: np.testing.assert_array_equal(ds[key],other[key])
            table=pd.read_csv(folder/'budget_audit.csv')
            z=ds.lev.values; edges=np.r_[0.,np.zeros(len(z))]
            for k in range(len(z)): edges[k+1]=2*z[k]-edges[k]
            T=ds.theta.values*(ds.pressure.values/1e5)**(287./1004.5)
            rho=ds.pressure.values/(T*(287.+461.6*ds.qv.values))
            rt=sum(ds[key].values for key in ['qv','qc','qr','qi','qs','qg'])
            storage=(rho*rt*np.diff(edges)).sum(axis=1)
            np.testing.assert_allclose(table.storage_before_mm,storage[:-1],atol=1e-11,rtol=1e-12)
            np.testing.assert_allclose(table.storage_after_mm,storage[1:],atol=1e-11,rtol=1e-12)
            self.assertLess(table.bookkeeping_residual_mm.abs().max(),1e-11)


if __name__=='__main__': unittest.main(verbosity=2)
