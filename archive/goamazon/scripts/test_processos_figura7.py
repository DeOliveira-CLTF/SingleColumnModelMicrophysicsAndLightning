import unittest
import tempfile
import pandas as pd
from scm_artigo import *
from executar_processos_figura7 import OUT,QS,CATEGORY
from executar_experimentos import write_nml,linux_command,linux_path
from instrumentar_processos import PROV,strip_hooks


class ProcessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=json.loads((OUT/'latest.json').read_text(encoding='utf-8'))
        cls.root=ROOT/cls.manifest['directory']

    def test_physical_statements_intact(self):
        baseline=json.loads((RESULTS/'ensemble_goamazon/manifest.json').read_text(encoding='utf-8'))
        for path in ['src/module_mp_thompson.F','src/wsm6/mp_wsm6.F90']:
            saved=PROV/Path(path).name
            self.assertEqual(strip_hooks((ROOT/path).read_text(encoding='utf-8')),saved.read_text(encoding='utf-8'))
            hashes={p.replace('\\','/'):v for p,v in baseline['identity']['sources'].items()}
            from historical_contracts import assert_audited_backup
            assert_audited_backup(self,saved)

    def test_all_original_outputs_and_species_closure(self):
        baseline=json.loads((RESULTS/'ensemble_goamazon/manifest.json').read_text(encoding='utf-8'))
        self.assertEqual(len(self.manifest['runs']),12)
        for eid,rec in self.manifest['runs'].items():
            folder=ROOT/rec['directory']; ds=load(folder/'profile_1.nc')
            old=load(ROOT/baseline['members'][eid]['directory']/'profile_1.nc')
            for key in old.variables: np.testing.assert_array_equal(ds[key],old[key])
            rates=load(folder/'process_tendencies.nc')
            names=[v.decode() if isinstance(v,bytes) else str(v) for v in rates.process_name.values]
            values=rates.delta_q.values
            residual=values[:,:,names.index('unattributed_residual'),:]
            self.assertLess(float(np.max(np.abs(residual))),1e-12)
            for iq,q in enumerate(QS[1:],1):
                # Exclude residual: test actual physical+numerical terms close independently.
                total=values[:,:,:,iq].sum(axis=2)-residual[:,:,iq]
                np.testing.assert_allclose(total,np.diff(ds[q].values,axis=0),atol=1e-12,rtol=1e-10)
            if eid.startswith('wsm6'):
                a=values[:,:,names.index('ws_paacw_snow'),4]
                b=values[:,:,names.index('ws_paacw_graupel'),5]
                np.testing.assert_array_equal(a,b)

    def test_diagnostics_off_regression(self):
        binary=self.root/'build/ensemble_driver'
        for eid in ['thompson_alta_w3','wsm6_media_w3']:
            rec=self.manifest['runs'][eid]
            with tempfile.TemporaryDirectory(dir=OUT) as tmp:
                folder=Path(tmp); cfg=dict(rec['config'],write_process_diagnostics=False)
                write_nml(folder/'ensemble.nml',cfg,'ensemble')
                linux_command([linux_path(binary)],folder)
                off=load(folder/'profile_1.nc'); on=load(ROOT/rec['directory']/'profile_1.nc')
                for key in on.variables: np.testing.assert_array_equal(off[key],on[key])


if __name__=='__main__': unittest.main(verbosity=2)
