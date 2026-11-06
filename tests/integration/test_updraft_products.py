"""Optional local regression of the executed 18-case experiment; no private data shipped."""
import json
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from scm_lightning.ensemble.updraft_sensitivity import OUT, RUNS, CONFIG
from scm_lightning.config import load_config
from scm_lightning.provenance import sha256

pytestmark=pytest.mark.skipif(not (OUT/'manifest.json').exists(), reason='18-case local products are not included in source distribution')


def test_eighteen_cases_frozen_design_and_precursor_identity():
    checks=json.loads((OUT/'run_checks.json').read_text())
    design=load_config(CONFIG)
    assert checks['run_count']==18 and checks['six_ERA5_controls_bit_identical']
    assert sha256(CONFIG)==checks['design_sha256']
    for item in checks['runs']:
        folder=RUNS/item['regime']/item['configuration']
        control=RUNS/'ERA5'/item['configuration']/'output.nc'
        with xr.open_dataset(folder/'output.nc',decode_times=False) as actual, xr.open_dataset(control,decode_times=False) as baseline:
            before=actual.time.values<=design['pulse']['start_s']
            np.testing.assert_array_equal(actual.pressure.values,baseline.pressure.values)
            for variable in ['theta','qv','qc','qr','qi','qs','qg','rainnc']:
                np.testing.assert_array_equal(actual[variable].values[before],baseline[variable].values[before])


def test_counterfactual_decomposition_and_available_zeros():
    f=pd.read_csv(OUT/'tables/response_decomposition_15min.csv.gz')
    assert len(f)==18*2*15*192
    valid=f[['X00','Xr0','X0r','Xrr']].notna().all(axis=1)
    np.testing.assert_allclose((f.state_response+f.velocity_on_control_state+f.interaction)[valid],f.total_response[valid],atol=1e-10,rtol=1e-10)
    np.testing.assert_allclose((f.state_response+f.velocity_on_perturbed_state)[valid],f.total_response[valid],atol=1e-10,rtol=1e-10)
    p=pd.read_csv(OUT/'tables/predictors_15min.csv.gz')
    assert (p.loc[p.predictor_integral==0,'status']=='available').all()
    assert (p.loc[p.predictor_integral.isna(),'status']=='unavailable').all()
    assert (p.groupby(['regime','configuration','velocity','method']).size()==192).all()


def test_water_stage_accounting_and_forcing_units():
    checks=json.loads((OUT/'run_checks.json').read_text())
    for item in checks['runs']:
        regime,configuration=item['regime'],item['configuration']
        stage=pd.read_csv(OUT/'tables'/(regime+'_'+configuration+'_water_budget.csv.gz'))
        assert len(stage)==5760 and np.isfinite(stage.values).all()
        assert stage.accounting_residual_ref.abs().max()<2e-11
        assert stage.adv_reconstruction_error_ref.abs().max()<2e-11
        with xr.open_dataset(OUT/'netcdf'/(regime+'_'+configuration+'.nc'),decode_times=False) as d:
            for name in ['w_LS','w_perturbation','w_effective','w_EDMF']:assert d[name].units=='m s-1'
            np.testing.assert_array_equal(d.w_effective,d.w_LS+d.w_perturbation)
            assert d.attrs['mynn_uedmf_mp']==0
            np.testing.assert_allclose(d.precipitation_rate.values*30.,d.precipitation_amount.values,atol=1e-12)
            np.testing.assert_allclose(d.precipitation_amount.sum(),stage.precip_kg_m2.sum(),atol=1e-9)
            assert d.time.attrs['bounds']=='time_bounds'
            assert d.time_bounds.values[-1,-1]==172800.
