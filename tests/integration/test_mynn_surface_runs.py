import json
import numpy as np
import pandas as pd
import pytest
from scm_lightning.evaluation.mynn_surface import OUT


def product(name):
    p=OUT/name
    if not p.exists():pytest.skip('Surface experiment not yet produced')
    return p


def test_surface_scenarios_and_reference_identity():
    r=json.loads(product('run_tests.json').read_text())
    assert len(r['runs'])==12
    reference=[v for v in r['runs'] if v['scenario'] in ['reference_off','corrected_bulk']]
    assert len(reference)==8 and all(v['baseline_all_variables_bit_identical'] for v in reference)
    assert r['ensemble_reruns']==0


def test_surface_contracts_and_no_automatic_reference_selection():
    r=json.loads(product('validation.json').read_text())
    assert len(r['checks'])==6 and r['protected_originals_unchanged']
    assert r['selected_reference'] is None
    for v in r['checks']:assert v['native_Beljaars_wstar_reconstruction']=='passed'
    for cfg in ['TH_CF','WS_CF']:
        d=pd.read_csv(product(f'bulk_memory_{cfg}_surface.csv'))
        np.testing.assert_array_equal(d.incoming_hfx.values[1:],d.applied_hfx.values[:-1])
        assert d.wstar.max()>0
        d=pd.read_csv(product(f'prescribed_era5_{cfg}_surface.csv'))
        np.testing.assert_array_equal(d.applied_hfx,d.forcing_hfx)
        np.testing.assert_array_equal(d.applied_lh,d.forcing_lh)


def test_land_sources_are_real_and_not_mapped_to_mavail():
    r=json.loads(product('land_context.json').read_text())
    assert {v['event'] for v in r}=={'argentina_20181213','argentina_20191028'}
    assert all(v['mavail_assigned'] is None for v in r)
    assert all(v['mapped_soil_land_fraction_min']>0 for v in r)


def test_electrical_diagnostics_are_uncalibrated_and_separate():
    d=pd.read_csv(product('raw_predictor_contrasts.csv'))
    assert len(d)==12
    assert json.loads(product('predictor_batch.json').read_text())['complete']
    assert set(d.method)=={'IxG','MUCAPExP_300hPa_SCM','McCaulF3'}
    assert not d.calibrated.any()
    assert len(list((OUT/'figures').glob('*.png')))==3
