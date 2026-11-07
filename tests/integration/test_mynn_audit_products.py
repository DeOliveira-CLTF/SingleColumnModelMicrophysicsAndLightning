"""Local run products are optional; tests never launch SCM ensembles."""
import json
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from scm_lightning.evaluation.mynn_audit import OUT,RUNS


def product(name):
    path=OUT/name
    if not path.exists():pytest.skip('Isolated MYNN audit has not been run')
    return path


def test_six_original_controls_bit_identical():
    record=json.loads(product('original_run_tests.json').read_text())
    assert len(record['runs'])==6
    assert all(r['all_native_variables_bit_identical'] for r in record['runs'])


def test_final_binary_reference_controls_bit_identical():
    record=json.loads(product('reference_v2_run_tests.json').read_text())
    assert len(record['runs'])==6
    assert all(r['all_native_variables_bit_identical'] for r in record['runs'])


def test_shadow_call_does_not_modify_state_and_missing_is_not_zero():
    record=json.loads(product('no_shadow_v3_run_tests.json').read_text())
    assert all(r['all_native_variables_bit_identical'] for r in record['runs'])
    for cfg in ['TH_CF','WS_CF']:
        with xr.open_dataset(RUNS/'no_shadow_v3'/cfg/'output.nc.mynn_audit.nc') as d:
            assert d.attrs['local_tendencies_available']==0
            assert np.isnan(d.qc_tendency_local).all() and np.isnan(d.qi_tendency_local).all()


def test_original_isotherms_and_active_counts():
    table=pd.read_csv(product('plume_audit.csv')).query("phase == 'original'").set_index('configuration')
    assert table.loc['TH_CF','active_samples']==1159
    assert table.loc['WS_CF','active_samples']==1213
    assert (table[['reaches_minus15_samples','reaches_260K_samples']]==0).all().all()


def test_corrected_controls_and_preservation_gate():
    record=json.loads(product('audit_tests.json').read_text())
    assert record['protected_files_unchanged'] and not record['calibration_changed']
    assert record['corrected_model_requires_2019_recalibration']
    for cfg in ['TH_CF','WS_CF']:
        check=record['checks']['corrected/'+cfg]
        assert check['surface_saturation_follows_current_Tskin_pressure']=='passed'
        assert check['stage_accounting_max_residual_kg_m2']<1e-10
        assert check['condensate_applied_exactly_once_and_clipping_accounted']=='passed'


def test_only_three_figures_and_raw_corrected_predictors():
    product('audit_tests.json')
    assert len(list((OUT/'figures').glob('*.png')))==3
    for cfg in ['TH_CF','WS_CF']:
        identity=json.loads(product(f'corrected_predictors/{cfg}/audit_identity.json').read_text())
        assert not identity['calibrated'] and not identity['observations_used']
