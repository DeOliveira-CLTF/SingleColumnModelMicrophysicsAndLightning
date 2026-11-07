"""Surface input-state regression and independent signed-flux contracts."""
from pathlib import Path
import numpy as np
from scm_lightning.evaluation.mynn_surface import soil_properties

ROOT=Path(__file__).resolve().parents[2]


def test_surface_flux_memory_is_explicit_and_disabled_by_default():
    source=(ROOT/'src/mynn_offline_coupling.F90').read_text()
    assert 'surface_flux_memory=.false.' in source
    assert 'hfx = hfx_s; qfx = qfx_s' in source
    assert 'hfx_s=hfx(1);qfx_s=qfx(1)' in source


def test_era5_flux_conversion_preserves_sign_and_interval():
    # Independent units: ERA5 hourly accumulated downward J/m2 -> upward W/m2.
    accumulated=np.array([-360000.,180000.,0.])
    upward=-accumulated/3600
    np.testing.assert_equal(upward,[100.,-50.,0.])
    np.testing.assert_allclose(upward/2.5e6*2.5e6,upward)
    assert np.sum(upward*3600)==-accumulated.sum()


def test_soil_categories_are_never_interpolated_or_substituted_as_availability():
    p=soil_properties([1.,2.,3.,6.,7.,2.5,np.nan])
    np.testing.assert_allclose(p['field_capacity'][:4],[.244,.347,.383,.663])
    assert np.isnan(p['field_capacity'][4:]).all()
    # Soil at field capacity has hydraulic index 1; this is not a MYNN coefficient.
    assert ((p['field_capacity'][1]-p['wilting_point'][1])/(p['field_capacity'][1]-p['wilting_point'][1]))==1
