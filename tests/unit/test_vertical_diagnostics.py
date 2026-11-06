import numpy as np
import pytest
from scm_lightning.lightning.adaptations import (
    interface_density, recover_velocity, crossing, lpi_1d, surface_cape,
)
from scm_lightning.lightning import ice_flux, mccaul
from scm_lightning.provenance import ROOT, resolve_path


def test_interface_density_uses_mynn_cross_weights_on_stretched_grid():
    # Interface 100 m above ground, centers 50 and 250 m: 3/4 lower density.
    np.testing.assert_allclose(interface_density([1.2, .8, .4], [100., 300., 500.]), [1.1, .65])


def test_density_weighted_speed_has_no_second_area_division():
    result = recover_velocity([2.4, 0.], [1.2, .9], [.1, 0.])
    assert result[0] == 2.
    assert np.isnan(result[1])
    with pytest.raises(ValueError):
        recover_velocity([1.], [0.], [.1])
    with pytest.raises(ValueError):
        recover_velocity([1.], [1.], [-.1])


def test_isotherm_chooses_lowest_downward_crossing_and_never_extrapolates():
    temperature = np.array([[270., 250., 270., 250.], [280., 275., 270., 265.]])
    result = crossing(temperature, [0., 1000., 2000., 3000.], 260.)
    assert result[0] == 500.
    assert np.isnan(result[1])


def test_flux_factors_are_colocated_before_multiplication():
    # At 260 K, interpolated concentration=2 and speed=3; flux must be 6,
    # whereas interpolating the products would spuriously give 7.
    assert ice_flux([270., 250.], [0., 1000.], [1., 1.], [1., 3.], [2., 4.]) == 6.
    assert ice_flux([270., 250.], [0., 1000.], [1., 1.], [1., 3.], [1., 1.]) == 0.
    f3 = mccaul([268.15, 248.15], [0., 1000.], [.001, .003], [2., 4.], 1.)
    assert f3 == pytest.approx(.95*.042*3.*2.+.05*.20)


def test_lpi_dry_and_equal_liquid_ice_limits():
    shape=(1, 2)
    one=np.ones(shape); zero=np.zeros(shape)
    # qg=qi=qs=1 -> effective Qi=1; liquid qc=1 -> epsilon=1.
    assert lpi_1d(one*260., [100., 300.], one, zero, one, one, one, one*2.)[0] == 4.
    assert lpi_1d(one*260., [100., 300.], zero, zero, zero, zero, zero, one)[0] == 0.
    assert np.isnan(lpi_1d(one*280., [100., 300.], one, zero, one, one, one, one)[0])


def test_cape_stable_and_conditionally_unstable_profiles():
    z=np.linspace(0., 15000., 31)
    p=1.e5*np.exp(-z/8000.)
    assert surface_cape(p, np.full_like(z, 280.), .001*np.exp(-z/2000.)) == pytest.approx(0., abs=1e-7)
    assert surface_cape(p, 300.-.0075*z, .018*np.exp(-z/2000.)) > 500.
    with pytest.raises(ValueError):
        surface_cape(p[::-1], 300.-.0075*z, .018*np.exp(-z/2000.))


def test_legacy_manifest_paths_resolve_to_existing_case_hierarchy():
    assert resolve_path('work/leona/ensemble40/TH_CF/21') == ROOT/'runs/argentina_20181213/thompson_CF/members/21'
    assert resolve_path(ROOT/'work/leona/controls/WS_CL/output.nc') == ROOT/'runs/argentina_20181213/wsm6_CL/control/output.nc'
