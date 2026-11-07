import numpy as np
import pytest
from metpy.calc import most_unstable_cape_cin,dewpoint,vapor_pressure
from metpy.units import units
from scm_lightning.lightning.thermodynamics import diagnose


def independent_profile(kind):
    p=np.array([1000,925,850,775,700,600,500,400,300,250,200,150,100,70,50.])*100
    z=8000*np.log(100000/p)
    if kind=='stable':t=300-2*z/1000
    else:t=303-7*z/1000+np.maximum(z-11000,0)*.006
    rv=.016*np.exp(-z/2200)
    if kind=='elevated':
        t[:2]-=8;rv[:2]*=.1;rv[2:5]*=1.4
    if kind=='multilayer':t[5:8]+=np.array([7.,-4.,7.])
    return p,t,rv,z


@pytest.mark.parametrize('kind',['stable','unstable','elevated','multilayer'])
def test_mucape_matches_official_metpy_with_si_inputs_and_sorting(kind):
    p,t,r,z=independent_profile(kind)
    td=dewpoint(vapor_pressure(p*units.Pa,r*units.dimensionless))
    a,b=most_unstable_cape_cin(p*units.Pa,t*units.K,td,depth=300*units.hPa)
    result=diagnose(p,t,r,z)
    assert result['cape_J_kg']==pytest.approx(a.to('J/kg').magnitude,rel=1e-12)
    assert result['cin_J_kg']==pytest.approx(b.to('J/kg').magnitude,rel=1e-12)
    reverse=diagnose(p[::-1],t[::-1],r[::-1],z[::-1])
    assert reverse['cape_J_kg']==pytest.approx(result['cape_J_kg'])
    assert result['parcel_pressure_Pa'] in p
    assert result['parcel_height_AGL_m']==z[np.flatnonzero(p==result['parcel_pressure_Pa'])[0]]
    assert result['status']==('valid_zero' if a.magnitude==0 else 'valid_positive')
    if kind=='stable':assert result['cape_J_kg']==0
    else:assert result['cape_J_kg']>0
    if kind=='multilayer':
        from metpy.calc import lfc,parcel_profile_with_lcl
        pp,tt,dd,path=parcel_profile_with_lcl(p*units.Pa,t*units.K,td)
        roots,_=lfc(pp,tt,dd,parcel_temperature_profile=path,which='all')
        assert np.size(roots.magnitude)>=2


def test_missing_incomplete_and_dry_profiles_are_not_filled_with_zero():
    p,t,r,z=independent_profile('unstable')
    short=diagnose(p[:6],t[:6],r[:6],z[:6])
    assert short['status']=='incomplete_profile' and np.isnan(short['cape_J_kg'])
    bad=t.copy();bad[5]=np.nan
    assert diagnose(p,bad,r,z)['status']=='invalid_profile'
    dry=r.copy();dry[-1]=0
    assert diagnose(p,t,dry,z)['status']=='undefined_dewpoint'
    duplicate=p.copy();duplicate[1]=duplicate[0]
    assert diagnose(duplicate,t,r,z)['status']=='invalid_profile'


@pytest.mark.parametrize('kind',['stable','unstable','multilayer'])
@pytest.mark.parametrize('selection',['bottom','most_cape'])
def test_virtual_temperature_limits_and_signed_areas_reconstruct_provider(kind,selection):
    from scm_lightning.lightning.parcel_levels import integration_limits
    from metpy.calc import parcel_profile_with_lcl,cape_cin
    p,t,r,z=independent_profile(kind)
    pressure=(p/100)*units.hPa;temperature=t*units.K
    td=dewpoint(vapor_pressure(pressure,r*units.dimensionless))
    prepared=parcel_profile_with_lcl(pressure,temperature,td)
    el_selection='top' if selection=='bottom' else 'most_cape'
    expected,_=cape_cin(*prepared,which_lfc=selection,which_el=el_selection)
    result=integration_limits(pressure,temperature,td,which_lfc=selection,which_el=el_selection)
    assert result['signed_integral_J_kg']==pytest.approx(expected.to('J/kg').magnitude,abs=1e-7,rel=1e-10)
    assert result['positive_buoyancy_contribution_J_kg']>=0
    assert result['negative_buoyancy_contribution_J_kg']<=0
    assert result['signed_integral_J_kg']==pytest.approx(result['positive_buoyancy_contribution_J_kg']+result['negative_buoyancy_contribution_J_kg'],abs=1e-7)
