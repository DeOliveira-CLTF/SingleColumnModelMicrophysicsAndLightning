"""Audit the *virtual-temperature* integration limits used by MetPy CAPE.

Diagnostic reconstruction only: these values never replace the official
most_unstable_cape_cin predictor. Private zero-crossing helpers are versioned
with the installed MetPy implementation, checked by provider-equality tests.
"""
import numpy as np
from scipy.integrate import trapezoid
from metpy.calc import (parcel_profile_with_lcl,lcl,saturation_mixing_ratio,
    virtual_temperature,virtual_temperature_from_dewpoint,lfc,el)
from metpy.calc.thermo import _find_append_zero_crossings,_less_or_close,_greater_or_close
from metpy.constants import Rd
from metpy.units import units


def integration_limits(pressure,temperature,dewpoints,which_lfc='bottom',which_el='top',prepared=None):
    p,t,td,path=prepared if prepared is not None else parcel_profile_with_lcl(pressure,temperature,dewpoints)
    # cape_cin removes nonfinite auxiliary LCL samples before integration.
    # Native inputs are validated upstream; this never fills missing humidity.
    finite=np.isfinite(p.magnitude)&np.isfinite(t.magnitude)&np.isfinite(td.magnitude)&np.isfinite(path.magnitude)
    p,t,td,path=(a[finite] for a in [p,t,td,path])
    pressure_lcl,_=lcl(p[0],t[0],td[0])
    r=np.where(p>pressure_lcl,saturation_mixing_ratio(p[0],td[0]),saturation_mixing_ratio(p,path))
    tv=virtual_temperature_from_dewpoint(p,t,td);parcel_tv=virtual_temperature(path,r)
    lp,_=lfc(p,tv,td,parcel_temperature_profile=parcel_tv,which=which_lfc)
    ep,_=el(p,tv,td,parcel_temperature_profile=parcel_tv,which=which_el)
    result={'LFC_Pa':float(lp.to('Pa').magnitude),'EL_Pa':float(ep.to('Pa').magnitude),
        'integration_top_Pa':float((p[-1] if np.isnan(ep) else ep).to('Pa').magnitude),
        'signed_integral_J_kg':0.,'positive_buoyancy_contribution_J_kg':0.,'negative_buoyancy_contribution_J_kg':0.}
    if np.isnan(lp):return result
    x,y=_find_append_zero_crossings(np.copy(p),(parcel_tv-tv).to(units.K))
    x_pa=x.to('Pa').magnitude
    mask=_less_or_close(x_pa,result['LFC_Pa'])&_greater_or_close(x_pa,result['integration_top_Pa'])
    xx,yy=x_pa[mask],y[mask].to('K').magnitude
    factor=float(Rd.to('J/kg/K').magnitude)
    result.update(signed_integral_J_kg=factor*float(trapezoid(yy,np.log(xx))),
        positive_buoyancy_contribution_J_kg=factor*float(trapezoid(np.maximum(yy,0),np.log(xx))),
        negative_buoyancy_contribution_J_kg=factor*float(trapezoid(np.minimum(yy,0),np.log(xx))))
    return result
