"""Published predictors; unsupported inputs return unavailable, never fabricated zero."""
import numpy as np
from scm_lightning.forcing.contracts import layer_integral, isotherm

def ice_product(rho_d, ri, rg, dz):
    """Romps 2019 eq.4: integral of mass concentrations, kg2 m-5.

    Hydrometeor r is kg per kg DRY air, therefore use dry-air density.
    """
    return layer_integral((np.asarray(rho_d)*ri)*(np.asarray(rho_d)*rg),dz)

def ice_flux(t,z,rho_d,ri,w_convective):
    speed = isotherm(t,z,260.,w_convective)
    flux = isotherm(t,z,260.,np.asarray(rho_d)*ri)*speed
    return flux if speed > 1 else (0. if np.isfinite(speed) else np.nan)

def mccaul(t,z,rg,w_convective,vii):
    # McCaul et al.2009: qg in g/kg for F1; VII in kg/m2 for F2.
    flux=isotherm(t,z,258.15,np.asarray(w_convective))*isotherm(t,z,258.15,np.asarray(rg)*1000.)
    if not np.isfinite(flux):return np.nan
    f1 = .042*max(0.,flux)
    f2 = .20*vii
    if not np.isfinite(isotherm(t,z,258.15,w_convective)):
        return np.nan
    return .95*f1 + .05*f2

def diagnose(fields, metadata):
    """One column/time, explicit semantics; returns method,value,status,reason.

    Original formulations require their original support. Exploratory column/plume
    combinations are implemented separately in lightning.adaptations.
    """
    rows=[]
    def add(name,value=np.nan,reason=''):
        if not np.isfinite(value) and not reason:reason='Missing isotherm or nonfinite required input'
        rows.append({'method':name,'value':float(value),'status':'available' if np.isfinite(value) else 'unavailable','reason':reason})
    required = {'rho_d','qi','qg','dz'}
    if required <= fields.keys() and metadata.get('mass_basis') == 'dry_air':
        add('IxG',ice_product(fields['rho_d'],fields['qi'],fields['qg'],fields['dz']))
    else:
        add('IxG',reason='Requires dry-air density and dry-air mixing ratios, layer thickness')
    if {'cape','precipitation'} <= fields.keys() and metadata.get('cape_definition') == 'near_surface_isentropic':
        add('CAPExP',fields['cape']*fields['precipitation'])
    else:
        add('CAPExP',reason='Requires near-surface isentropic CAPE and interval-aligned total precipitation')
    resolved = metadata.get('velocity_support') == 'resolved_convective' and metadata.get('hydrometeor_support') == 'same_resolved_column'
    if resolved and 'w_convective' in fields:
        w = fields['w_convective']
        add('PW10',fields['precipitation']*(np.max(w)>10) if 'precipitation' in fields else np.nan,'Requires total precipitation')
        if {'temperature','z','rho_d','qi'} <= fields.keys():
            add('IFluxT',ice_flux(fields['temperature'],fields['z'],fields['rho_d'],fields['qi'],w))
        else:
            add('IFluxT',reason='Missing temperature/cloud ice/density/height')
        if {'temperature','z','qg','vii'} <= fields.keys():
            add('McCaul2009',mccaul(fields['temperature'],fields['z'],fields['qg'],w,fields['vii']))
        else:
            add('McCaul2009',reason='Missing temperature/graupel/VII/height')
        add('PR92_wmax',5.e-6*max(0.,np.max(w))**4.5)
    else:
        for method in ['PW10','IFluxT','McCaul2009','PR92_wmax']:
            add(method,reason='Original CRM support unavailable; separate LS/EDMF 1D adaptations are operational in lightning.adaptations')
    for method in ['precipitation_ice_mass','IWP18','ice_mass_flux_product','updraft_volume','PR92_echo_top']:
        add(method,reason='Requires resolved storm geometry/reflectivity/velocity or fall speed; no SCM substitution')
    return rows
