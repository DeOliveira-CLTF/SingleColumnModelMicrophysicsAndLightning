"""Official MetPy MUCAPE, dry-air mixing ratios, native levels only.

No parcel is selected from lightning observations. The search layer is the
lowest 300 hPa; maximum equivalent potential temperature defines the parcel.
"""
import warnings
import numpy as np

DEPTH_PA=30000.
DEFINITION='MUCAPE_MetPy_lowest_300hPa_max_theta_e_liquid_pseudoadiabatic'


def profile(p,t,rv,z):
    arrays=[np.asarray(x,float) for x in [p,t,rv,z]]
    if any(a.ndim!=1 for a in arrays) or len({len(a) for a in arrays})!=1:raise ValueError('Unequal 1D thermodynamic profiles')
    if not np.isfinite(arrays).all():raise ValueError('Missing thermodynamic sample; no imputation')
    p,t,rv,z=arrays
    order=np.argsort(-p);p,t,rv,z=(a[order] for a in arrays)
    if np.any(p<=0) or np.any(t<=0) or np.any(rv<0) or np.any(np.diff(p)>=0) or np.any(np.diff(z)<=0):
        raise ValueError('Invalid pressure, height, temperature or mixing ratio')
    return p,t,rv,z


def diagnose(p,t,rv,z,depth_pa=DEPTH_PA,comparison=False):
    from metpy.calc import (dewpoint,vapor_pressure,most_unstable_parcel,most_unstable_cape_cin,
        surface_based_cape_cin,parcel_profile_with_lcl,cape_cin,lfc,el)
    from metpy.units import units
    names=['cape_J_kg','cin_J_kg','parcel_pressure_Pa','parcel_height_AGL_m','parcel_temperature_K',
        'parcel_mixing_ratio_kg_kg','parcel_index','search_layer_top_height_AGL_m','raw_cape_J_kg']
    result={name:np.nan for name in names};result.update(status='invalid_profile',definition=DEFINITION,search_depth_Pa=depth_pa)
    try:p,t,rv,z=profile(p,t,rv,z)
    except ValueError as e:result['reason']=str(e);return result
    # Do not silently shorten the parcel-search layer or fill a dry level.
    if len(p)<4 or p[-1]>p[0]-depth_pa or p[-1]>10000 or z[-1]-z[0]<15000:
        result.update(status='incomplete_profile',reason='Insufficient search/ascent depth (requires top <=100hPa and 15km vertical span)');return result
    if np.any(rv==0):
        result.update(status='undefined_dewpoint',reason='Exactly dry native level; no artificial humidity');return result
    pressure=p*units.Pa;temperature=t*units.K
    td=dewpoint(vapor_pressure(pressure,rv*units.dimensionless))
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always',UserWarning)
            _,_,_,idx=most_unstable_parcel(pressure,temperature,td,depth=depth_pa*units.Pa)
            cape,cin=most_unstable_cape_cin(pressure,temperature,td,depth=depth_pa*units.Pa)
            value=float(cape.to('J/kg').magnitude);inhibition=float(cin.to('J/kg').magnitude)
            result.update(raw_cape_J_kg=value,cin_J_kg=inhibition,parcel_index=int(idx),
                parcel_pressure_Pa=float(p[idx]),parcel_height_AGL_m=float(z[idx]),
                parcel_temperature_K=float(t[idx]),parcel_mixing_ratio_kg_kg=float(rv[idx]),
                search_layer_top_height_AGL_m=float(np.interp(np.log(p[0]-depth_pa),np.log(p[::-1]),z[::-1])),
                warnings=' | '.join(sorted({str(w.message) for w in caught})))
            if np.isfinite(value) and np.isfinite(inhibition) and value>=0:
                result.update(cape_J_kg=value,status='valid_zero' if value==0 else 'valid_positive',reason='')
            else:result.update(status='invalid_calculation',reason='Nonfinite or negative official CAPE/CIN')
            if comparison:
                sb,sb_cin=surface_based_cape_cin(pressure,temperature,td)
                ps,env,dew,parcel=parcel_profile_with_lcl(pressure,temperature,td)
                old,old_cin=cape_cin(ps,env,dew,parcel,which_lfc='most_cape',which_el='most_cape')
                result.update(SBCAPE_J_kg=float(sb.to('J/kg').magnitude),SBCIN_J_kg=float(sb_cin.to('J/kg').magnitude),
                    previous_most_cape_J_kg=float(old.to('J/kg').magnitude),previous_most_cin_J_kg=float(old_cin.to('J/kg').magnitude))
                for label,i,selection in [('SBCAPE',0,'bottom'),('previous_most_cape',0,'most_cape'),('MUCAPE',idx,'bottom')]:
                    pp,tt,dd,path=parcel_profile_with_lcl(pressure[i:],temperature[i:],td[i:])
                    lp,_=lfc(pp,tt,dd,parcel_temperature_profile=path,which=selection)
                    ep,_=el(pp,tt,dd,parcel_temperature_profile=path,which='most_cape' if selection=='most_cape' else 'top')
                    result[label+'_LFC_Pa']=float(lp.to('Pa').magnitude)
                    result[label+'_EL_Pa']=float(ep.to('Pa').magnitude)
                result['MUCAPE_top_limited']=bool(value>0 and not np.isfinite(result['MUCAPE_EL_Pa']))
    except (ValueError,IndexError,RuntimeError) as e:
        result.update(cape_J_kg=np.nan,status='invalid_calculation',reason=type(e).__name__+': '+str(e))
    return result
