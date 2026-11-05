import json
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
p=Path('build/output_mp_offline.nc')
with Dataset(p) as ds:
    a={k:np.asarray(v[:]) for k,v in ds.variables.items()}
    stats={k:{'min':float(v.min()),'max':float(v.max()),'nonfinite':int((~np.isfinite(v)).sum()),'masked':int(np.ma.getmaskarray(ds[k][:]).sum())} for k,v in a.items()}
    T=a['theta']*(a['pressure']/100000)**(287.0/1004.5)
    es=611.2*np.exp(17.67*(T-273.15)/(T-29.65))
    e=a['qv']*a['pressure']/(0.622+a['qv'])
    report={'variables':stats,'time_hours':(a['time'][[0,-1]]/3600).tolist(),'time_step_seconds':np.unique(np.diff(a['time'])).tolist(),'rain_final_mm':float(a['rainnc'][-1]),'rain_increment_sum_mm':float(a['rainncv'].sum()),'rain_monotonic':bool(np.all(np.diff(a['rainnc'])>=-1e-12)),'temperature_K':[float(T.min()),float(T.max())],'warm_RH_max_percent':float((100*e/es)[T>273.15].max()),'pressure_decreases_with_height':bool(np.all(np.diff(a['pressure'],axis=1)<0))}
Path('resultados/TH_CF_validation.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
