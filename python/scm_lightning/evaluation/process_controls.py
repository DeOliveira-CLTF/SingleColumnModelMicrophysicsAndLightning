from scm_lightning.ensemble.experiments import control_directory
"""Control microphysics paths and applied process increments (technical products)."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import xarray as xr
from scm_lightning.ensemble.experiments import CONFIGURATIONS
from scm_lightning.io import numeric_frame
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.evaluation.process_categories import CATEGORY

def main():
    out=ROOT/'results/argentina_20181213/control_diagnostics';out.mkdir(parents=True,exist_ok=True)
    for configuration in CONFIGURATIONS:
        folder=control_directory(configuration)
        ds=xr.load_dataset(folder/'output.nc',decode_times=False)
        z=ds.lev.values;zw=np.zeros(len(z)+1)
        for k,center in enumerate(z):zw[k+1]=2*center-zw[k]
        dz=np.diff(zw)
        t=ds.theta.values*(ds.pressure.values/1e5)**(2/7)
        rho_d=ds.pressure.values/(287*t*(1+ds.qv.values/(287/461.6)))
        initial_mass=rho_d[0]*dz
        paths=pd.DataFrame({'time_s':ds.time.values,'configuration':configuration})
        for q in ['qc','qr','qi','qs','qg']:
            paths[q]=np.sum(ds[q].values*rho_d*dz,axis=1)
        paths.to_csv(out/(configuration+'_paths.csv'),index=False)
        frames=[]
        cols=['time_s','k','dqv','dqc','dqr','dqi','dqs','dqg']
        source=folder/'output.nc.processes.csv'
        for raw in pd.read_csv(source,chunksize=150000):
            raw=numeric_frame(raw,cols)
            raw['category']=raw.process.map(CATEGORY).fillna('Unmapped process')
            mass=initial_mass[raw.k.astype(int).values-1]
            for q in ['qv','qc','qr','qi','qs','qg']:raw['d'+q]=raw['d'+q]*mass/30.
            frames.append(raw.groupby(['time_s','category'])[['d'+q for q in ['qv','qc','qr','qi','qs','qg']]].sum())
        rates=pd.concat(frames).groupby(level=[0,1]).sum().reset_index()
        rates['configuration']=configuration
        rates.to_csv(out/(configuration+'_processes.csv.gz'),index=False)
        write_json(out/(configuration+'_manifest.json'),{'simulation_sha256':sha256(folder/'output.nc'),
            'process_csv_sha256':sha256(source),'paths_units':'kg m-2 using instantaneous dry-air density',
            'process_rate_units':'kg m-2 s-1 per INITIAL dry-air cell mass; not instantaneous conservative mass flux',
            'density_contract':'MPAS Rd=287, Rv=461.6, kappa=2/7; fixed dz inferred from layer centers',
            'dt_s':30,'unmapped_processes':sorted(set(rates.category)-set(CATEGORY.values()))})
        print(configuration,flush=True)

if __name__=='__main__':main()
