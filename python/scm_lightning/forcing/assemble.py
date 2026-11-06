"""Assemble checked daily cache into complete forcing inputs, preserving raw files."""
from pathlib import Path
import argparse
import json
import sys
import xarray as xr
from scm_lightning.forcing import normalize
from scm_lightning.provenance import ROOT,write_json,sha256
from scm_lightning.config import load_config
import pandas as pd

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',default=str(ROOT/'configs/cases/argentina_20181213.yaml'))
    args=ap.parse_args();c=load_config(args.config)
    cache=ROOT/c.get('paths',{}).get('era5','data/raw/era5')
    out=ROOT/'data/intermediate/era5'/c.get('event_id','');out.mkdir(parents=True,exist_ok=True)
    acquisition=c.get('acquisition',{})
    days=pd.date_range(acquisition.get('start','2018-12-12'),acquisition.get('end','2018-12-15'),freq='D')
    records=[]
    for kind in ['pl','sl']:
        paths=[cache/f'reanalysis_{kind}_{day:%Y%m%d}.nc' for day in days]
        target=out/f'reanalysis_{kind}.nc'
        if target.exists():raise FileExistsError(target)
        data=xr.concat([normalize(xr.load_dataset(p)) for p in paths],dim='time').sortby('time')
        data.to_netcdf(target)
        records.append({'output':str(target),'sha256':sha256(target),'inputs':[{'path':str(p),'sha256':sha256(p)} for p in paths]})
    write_json(out/'manifest.json',records)

if __name__=='__main__':main()
