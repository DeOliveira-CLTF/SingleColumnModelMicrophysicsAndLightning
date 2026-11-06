"""Assemble checked daily cache into complete forcing inputs, preserving raw files."""
from pathlib import Path
import json
import sys
import xarray as xr
from scm_lightning.forcing import normalize
from scm_lightning.provenance import ROOT,write_json,sha256

def main():
    cache=ROOT/'data/raw/era5'
    out=ROOT/'data/intermediate/era5';out.mkdir(parents=True,exist_ok=True)
    records=[]
    for kind in ['pl','sl']:
        paths=[cache/f'reanalysis_{kind}_201812{day}.nc' for day in ['12','13','14','15']]
        target=out/f'reanalysis_{kind}.nc'
        if target.exists():raise FileExistsError(target)
        data=xr.concat([normalize(xr.load_dataset(p)) for p in paths],dim='time').sortby('time')
        data.to_netcdf(target)
        records.append({'output':str(target),'sha256':sha256(target),'inputs':[{'path':str(p),'sha256':sha256(p)} for p in paths]})
    write_json(out/'manifest.json',records)

if __name__=='__main__':main()
