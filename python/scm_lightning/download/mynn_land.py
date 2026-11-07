"""Official ERA5 soil/vegetation context; never equates volumetric water to mavail."""
import argparse
import json
from pathlib import Path
import zipfile
import pandas as pd
import xarray as xr
from scm_lightning.provenance import ROOT,sha256,write_json

VARIABLES=['volumetric_soil_water_layer_1','volumetric_soil_water_layer_2',
    'volumetric_soil_water_layer_3','volumetric_soil_water_layer_4',
    'high_vegetation_cover','low_vegetation_cover','soil_type',
    'type_of_high_vegetation','type_of_low_vegetation','land_sea_mask',
    'leaf_area_index_high_vegetation','leaf_area_index_low_vegetation',
    'total_precipitation','surface_latent_heat_flux','surface_sensible_heat_flux']


def acquire():
    import cdsapi
    client=cdsapi.Client(quiet=True,progress=False)
    outputs=[]
    for event,year,month,first,last in [('argentina_20181213','2018','12',6,15),('argentina_20191028','2019','10',21,30)]:
        from scm_lightning.config import load_config
        c=load_config(ROOT/f'configs/cases/{event}.yaml');w,s,e,n=c['area']
        folder=ROOT/'data/raw/era5/mynn_land'/event;folder.mkdir(parents=True,exist_ok=True)
        target=folder/'land_context.download';meta=folder/'source.json'
        request={'product_type':['reanalysis'],'variable':VARIABLES,'year':[year],'month':[month],
            'day':[f'{i:02}' for i in range(first,last+1)],'time':[f'{h:02}:00' for h in range(24)],
            'area':[n+.5,w-.5,s-.5,e+.5],'data_format':'netcdf','download_format':'unarchived'}
        write_json(folder/'request.json',{'dataset':'reanalysis-era5-single-levels','request':request,'SCM_box_wsen':c['area']})
        if meta.exists():
            record=json.loads(meta.read_text());assert sha256(target)==record['download_sha256']
            outputs.append(record);continue
        pending=folder/'land_context.pending'
        if not target.exists():
            client.retrieve('reanalysis-era5-single-levels',request,str(pending));pending.replace(target)
        files=[]
        if zipfile.is_zipfile(target):
            with zipfile.ZipFile(target) as archive:
                for member in archive.namelist():
                    if not member.endswith('.nc'):continue
                    path=folder/Path(member).name
                    if not path.exists():path.write_bytes(archive.read(member))
                    files.append(path)
        else:
            path=folder/'land_context.nc'
            if not path.exists():path.write_bytes(target.read_bytes())
            files=[path]
        variables={}
        for p in files:
            with xr.open_dataset(p) as d:
                for name in d.data_vars:variables[name]={'units':d[name].attrs.get('units'),'shape':list(d[name].shape)}
        record={'source':'ERA5 CDS official','dataset':'reanalysis-era5-single-levels','request':request,
            'download_sha256':sha256(target),'files':[{'path':str(p),'sha256':sha256(p)} for p in files],
            'variables':variables,'SCM_box_wsen':c['area'],'mavail_from_volumetric_water':False,
            'flux_convention':'downward positive; accumulated J m-2, ERA5 hourly interval ending at valid_time',
            'reference':'https://confluence.ecmwf.int/pages/viewpage.action?pageId=76414402'}
        write_json(meta,record);outputs.append(record);print('ERA5 land acquired: '+event,flush=True)
    return outputs


if __name__=='__main__':acquire()
