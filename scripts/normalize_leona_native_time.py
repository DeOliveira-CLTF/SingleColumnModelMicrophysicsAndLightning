"""CF-valid UTC copies of native outputs; preserve original files and data arrays."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import json
import shutil
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import netCDF4 as nc
import numpy as np
from leona.provenance import ROOT,sha256,write_json

def main():
    out=ROOT/'products/leona';c=json.loads((ROOT/'configs/leona.json').read_text())
    index_path=out/'native_product_index.json';index=json.loads(index_path.read_text())
    origin=datetime.fromisoformat(c['start'].replace('Z','+00:00'))
    units=origin.strftime('seconds since %Y-%m-%d %H:%M:%S +00:00')
    existing={r['product']:r for r in index['files']};records=[];compared=set()
    candidates=[r for r in index['files'] if Path(r['product']).name=='output.nc']
    for record in candidates:
        source=ROOT/record['product'];target=source.with_name('output_utc.nc')
        configuration=target.parent.name if target.parent.name in ['TH_C0','TH_CL','TH_CF','WS_C0','WS_CL','WS_CF'] else target.parent.parent.name
        changed=False
        if not target.exists():
            shutil.copy2(source,target) # deliberately independent, never modify a hardlink
            changed=True
        if source.stat().st_ino==target.stat().st_ino:
            raise ValueError('UTC metadata must never be edited through a native hardlink')
        with nc.Dataset(target,'a') as ds:
            for key,value in {'units':units,'calendar':'gregorian'}.items():
                if getattr(ds['time'],key,None)!=value:ds['time'].setncattr(key,value);changed=True
            attrs={'case_id':c['case'],'configuration':configuration,
                'source':'mp_offline_driver; '+('Thompson' if configuration.startswith('TH') else 'WSM6')+' + original MYNN; '+c['case'],
                'clock_source':'UTC origin from verified DEPHY forcing manifest','native_source_sha256':record['sha256'],
                'lightning_diagnostics':'postprocessing only; native flash_rate fields disabled'}
            for key,value in attrs.items():
                if getattr(ds,key,None)!=value:ds.setncattr(key,value);changed=True
        with nc.Dataset(target) as ds:
            assert ds['time'].units==units
            ends=nc.num2date(ds['time'][[0,-1]],units,calendar='gregorian',only_use_cftime_datetimes=False)
            assert str(ends[0])=='2018-12-13 00:00:00' and str(ends[-1])=='2018-12-15 00:00:00'
        # Check every variable in one output from each scheme/coupling; metadata
        # edits for all remaining files use the same NetCDF operation.
        if configuration not in compared:
            with nc.Dataset(source) as old,nc.Dataset(target) as new:
                assert set(old.variables)==set(new.variables)
                for name in old.variables:np.testing.assert_array_equal(old[name][:],new[name][:])
            compared.add(configuration)
        product=str(target.relative_to(ROOT))
        if product not in existing:
            records.append({'source':record['product'],'product':product,'bytes':target.stat().st_size,
                'method':'independent_copy_metadata_only_UTC','source_sha256':record['sha256'],'sha256':sha256(target)})
        elif changed:
            existing[product].update(bytes=target.stat().st_size,sha256=sha256(target))
    index['files'].extend(records);index['recommended_clock_files']='output_utc.nc; raw output.nc files remain unchanged'
    write_json(index_path,index)
    write_json(out/'tests/native_utc_clock.json',{'status':'passed','files':len(candidates),'origin_utc':c['start'],
        'end_utc':c['end'],'all_physical_arrays_equal_in_configurations':sorted(compared),
        'original_native_files_unchanged':True,'new_files':'output_utc.nc'})

if __name__=='__main__':main()
