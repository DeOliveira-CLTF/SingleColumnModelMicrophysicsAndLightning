"""Official CPTEC MERGE hourly archive, raw GRIB and descriptor provenance."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
import requests
from scm_lightning.provenance import ROOT,sha256,write_json
from scm_lightning.config import load_config

BASE='https://ftp.cptec.inpe.br/modelos/tempo/MERGE/GPM/HOURLY/'
CACHE=ROOT/'data/raw/precipitation/merge'
STAC='https://data.inpe.br/bdc/stac/v1/collections/prec_merge_hourly-1'


def acquire(event):
    c=load_config(ROOT/'configs/cases'/(event+'.yaml'));folder=CACHE/event;folder.mkdir(parents=True,exist_ok=True)
    # File validity labels the END of the preceding one-hour accumulation.
    times=pd.date_range(c['start'],c['end'],freq='h',inclusive='right')
    def get(stamp):
        stem=f'MERGE_CPTEC_{stamp:%Y%m%d%H}';base=BASE+f'{stamp:%Y/%m/%d}/'+stem
        entries=[]
        for extension in ['.grib2','.ctl']:
            target=folder/(stem+extension);url=base+extension
            if not target.exists():
                r=requests.get(url,timeout=60);r.raise_for_status();body=r.content
                if extension=='.grib2' and body[:4]!=b'GRIB':raise ValueError('Not a GRIB payload: '+url)
                temp=target.with_suffix(extension+'.pending');temp.write_bytes(body);temp.replace(target)
            entries.append({'path':str(target),'url':url,'sha256':sha256(target),'bytes':target.stat().st_size})
        return {'valid_time_UTC':str(stamp),'interval_start_UTC':str(stamp-pd.Timedelta(hours=1)),'sources':entries}
    with ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(get,times))
    documentation=CACHE/'prec_merge_hourly-1.json'
    if not documentation.exists():
        r=requests.get(STAC,timeout=30);r.raise_for_status();documentation.write_text(r.text,encoding='utf-8')
    write_json(folder/'download_manifest.json',{'event':event,'files':records,'area_wsen':c['area'],
        'product':'CPTEC/INPE MERGE GPM HOURLY','STAC_collection_version':'1',
        'archive_algorithm_version':'not encoded in legacy GRIB; do not assert current V07B for historical files',
        'documentation':{'url':STAC,'sha256':sha256(documentation)},'time_convention':'one-hour accumulation ending at filename/GRIB valid_time',
        'resolution_degrees':.1,'interval_count':len(records),'independent_of_SCM_ERA5_forcing':True})
    print(f'MERGE {event}: {len(records)} actual hourly files.',flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--event',choices=['argentina_20181213','argentina_20191028','both'],default='both');a=p.parse_args()
    for event in ['argentina_20181213','argentina_20191028'] if a.event=='both' else [a.event]:acquire(event)


if __name__=='__main__':main()
