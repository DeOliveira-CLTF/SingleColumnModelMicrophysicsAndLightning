from scm_lightning.config import load_config
"""Local availability tables and sounding QA figure; no result interpretation."""
from pathlib import Path
import importlib.metadata
import json
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scm_lightning.provenance import ROOT,write_json,sha256,git

def main():
    out=ROOT/'results/argentina_20181213';(out/'tables').mkdir(parents=True,exist_ok=True)
    quality=json.loads((ROOT/'data/processed/profiles/quality.json').read_text())
    pd.DataFrame(quality).to_csv(out/'tables/sounding_quality.csv',index=False)
    reg=load_config(ROOT/'configs/lightning/schemes.yaml')
    pd.DataFrame(reg['methods']).to_csv(out/'tables/parameterization_requirements.csv',index=False)
    availability=out/'soundings/availability.csv'
    if availability.exists():
        f=pd.read_csv(availability)
        f.to_csv(out/'tables/data_availability_wyoming.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(13,6))
    for p in sorted((ROOT/'data/processed/profiles').glob('*.csv')):
        f=pd.read_csv(p);name=p.stem
        label=name if 'time' not in f else f'{name}: {f.time.iloc[0]} to {f.time.iloc[-1]}'
        axes[0].plot(f.temperature_C,f['geopotential height_m']/1000,label=label)
        axes[1].plot(f['relative humidity_%'],f['geopotential height_m']/1000,label=name)
        if {'latitude','longitude'}<=set(f):axes[2].plot(f.longitude,f.latitude,label=name)
    axes[0].set(xlabel='Temperature (C)',ylabel='Geopotential height (km)')
    axes[1].set(xlabel='Relative humidity (%)',ylabel='Geopotential height (km)')
    axes[2].set(xlabel='Longitude',ylabel='Latitude')
    for ax in axes:ax.grid(alpha=.2)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.01),fontsize=7)
    fig.tight_layout(rect=[0,.16,1,1]);(out/'supplement').mkdir(exist_ok=True)
    fig.savefig(out/'supplement/S01_sounding_QA.png',dpi=180,bbox_inches='tight');plt.close(fig)
    names=['numpy','pandas','xarray','netCDF4','scipy','matplotlib','metpy','requests','cdsapi','PyYAML','pytest']
    versions={n:importlib.metadata.version(n) for n in names}
    write_json(out/'tests/environment.json',{'python':sys.version,'packages':versions,
        'head':git('rev-parse','HEAD'),'dephy_head':git('rev-parse','HEAD',cwd=ROOT/'external/DEPHY-SCM')})
    (out/'tests/requirements-executed.txt').write_text('\n'.join(n+'=='+versions[n] for n in names)+'\n')
    files=[p for p in out.rglob('*') if p.is_file() and p.name!='product_manifest.json']
    archive=out/'native_product_index.json'
    known={r['product']:r['sha256'] for r in json.loads(archive.read_text())['files']} if archive.exists() else {}
    write_json(out/'product_manifest.json',[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,
        'sha256':known.get(str(p.relative_to(ROOT))) or sha256(p),'redistribution':'local_only'} for p in files])

if __name__=='__main__':main()
