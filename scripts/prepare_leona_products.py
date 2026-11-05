"""Local availability tables and sounding QA figure; no result interpretation."""
from pathlib import Path
import importlib.metadata
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from leona.provenance import ROOT,write_json,sha256,git

def main():
    out=ROOT/'products/leona';(out/'tables').mkdir(parents=True,exist_ok=True)
    quality=json.loads((ROOT/'data_private/leona/processed/soundings/quality.json').read_text())
    pd.DataFrame(quality).to_csv(out/'tables/sounding_quality.csv',index=False)
    reg=json.loads((ROOT/'configs/parameterizations.json').read_text())
    pd.DataFrame(reg['methods']).to_csv(out/'tables/parameterization_requirements.csv',index=False)
    availability=out/'soundings/availability.csv'
    if availability.exists():
        f=pd.read_csv(availability)
        f.to_csv(out/'tables/data_availability_wyoming.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(13,6))
    for p in sorted((ROOT/'data_private/leona/processed/soundings').glob('*.csv')):
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
    names=['numpy','pandas','xarray','netCDF4','scipy','matplotlib','requests','cdsapi','pytest']
    write_json(out/'tests/environment.json',{'python':sys.version,'packages':{n:importlib.metadata.version(n) for n in names},
        'head':git('rev-parse','HEAD'),'dephy_head':git('rev-parse','HEAD',cwd=ROOT/'era5_forcing/DEPHY-SCM')})
    files=[p for p in out.rglob('*') if p.is_file() and p.name!='product_manifest.json']
    write_json(out/'product_manifest.json',[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha256(p),'redistribution':'local_only'} for p in files])

if __name__=='__main__':main()
