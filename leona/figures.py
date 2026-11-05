"""Four technical figures, produced only from supplied valid case products."""
from pathlib import Path
import numpy as np

def save(fig,path):
    import matplotlib.pyplot as plt
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=180,bbox_inches='tight');plt.close(fig)

def case_figure(meteorology,box,soundings,path):
    """meteorology: regional xarray t2m/mean sea level pressure at selected time."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    if box is None: raise ValueError('Fixed meteorological area required')
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    meteorology.t2m.plot(ax=axes[0],cbar_kwargs={'label':'2-m temperature (K)'})
    w,s,e,n=box;axes[0].add_patch(Rectangle((w,s),e-w,n-s,fill=False,color='red'))
    for label,f in soundings.items():
        axes[1].plot(f.temperature_C,f['geopotential height_m']/1000,label=label)
    axes[1].set(xlabel='Temperature (C)',ylabel='Height (km)');axes[1].legend()
    save(fig,path)

def prediction_figure(heldout,path):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(12,5))
    obs=heldout.groupby('time')['count'].first();ax.plot(obs.index,obs,label='Observed',color='black')
    for (method,configuration),f in heldout.groupby(['method','configuration']):
        ensemble=f[f.member!=0] if (f.member!=0).any() else f
        g=ensemble.groupby('time').prediction_A
        ax.plot(g.mean().index,g.mean(),label=f'{method}/{configuration}')
        ax.fill_between(g.mean().index,g.quantile(.05),g.quantile(.95),alpha=.12)
    ax.set(ylabel='Flashes / 15 min',xlabel='UTC (held-out intervals)');ax.legend(fontsize=7)
    save(fig,path)

def comparison_figure(metrics,contrasts,path):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    m=metrics[metrics.comparison=='A'];labels=m.method+'/'+m.configuration
    axes[0].barh(labels,m.rmse);axes[0].set(xlabel='Held-out count RMSE')
    axes[1].barh(contrasts.method+'/'+contrasts.configuration,contrasts.rmse_change)
    axes[1].axvline(0,color='black');axes[1].set(xlabel='RMSE change with fixed reference scale')
    fig.tight_layout();save(fig,path)

def microphysics_figure(time,ice_paths,processes,path):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(12,7),sharex=True)
    for name,value in ice_paths.items():axes[0].plot(time,value,label=name)
    for name,value in processes.items():axes[1].plot(time,value,label=name)
    axes[0].set(ylabel='Species path (kg m-2)');axes[1].set(ylabel='Process transfer (kg m-2 s-1)',xlabel='UTC')
    for ax in axes:ax.legend(fontsize=7)
    fig.tight_layout();save(fig,path)
