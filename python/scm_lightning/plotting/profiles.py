"""Plot each launch without combining different stations into a synthetic column."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patheffects


def draw(axes,records,map_axes=()):
    keys=[];groups={};frames=[];metadata=[]
    styles=['-','--','-.',':',(0,(5,1,1,1)),(0,(3,1)),(0,(1,1))]
    for record in sorted(records,key=lambda r:r['launch_UTC']):
        key=(record['station'],round(record['latitude'],2),round(record['longitude'],2))
        if key not in keys:keys.append(key)
        index=keys.index(key);number=groups.get(key,0);groups[key]=number+1
        color=plt.get_cmap('tab10')(index);style=styles[number%len(styles)]
        alias=('CSU'+str(index+1)) if record['source']=='CSU_RELAMPAGO_in_projetoCHUVA' else record['station']
        label=alias+' '+pd.Timestamp(record['launch_UTC']).strftime('%m/%d %H:%M')
        full_depth=any(v.get('covers_20km_AGL',False) for v in record['cases'].values())
        if not full_depth:label+=' *'
        f=pd.read_csv(record['normalized']);f=f[f.thermodynamic_valid & f.ascending_segment].copy()
        f['launch_label']=label;frames.append(f)
        axes[0].plot(f.temperature_C,f.height_geopotential_m/1000,c=color,ls=style,lw=1.1,label=label)
        axes[1].plot(1000*f.rv_kg_kg_dry,f.height_geopotential_m/1000,c=color,ls=style,lw=1.1,label=label)
        if number==0:
            for ax in map_axes:
                ax.scatter(record['longitude'],record['latitude'],marker='^',s=28,c=[color],zorder=7)
                ax.annotate(alias,(record['longitude'],record['latitude']),xytext=(4,5),textcoords='offset points',fontsize=7,zorder=8,
                    path_effects=[patheffects.withStroke(linewidth=2,foreground='white')])
        metadata.append({'station':record['station'],'display_site':alias,'latitude':record['latitude'],'longitude':record['longitude'],
            'actual_launch_UTC':record['launch_UTC'],'source':record['file'],'source_sha256':record['sha256'],
            'label':label,'partial_for_20km_AGL_column':not full_depth,'asterisk_meaning':'incomplete depth relative to 20km AGL SCM top',
            'case_representativeness':record['cases']})
    for ax in axes:
        ax.set_ylabel('Geopotential height (km MSL)');ax.grid(alpha=.2)
    axes[0].set_xlabel('Temperature (°C)');axes[1].set_xlabel('Mixing ratio (g kg$^{-1}$)')
    axes[0].legend(loc='upper right',fontsize=6.5,frameon=False,ncol=2 if len(records)>5 else 1)
    return pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(),metadata
