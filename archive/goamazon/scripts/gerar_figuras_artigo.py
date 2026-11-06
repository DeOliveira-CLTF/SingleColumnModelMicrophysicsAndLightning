"""Publication figures, identical species scales between schemes; no fabricated data."""
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm, TwoSlopeNorm
import pandas as pd
from scm_artigo import *
from validar_resultados import validate_file

COLORS={'TH':'#0072B2','WS':'#D55E00'}
STYLES={'C0':':','CL':'--','CF':'-'}
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,'legend.fontsize':9,
                     'savefig.dpi':300,'pdf.fonttype':42,'font.family':'DejaVu Sans'})

def save(fig,name,out,tight=True):
    bbox='tight' if tight else None
    fig.savefig(out/(name+'.png'),dpi=300,bbox_inches=bbox)
    fig.savefig(out/(name+'.pdf'),bbox_inches=bbox)
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--results',type=Path,default=RESULTS); ap.add_argument('--output',type=Path,default=ROOT/'figuras_artigo'); args=ap.parse_args()
    out=args.output; out.mkdir(exist_ok=True,parents=True)
    records={}; ref=parse_nml(nml_path('TH_CF',args.results))
    for eid in IDS:
        path=args.results/f'{eid}.nc'
        if not path.exists(): continue
        cfg=parse_nml(nml_path(eid,args.results)); qa=validate_file(path,cfg,eid,ref)
        if not qa['valid']: raise ValueError(f'{eid}: {qa["errors"]}')
        ds=load(path); m,s=metrics(ds,cfg); records[eid]=(ds,m,s)
    forcing=load(INPUT); ds=next(iter(records.values()))[0]
    fig,ax=plt.subplots(2,3,figsize=(14,8))
    z=ds.lev.values/1000; T=temperature(ds)
    ax[0,0].plot(T[0]-273.15,z,color='#333333'); ax[0,0].set(xlabel='Initial temperature (°C)',ylabel='Height (km)',title='A  Thermal profile')
    ax[0,0].axhspan(16.5,20,color='grey',alpha=.15,label='Artificial upper layer'); ax[0,0].legend(loc='lower left',fontsize=8)
    ax[0,1].plot(ds.qv[0]*1000,z,color='#009E73'); ax[0,1].set(xlabel='Water vapor mixing ratio (g/kg)',ylabel='Height (km)',title='B  Initial humidity')
    h=forcing.time.values/3600
    ax[0,2].plot(h,forcing.ts_forc.values-273.15,label='Skin temperature',color='#CC79A7'); ax[0,2].set(xlabel='Time since t₀ (h)',ylabel='Temperature (°C)',title='C  Surface condition')
    ax[1,0].plot(h,forcing.hfss.values,label='Sensible',color='#E69F00'); ax[1,0].plot(h,forcing.hfls.values,label='Latent',color='#009E73'); ax[1,0].legend()
    ax[1,0].set(xlabel='Time since t₀ (h)',ylabel='Flux in input file (W/m²)',title='D  Available surface fluxes')
    for axis,var,title,scale,unit in [(ax[1,1],'wa','E  Prescribed vertical velocity',1,'m/s'),(ax[1,2],'tnta_adv','F  Prescribed advective tendency',3600,'K/h')]:
        a=forcing[var].values*scale; lim=float(np.max(np.abs(a)))
        im=axis.pcolormesh(h,forcing.lev.values/1000,a.T,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto',rasterized=True)
        colorbar=fig.colorbar(im,ax=axis,label=unit,pad=.08)
        colorbar.ax.tick_params(labelsize=8)
        axis.set(xlabel='Time since t₀ (h)',ylabel='Height (km)',title=title)
        axis.yaxis.labelpad=12
    for axis in ax.flat: axis.grid(alpha=.15)
    forcing_duration_h=(forcing.time.values[-1]-forcing.time.values[0])/3600
    fig.suptitle('GOAMAZON/SINGLE — Initial atmospheric state and prescribed forcing',fontsize=13,y=.98)
    fig.text(.5,.945,f'5–6 October 2014 | 8°S, 63°W | {forcing_duration_h:.2f} h of forcing data',
             ha='center',va='top',fontsize=11)
    fig.subplots_adjust(top=.84,wspace=.62,hspace=.3)
    save(fig,'Figura1_caracterizacao',out)

    if all(e in records for e in ['TH_CF','WS_CF']):
        fig,axes=plt.subplots(4,2,figsize=(16,18),layout='constrained',sharex=True,sharey=True)
        for row,q in enumerate(['qc','qi','qs','qg']):
            vmax=max(float(records[e][0][q].max())*1000 for e in ['TH_CF','WS_CF'])
            for col,eid in enumerate(['TH_CF','WS_CF']):
                d=records[eid][0]; axis=axes[row,col]; a=d[q].values*1000; h=d.time.values/3600
                if a.max()==0:
                    axis.set_facecolor('#f3f3f3'); axis.text(.5,.72,'Zero field throughout integration',ha='center',transform=axis.transAxes,fontsize=22)
                else:
                    im=axis.pcolormesh(h,d.lev.values/1000,a.T,shading='auto',cmap='viridis',norm=PowerNorm(.5,vmin=0,vmax=vmax),rasterized=True)
                    colorbar=fig.colorbar(im,ax=axis,label=f'{q} (g/kg)',shrink=.85)
                    colorbar.ax.tick_params(labelsize=18)
                    colorbar.ax.yaxis.label.set_size(20)
                cs=axis.contour(h,d.lev.values/1000,(temperature(d)-273.15).T,levels=[-40,-15,0],colors='#666666' if a.max()==0 else '#eeeeee',linewidths=.7)
                axis.clabel(cs,fmt='%d °C',fontsize=16)
                if q=='qi' and eid=='TH_CF': axis.text(.02,.93,f'Maximum: {a.max():.3g} g/kg',color='white',transform=axis.transAxes,fontsize=18)
                axis.set(title=f'{eid.replace("_","-")} · {q}',ylabel='Height (km)',ylim=(0,20),xlim=(0,h[-1]))
                axis.set_title(f'{eid.replace("_","-")} · {q}',fontsize=22)
                axis.set_ylabel('Height (km)',fontsize=20)
                if row==3: axis.set_xlabel('Time since t₀ (h)',fontsize=20)
                axis.tick_params(axis='both',labelsize=18)
        fig.suptitle('Hydrometeor partitioning under full coupling\nSame scale between schemes for each species; square-root normalization',fontsize=28)
        save(fig,'Figura2_hidrometeoros',out,tight=False)

    def temporal(eids,name):
        fig,axes=plt.subplots(5,1,figsize=(18,21),layout='constrained',sharex=True)
        for eid in eids:
            if eid not in records: continue
            d,m,s=records[eid]; h=s['time_h']; color=COLORS[eid[:2]]; ls=STYLES[eid[3:]]
            for axis,key in zip([axes[0],axes[2],axes[3],axes[4]],['qg_mixed_kg_m2','flash_min','rain_rate_mm_h','rain_mm']):
                axis.plot(h,s[key],label=eid.replace('_','-'),color=color,ls=ls,lw=1.3)
            if len(eids)==2:
                for q,qstyle in [('qc','-'),('qs','--'),('qi',':')]:
                    axes[1].plot(h,s[f'{q}_mixed_kg_m2'],color=color,ls=qstyle,label=f'{eid}: {q}',lw=1)
                for axis,key in [(axes[0],'qg_mixed_mass_peak_h'),(axes[2],'flash_peak_h'),(axes[3],'rain_peak_h')]:
                    if m[key] is not None: axis.axvline(m[key],color=color,ls=':',alpha=.6); axis.text(m[key],.96 if eid.startswith('TH') else .83,f'{m[key]:.2f} h',color=color,transform=axis.get_xaxis_transform(),fontsize=16,ha='right')
            else:
                axes[1].plot(h,s['qs_mixed_kg_m2'],color=color,ls=ls,label=eid.replace('_','-'))
        labels=['Mixed-phase graupel\n(kg/m²)','Other mixed-phase species\n(kg/m²)' if len(eids)==2 else 'Mixed-phase snow\n(kg/m²)','Nominal electrical rate\n(flashes/min)','Precipitation rate\n(mm/h)','Accumulated precipitation\n(mm)']
        for i,(axis,label) in enumerate(zip(axes,labels)):
            axis.set_ylabel(label,fontsize=20); axis.set_ylim(bottom=0); axis.grid(alpha=.2); axis.legend(ncol=3,loc='upper left',fontsize=18); axis.set_title('ABCDE'[i],loc='left',fontweight='bold',fontsize=22); axis.tick_params(axis='both',labelsize=18)
        axes[-1].set_xlabel('Time since t₀ (h)',fontsize=20)
        fig.suptitle('Temporal evolution · uncalibrated lightning diagnostic',fontsize=28)
        save(fig,name,out)
    if all(e in records for e in ['TH_CF','WS_CF']): temporal(['TH_CF','WS_CF'],'Figura3_evolucao_temporal')
    temporal(IDS,'Figura3_suplementar_seis_configuracoes')

    fig,axes=plt.subplots(2,2,figsize=(18,14),layout='constrained')
    for axis,key,label,letter in zip(axes.flat,['rain_final_mm','qg_mixed_mass_max_kg_m2','flash_max_min','flash_integral_nominal'],['Accumulated precipitation (mm)','Maximum mixed-phase graupel mass (kg/m²)','Maximum electrical rate (flashes/min)','Nominal electrical integral (flashes)'],'ABCD'):
        for i,scheme in enumerate(['TH','WS']):
            for j,c in enumerate(['C0','CL','CF']):
                eid=f'{scheme}_{c}'; x=j+(i-.5)*.36
                if eid not in records:
                    axis.text(x,0,'Missing',rotation=90,va='bottom',ha='center',fontsize=18); continue
                v=records[eid][1][key]
                axis.bar(x,v,width=.34,color=COLORS[scheme],hatch='' if scheme=='TH' else '//',label={'TH':'Thompson','WS':'WSM6'}[scheme] if j==0 else None)
                axis.annotate(f'{v:.4g}',(x,v),xytext=(0,6),textcoords='offset points',ha='center',fontsize=18)
        axis.set_xticks(range(3),['C0','CL','CF']); axis.set_ylabel(label,fontsize=20); axis.set_title(letter,loc='left',fontweight='bold',fontsize=22); axis.set_ylim(0,max([records[e][1][key] for e in records]+[1e-10])*1.25); axis.grid(axis='y',alpha=.2); axis.legend(loc='upper left',ncol=2,fontsize=18); axis.tick_params(axis='both',labelsize=18)
    fig.suptitle('Sensitivity to microphysics and MYNN condensate delivery\nC0: not applied · CL: local only · CF: local + plumes',fontsize=28)
    save(fig,'Figura4_comparacao_experimentos',out)

if __name__=='__main__': main()
