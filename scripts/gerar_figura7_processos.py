"""Figure 7: applied microphysical sources of total ice and graupel.

Static publication figure; all 12 members and individual rates remain available.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scm_artigo import *
from executar_processos_figura7 import OUT,CATEGORY,QS

PLOT_GROUPS=['Deposição/sublimação','Nucleação de gelo','Congelamento','Riming (água de nuvem)',
             'Coleta chuva–gelo','Coleta gelo–gelo','Autoconversão/agregação','Gelo secundário']
LABELS=['Deposition','Ice nucleation','Freezing','Riming: cloud water',
        'Rain–ice collection','Ice–ice collection','Conversion / aggregation','Secondary ice']


def table_md(df):
    def f(v): return f'{v:.6g}' if isinstance(v,(float,np.floating)) else str(v)
    return '\n'.join(['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']+
                     ['| '+' | '.join(f(v) for v in row)+' |' for row in df.itertuples(index=False,name=None)])


def main():
    manifest=json.loads((OUT/'latest.json').read_text(encoding='utf-8')); root=ROOT/manifest['directory']
    rows=[]; balances=[]
    for eid,rec in manifest['runs'].items():
        folder=ROOT/rec['directory']; rates=load(folder/'process_tendencies.nc'); ds=load(folder/'profile_1.nc')
        assert sha(folder/'process_tendencies.nc')==rec['validation']['diagnostics_sha256']
        names=[x.decode() if isinstance(x,bytes) else str(x) for x in rates.process_name.values]
        mass=rates.reference_dry_mass.values
        vals=rates.delta_q.values*mass[None,:,None,None]
        for name,j in zip(names,range(len(names))):
            for target,data in [('gelo_total',vals[:,:,j,3:].sum(axis=2)),('graupel',vals[:,:,j,5]),('agua_nuvem',vals[:,:,j,1])]:
                rows.append(dict(experiment=eid,target=target,process=name,category=CATEGORY[name],
                                 production_kg_m2=float(np.maximum(data,0).sum()),
                                 loss_kg_m2=float(-np.minimum(data,0).sum()),net_kg_m2=float(data.sum())))
        for j,q in enumerate(QS[1:],1):
            initial=float(np.dot(ds[q].values[0],mass)); final=float(np.dot(ds[q].values[-1],mass))
            tendency=float(vals[:,:,:,j].sum())
            assert abs(final-initial-tendency)<1e-10
            balances.append(dict(experiment=eid,species=q,initial_kg_m2=initial,final_kg_m2=final,
                                 applied_net_kg_m2=tendency,closure_kg_m2=final-initial-tendency))
    detail=pd.DataFrame(rows); detail.to_csv(root/'figure7_individual_processes.csv',index=False,float_format='%.17g')
    groups=detail.groupby(['experiment','target','category'],as_index=False)[['production_kg_m2','loss_kg_m2','net_kg_m2']].sum()
    groups.to_csv(root/'figure7_grouped_processes.csv',index=False,float_format='%.17g')
    pd.DataFrame(balances).to_csv(root/'species_balances.csv',index=False,float_format='%.17g')
    # Gross generation excludes transport between cells and numerical corrections.
    physical=groups[groups.category.isin(PLOT_GROUPS)]
    graupel=physical[physical.target=='graupel'].copy()
    graupel['share_percent']=graupel.production_kg_m2/graupel.groupby('experiment').production_kg_m2.transform('sum').replace(0,np.nan)*100
    graupel.to_csv(root/'graupel_source_shares.csv',index=False,float_format='%.17g')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':20,'pdf.fonttype':42,'axes.titlesize':22,'axes.labelsize':20,'xtick.labelsize':18,'ytick.labelsize':18,'legend.fontsize':18,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(22,16))
    colors=['#0072B2','#D55E00']; letters=['a','b','c','d']
    plotted=[]
    for row,layer in enumerate(['media','alta']):
        for col,target in enumerate(['gelo_total','graupel']):
            ax=axes[row,col]; ymax=0.
            for si,scheme in enumerate(['thompson','wsm6']):
                eid=f'{scheme}_{layer}_w3'
                series=groups[(groups.experiment==eid)&(groups.target==target)].set_index('category').production_kg_m2
                values=series.reindex(PLOT_GROUPS,fill_value=0).values
                ypos=np.arange(len(values))+(si-.5)*.34
                ax.barh(ypos,values,height=.31,color=colors[si],label=['Thompson','WSM6'][si])
                ymax=max(ymax,float(values.max()))
                for yy,v,cat in zip(ypos,values,PLOT_GROUPS):
                    plotted.append(dict(panel=letters[row*2+col],experiment=eid,target=target,category=cat,production_kg_m2=v))
                    if v>1e-8:
                        label=f'{v:.3f}' if v>=.001 else f'{v:.1e}'
                        ax.annotate(label,(v,yy),xytext=(8,0),textcoords='offset points',va='center',fontsize=16,color=colors[si])
            ax.set_yticks(np.arange(len(LABELS)),LABELS,fontsize=18)
            ax.invert_yaxis(); ax.set_xlim(0,max(ymax*1.33,.01))
            ax.grid(axis='x',alpha=.18); ax.set_axisbelow(True)
            ax.set_xlabel('Cumulative gross production in column (kg m$^{-2}$)',fontsize=20)
            region='Mid-level · 3.5–5.5 km' if row==0 else 'Upper level · 7.5–9.5 km'
            species='total ice ($q_i+q_s+q_g$)' if col==0 else 'graupel ($q_g$)'
            ax.set_title(f'({letters[row*2+col]}) {region}\nSources of {species}',loc='left',fontsize=22,pad=16)
            if row==0: ax.text(.98,.02,'Thompson: all sources = 0',transform=ax.transAxes,ha='right',fontsize=18,color=colors[0])
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.935),ncol=2,frameon=False,fontsize=22)
    fig.suptitle('Microphysical contributions to ice and graupel production',fontsize=32,y=.985)
    fig.text(.5,.93,'Same members as Figure 6 · idealized GOAMAZON · 3 m/s pulse · 60 levels · Δt = 20 s · 2 h',ha='center',fontsize=20)
    fig.subplots_adjust(left=.19,right=.98,top=.84,bottom=.10,wspace=.9,hspace=.5)
    out=ROOT/'figuras_artigo'; out.mkdir(exist_ok=True)
    for ext in ['png','pdf']:
        fig.savefig(out/f'Figura7_processos_microfisicos.{ext}',dpi=240,
                    facecolor='white',bbox_inches='tight',pad_inches=.05)
    plt.close(fig)
    pd.DataFrame(plotted).to_csv(root/'figure7_plotted_values.csv',index=False,float_format='%.17g')
    # Explicitly verify the zero-ice cases stated on the figure.
    zeros=physical[((physical.experiment.str.contains('_baixa_'))|(physical.experiment.str.endswith('_w0')))&(physical.target=='graupel')]
    assert zeros.production_kg_m2.max()<1e-12
    top=graupel[graupel.share_percent>0.001].sort_values(['experiment','share_percent'],ascending=[True,False])
    report='''# Figura 7 — tendências microfísicas aplicadas

## Escopo

Foram reproduzidos os 12 membros da Figura 6, usando exclusivamente GOAMAZON original. A Figura 7 destaca as camadas média e alta com pulso de 3 m/s, nas quais aparecem as diferenças de gelo e eletrificação. Os demais membros estão nos arquivos de diagnóstico e não produzem graupel. Há, porém, nucleação de gelo no WSM6 baixo e nos controles WSM6 sem pulso; nos controles médio e alto há também deposição. Portanto, ausência de graupel/EP não significa ausência de processos de gelo. Esta comparação preserva a perturbação histórica da Figura 6; não mistura o ensemble posterior de aporte controlado.

## Método e validação

As inserções dentro de Thompson e WSM6 são delimitadas e somente observacionais. Remover os blocos de diagnóstico recupera exatamente o texto anterior de cada parametrização, cujos arquivos e hashes foram preservados em `proveniencia`. Todas as variáveis dos 12 perfis são idênticas bit a bit às saídas originais da Figura 6. Testes adicionais com o diagnóstico desligado também passaram para ambos os esquemas.

O registro é feito após a limitação das taxas. As tendências volumétricas do Thompson são divididas pela densidade usada naquele ponto e multiplicadas pelo passo efetivo; condensação e evaporação já expressas como razão de mistura não recebem outra divisão. No WSM6, as taxas finais são multiplicadas por dtcld em cada subpasso. Congelamentos/fusões diretos registram o incremento realmente aplicado. Sedimentação é obtida pela diferença de estado/tendência antes e depois de sua rotina, incluindo subpassos. Clipping, limiares e substituições de tendências ficam separados dos processos físicos.

Os arquivos NetCDF contêm Δq(time,level,process,species), em kg/kg por passo, para vapor, água de nuvem, chuva, gelo de nuvem, neve e graupel; dividir pelo passo fornece kg/kg/s. As tabelas incluem produção, perda e saldo por processo. A figura integra a parte positiva de cada processo em tempo e altura com a massa seca inicial de referência. Para gelo total, soma-se primeiro Δqi+Δqs+Δqg por processo/célula/passo, evitando contar a simples conversão de gelo em neve como geração de gelo. Para graupel, a conversão de neve em graupel é uma fonte legítima de qg.

Sedimentação é redistribuição/saída e não é incluída na produção microfísica bruta. A ponderação fixa pode produzir saldo de sedimentação que não coincide com a precipitação; o orçamento de densidade já documentado continua aplicável. A mesma água pode atravessar vários processos, portanto fluxos brutos não são estoques, chuva acumulada ou frações causais da precipitação.

Riming significa aqui coleta de água de nuvem por gelo/neve/graupel. A coleta de chuva por partículas congeladas é mostrada separadamente; não se soma novamente riming dentro de uma categoria genérica de coleta. Nucleação de gelo é separada da deposição sobre gelo existente. As tendências de número de partículas do Thompson não fazem parte deste orçamento de massa; podem influenciar as taxas, mas não foram isoladas causalmente.

## Fontes de graupel

'''
    report+=table_md(top[['experiment','category','production_kg_m2','share_percent']])
    report+='''

O Thompson médio não apresenta fontes de massa de gelo ou graupel nesta execução. No Thompson alto, a produção de graupel vem da interação chuva–gelo de nuvem (`prg_rci`, com as perdas correspondentes `prr_rci` e `pri_rci`). No WSM6 médio, coleta chuva–gelo e riming são as contribuições dominantes; no WSM6 alto domina a conversão de neve em graupel (`pgaut`). Portanto, “congelamento” isolado não descreve a principal fonte direta de graupel nesses membros.

Uma particularidade desta versão do WSM6 é aplicar `paacw` tanto à neve quanto ao graupel, com duas retiradas correspondentes da água de nuvem. Foram registradas essas duas transferências exatamente como executadas, sem corrigi-las ou substituir por taxas potenciais `psacw/pgacw`. Elas não são duplicação do diagnóstico: aparecem nas próprias equações aplicadas. A relevância física dessa implementação merece revisão separada, sem alterar silenciosamente o experimento.

Esta é uma atribuição das transferências realizadas, não um experimento causal de desligamento de processos. Deposição e congelamento podem criar precursores que depois alimentam coleta ou conversão; a maior fonte direta de qg não identifica sozinha toda a cadeia causal. Tampouco a produção de qg determina isoladamente os flashes: o EP depende da coexistência com neve, temperatura, água líquida e w. A sensibilidade à resolução documentada anteriormente permanece.

## Arquivos

'''
    report+=f"Execução: `{manifest['directory']}`. `process_integrals.csv` guarda os processos individuais e seis espécies; `figure7_individual_processes.csv` e `figure7_grouped_processes.csv` guardam fontes/perdas de gelo total, graupel e água de nuvem; `figure7_plotted_values.csv` contém exatamente as barras. Cada membro tem `process_tendencies.nc`, CSV de taxas aplicadas, séries e configuração. `species_balances.csv` contém o fechamento integrado.\n\n"
    maxres=max(r['validation']['absolute_unattributed_mm'] for r in manifest['runs'].values())
    report+=f"Maior soma absoluta do resíduo não atribuído, sobre todas as espécies/células/passos de um membro: {maxres:.4g} kg/m². Esse resíduo foi verificado antes da geração da figura, não usado para ajustar os processos.\n"
    (OUT/'RELATORIO_FIGURA7.md').write_text(report,encoding='utf-8')
    caption='Figura 7. Fontes microfísicas de gelo total e graupel nos membros médio e alto com pulso de 3 m/s da Figura 6. Barras integram as contribuições positivas aplicadas após limitadores durante duas horas, ponderadas pela massa seca inicial. O gelo total reúne qi, qs e qg antes de separar contribuições positivas. Riming (coleta de água de nuvem) é distinto das demais coletas. Perdas, sedimentação e ajustes estão nas tabelas de balanço. As execuções são idênticas bit a bit às originais; não foi alterada a física. A comparação é diagnóstica e não implica convergência numérica ou calibração observacional do EP.\n'
    (out/'Figura7_legenda.txt').write_text(caption,encoding='utf-8')
    save_json(out/'Figura7_proveniencia.json',dict(manifest_sha256=sha(root/'manifest.json'),table_sha256=sha(root/'figure7_plotted_values.csv'),
              png_sha256=sha(out/'Figura7_processos_microfisicos.png'),pdf_sha256=sha(out/'Figura7_processos_microfisicos.pdf'),script_sha256=sha(Path(__file__))))
    print(top[['experiment','category','production_kg_m2','share_percent']].to_string(index=False))


if __name__=='__main__': main()
