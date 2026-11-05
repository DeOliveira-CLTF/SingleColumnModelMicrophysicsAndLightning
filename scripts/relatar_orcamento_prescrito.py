"""Write traceable budget and controlled-resolution tables, without figures."""
import pandas as pd
from scm_artigo import *
from orcamento_coluna_prescrita import OUT as BUDGET
from resolucao_controlada import OUT as RESOLUTION


def markdown(frame):
    cols=list(frame.columns)
    def fmt(v):
        if pd.isna(v): return '—'
        return f'{v:.8g}' if isinstance(v,(float,np.floating)) else str(v)
    return '\n'.join(['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']+
                     ['| '+' | '.join(map(fmt,row))+' |' for row in frame.itertuples(index=False,name=None)])


def main():
    budget=json.loads((BUDGET/'latest.json').read_text(encoding='utf-8'))
    resolution=json.loads((RESOLUTION/'latest.json').read_text(encoding='utf-8'))
    b=pd.read_csv(ROOT/budget['directory']/'summary.csv')
    r=pd.read_csv(ROOT/resolution['directory']/'summary.csv')
    contrasts=[]
    fields=['rain_final_mm','qg_max_g_kg','flash_integral_nominal']
    for eid,g in r.groupby('experiment'):
        for label,a,c in [('tempo_240_10_para_5',g[(g.nz==240)&(g.dt_s==10)],g[(g.nz==240)&(g.dt_s==5)]),
                          ('vertical_120_para_240_dt5',g[(g.nz==120)&(g.dt_s==5)],g[(g.nz==240)&(g.dt_s==5)])]:
            a=a.iloc[0]; c=c.iloc[0]
            for field in fields:
                contrasts.append(dict(experiment=eid,comparison=label,metric=field,before=a[field],after=c[field],
                                      change_percent=100*(c[field]/a[field]-1) if a[field]!=0 else None))
    contrasts=pd.DataFrame(contrasts)
    contrasts.to_csv(ROOT/resolution['directory']/'contrasts.csv',index=False,float_format='%.17g')
    report='''# Orçamento da coluna prescrita e resolução com aporte controlado

Foi adotada a escolha do usuário: manter o experimento prescrito e explicitar fontes, sumidouros e fluxos. Apenas `data_goamazon/GOAMAZON_SINGLE_SCM_driver.nc` foi utilizado. As saídas anteriores foram preservadas. Não foram produzidas novas figuras.

## Definições e limites

As razões de mistura são ponderadas por massa seca anterior à perturbação (m_ref), fixa durante a integração. O inventário Euleriano W usa rho_d = p/[T (Rd + Rv rv)] a cada saída. A fonte inicial é externa, a temperatura fica fixa e Lv × dose é um diagnóstico latente, não um orçamento de energia fechado. Os valores de dose não são frações causais da chuva.

O driver observa o estado antes e depois do transporte e da microfísica. O transporte atua somente em theta e vapor; não transporta os hidrometeoros com w. Sedimentação pertence à etapa microfísica. `rainnc` é precipitação total em equivalente de água: neve e graupel são componentes, não parcelas adicionais a somar novamente.

A auditoria reconstrói cada subpasso do upwind. Define fluxo de massa de referência nas interfaces por interpolação linear de rho_ref × w, com extensão constante nos extremos, e doador upwind de vapor. O fluxo de água inferior entra com sinal positivo, e o superior sai. A decomposição é:

`transporte = fluxo_inferior − fluxo_superior + Σ rv × divergência(fluxo_de_ar) × dt + correção_do_estêncil`.

A correção mede a diferença entre o operador advectivo efetivamente usado e essa representação em fluxos. A separação depende da reconstrução de interfaces documentada; ela não é uma medição de circulação horizontal. O termo de divergência é a contribuição da continuidade de ar ausente no experimento prescrito, não uma forçante horizontal observada. O contorno de entrada usa gradiente zero; um nível extremo congelado pelo operador não implica fluxo físico zero.

O segundo balanço usa a densidade nativa de cada esquema antes da chamada: WSM6 recebe p/[Rd T (1 + 0,61 rv)]; Thompson calcula 0,622 p/[287,04 T (rv + 0,622)]. Nenhuma delas é substituída silenciosamente pela massa seca de referência.

`ΔW + precipitação = transporte_ref + resíduo_micro_nativo + conversão_de_pesos + ajuste_de_inventário_por_densidade`.

O resíduo nativo é `Σ m_nativa Δrt_micro + precipitação`. Conversão de pesos é `Σ (m_ref − m_nativa) Δrt_micro`. O ajuste de densidade é a diferença entre a mudança Euleriana de inventário e a mudança com pesos fixos. Esses dois últimos termos são diagnósticos de representação, não fontes físicas independentes. O resíduo microfísico permanece sem atribuição a processos: pode conter limitadores, limiares e particularidades da sedimentação. Fechamento contábil não demonstra conservação física nem identifica sua causa.

## Quatro integrações controladas de referência, 60 níveis e 20 s

'''
    report+=markdown(b[['experiment','precip_mm','transport_ref_mm','prescribed_divergence_mm','stencil_correction_mm','micro_native_residual_mm','weight_conversion_mm','density_remap_mm']])
    report+='\n\nTodos os quatro testes com diagnóstico ligado/desligado produziram saídas idênticas bit a bit. O máximo resíduo contábil por passo foi '+f'{b.max_step_bookkeeping_residual_mm.max():.3g}'+' mm. Isso verifica a instrumentação e a identidade de orçamento, mantendo os resíduos nativos acima explícitos.\n'
    report+='\n## Matriz de resolução controlada\n\nSão 36 combinações: quatro casos, 60/120/240 níveis e passos de 20/10/5 s; quatro referências reutilizadas e 32 novas integrações. Todas têm campos finitos, razões de mistura não negativas, cronologia validada, reconstrução do EP compatível e orçamento auditado. A dose inicial permanece igual entre grades, com diferença menor que 10⁻¹¹ mm. A projeção das bordas é por sobreposição física das células; os perfis de fundo ainda diferem por regradeamento.\n\n'
    report+=markdown(r[((r.nz==60)&(r.dt_s==20))|((r.nz==240)&(r.dt_s==5))][['experiment','nz','dt_s','initial_added_mm']+fields])
    report+='\n\nOs últimos refinamentos separados são:\n\n'+markdown(contrasts)
    report+='''

Estes contrastes medem sensibilidade; não constituem prova de convergência assintótica nem permitem estimar ordem formal com segurança. A fonte inicial agora é equivalente, mas o fundo regradeado, a difusão numérica, os passos microfísicos e a amostragem do EP ainda participam da resposta. Integrais de flashes são nominais, não calibração observacional.

Persistem diferenças relevantes no último refinamento vertical a 5 s: a chuva do Thompson alto aumenta 26,61%, e o máximo de graupel do WSM6 alto aumenta 28,85%. No refinamento temporal de 10 para 5 s em 240 níveis, o máximo de graupel do WSM6 médio diminui 12,46%. Portanto, o controle da dose removeu uma confusão importante, mas ainda não sustenta convergência dos diagnósticos de interesse. Nos campos nulos do Thompson médio, a variação percentual é indefinida e aparece como travessão.

## Próximas etapas ainda não concluídas

O orçamento registra a mudança microfísica líquida total. Ainda faltam as taxas efetivamente aplicadas por condensação, deposição, congelamento e coleta/riming, com limitadores, subpassos e atualizações diretas. Não se pode atribuir a diferença Thompson–WSM6 a uma dessas taxas apenas a partir deste balanço. Riming deverá ser tratado como subconjunto da coleta.

A reconstrução do EP validou a fórmula atual; os testes contrafactuais separados de w, janela térmica e conversão EP→flashes ainda estão pendentes. Sem observações compatíveis, qualquer variação de coeficientes será sensibilidade, não calibração.

## Proveniência e reprodução

'''
    report+=f"Orçamento: `{budget['directory']}`. Resolução: `{resolution['directory']}`. Os manifests registram configurações e hashes. Scripts: `orcamento_coluna_prescrita.py`, `resolucao_controlada.py`, `relatar_orcamento_prescrito.py`; testes em `test_orcamento_prescrito.py`.\n"
    (BUDGET/'RELATORIO_ORCAMENTO_E_RESOLUCAO.md').write_text(report,encoding='utf-8')
    print(contrasts.to_string(index=False))


if __name__=='__main__': main()
