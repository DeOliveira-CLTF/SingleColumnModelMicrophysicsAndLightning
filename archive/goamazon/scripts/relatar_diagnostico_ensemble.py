"""Summarize the completed sensitivity study; no figures or new integrations."""
import pandas as pd
from diagnosticar_ensemble_goamazon import *


def tab(headers,rows):
    def val(x):
        if isinstance(x,(float,np.floating)): return 'indefinido' if np.isnan(x) else f'{x:.5g}'
        return str(x)
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(val(x) for x in row)+' |' for row in rows])


def main():
    manifest=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
    df=pd.read_csv(OUT/'sensibilidade_resolucao.csv')
    if len(df)!=36 or not (df.status=='valid').all(): raise ValueError('Incomplete/invalid sensitivity matrix')
    original=json.loads((BASE/'manifest.json').read_text(encoding='utf-8'))
    details=[]; events=[]
    for eid in CASES:
        r=original['members'][eid]; ds=load(ROOT/r['directory']/'profile_1.nc')
        d,e=detailed(eid,ds,r['config']); details.append(d); events+=e
    pd.DataFrame(details).to_csv(OUT/'diagnosticos_quatro_membros.csv',index=False,float_format='%.17g')
    pd.DataFrame(events).to_csv(OUT/'eventos_EP.csv',index=False,float_format='%.17g')
    detail=pd.DataFrame(details).set_index('experiment'); event=pd.DataFrame(events)
    keys=['rain_final_mm','qg_max_g_kg','qg_mixed_mass_max_kg_m2','flash_max_min','flash_integral_nominal','water_added_fixed_dry_air_mm']
    contrasts=[]; field_errors=[]
    for case in CASES:
        frame=df[df.case==case].set_index(['nz','dt_s'])
        for kind,a,b in [('temporal_60_20to10',(60,10),(60,20)),('temporal_60_10to5',(60,5),(60,10)),
                         ('temporal_240_20to10',(240,10),(240,20)),('temporal_240_10to5',(240,5),(240,10)),
                         ('vertical_dt5_60to120',(120,5),(60,5)),('vertical_dt5_120to240',(240,5),(120,5)),
                         ('baseline_to_finest',(240,5),(60,20))]:
            for key in keys:
                av=float(frame.loc[a,key]); bv=float(frame.loc[b,key])
                contrasts.append(dict(case=case,comparison=kind,metric=key,coarse=bv,fine=av,
                                      difference=av-bv,relative_percent=100*(av-bv)/bv if bv else None))
        refrec=manifest['runs'][f'{case}_nz240_dt5']; ref=load(ROOT/refrec['directory']/'profile_1.nc')
        refT=temperature(ref)
        for nz in [60,120,240]:
            for dt in [20,10,5]:
                if nz==240 and dt==5: continue
                rec=manifest['runs'][f'{case}_nz{nz}_dt{dt}']; path=ROOT/rec['directory']/'profile_1.nc'
                if sha(path)!=rec['output_sha256']: raise ValueError('Changed output')
                ds=load(path); temp=temperature(ds)
                # Common timestamps every 20 s; reference centers inside both domains.
                ir=np.arange(0,len(ref.time),4); ii=np.arange(0,len(ds.time),int(20/dt))
                valid=(ref.lev.values>=float(ds.lev.min()))&(ref.lev.values<=float(ds.lev.max()))
                zr=ref.lev.values[valid]; dz=np.diff(interfaces(ref.lev.values,20000))[valid]
                for key in ['qc','qi','qs','qg','temperature']:
                    a=temp[ii] if key=='temperature' else ds[key].values[ii]
                    b=refT[ir][:,valid] if key=='temperature' else ref[key].values[ir][:,valid]
                    mapped=np.stack([np.interp(zr,ds.lev.values,row) for row in a])
                    absolute=float(np.sum(np.abs(mapped-b)*dz)/(len(ii)*dz.sum()))
                    denominator=float(np.sum(np.abs(b)*dz))
                    field_errors.append(dict(case=case,nz=nz,dt_s=dt,variable=key,
                                             mean_abs_error_native_units=absolute,
                                             relative_L1_percent=100*np.sum(np.abs(mapped-b)*dz)/denominator if denominator and key!='temperature' else None))
    contrast=pd.DataFrame(contrasts)
    contrast.to_csv(OUT/'contrastes_resolucao.csv',index=False,float_format='%.17g')
    pd.DataFrame(field_errors).to_csv(OUT/'diferencas_perfis_vs_240_5.csv',index=False,float_format='%.17g')
    w=pd.read_csv(OUT/'agua_introduzida_12_membros.csv').set_index('experiment')
    water_table=tab(['Caso','Água acrescentada (mm)','Chuva (mm)','Aporte / chuva (%)','Δ proxy latente (MJ/m²)'],
                    [[e,w.loc[e,'water_added_fixed_dry_air_mm'],w.loc[e,'rain_final_mm'],w.loc[e,'added_over_rain_percent'],w.loc[e,'latent_energy_proxy_delta_MJ_m2']] for e in CASES])
    transport_table=tab(['Caso','ΔW + P (mm)','Advecção reconstruída (mm)','Restante (mm)'],
                        [[e,detail.loc[e,'apparent_water_residual_final_mm'],detail.loc[e,'water_advection_cumulative_final_mm'],detail.loc[e,'water_remainder_after_advection_final_mm']] for e in CASES])
    hydro_table=tab(['Caso','qc máx. (g/kg)','qi máx. (g/kg)','qs máx. (g/kg)','qg máx. (g/kg)','Início EP (min)','Pico EP (min)'],
                    [[e]+[detail.loc[e,k+'_max_g_kg'] for k in ['qc','qi','qs','qg']]+[detail.loc[e,'flash_onset_h']*60,detail.loc[e,'flash_peak_h']*60] for e in CASES])
    peak=event[event.event=='pico_EP']
    peak_table=tab(['Caso','Hora (min)','z dominante (km)','T (°C)','w utilizado (m/s)','qc (g/kg)','qs (g/kg)','qg (g/kg)','EP de centros fora da faixa (%)'],
                   [[r.experiment,r.time_s/60,r.ep_dominant_z_m/1000,r.ep_dominant_temperature_C,r.ep_dominant_w_m_s,
                     r.qc_at_ep_dominant_g_kg,r.qs_at_ep_dominant_g_kg,r.qg_at_ep_dominant_g_kg,r.ep_outside_centers_percent] for r in peak.itertuples()])
    resolution_rows=[]
    for e in CASES:
        f=df[df.case==e].set_index(['nz','dt_s']); coarse=f.loc[(60,20)]; fine=f.loc[(240,5)]
        resolution_rows.append([e,coarse.rain_final_mm,fine.rain_final_mm,coarse.qg_max_g_kg,fine.qg_max_g_kg,coarse.flash_integral_nominal,fine.flash_integral_nominal])
    resolution_table=tab(['Caso','P 60/20','P 240/5','qg máx. 60/20','qg máx. 240/5','Integral 60/20','Integral 240/5'],resolution_rows)
    finest_rows=[]
    for e in CASES:
        for kind in ['temporal_240_10to5','vertical_dt5_120to240']:
            subset=contrast[(contrast.case==e)&(contrast.comparison==kind)].set_index('metric')
            finest_rows.append([e,kind]+[subset.loc[k,'relative_percent'] for k in keys])
    finest_table=tab(['Caso','Comparação','Δ P (%)','Δ qg máx. (%)','Δ Mqg máx. (%)','Δ EP máx. (%)','Δ integral (%)','Δ aporte (%)'],finest_rows)
    dose_rows=[]
    for layer,e in [('média','wsm6_media_w3'),('alta','wsm6_alta_w3')]:
        for nz in [60,120,240]:
            r=df[(df.case==e)&(df.nz==nz)&(df.dt_s==5)].iloc[0]
            dose_rows.append([layer,nz,r.water_added_fixed_dry_air_mm,r.saturated_cells,r.saturated_thickness_m,r.dz_in_saturation_min_m,r.dz_in_saturation_max_m])
    dose_table=tab(['Camada','Níveis','Aporte (mm)','Células','Espessura efetiva (m)','Δz mín. (m)','Δz máx. (m)'],dose_rows)
    report=f'''# Diagnóstico científico do ensemble idealizado GOAMAZON

## Resultado principal

As três limitações levantadas são relevantes. O aporte inicial de água é expressivo, a coluna recebe contribuição do transporte prescrito durante a integração, e a convergência quantitativa de graupel/EP não está demonstrada. A ocorrência qualitativa de EP permaneceu igual nas nove resoluções de cada caso: WSM6 médio e ambos os altos ativos; Thompson médio inativo. Isso não torna as magnitudes independentes da grade ou do passo.

Foram reutilizados os quatro casos de 60 níveis/20 s e executadas 32 integrações adicionais em pastas novas, totalizando 36 configurações avaliadas. Todas passaram nos controles implementados. Nenhuma figura nova foi criada e nenhuma parametrização ou entrada original foi alterada. Aprovação de integridade é distinta de convergência e de validação observacional.

## 1. Água e energia introduzidas

Reconstruí o perfil original na grade de cada membro e comparei-o ao registro inicial perturbado. Como rv, rl e ri são razões por massa de ar seco, a medida principal usa a densidade de ar seco anterior à perturbação, `rho_d0 = p/[T(Rd + Rv rv0)]`, e soma `rho_d0 × Δ(rv+rl+ri) × dz`. O resultado em kg/m² equivale numericamente a mm de água. É um aporte equivalente com massa de ar seco congelada, não um orçamento Euleriano fechado.

{water_table}

Nos três membros eletrificados, o aporte equivale a 30,89%, 336,56% e 103,27% da chuva, respectivamente para WSM6 médio, Thompson alto e WSM6 alto. **Esses percentuais não são frações da chuva originadas na perturbação.** Sem traçadores e um orçamento completo não é possível identificar a origem da água precipitada. Na camada alta, a água acrescentada é maior que o acumulado final de ambos os esquemas. Na média, a chuva WSM6 excede o aporte; isso demonstra disponibilidade de outras fontes/reservatórios, não a eficiência de conversão do vapor acrescentado.

A coluna original tem cerca de 57,23 mm de água sob essa ponderação. A saturação acrescenta 2,633 mm na camada baixa, 2,080 mm na média e 0,819 mm na alta. Todos os acréscimos líquidos destes casos são de vapor; qc e qi não mudam imediatamente na perturbação, e a temperatura inicial é mantida. Os pares sem pulso têm o mesmo aporte e chuva nula em 2 h: acrescentar água, sozinho, não foi suficiente para chover nessa janela.

O CSV também fornece a estimativa com a densidade aproximada usada anteriormente (2,092 mm na média e 0,820 mm na alta), preservando a rastreabilidade, e a diferença de inventário Euleriano com densidade recalculada após a perturbação (2,055 e 0,816 mm). As diferenças decorrem das convenções de massa/densidade; não devem ser confundidas com precisão observacional.

O proxy `Δ∫rho_d0 [Lv rv − Lf qi] dz`, com Lv=2,50×10⁶ e Lf=3,50×10⁵ J/kg do código, cresce cerca de 5,20 MJ/m² na média e 2,05 MJ/m² na alta. Não houve compensação térmica inicial. Esse é um diagnóstico parcial de conteúdo energético latente, não o orçamento completo de entalpia/energia do sistema ou uma estimativa do aquecimento posterior.

### Balanço aparente durante a integração

Para separar o aporte inicial da evolução, calculei `W(t)=Σrho_d0 (qv+qc+qr+qi+qs+qg) dz`. Depois reconstruí, passo a passo, exatamente a atualização upwind do vapor a partir da saída anterior e do w usado no passo, incluindo os subpassos CFL. A coluna seguinte é a soma dessas tendências de transporte com os mesmos pesos fixos.

{transport_table}

Aqui ΔW parte do estado **já perturbado**, portanto não inclui o aporte inicial. O restante é `ΔW + precipitação − advecção reconstruída`. Não deve ser apresentado automaticamente como erro de conservação da microfísica: mistura pesos fixos de ar seco com a conversão de precipitação do driver, enquanto temperatura/densidade evoluem e não há prognóstico de continuidade de massa de ar. É um diagnóstico aparente, com convenção declarada.

A advecção implementada é `−w ∂r/∂z`, não uma divergência de fluxo de água conservativa. Com densidade fixa de referência, sua integral contém termos de contorno e `r ∂(rho_d0 w)/∂z`. A combinação de w prescrito, ausência de continuidade da coluna e contornos com gradiente nulo permite ganho integrado de água. Consequentemente, não é defensável atribuir os totais de chuva apenas ao reservatório inicial nem somente à saturação artificial. Um traçador passivo e um orçamento de fluxo coerente seriam necessários para atribuição causal.

## 2. Evolução vertical dos quatro membros

O controle Thompson médio é necessário para comparar o caso WSM6 médio. Os quatro usam a mesma entrada; os pares de camada têm condições iniciais e perturbações idênticas. Os arquivos derivados preservam todas as amostras de qc, qi, qs, qg, temperatura e w reconstruído, juntamente com EP e precipitação. Há CSVs de perfis em t=0, 15 e 30 min, início/pico do EP, pico de graupel, pico de chuva e fim da integração.

{hydro_table}

**Camada média.** Thompson não forma qi, qs nem qg durante a janela na resolução original; o WSM6 forma as três espécies e ativa EP a 9,33 min. Portanto, a ausência elétrica Thompson tem uma explicação imediata no diagnóstico: falta graupel e neve, que entram no integrando. As saídas não identificam qual taxa microfísica individual causou essa partição; não atribuo a diferença especificamente a nucleação, coleta ou riming sem orçamentos de processos. A água de nuvem máxima é maior no Thompson, apesar da ausência de gelo.

**Camada alta.** Ambos formam neve e graupel, mas Thompson retém muito mais qc e muito menos qi e qg que WSM6. O EP começa a 21,67 min no Thompson e 10,67 min no WSM6. A região de graupel se estende para níveis mais baixos ao longo do tempo; a sequência é compatível com redistribuição/sedimentação e evolução microfísica, sem permitir separar quantitativamente esses processos. Os máximos elétricos ocorrem depois do término do pulso adicional: 65,0 min no Thompson e 37,67 min no WSM6. O w de fundo continua presente.

### Condições no máximo elétrico

{peak_table}

O nível dominante é o centro com maior densidade local do integrando EP, antes da quadratura vertical. No Thompson alto, o máximo elétrico acontece com qc e qi nulos nesse nível, porém qs e qg positivos. Isso é permitido pela fórmula: qi não entra diretamente; qs fornece a aproximação de concentração de cristais, qc desloca a temperatura de reversão, e o fator de ascensão é `1 + max(w,0)/5`, que não zera quando w=0. Portanto, coincidência de qc, qi, qs e qg ou ascensão forte não é um requisito lógico para EP positivo neste código. Trata-se de interpretação da formulação, não de validação do mecanismo atmosférico.

Reconstruí integralmente o EP em Python e comparei cada passo posterior ao inicial com o NetCDF (rtol 5×10⁻¹¹, atol 10⁻¹²). O w associado ao registro tₙ foi calculado em tₙ₋₁, antes do avanço do passo; ignorar esse deslocamento introduziria erro na comparação temporal. O campo w é reconstruído a partir de uma função conhecida, não uma variável prognóstica salva.

Há ainda uma limitação espacial do EP: o código integra um intervalo quando **pelo menos um** de seus extremos está entre −40 e 0 °C. Assim, o outro extremo pode contribuir fora dessa faixa. No pico do WSM6 médio, cerca de 7,70% do EP integrado vem de centros fora da faixa, sobretudo junto à fronteira quente. Isso não equivale a 7,70% da taxa elétrica, pois a taxa é proporcional à raiz quadrada do EP. O comportamento foi reproduzido, não corrigido, e pode participar da sensibilidade vertical.

## 3. Sensibilidade a Δz e Δt

Matriz cruzada: N=60,120,240; Δt=20,10,5 s; topo 20 km; primeira espessura 50,25,12,5 m, respectivamente. Cada integração dura 2 h. A malha continua esticada, por isso o Δz relevante na camada de interesse é muito maior que a primeira espessura. Temperatura, entrada, pulso e rotinas físicas seguem os mesmos parâmetros. A saturação é reaplicada aos centros de cada grade; não há remapeamento conservativo da perturbação.

{dose_table}

Em particular, a camada alta original abrange 2.435,8 m de células inteiras, embora a faixa nominal tenha 2.000 m. O refino modifica a espessura efetiva e a água acrescentada. Isso mistura erro de discretização do transporte, amostragem da microfísica/EP e mudança da condição inicial efetivamente aplicada. **Este ensaio é uma análise de sensibilidade de resolução; não é uma demonstração limpa de ordem de convergência.**

O passo Δt controla também as chamadas da microfísica, o acoplamento sequencial e a amostragem do EP. Logo, sua sensibilidade não isola a difusão do upwind. Da mesma forma, os testes verticais alteram a quadratura do EP e o suporte da saturação, além do transporte.

Comparação da configuração original com a mais refinada testada (chuva em mm, qg em g/kg, integral elétrica nominal):

{resolution_table}

Mudanças nos últimos refinamentos, percentuais com denominador na configuração anterior:

{finest_table}

A chuva pode parecer estável enquanto o máximo de graupel ainda muda. No WSM6 médio, reduzir o passo de 10 para 5 s com 240 níveis altera pouco a chuva, mas reduz o máximo de graupel em mais de 10%. No Thompson alto, dobrar de 120 para 240 níveis a 5 s muda a chuva e o máximo de graupel em dezenas de porcento. Não é adequado declarar esses diagnósticos quantitativamente convergidos.

A ocorrência qualitativa de EP é robusta dentro desta matriz limitada: os três casos continuam ativos e Thompson médio permanece inativo. Os totais elétricos e os máximos, contudo, dependem da resolução. Não fiz extrapolação de Richardson nem estimei ordem assintótica: há respostas não monotônicas, mudanças de amostragem da saturação e efeitos dos limiares físicos. O nível 240/5 é uma referência de comparação, não uma solução exata.

O arquivo `diferencas_perfis_vs_240_5.csv` compara qc, qi, qs, qg e T a cada 20 s, interpolando linearmente para centros da grade refinada dentro do domínio comum e ponderando pelas espessuras dessa grade. Fornece erro absoluto médio (kg/kg ou K) e diferença L1 relativa dos hidrometeoros; séries de referência nulas têm percentual indefinido. A interpolação é somente diagnóstica, não remapeamento conservativo nem medida isolada de erro de truncamento.

## 4. Implicações para o artigo

Os resultados devem ser chamados de **experimentos idealizados baseados no perfil inicial GOAMAZON**. Não reproduzem a evolução do evento: usam pressão e wa inicial fixas, sem tendências horizontais variáveis, fluxos de superfície e MYNN. A classificação anterior de “aprovado” é apenas de integridade/consistência dos controles implementados; não implica conservação, convergência ou habilidade observacional.

Antes de interpretar quantitativamente totais de chuva ou eletricidade, são prioritários: (1) definir uma perturbação de água e energia controlada, com dose e suporte físico equivalentes entre grades; (2) obter orçamento de massa coerente com o transporte e, para atribuição, traçar a água acrescentada; (3) continuar os testes temporais/verticais depois de controlar a mudança de entrada entre grades; (4) registrar tendências microfísicas se o objetivo for atribuir diferenças a processos específicos. Nenhuma dessas mudanças físicas foi feita silenciosamente neste diagnóstico.

Não produzi novas figuras nem substituí os resultados anteriores. A Figura 6 existente deve ser tratada como síntese exploratória de 60 níveis/20 s, acompanhada destas ressalvas. O contraste qualitativo entre esquemas é sustentado nesta matriz; as magnitudes e as explicações causais permanecem limitadas.

## Arquivos reproduzíveis

- `agua_introduzida_12_membros.csv`: aporte inicial, convenções de massa, proxy latente e relação com chuva.
- `diagnosticos_quatro_membros.csv`, `eventos_EP.csv`: métricas, eventos e balanço aparente.
- `perfis/*_diagnostico.nc`: qc, qi, qs, qg, T, w utilizado e EP local em todos os tempos; `*_series.csv` e `*_eventos_verticais.csv` para inspeção.
- `sensibilidade_resolucao.csv`: 36 configurações; `contrastes_resolucao.csv` e `diferencas_perfis_vs_240_5.csv`: comparações automáticas.
- `manifest.json` e manifests datados: configuração, hashes, validações e diretórios de execução. Os quatro baselines foram reutilizados; 32 casos novos foram executados.

Comandos na raiz: `python scripts/diagnosticar_ensemble_goamazon.py` e `python scripts/relatar_diagnostico_ensemble.py`. O primeiro reutiliza arquivos compatíveis; o segundo apenas analisa e relata. Não geram figuras.
'''
    (OUT/'RELATORIO_DIAGNOSTICO_CIENTIFICO.md').write_text(report,encoding='utf-8')
    save_json(OUT/'auditoria_diagnostico.json',{'input_sha256':sha(INPUT),
        'all_baseline_outputs_unchanged':all(sha(ROOT/r['directory']/'profile_1.nc')==r['output_sha256'] for r in original['members'].values()),
        'all_sources_unchanged':all(sha(ROOT/p)==h for p,h in original['identity']['sources'].items()),
        'valid_cases':len(df),'new_runs':sum(r['action']=='executed' for r in manifest['runs'].values()),
        'scripts_sha256':{p.name:sha(p) for p in [ROOT/'scripts/diagnosticar_ensemble_goamazon.py',ROOT/'scripts/relatar_diagnostico_ensemble.py']}})
    print(resolution_table); print(finest_table); print(transport_table)


if __name__=='__main__': main()
