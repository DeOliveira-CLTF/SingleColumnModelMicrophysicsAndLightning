# Orçamento da coluna prescrita e resolução com aporte controlado

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

| experiment | precip_mm | transport_ref_mm | prescribed_divergence_mm | stencil_correction_mm | micro_native_residual_mm | weight_conversion_mm | density_remap_mm |
| --- | --- | --- | --- | --- | --- | --- | --- |
| wsm6_media_w3 | 6.7487141 | 9.0674779 | 8.9349788 | 0.12970174 | -0.0023353575 | 0.070454894 | 0.17492421 |
| thompson_media_w3 | 2.4154881 | 9.0740297 | 8.9412173 | 0.13001513 | -0.00015840228 | 0.01501712 | 0.18400565 |
| thompson_alta_w3 | 0.20243395 | 3.0293631 | 2.9314057 | 0.097957358 | -0.00011253385 | 0.0067955577 | 0.12222026 |
| wsm6_alta_w3 | 0.69155632 | 3.063775 | 2.9667329 | 0.097042125 | -0.00053538568 | 0.0019208102 | 0.16140613 |

Todos os quatro testes com diagnóstico ligado/desligado produziram saídas idênticas bit a bit. O máximo resíduo contábil por passo foi 1.43e-14 mm. Isso verifica a instrumentação e a identidade de orçamento, mantendo os resíduos nativos acima explícitos.

## Matriz de resolução controlada

São 36 combinações: quatro casos, 60/120/240 níveis e passos de 20/10/5 s; quatro referências reutilizadas e 32 novas integrações. Todas têm campos finitos, razões de mistura não negativas, cronologia validada, reconstrução do EP compatível e orçamento auditado. A dose inicial permanece igual entre grades, com diferença menor que 10⁻¹¹ mm. A projeção das bordas é por sobreposição física das células; os perfis de fundo ainda diferem por regradeamento.

| experiment | nz | dt_s | initial_added_mm | rain_final_mm | qg_max_g_kg | flash_integral_nominal |
| --- | --- | --- | --- | --- | --- | --- |
| wsm6_media_w3 | 60 | 20 | 2.1232163 | 6.7487141 | 0.68652144 | 4.0283081 |
| wsm6_media_w3 | 240 | 5 | 2.1232163 | 6.7403734 | 0.67787459 | 3.5046546 |
| thompson_media_w3 | 60 | 20 | 2.1232163 | 2.4154881 | 0 | 0 |
| thompson_media_w3 | 240 | 5 | 2.1232163 | 2.3925367 | 0 | 0 |
| thompson_alta_w3 | 60 | 20 | 0.66013065 | 0.20243395 | 0.0097300341 | 8.8066685 |
| thompson_alta_w3 | 240 | 5 | 0.66013065 | 0.3167967 | 0.015107228 | 8.1051561 |
| wsm6_alta_w3 | 60 | 20 | 0.66013065 | 0.69155632 | 0.3720698 | 12.481172 |
| wsm6_alta_w3 | 240 | 5 | 0.66013065 | 0.6708591 | 0.59284626 | 11.9774 |

Os últimos refinamentos separados são:

| experiment | comparison | metric | before | after | change_percent |
| --- | --- | --- | --- | --- | --- |
| thompson_alta_w3 | tempo_240_10_para_5 | rain_final_mm | 0.32825287 | 0.3167967 | -3.4900441 |
| thompson_alta_w3 | tempo_240_10_para_5 | qg_max_g_kg | 0.01506355 | 0.015107228 | 0.28995873 |
| thompson_alta_w3 | tempo_240_10_para_5 | flash_integral_nominal | 8.0812651 | 8.1051561 | 0.29563372 |
| thompson_alta_w3 | vertical_120_para_240_dt5 | rain_final_mm | 0.25021441 | 0.3167967 | 26.610094 |
| thompson_alta_w3 | vertical_120_para_240_dt5 | qg_max_g_kg | 0.014773148 | 0.015107228 | 2.2614009 |
| thompson_alta_w3 | vertical_120_para_240_dt5 | flash_integral_nominal | 8.7972715 | 8.1051561 | -7.8673868 |
| thompson_media_w3 | tempo_240_10_para_5 | rain_final_mm | 2.3889957 | 2.3925367 | 0.14822348 |
| thompson_media_w3 | tempo_240_10_para_5 | qg_max_g_kg | 0 | 0 | — |
| thompson_media_w3 | tempo_240_10_para_5 | flash_integral_nominal | 0 | 0 | — |
| thompson_media_w3 | vertical_120_para_240_dt5 | rain_final_mm | 2.4046209 | 2.3925367 | -0.50254157 |
| thompson_media_w3 | vertical_120_para_240_dt5 | qg_max_g_kg | 0 | 0 | — |
| thompson_media_w3 | vertical_120_para_240_dt5 | flash_integral_nominal | 0 | 0 | — |
| wsm6_alta_w3 | tempo_240_10_para_5 | rain_final_mm | 0.67262789 | 0.6708591 | -0.26296641 |
| wsm6_alta_w3 | tempo_240_10_para_5 | qg_max_g_kg | 0.57272823 | 0.59284626 | 3.5126668 |
| wsm6_alta_w3 | tempo_240_10_para_5 | flash_integral_nominal | 11.969631 | 11.9774 | 0.064903332 |
| wsm6_alta_w3 | vertical_120_para_240_dt5 | rain_final_mm | 0.68417175 | 0.6708591 | -1.9458047 |
| wsm6_alta_w3 | vertical_120_para_240_dt5 | qg_max_g_kg | 0.46008921 | 0.59284626 | 28.854634 |
| wsm6_alta_w3 | vertical_120_para_240_dt5 | flash_integral_nominal | 12.203938 | 11.9774 | -1.8562697 |
| wsm6_media_w3 | tempo_240_10_para_5 | rain_final_mm | 6.7366312 | 6.7403734 | 0.055549834 |
| wsm6_media_w3 | tempo_240_10_para_5 | qg_max_g_kg | 0.7743897 | 0.67787459 | -12.463378 |
| wsm6_media_w3 | tempo_240_10_para_5 | flash_integral_nominal | 3.6407249 | 3.5046546 | -3.737451 |
| wsm6_media_w3 | vertical_120_para_240_dt5 | rain_final_mm | 6.7439656 | 6.7403734 | -0.053265433 |
| wsm6_media_w3 | vertical_120_para_240_dt5 | qg_max_g_kg | 0.63790989 | 0.67787459 | 6.2649444 |
| wsm6_media_w3 | vertical_120_para_240_dt5 | flash_integral_nominal | 3.5234187 | 3.5046546 | -0.53255553 |

Estes contrastes medem sensibilidade; não constituem prova de convergência assintótica nem permitem estimar ordem formal com segurança. A fonte inicial agora é equivalente, mas o fundo regradeado, a difusão numérica, os passos microfísicos e a amostragem do EP ainda participam da resposta. Integrais de flashes são nominais, não calibração observacional.

Persistem diferenças relevantes no último refinamento vertical a 5 s: a chuva do Thompson alto aumenta 26,61%, e o máximo de graupel do WSM6 alto aumenta 28,85%. No refinamento temporal de 10 para 5 s em 240 níveis, o máximo de graupel do WSM6 médio diminui 12,46%. Portanto, o controle da dose removeu uma confusão importante, mas ainda não sustenta convergência dos diagnósticos de interesse. Nos campos nulos do Thompson médio, a variação percentual é indefinida e aparece como travessão.

## Próximas etapas ainda não concluídas

O orçamento registra a mudança microfísica líquida total. Ainda faltam as taxas efetivamente aplicadas por condensação, deposição, congelamento e coleta/riming, com limitadores, subpassos e atualizações diretas. Não se pode atribuir a diferença Thompson–WSM6 a uma dessas taxas apenas a partir deste balanço. Riming deverá ser tratado como subconjunto da coleta.

A reconstrução do EP validou a fórmula atual; os testes contrafactuais separados de w, janela térmica e conversão EP→flashes ainda estão pendentes. Sem observações compatíveis, qualquer variação de coeficientes será sensibilidade, não calibração.

## Proveniência e reprodução

Orçamento: `resultados\ensemble_goamazon\orcamento_prescrito\20261105T012206982025Z`. Resolução: `resultados\ensemble_goamazon\resolucao_controlada\20261105T012316616253Z`. Os manifests registram configurações e hashes. Scripts: `orcamento_coluna_prescrita.py`, `resolucao_controlada.py`, `relatar_orcamento_prescrito.py`; testes em `test_orcamento_prescrito.py`.
