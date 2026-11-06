# Figura 7 — tendências microfísicas aplicadas

## Escopo

Foram reproduzidos os 12 membros da Figura 6, usando exclusivamente GOAMAZON original. A Figura 7 destaca as camadas média e alta com pulso de 3 m/s, nas quais aparecem as diferenças de gelo e eletrificação. Os demais membros estão nos arquivos de diagnóstico e não produzem graupel. Há, porém, nucleação de gelo no WSM6 baixo e nos controles WSM6 sem pulso; nos controles médio e alto há também deposição. Portanto, ausência de graupel/EP não significa ausência de processos de gelo. Esta comparação preserva a perturbação histórica da Figura 6; não mistura o ensemble posterior de aporte controlado.

## Método e validação

As inserções dentro de Thompson e WSM6 são delimitadas e somente observacionais. Remover os blocos de diagnóstico recupera exatamente o texto anterior de cada parametrização, cujos arquivos e hashes foram preservados em `proveniencia`. Todas as variáveis dos 12 perfis são idênticas bit a bit às saídas originais da Figura 6. Testes adicionais com o diagnóstico desligado também passaram para ambos os esquemas.

O registro é feito após a limitação das taxas. As tendências volumétricas do Thompson são divididas pela densidade usada naquele ponto e multiplicadas pelo passo efetivo; condensação e evaporação já expressas como razão de mistura não recebem outra divisão. No WSM6, as taxas finais são multiplicadas por dtcld em cada subpasso. Congelamentos/fusões diretos registram o incremento realmente aplicado. Sedimentação é obtida pela diferença de estado/tendência antes e depois de sua rotina, incluindo subpassos. Clipping, limiares e substituições de tendências ficam separados dos processos físicos.

Os arquivos NetCDF contêm Δq(time,level,process,species), em kg/kg por passo, para vapor, água de nuvem, chuva, gelo de nuvem, neve e graupel; dividir pelo passo fornece kg/kg/s. As tabelas incluem produção, perda e saldo por processo. A figura integra a parte positiva de cada processo em tempo e altura com a massa seca inicial de referência. Para gelo total, soma-se primeiro Δqi+Δqs+Δqg por processo/célula/passo, evitando contar a simples conversão de gelo em neve como geração de gelo. Para graupel, a conversão de neve em graupel é uma fonte legítima de qg.

Sedimentação é redistribuição/saída e não é incluída na produção microfísica bruta. A ponderação fixa pode produzir saldo de sedimentação que não coincide com a precipitação; o orçamento de densidade já documentado continua aplicável. A mesma água pode atravessar vários processos, portanto fluxos brutos não são estoques, chuva acumulada ou frações causais da precipitação.

Riming significa aqui coleta de água de nuvem por gelo/neve/graupel. A coleta de chuva por partículas congeladas é mostrada separadamente; não se soma novamente riming dentro de uma categoria genérica de coleta. Nucleação de gelo é separada da deposição sobre gelo existente. As tendências de número de partículas do Thompson não fazem parte deste orçamento de massa; podem influenciar as taxas, mas não foram isoladas causalmente.

## Fontes de graupel

| experiment | category | production_kg_m2 | share_percent |
| --- | --- | --- | --- |
| thompson_alta_w3 | Coleta chuva–gelo | 0.0157129 | 100 |
| wsm6_alta_w3 | Autoconversão/agregação | 0.461814 | 87.7712 |
| wsm6_alta_w3 | Coleta gelo–gelo | 0.0244567 | 4.64817 |
| wsm6_alta_w3 | Deposição/sublimação | 0.0227138 | 4.31693 |
| wsm6_alta_w3 | Riming (água de nuvem) | 0.0171398 | 3.25755 |
| wsm6_alta_w3 | Coleta chuva–gelo | 3.21042e-05 | 0.00610165 |
| wsm6_media_w3 | Coleta chuva–gelo | 0.601635 | 54.4746 |
| wsm6_media_w3 | Riming (água de nuvem) | 0.382944 | 34.6734 |
| wsm6_media_w3 | Coleta gelo–gelo | 0.0619295 | 5.60736 |
| wsm6_media_w3 | Deposição/sublimação | 0.0338672 | 3.06648 |
| wsm6_media_w3 | Autoconversão/agregação | 0.0240211 | 2.17498 |
| wsm6_media_w3 | Congelamento | 3.4811e-05 | 0.00315194 |

O Thompson médio não apresenta fontes de massa de gelo ou graupel nesta execução. No Thompson alto, a produção de graupel vem da interação chuva–gelo de nuvem (`prg_rci`, com as perdas correspondentes `prr_rci` e `pri_rci`). No WSM6 médio, coleta chuva–gelo e riming são as contribuições dominantes; no WSM6 alto domina a conversão de neve em graupel (`pgaut`). Portanto, “congelamento” isolado não descreve a principal fonte direta de graupel nesses membros.

Uma particularidade desta versão do WSM6 é aplicar `paacw` tanto à neve quanto ao graupel, com duas retiradas correspondentes da água de nuvem. Foram registradas essas duas transferências exatamente como executadas, sem corrigi-las ou substituir por taxas potenciais `psacw/pgacw`. Elas não são duplicação do diagnóstico: aparecem nas próprias equações aplicadas. A relevância física dessa implementação merece revisão separada, sem alterar silenciosamente o experimento.

Esta é uma atribuição das transferências realizadas, não um experimento causal de desligamento de processos. Deposição e congelamento podem criar precursores que depois alimentam coleta ou conversão; a maior fonte direta de qg não identifica sozinha toda a cadeia causal. Tampouco a produção de qg determina isoladamente os flashes: o EP depende da coexistência com neve, temperatura, água líquida e w. A sensibilidade à resolução documentada anteriormente permanece.

## Arquivos

Execução: `resultados\ensemble_goamazon\processos\20261105T014954121326Z`. `process_integrals.csv` guarda os processos individuais e seis espécies; `figure7_individual_processes.csv` e `figure7_grouped_processes.csv` guardam fontes/perdas de gelo total, graupel e água de nuvem; `figure7_plotted_values.csv` contém exatamente as barras. Cada membro tem `process_tendencies.nc`, CSV de taxas aplicadas, séries e configuração. `species_balances.csv` contém o fechamento integrado.

Maior soma absoluta do resíduo não atribuído, sobre todas as espécies/células/passos de um membro: 7.629e-13 kg/m². Esse resíduo foi verificado antes da geração da figura, não usado para ajustar os processos.
