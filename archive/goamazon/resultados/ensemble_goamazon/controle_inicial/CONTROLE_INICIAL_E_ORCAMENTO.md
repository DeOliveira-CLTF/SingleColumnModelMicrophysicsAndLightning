# Controle inicial e sequência do orçamento hídrico

## Etapa 1 implementada e verificada

O programa `scripts/controlar_perturbacao_inicial.py` define uma fonte volumétrica de vapor no perfil nativo GOAMAZON de 10 m:

`s(z) = rho_d0(z) × max[r_sat(p,T) − rv0(z), 0]`.

As fórmulas de saturação são as mesmas da perturbação original: água acima de −20 °C e gelo abaixo. A função s é interpolada linearmente entre os níveis nativos. Seu integral é calculado apenas no intervalo físico especificado, antes de escolher a grade do modelo. Cada célula recebe o integral da mesma função na interseção entre suas interfaces e a camada. A razão de mistura acrescentada é esse aporte dividido pela massa de ar seco de referência daquela célula, calculada do estado inicial do driver.

| Camada física | Aporte fixado pelo perfil nativo | Fonte latente declarada Lv × aporte |
| --- | ---: | ---: |
| 500–1500 m | 2,744294626875 mm | 6,860736567188 MJ/m² |
| 3500–5500 m | 2,123216322894 mm | 5,308040807236 MJ/m² |
| 7500–9500 m | 0,660130647472 mm | 1,650326618680 MJ/m² |

Foram verificadas as três camadas em 60, 120 e 240 níveis: nove testes de inicialização, todos aprovados. A diferença do aporte aplicado entre grades ficou abaixo de 10⁻¹¹ mm. A soma das sobreposições corresponde a 1.000, 2.000 e 2.000 m, respectivamente. Esses valores não foram ajustados aos resultados da grade de 60 níveis; vêm de uma referência comum anterior à discretização do modelo.

Isso controla o suporte físico da fonte em sentido de volume finito. A célula que cruza uma borda recebe uma média da contribuição parcial; o modelo não resolve a descontinuidade subcelular. Portanto, não se afirma identidade ponto a ponto nem RH final exatamente igual a 100% em todos os centros. Foi priorizada a mesma distribuição física projetada e a mesma massa equivalente acrescentada. O perfil de fundo ainda tem erros de regradeamento próprios, que deverão ser medidos na etapa 3.

Temperatura, qc e qi iniciais permanecem iguais aos perfis originais regradeados. O vapor é uma fonte externa, acompanhada de Lv × ΔM como diagnóstico energético latente. Não houve compensação de temperatura e não se declara conservação de energia. Esse termo não é a entalpia completa da fonte: uma formulação energeticamente fechada exigiria especificar o estado do reservatório de água e o trabalho/troca com a coluna. O orçamento adota a massa de ar seco anterior à perturbação como referência; a mudança de inventário Euleriano por ajuste de densidade deve ser contabilizada separadamente.

O novo campo opcional `path_initial_perturbation` lê o arquivo de projeção no driver. Quando vazio, permanece o caminho histórico de `saturate_layer`. Uma integração de regressão de Thompson médio confirmou igualdade bit a bit de todas as variáveis com a execução anterior. Não foram executadas integrações científicas com a nova perturbação: os testes controlados têm duração zero. A regressão usa o caminho antigo por duas horas. Microfísica, EP, rotinas de saturação e pulso não foram alterados.

Os arquivos imutáveis desta execução, logs de compilação, projeções e resultados estão na pasta datada indicada em `latest.json`. Reproduzir com `python scripts/controlar_perturbacao_inicial.py`; o programa cria outra pasta e preserva a anterior.

## Decisão adotada para a etapa 2

O pulso varia com a altura. Com pressão/densidade prescritas, uma coluna sem continuidade do ar não fornece, por si só, um orçamento fechado de água. Existem dois desenhos distintos:

1. **Experimento prescrito com fontes explícitas.** Manter o desenho idealizado e registrar armazenamento, advecção, fluxos de contorno, termo associado à divergência do fluxo de ar, sedimentação/precipitação e mudanças de inventário por densidade. O orçamento fecha com fontes explicitamente reconhecidas, mas a coluna não é isolada. Isso preserva a comparação atual; não se deve chamar o termo de compensação horizontal de observação GOAMAZON.
2. **Novo experimento com continuidade.** Introduzir massa de ar consistente, transporte em forma de fluxo e contornos definidos, com uma hipótese coerente para pressão/temperatura e compensação horizontal do pulso. Apenas trocar `−w ∂r/∂z` por `−∂(rho w r)/∂z` mantendo todas as outras hipóteses não resolve automaticamente a consistência do modelo. Essa escolha altera o problema dinâmico e exige novos testes de positividade, massa e energia.

O usuário escolheu explicitamente a opção 1: manter o experimento prescrito e explicitar fontes, sumidouros e fluxos de contorno. O orçamento foi instrumentado sem mudar o transporte ou as parametrizações. Quatro comparações com diagnóstico ligado/desligado foram idênticas bit a bit. A matriz controlada de 36 combinações foi executada após essa validação; os resíduos microfísicos permanecem explícitos. Consulte `../orcamento_prescrito/RELATORIO_ORCAMENTO_E_RESOLUCAO.md` e os manifests `latest.json` de `orcamento_prescrito` e `resolucao_controlada`. Os manifests históricos de inicialização registram o estado anterior à decisão e foram preservados.

## Critérios para as etapas seguintes

**2 — Orçamento:** registrar estados antes/depois de transporte, antes/depois de microfísica e os fluxos de precipitação no mesmo passo. Identificar densidades e unidades: Thompson calcula densidade a partir de p, T e qv dentro de sua rotina; WSM6 recebe uma densidade do driver. Não misturar indiscriminadamente ponderações de ar seco com densidade aproximada. Mostrar explicitamente qualquer termo não explicado; um fechamento por definição algébrica não prova conservação física.

**3 — Resolução:** repetir a matriz somente após a escolha do orçamento e sua validação. Reutilizar a mesma fonte física projetada e sua dose. Separar sensibilidade temporal de vertical, medir mudanças do fundo regradeado e reconhecer que Δt afeta transporte, microfísica e amostragem do EP. Não declarar convergência apenas pela ocorrência de eletrificação ou pela estabilidade da chuva.

**4 — Tendências microfísicas:** registrar taxas efetivamente aplicadas após limitadores, incluindo subpassos internos, mudanças diretas de estado e sedimentação. O levantamento inicial identificou `prw_vcd`, `pri_ide`, `prs_sde`, `prg_gde`, `pri_wfz`, `pri_rfz`, `prs_scw`, `prg_scw` e `prg_gcw` no Thompson; no WSM6, `pcond`, `pidep`, `psdep`, `pgdep`, `psacw`, `pgacw` e congelamentos aplicados diretamente por `pfrzdtc/pfrzdtr`. Isso é um mapa preliminar, não uma equivalência validada entre categorias. Unidades, sinais, fatores de densidade e limitadores precisam ser auditados antes de somar taxas. Riming é uma modalidade de coleta: não somar ambos como categorias independentes e duplicar transferências. Exigir fechamento por espécie e regressão bit a bit com diagnóstico ligado/desligado, sem alterar as equações físicas.

**5 — EP:** testar a fórmula separadamente sobre estados congelados, sem realimentação na microfísica: variação de w, máscara térmica por centros/intervalos e parâmetros de conversão. A máscara atual inclui um intervalo quando pelo menos um extremo pertence à faixa térmica. A lei de conversão EP→flashes continua nominal; sem observações compatíveis, a etapa pode avaliar sensibilidade e estrutura, mas não realizar calibração observacional defensável.
