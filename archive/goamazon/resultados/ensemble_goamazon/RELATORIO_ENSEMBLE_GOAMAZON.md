# Ensemble GOAMAZON — resultados

**Atualização científica:** o [diagnóstico de água, perfis e resolução](diagnostico_cientifico/RELATORIO_DIAGNOSTICO_CIENTIFICO.md) mostra aporte artificial relevante e convergência quantitativa ainda não demonstrada. A aprovação abaixo significa integridade numérica, não conservação ou convergência. A Figura 6 é exploratória.

Doze membros aprovados nos controles implementados. Dados originais em `data_goamazon/`; nenhuma entrada ERA5 é lida. As rotinas de saturação, pulso, microfísica e EP permanecem inalteradas.

## Desenho e interpretação

Perfil inicial (pa, theta, rv, rl, ri) e wa do primeiro registro do GOAMAZON. Pressão e wa de fundo fixas durante 2 h; theta e vapor evoluem pela advecção vertical e microfísica do ensemble original. Sem MYNN, tendências horizontais, atualização temporal das forçantes ou fluxos de superfície. Portanto, este é um experimento idealizado com o perfil GOAMAZON, distinto da matriz C0/CL/CF e da reprodução de 35,26 h.

Grade de 60 camadas, topo 20 km, primeira espessura 50 m, passo 20 s. Saturação aplicada uma vez nos centros entre 500–1500, 3500–5500 e 7500–9500 m, com RH alvo 1. A rotina original usa saturação sobre água acima de −20 °C e sobre gelo abaixo. Acrescenta vapor onde necessário, sem compensação energética: uma perturbação idealizada, não uma evolução fechada conservativa. Excesso preexistente vira condensado conforme a rotina original. As massas de água acrescentada e alturas efetivas estão nas validações individuais.

Pulso original: 0,5[1−cos(2πt/1800)] × exp[−0,5((z−zc)/1000)²] × amplitude, durante 30 min, centrado no meio de cada camada. Amplitude 0 ou 3 m/s. Pulso zero mantém wa inicial. A discretização vertical pode não amostrar exatamente o máximo de 3 m/s.

O transporte vertical no caminho GOAMAZON foi corrigido: upwind de primeira ordem com distâncias reais entre centros e subpassos CFL ≤0,8, mantendo a equação −w ∂φ/∂z. Contorno de entrada com gradiente nulo. O método centrado explícito anterior gerou temperaturas de 885 K antes de NaNs no teste de pulso baixo, em 1.260 s; essa tentativa foi preservada e excluída das tabelas. A mudança é numérica e introduz difusão de primeira ordem; não se afirma equivalência dos resultados com o transporte antigo nem convergência de resolução. As rotinas físicas e o pulso não foram alterados.

## Síntese quantitativa

| Experimento | Chuva (mm) | qg máx. (g/kg) | Taxa elétrica máx. (min⁻¹) | Integral nominal |
| --- | --- | --- | --- | --- |
| thompson_baixa_w0 | 0 | 0 | 0 | 0 |
| thompson_baixa_w3 | 8.35964 | 0 | 0 | 0 |
| thompson_media_w0 | 0 | 0 | 0 | 0 |
| thompson_media_w3 | 2.41077 | 0 | 0 | 0 |
| thompson_alta_w0 | 0 | 0 | 0 | 0 |
| thompson_alta_w3 | 0.243448 | 0.0166899 | 0.150707 | 10.2492 |
| wsm6_baixa_w0 | 0 | 0 | 0 | 0 |
| wsm6_baixa_w3 | 11.8454 | 0 | 0 | 0 |
| wsm6_media_w0 | 0 | 0 | 0 | 0 |
| wsm6_media_w3 | 6.73227 | 0.670557 | 0.136833 | 3.96488 |
| wsm6_alta_w0 | 0 | 0 | 0 | 0 |
| wsm6_alta_w3 | 0.793367 | 0.398394 | 0.265208 | 13.014 |

3 de 12 membros apresentam EP positivo. Taxas e integrais são nominais, sem calibração. Horários de séries nulas permanecem ausentes no CSV. Indicador de censura identifica atividade no último passo.

Máximos globais da Tabela 4 incluem o estado inicial perturbado; o resumo nativo do driver considera somente os passos posteriores. A validação compara esse resumo com o mesmo subconjunto, evitando uma falsa discrepância.

## Validação e proveniência

Verificados pressão estritamente decrescente, finitude, não negatividade, fechamento da chuva incremental, 361 registros em 2 h, estado inicial regradeado, saturação contra os polinômios RSLF/RSIF, identidade inicial entre esquemas/pulsos pareados e integral elétrica contra o resumo Fortran. Todos os membros usam exatamente o mesmo arquivo fonte. A camada superior artificial do GOAMAZON e a incerteza geográfica/temporal original permanecem limitações. Não há validação observacional.

Execuções e logs em subpastas individuais. `manifest.json` registra fontes, entrada, binário, compilação e hashes das saídas. O teste ERA5 antigo permanece apenas como histórico, fora destes resultados.

**Figura 6.** Chuva final, máximo da massa de graupel entre −40 e 0 °C, taxa elétrica máxima e integral nominal nos 12 membros. Cores distinguem microfísicas; cada camada aparece com pulso adicional 0/3 m/s. Escalas lineares e zeros explícitos; sem inferência de habilidade observacional.

Reprodução: `python scripts/executar_ensemble_goamazon.py`. Reutiliza somente execuções com hashes compatíveis; `--rerun` cria novas pastas e preserva as anteriores.
