# Resultados — microfísica e diagnóstico elétrico no SCM

Relatório gerado dos seis NetCDF locais. Commit: `f4096d1b4526ffac9a2167fc852d5201120245fa`. EP original, sem recalibração; nenhuma alteração nas parametrizações físicas ou nos arquivos de entrada. TH-CF e WS-CF foram reutilizados; TH-C0, TH-CL, WS-C0 e WS-CL foram executados, nessa ordem. Os seis passaram nos testes de integridade e consistência, com ressalvas físicas documentadas. Os valores abaixo são resultados computados, salvo identificação explícita de referência externa.

| Experimento | Chuva final (mm) | Chuva 35 h (mm) | qg máximo (g/kg) | Massa qg mista máx. (kg/m²) | Taxa elétrica máx. (min⁻¹) | Integral nominal |
| --- | --- | --- | --- | --- | --- | --- |
| TH-C0 | 1,6884 | 1,6440 | 0,000000 | 0,000000 | 0,000000 | 0,000000 |
| TH-CL | 0,4764 | 0,4459 | 0,000000 | 0,000000 | 0,000000 | 0,000000 |
| TH-CF | 0,9226 | 0,8869 | 0,000000 | 0,000000 | 0,000000 | 0,000000 |
| WS-C0 | 3,5270 | 3,5240 | 0,023385 | 0,005090 | 0,017634 | 3,910802 |
| WS-CL | 2,1320 | 2,1270 | 0,022514 | 0,005250 | 0,021117 | 3,577015 |
| WS-CF | 2,7068 | 2,7000 | 0,025594 | 0,005656 | 0,021313 | 2,939459 |

## 3.1. Reprodução do caso amazônico

Cada arquivo possui 4.232 instantes e 60 camadas, com passo de 30 s. O último registro é 126.930 s (35,258333 h), e não 35 h. Os estados iniciais e os eixos são idênticos entre os seis experimentos. Os CF reproduzem as referências fornecidas pelo usuário dentro de meia unidade da última casa decimal informada: TH-CF 0,922617 mm; WS-CF 2,706833 mm, graupel máximo 2,559396 × 10⁻⁵ kg/kg, taxa máxima 0,02131269 min⁻¹ e integral 2,939459.

A comparação seguinte usa a tabela da seção 7.4.2 do relatório técnico local de Kubota (`RPQ_driver_offline_microfisica_MONAN_MPAS.pdf`). São referências arredondadas de outra execução, não observações. A referência de aproximadamente 35 h não garante um instante final exatamente comum; por isso as duas janelas atuais são mostradas.

| Caso | Relatório (mm, aproximado) | Atual final (mm) | Atual 35 h (mm) | Diferença final (%) |
| --- | --- | --- | --- | --- |
| TH-C0 | 1,69 | 1,688433 | 1,644034 | -0,09 |
| TH-CL | 0,49 | 0,476352 | 0,445877 | -2,79 |
| TH-CF | 0,90 | 0,922617 | 0,886937 | 2,51 |
| WS-C0 | 3,19 | 3,527034 | 3,524015 | 10,57 |
| WS-CL | 2,29 | 2,131987 | 2,126974 | -6,90 |
| WS-CF | 2,72 | 2,706833 | 2,700031 | -0,48 |

A maior discrepância de chuva é WS-C0 (10,57%). Em 35 h ainda há 3,524015 mm: a diferença não é explicada apenas pelos 15,5 minutos adicionais. O máximo de graupel WS-C0 atual é 2,338478 × 10⁻⁵ kg/kg, contra aproximadamente 1,1 × 10⁻⁵ no relatório; WS-CL e WS-CF atuais são 2,251424 × 10⁻⁵ e 2,559396 × 10⁻⁵, ante 2,2 × 10⁻⁵ e 2,5 × 10⁻⁵. Não houve ajuste para forçar concordância. Versão exata das fontes, entradas, ambiente e detalhes dos namelists históricos são possíveis explicações a investigar, não causas demonstradas. Não há evidência suficiente para atribuir a diferença somente a arredondamento ou ao compilador.

Os metadados de SCM e DEF declaram latitude −8°, longitude −63°, início em 05/10/2014 às 12:00 e fim em 06/10/2014 às 00:00. Contudo, o eixo numérico do SCM termina em 126.921,09375 s (35,255859 h), incompatível com as 12 h desses atributos. O driver avança até o primeiro passo que ultrapassa o último instante da forçante, 8,90625 s além dele. A referência temporal textual, se aplicada literalmente, levaria o último registro do modelo a 06/10 às 23:15:30; essa conversão não é adotada como data observacional comprovada.

O gerador local `era5_forcing/DEPHY-SCM/GOAMAZON/SINGLE/driver_DEF.py` contém as coordenadas e datas declaradas e converte o tempo por `(timeForc − 278.5) × 86400`; isso identifica a transformação, mas não comprova a geolocalização do perfil original. O nome GOAMAZON não demonstra que a coluna seja Manaus ou Manacapuru. A documentação pública do [produto de forçantes ARM GOAmazon](https://portal.nersc.gov/cfs/capt/ARMVAP/GOAmazon.html) descreve uma média regional em torno do radar de Manaus. É necessário reconciliar esse domínio com os arquivos locais antes de uma comparação espacial. A Figura 1 apresenta os perfis e as forçantes efetivamente disponíveis, sem mapa de localização presumida.

Não há NaN, infinito, razão de mistura negativa ou inconsistência entre chuva acumulada e incremental. A pressão decresce com a altitude. A RH inicial sobre água, diagnosticada pela fórmula de Bolton nas células acima de 0 °C, alcança 101,12655% em uma célula; depois, o maior valor entre os seis casos é 100,14484%. A supersaturação inicial não foi corrigida. A camada artificial acima de aproximadamente 16,5 km produz temperaturas próximas de 161,68 K; esse extremo limita a interpretação da alta coluna e não constitui validação atmosférica. Detalhes em [relatório de qualidade](validacao/RELATORIO_QUALIDADE.md).

## 3.2. Diferenças na partição dos hidrometeoros

| Espécie | TH-CF máx. (g/kg) | WS-CF máx. (g/kg) | TH-CF hora / altura (km) | WS-CF hora / altura (km) |
| --- | --- | --- | --- | --- |
| qc | 1,886987 | 0,967120 | 27,558 / 3,171 | 30,692 / 4,383 |
| qr | 0,055632 | 0,136655 | 33,942 / 3,618 | 31,083 / 1,778 |
| qi | 0,000964 | 0,174107 | 30,217 / 14,030 | 31,225 / 8,960 |
| qs | 0,189586 | 0,188451 | 33,775 / 4,963 | 30,808 / 4,383 |
| qg | 0,000000 | 0,025594 | ausente | 30,875 / 4,383 |

No CF, Thompson acumula mais água de nuvem e muito menos gelo de nuvem que WSM6 nos máximos globais. Os máximos de neve são próximos, mas essa semelhança pontual não implica igual massa vertical nem evolução temporal. A Figura 2 mostra a distribuição tempo–altura e as isotermas; o gelo de Thompson é pequeno, porém não nulo. O graupel é exatamente nulo em todos os TH, enquanto os três WS têm graupel na região de −40 a 0 °C.

O máximo global de graupel WS-CF (0,025594 g/kg) difere do máximo restrito à fase mista (0,025153 g/kg); não são medidas intercambiáveis. O máximo global ocorre a 4,383 km. A massa mista máxima cresce de 0,005090 para 0,005250 e 0,005656 kg/m² entre C0, CL e CF. Nenhuma dessas medidas é fluxo de gelo.

Método: `T = theta × (p/100000)^(287/1004.5)`, com constantes do driver. A densidade é `p/[287 T (1 + 0.61 qv)]`. Como o arquivo fornece centros aritméticos de camada, as interfaces são recuperadas por `zw[0]=0; zw[k+1]=2 lev[k]−zw[k]`. Trata-se da inversão exata da grade usada, com topo verificado a 10⁻⁷ m; as espessuras variam de 50 a 1.054,82 m. O conteúdo é `Σ rho q dz`, em kg/m², selecionando a camada inteira quando a temperatura de seu centro está entre 233,15 e 273,15 K. A aproximação está na quadratura por células, especialmente nas fronteiras térmicas, e na densidade calculada do estado salvo. A densidade interna usada pelo EP é calculada em outra etapa do passo; não se afirma reprodução exata de seu integrando a partir da saída.

## 3.3. Resposta da atividade elétrica

Todos os TH apresentam taxa e integral nulas; seus horários de início e máximo e suas correlações elétricas são ausentes, não zero. Isso, isoladamente, não provaria funcionamento da microfísica. A evolução de hidrometeoros, a chuva não nula, a inspeção do caminho de chamada no código e os diagnósticos positivos dos WS sustentam que as rotinas foram exercitadas. Não foram adicionados contadores ou orçamentos de processos.

| Caso | Início (h) | Máximo (h) | Duração ativa na janela (h) | Último positivo (h) | Correlação EP × massa qg |
| --- | --- | --- | --- | --- | --- |
| WS-C0 | 30,117 | 32,025 | 5,150 | 35,258 | 0,959 |
| WS-CL | 30,383 | 31,892 | 4,883 | 35,258 | 0,948 |
| WS-CF | 30,233 | 30,808 | 5,033 | 35,258 | 0,877 |

Os três WS permanecem eletricamente ativos no último instante. Suas durações são tempos ativos dentro da janela e estão censuradas pelo término da integração; o encerramento completo dos eventos é desconhecido. A integral segue exatamente `Σ flash_rate[1:] × diff(time)/60`, por retângulos à direita, e representa descargas nominais, sem calibração observacional. A divisão IC/CG é fixa e não é previsão independente.

Em WS-CF, o máximo elétrico ocorre em 30,808 h, o da massa de graupel em 30,850 h e o da taxa de chuva em 31,000 h. O máximo elétrico antecede o de massa de graupel por 2,50 min e o de chuva por 11,50 min. São diferenças de horários de máximos, não evidência causal nem estimativas otimizadas por correlação cruzada.

A correlação de Pearson entre EP e massa mista de graupel é 0,959, 0,948 e 0,877 para WS-C0, CL e CF; com o máximo pontual misto, é 0,960, 0,942 e 0,860. Foram usados todos os instantes da janela comum, inclusive zeros, sem testes de significância ou independência temporal. O EP depende diretamente de hidrometeoros, temperatura e velocidade vertical prescrita; essas correlações internas não validam o diagnóstico.

Há coexistência positiva de qc, qs e qg na mesma célula de fase mista durante 5,150 h (WS-C0), 3,458 h (WS-CL) e 2,983 h (WS-CF). Isso corresponde a 100%, 70,8% e 59,3% do tempo eletricamente ativo, respectivamente. Positividade estrita inclui traços pequenos. A taxa elétrica não exige que esse critério diagnóstico simplificado de coexistência seja verdadeiro em todos os passos. Não se identificam taxas individuais de riming, deposição ou coleta a partir dessas saídas.

## 3.4. Sensibilidade ao acoplamento com o MYNN

C0 mantém MYNN ativo: apenas suas tendências de condensado deixam de ser aplicadas. CL aplica a parcela local; CF inclui também a contribuição das plumas. As mudanças seguintes são respostas da simulação a esses controles, não decomposições exatas das taxas dos processos microfísicos.

| Contraste (A − B) | Δ chuva (mm) | Δ chuva (%) | Δ taxa elétrica máx. (min⁻¹) | Δ integral nominal |
| --- | --- | --- | --- | --- |
| WS_C0 - TH_C0 | 1,838602 | 108,89 | 0,017634 | 3,910802 |
| WS_CL - TH_CL | 1,655635 | 347,57 | 0,021117 | 3,577015 |
| WS_CF - TH_CF | 1,784216 | 193,39 | 0,021313 | 2,939459 |
| TH_CL - TH_C0 | -1,212081 | -71,79 | 0,000000 | 0,000000 |
| WS_CL - WS_C0 | -1,395047 | -39,55 | 0,003483 | -0,333787 |
| TH_CF - TH_CL | 0,446265 | 93,68 | 0,000000 | 0,000000 |
| WS_CF - WS_CL | 0,574846 | 26,96 | 0,000196 | -0,637556 |
| TH_CF - TH_C0 | -0,765816 | -45,36 | 0,000000 | 0,000000 |
| WS_CF - WS_C0 | -0,820201 | -23,25 | 0,003678 | -0,971343 |

A mistura local reduz a chuva em 71,79% no Thompson e 39,55% no WSM6. O acréscimo das plumas aumenta a chuva em relação a CL em 93,68% e 26,96%, respectivamente, sem recuperar o acumulado C0. O efeito total CF − C0 é -45,36% no Thompson e -23,25% no WSM6.

No WSM6, CF aumenta o pico elétrico em 20,86% em relação a C0, mas diminui a integral em 24,84%. CL também eleva o pico e reduz a integral. Portanto, “aumentar a atividade elétrica” depende da métrica: pico e atividade acumulada respondem de modo distinto. A passagem CL → CF aumenta a chuva, mas reduz a integral elétrica em 17,82%. Chuva e eletrificação não variam sempre no mesmo sentido.

O contraste de diferenças `(WS-CL − WS-C0) − (TH-CL − TH-C0)` para chuva é -0,182967 mm; o correspondente com CF é -0,054385 mm. A diferença absoluta entre esquemas muda com o acoplamento. Isso evidencia interação na resposta deste caso; não autoriza assumir linearidade nem extrapolar uma sensibilidade universal. Percentuais elétricos com base Thompson nula são indefinidos no CSV.

## 3.5. Avaliação observacional

Não foi realizada validação observacional e não foi produzida Figura 5. Foram identificadas fontes candidatas de precipitação e eletricidade de 2014, mas não foi estabelecida correspondência espacial e temporal com esta entrada. A pendência é a compatibilidade do caso, e não uma afirmação de inexistência de observações. O [inventário de fontes e critérios](observacoes_pendentes.md) registra cobertura, acesso e diferenças entre strokes CG e flashes totais. Nenhum valor de outro dia ou região foi usado para validar estas simulações.

## 4. Discussão

Neste caso prescrito, a escolha da microfísica altera a partição condensada e a ativação do EP: só WSM6 gera graupel e diagnóstico elétrico positivo. A conclusão é sobre a configuração testada, não sobre incapacidade geral de Thompson de produzir graupel, nem superioridade observacional de WSM6. Sem orçamentos microfísicos não é possível atribuir as diferenças especificamente a riming ou autoconversão; explicações desse tipo permanecem hipóteses.

O maior pico de graupel ou de EP não determina sozinho a integral elétrica: WS-CF tem os maiores picos, mas WS-C0 a maior integral. A persistência, a distribuição vertical e a combinação dos campos entram na resposta. Os máximos e as correlações disponíveis não permitem concluir que um desses fatores predomina causalmente. A Figura 3 e as séries completas permitem examinar essa distinção sem reduzir a análise a um único máximo.

O modelo de coluna única não representa a organização horizontal, a propagação e a advecção lateral interativa de tempestades. A ascensão de grande escala é prescrita. Não há ensaio de sensibilidade às forçantes nem estimativa de incerteza entre casos. Os fluxos hfss/hfls exibidos na Figura 1 são campos disponíveis no arquivo; neste caminho, o MYNN calcula fluxos de superfície a partir da temperatura de pele e das variáveis da coluna. Sua exibição não significa imposição direta daqueles dois fluxos. Também não se inferiu uma solução radiativa apenas de atributos da entrada.

As magnitudes elétricas são nominais. Os coeficientes e a partição IC/CG não foram recalibrados; a área observacional de uma tempestade não está automaticamente associada à coluna. A discordância dos metadados, a camada superior artificial e a ausência de observações pareadas limitam a representatividade amazônica que se pode defender. Os resultados sustentam uma comparação controlada sob a mesma entrada, não uma previsão validada de descargas em Manaus.

A extensão ERA5 foi mantida separada. Um teste de duas horas do `ensemble_driver` terminou com código zero, mas reprovou a consistência vertical da pressão e a semântica da umidade: umidade específica foi passada como razão de mistura, com diferença de aproximadamente −1,735% no nível examinado. Os 12 namelists amazônicos foram preparados, mas os membros não foram promovidos a resultados científicos. Não foram produzidas Tabela 4 ou Figura 6. Consulte [pendência ERA5](ensemble_amazonia/PENDENCIAS.md). Não se alterou Fortran para contornar o problema.

## 5. Conclusões preliminares

1. WSM6 produz mais chuva que Thompson em C0, CL e CF: diferenças de 1,838602, 1,655635 e 1,784216 mm, respectivamente.
2. Os três WS têm graupel misto e EP positivo; os três TH têm graupel e EP nulos, apesar de chuva e evolução condensada não nulas.
3. A entrega completa de condensado reduz a chuva ante C0 em ambos os esquemas. No WSM6, eleva o pico elétrico, mas reduz a integral nominal de 3,910802 para 2,939459.
4. Chuva e eletrificação podem divergir: CL → CF aumenta a chuva WSM6 e diminui sua integral elétrica. O efeito do acoplamento depende da microfísica e da métrica.
5. Picos pontuais não bastam para explicar a atividade acumulada. Os dados são compatíveis com influência da persistência e estrutura de fase mista, mas não identificam um mecanismo causal predominante.
6. As seis integrações estão numericamente consistentes com as verificações implementadas. Validação observacional, reconciliação da procedência e generalização para outros casos continuam pendentes; não se declara um esquema superior.

## Legendas científicas das figuras

**Figura 1.** Perfis iniciais e forçantes do arquivo GOAMAZON/SINGLE usado. Temperatura e vapor são apresentados na grade de 60 camadas; a faixa cinza indica a camada superior artificial. Temperatura de pele, fluxos disponíveis, ascensão e tendência advectiva são mostrados nos tempos nativos da entrada. hfss/hfls são referências do arquivo, não fluxos diretamente impostos neste caminho MYNN. Tempo relativo; coordenadas declaradas sem validação geográfica.

**Figura 2.** Razões de mistura de água de nuvem, gelo de nuvem, neve e graupel em TH-CF e WS-CF (g/kg), com isotermas de −40, −15 e 0 °C. Cada espécie tem a mesma escala nos dois esquemas e normalização por raiz quadrada. A ausência de graupel Thompson é explicitada; o pequeno gelo Thompson não é zero.

**Figura 3.** Séries TH-CF e WS-CF de massas de hidrometeoros na fase mista, taxa elétrica nominal, taxa de precipitação e acumulado. Marcas verticais indicam máximos de massa de graupel, taxa elétrica e taxa de chuva quando definidos. Massa em kg/m² calculada por células entre −40 e 0 °C. A versão suplementar inclui os seis casos; seu segundo painel mostra apenas neve para preservar legibilidade. Atividade elétrica no fim da janela é censurada.

**Figura 4.** Sensibilidade de chuva acumulada, massa mista máxima de graupel, taxa elétrica máxima e integral elétrica nominal à microfísica e ao condensado MYNN. C0 retém MYNN, mas não aplica condensado; CL aplica a contribuição local; CF acrescenta plumas. Barras em escala linear, zeros explícitos; integração comum de 35,258333 h. Valores elétricos não calibrados, sem barras de incerteza porque se trata de um único caso determinístico por configuração.

## Arquivos e rastreabilidade

[Tabela 1](tabelas/Tabela1_configuracoes.csv), [Tabela 2](tabelas/Tabela2_resultados.csv), [Tabela 3](tabelas/Tabela3_contrastes.csv) e [definições das métricas](tabelas/metadados.json). Séries de precisão completa em `series/`; figuras PNG 300 dpi e PDF em `../figuras_artigo/`. Os mapas de calor são rasterizados dentro dos PDFs; eixos, texto e curvas permanecem vetoriais.

O atributo global `source` não identifica o esquema: o escritor menciona Thompson também em WSM6. Identificação baseada nos namelists e manifests. Os arquivos CF têm limite de proveniência retrospectiva: os hashes foram registrados após sua execução original, não são prova contemporânea do ambiente daquela execução. As quatro novas execuções possuem registros de saída, código de retorno e hashes. O [manifest de auditoria](proveniencia/auditoria_final.json) verifica preservação dos seis NetCDF, fontes e entrada contra o registro desta tarefa; as parametrizações e o driver principal permanecem preservados; o ensemble recebeu adaptação autorizada de entrada e transporte numérico. A pequena correção anterior do Makefile registra dependências de módulos para compilação paralela, sem alteração física. Os logs vazios de execuções bem-sucedidas foram preservados; código de retorno e validação não dependem de mensagens de progresso.

O procedimento reproduzível e suas dependências estão em `../README_ANALISE_ARTIGO.md`. Este relatório é regenerado por `scripts/gerar_relatorio_artigo.py`; tabelas e figuras derivam dos NetCDF, e os números externos de comparação são identificados explicitamente.

## Atualização: entrada original e ensemble GOAMAZON

A entrada ativa está em `data_goamazon/GOAMAZON_SINGLE_SCM_driver.nc`, com o mesmo SHA-256 do arquivo utilizado na matriz principal. A pasta antiga `data` foi removida por solicitação do usuário. Namelists e manifests históricos preservam os caminhos da época; novas execuções usam somente o caminho atual.

A extensão agora usa exclusivamente o perfil inicial GOAMAZON e wa inicial fixa, com as mesmas perturbações de saturação e pulso do ensemble idealizado anterior. Os 12 membros foram executados e validados; resultados separados na [Tabela 4](tabelas/Tabela4_ensemble_goamazon.csv), Figura 6 e [relatório do ensemble](ensemble_goamazon/RELATORIO_ENSEMBLE_GOAMAZON.md). O teste ERA5 reprovado é apenas histórico. Para o caminho GOAMAZON foi necessário substituir a advecção centrada explícita instável por transporte upwind com subpassos CFL, sem alterar microfísica, EP ou perturbações. A pressão/wa fixas e a ausência de MYNN e forçantes horizontais distinguem esta extensão da matriz principal.

### Ressalva posterior sobre a extensão idealizada

O [diagnóstico científico posterior](ensemble_goamazon/diagnostico_cientifico/RELATORIO_DIAGNOSTICO_CIENTIFICO.md) quantifica a água acrescentada e o ganho aparente pelo transporte, compara os perfis e avalia 36 configurações de resolução. A convergência quantitativa de graupel/EP não está demonstrada; a Figura 6 permanece exploratória. A ocorrência qualitativa de eletrificação é robusta nas resoluções ensaiadas, mas aprovação numérica não implica conservação, convergência ou reprodução do evento.
