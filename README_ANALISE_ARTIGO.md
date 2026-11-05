# Análise dos seis experimentos GOAMAZON

## Figura 7: processos microfísicos aplicados

A instrumentação de massa por processo foi concluída para os 12 membros originais da Figura 6. Todas as variáveis permaneceram idênticas bit a bit às saídas arquivadas, com fechamento por espécie e testes de diagnóstico ligado/desligado. Figura: `figuras_artigo/Figura7_processos_microfisicos.png` e `.pdf`; legenda em `Figura7_legenda.txt`. Relatório, taxas por nível/tempo/processo, tabelas de fontes/perdas e proveniência: `resultados/ensemble_goamazon/processos/RELATORIO_FIGURA7.md` e pasta indicada em `processos/latest.json`.

Reprodução: `python scripts/executar_processos_figura7.py`; `python scripts/gerar_figura7_processos.py`. Validação: `python -m unittest discover -s scripts -p test_processos_figura7.py`. As equações originais são preservadas por blocos de instrumentação delimitados. A Figura 7 compara fontes diretas de gelo e graupel, não causalidade por desligamento de processos; não mistura a inicialização histórica com o estudo posterior de aporte controlado. A análise contrafactual separada de EP permanece pendente.

## Controle inicial de aporte entre grades

Foi implementado um caminho opcional de perturbação com aporte definido no perfil nativo GOAMAZON e projeção pelas sobreposições das células. Nove inicializações verificam o mesmo aporte em 60/120/240 níveis, e uma regressão confirma que o caminho antigo permanece idêntico bit a bit. O usuário escolheu manter a coluna prescrita com fontes e fluxos explícitos. O orçamento por etapas foi validado e a matriz controlada de 36 combinações executada, preservando as saídas anteriores. Consulte `resultados/ensemble_goamazon/orcamento_prescrito/RELATORIO_ORCAMENTO_E_RESOLUCAO.md`. O resíduo nativo da microfísica permanece explícito; fechamento contábil não prova conservação física. As taxas internas por processo foram instrumentadas e analisadas nos membros originais da Figura 6 (ver Figura 7 acima); a extensão ao ensemble de aporte controlado e a análise contrafactual separada de EP permanecem pendentes. Comandos: `python scripts/orcamento_coluna_prescrita.py`, `python scripts/resolucao_controlada.py` e `python scripts/relatar_orcamento_prescrito.py`. Verificação inicial: `python scripts/controlar_perturbacao_inicial.py`.

## Atualização científica: água, perfis e resolução do ensemble

O diagnóstico posterior está em `resultados/ensemble_goamazon/diagnostico_cientifico/RELATORIO_DIAGNOSTICO_CIENTIFICO.md`. Foram analisados os três membros eletrificados e Thompson médio como controle, com 32 novas integrações e quatro baselines reutilizados (60/120/240 níveis × 20/10/5 s). Nenhuma figura nova foi produzida. O aporte artificial é relevante e a convergência quantitativa de graupel/EP ainda não foi demonstrada; “aprovado” abaixo se refere apenas aos controles de integridade implementados.

Reprodução: `python scripts/diagnosticar_ensemble_goamazon.py`, `python scripts/relatar_diagnostico_ensemble.py`, `python scripts/test_diagnostico_ensemble.py`. CSVs e NetCDF de perfis documentam água/energia introduzidas, advecção reconstruída, evolução vertical e sensibilidade. O ensaio é idealizado com perfil inicial e wa fixa, não reprodução do evento GOAMAZON.

Os resultados científicos principais estão em `resultados/RELATORIO_RESULTADOS_ARTIGO.md`; os seis NetCDF estão em `resultados/`. TH-CF e WS-CF foram reutilizados sem alteração; TH-C0, TH-CL, WS-C0 e WS-CL foram executados nessa ordem. As Tabelas 1–3, Figuras 1–4 e Figura 3 suplementar foram produzidas dos arquivos locais. Todas as figuras têm PNG 300 dpi e PDF; os mapas de calor são rasterizados no PDF, com texto/eixos vetoriais.

As seis saídas passaram nos controles implementados, com ressalvas de RH inicial, camada superior artificial e procedência geográfica/temporal. Isso não significa validação observacional. EP e IC/CG não foram recalibrados. Parametrizações físicas, constantes e entradas permanecem inalteradas. O driver do ensemble recebeu suporte autorizado ao GOAMAZON e uma correção de transporte numérico; o driver principal não foi modificado. Não houve commit automático. O Makefile registra as dependências necessárias à compilação paralela e ao leitor GOAMAZON do ensemble.

## Ambiente e instalação

Execução realizada em Windows, Python de `C:\Users\teixe\anaconda3\python.exe`, com compilação e executável Linux no WSL `Ubuntu-24.04`. O script usa WSL automaticamente no Windows e execução direta no Linux. Dependências nativas: GNU Fortran, make, NetCDF C/Fortran e utilitários NetCDF. No Ubuntu, quando ausentes: `sudo apt install gfortran make libnetcdf-dev libnetcdff-dev netcdf-bin`. Os detalhes reais de compilador e `nf-config` estão em `resultados/proveniencia/execution_20261102T214208919408Z.json`.

Instale as dependências Python em ambiente próprio, a partir da raiz:

```bash
python -m venv .venv-analise
# Linux: source .venv-analise/bin/activate
# PowerShell: .\.venv-analise\Scripts\Activate.ps1
python -m pip install -r requirements-analysis.txt
```

O ambiente Anaconda utilizado tinha extensões opcionais antigas numexpr/bottleneck; pandas emitiu avisos e utilizou seus caminhos disponíveis. Elas não são dependências deste procedimento. As versões efetivamente usadas das cinco bibliotecas principais foram fixadas no arquivo de requisitos.

## Reproduzir o procedimento

Na raiz do repositório, com Python configurado:

```bash
python scripts/test_artigo.py
python scripts/executar_experimentos.py
python scripts/validar_resultados.py
python scripts/analisar_experimentos.py
python scripts/gerar_figuras_artigo.py
python scripts/gerar_relatorio_artigo.py
python scripts/executar_ensemble_goamazon.py
python scripts/test_ensemble_goamazon.py
```

No PowerShell, se `python` não apontar ao ambiente desejado, substitua-o por `& 'C:\Users\teixe\anaconda3\python.exe'`. As chamadas WSL podem precisar da permissão de execução do aplicativo hospedeiro.

O executor valida e reutiliza arquivos existentes. Não repete simulações válidas por conveniência. Resultados novos são executados em `resultados/runs/ID_DATA/`, validados e só então copiados para o nome canônico. Os logs ficam em `resultados/logs/`; namelists em `resultados/namelists/`; hashes, código de retorno, fontes e ambiente em `resultados/proveniencia/`. Código de retorno zero é necessário, mas não suficiente. Logs podem estar vazios em execuções normais.

```bash
# Reexecução explícita de um caso; o anterior é copiado para versions/
python scripts/executar_experimentos.py --experiments TH_CL --rerun
# Compilação limpa explícita; primeiro preserva TODO o build/ em versions/
python scripts/executar_experimentos.py --rebuild
```

Antes de novas simulações, um executável só é reutilizado quando seu hash e os hashes das fontes coincidem com um manifest conhecido; caso contrário ocorre compilação limpa. `--rebuild` usa o Makefile existente: `make clean` seguido de `make -j4`. Nunca execute manualmente `make clean` antes de preservar saídas próprias de `build/`. Objetos antigos não são evidência de compilação válida. O baseline deste trabalho é o manifest de 02/11/2026, preservado para auditoria. Após reexecução deliberada com mudança de dados/código, a auditoria final rejeita diferenças desse baseline: reveja cientificamente a mudança e estabeleça um novo baseline explicitamente, sem apagar o anterior.

## Scripts e métricas

- `scm_artigo.py`: leitura `decode_times=False`, constantes, grade, métricas, configurações e hashes.
- `validar_resultados.py`: variáveis/dimensões/unidades, finitude, sinais, pressão, chuva incremental, partição elétrica, grade, duração, estado inicial e identidade entre experimentos. Se há falha, não corrige dados nem a mascara com zeros.
- `analisar_experimentos.py`: Tabelas 1–3 e séries de alta precisão. Funciona com subconjunto de saídas, mantendo o namelist TH-CF de referência. Ausências têm status próprio e campos vazios; percentuais com base zero são indefinidos.
- `gerar_figuras_artigo.py`: figuras dos dados válidos disponíveis; painéis CF somente quando ambos existem. Ausência não é barra zero. Requer pelo menos uma saída válida.
- `gerar_relatorio_artigo.py`: narrativa e auditoria para os seis casos completos; recusa tabelas desatualizadas. A discussão é específica deste caso, não um classificador genérico de hipóteses para qualquer simulação futura.
- `executar_ensemble_goamazon.py`: executa/valida os 12 membros com dados originais GOAMAZON; gera Tabela 4, Figura 6 e relatório separado. Compila em pasta nova e preserva as execuções anteriores. `--rerun` cria outra versão.
- `avaliar_ensemble_amazonia.py`: diagnóstico ERA5 histórico, retirado do fluxo ativo; não é necessário executá-lo.

O parser de namelists foi feito para as atribuições escalares desta matriz. Não é um parser Fortran genérico de arrays ou strings com comentários embutidos. Os namelists multicoluna do ensemble de exemplo não são lidos por ele.

A chuva incremental é dividida pelo intervalo temporal correspondente. A taxa no primeiro instante é ausente. A integral elétrica é a soma retangular à direita, em descargas nominais. Critérios de início e atividade usam positividade estrita e podem incluir traços numéricos; início de chuva não é necessariamente início de chuva mensurável. Se a atividade persiste no último passo, a duração do evento é censurada.

Temperatura usa Rd=287 e cp=1004,5; densidade usa o estado salvo e a aproximação virtual do driver. A inversão das interfaces é exata para centros aritméticos com base no solo, mas a máscara de fase mista seleciona células inteiras pelos centros. Massa integrada é kg/m², não fluxo de gelo. Densidade recalculada da saída não é necessariamente a densidade interna do EP no mesmo passo. Correlações só existem para séries não constantes, sobre toda a janela; sem inferência causal ou teste de significância com amostras autocorrelacionadas. Unidades e definições completas em `resultados/tabelas/metadados.json`.

## Qualidade, testes e limitações

Os testes reproduzem os CF dentro do arredondamento informado pelo usuário, verificam zeros, fim censurado, não finitos, intervalos irregulares, arquivos ausentes, configurações incompatíveis e condições iniciais. Testes adicionais verificam a saída parcial e integram eletricidade de forma independente. Logs em `resultados/validacao/`.

`resultados/validacao/RELATORIO_QUALIDADE.md` separa RH inicial/evolução. `resultados/proveniencia/auditoria_final.json` registra as fontes, o PDF técnico local e a preservação contra o manifest da execução. Namelists e manifests identificam o esquema; o atributo global `source`, que pode citar Thompson em WSM6, não é usado para isso. Os CF anteriores têm proveniência parcialmente retrospectiva, explicitada no relatório.

O arquivo SCM declara 8°S, 63°W e uma janela textual de 12 h, mas contém 35,255859 h de forçantes. O modelo principal termina em 35,258333 h. Não houve mapa nem associação presumida com Manaus. `resultados/observacoes_pendentes.md` documenta fontes candidatas e a falta de pareamento; não existe Figura 5. O smoke ERA5 reprovado foi preservado como histórico; foi substituído no escopo pelo ensemble GOAMAZON abaixo.

## Ensemble GOAMAZON original — extensão concluída

A única entrada ativa é `data_goamazon/GOAMAZON_SINGLE_SCM_driver.nc`. Seu SHA-256 coincide com o registrado para a antiga cópia `data/`, removida a pedido do usuário. Os scripts usam o novo caminho; namelists históricos permanecem intactos como evidência da época. O validador aceita o caminho histórico somente quando o hash da entrada atual confirma equivalência.

Foram executados 12 membros: Thompson/WSM6 × saturação em 500–1500, 3500–5500 e 7500–9500 m × pulso adicional 0/3 m/s. Todos têm 2 h, passo 20 s e 60 níveis. Usa-se o perfil inicial GOAMAZON e wa inicial fixa, sem MYNN nem a sequência de tendências horizontais/superficiais. São as mesmas rotinas e amplitudes de perturbação do ensemble ERA5, aplicadas ao novo perfil; não se misturam esses resultados com C0/CL/CF.

O método centrado explícito de advecção do ensemble tornou-se instável no pulso forte. No caminho GOAMAZON, foi substituído por upwind de primeira ordem com subpassos CFL ≤0,8 e distâncias entre centros. A equação de advecção e as rotinas físicas não foram alteradas. Esse método é difusivo; não foi feito estudo de convergência de grade. O teste instável está preservado e excluído dos resultados aceitos. O caminho legado ERA5 não foi reparado nem validado por esta extensão.

Os 12 membros passaram em verificação individual, pareamento de condições iniciais, teste da perturbação e fechamento das integrais com o resumo Fortran. Dois testes adicionais verificam diretamente transporte/pulso Fortran e preservação das parametrizações. Resultados: `resultados/tabelas/Tabela4_ensemble_goamazon.csv`, `figuras_artigo/Figura6_ensemble_goamazon.png` e `.pdf`, `resultados/ensemble_goamazon/RELATORIO_ENSEMBLE_GOAMAZON.md`. Logs e manifests ficam nessa última pasta. Os seis casos principais não foram reexecutados nesta extensão.

Arquivos derivados (CSV, gráficos, relatório e qualidade) são regeneráveis e atualizados pelos scripts. Os NetCDF originais, execuções brutas e manifests históricos são preservados. Documentos e alterações preexistentes fora desse fluxo não foram removidos.
