# Single-column microphysics and lightning — LEONA

Fluxo principal: tempestade de 13–14/12/2018. Duas perguntas operacionais:
1. Desempenho fora do ajuste das parametrizações de raios totais calibradas.
2. Sensibilidade às configurações Thompson/WSM6 × C0/CL/CF, com coeficiente de referência fixo.

## Estrutura

- `leona/`: aquisição, contratos, DEF/SCM, eventos/flashes, proxies, ensemble e avaliação.
- `configs/leona.json`: período UTC, referência TH_C0, divisão contígua, sementes e pendências.
- `configs/parameterizations.json`: equações, unidades, coeficientes e requisitos publicados.
- `tests/`: testes sintéticos de contratos e integração DEPHY; nunca representam observações.
- `data_private/leona/`: originais copiados com SHA256, cache e sondagens processadas; ignorado pelo Git.
- `products/leona/`: tabelas, figuras, manifests e resultados de testes locais; ignorado pelo Git.
- `work/leona/`: build isolado, logs e execuções de verificação; ignorado pelo Git.
- `src/`, `era5_forcing/DEPHY-SCM/`: física e biblioteca existente preservadas.
- `docs/historical/README_pre_leona.md`, `scripts/`, `resultados/`, `figuras_artigo/`: reprodução histórica. Não são o fluxo LEONA.

## Ambiente e comandos

Python 3.11+, dependências em `requirements-leona.txt`; DEPHY importado de
`era5_forcing/DEPHY-SCM` (commit e caminho importado registrados no manifest).
Fortran: gfortran + NetCDF-Fortran (`nf-config`) + make. No Windows, executar
compilação e simulações no WSL; scripts Python podem rodar no Windows.

```text
python -m leona.cli --help
python -m leona.cli audit
python -m leona.cli cache-inputs
python -m leona.cli soundings
python -m leona.cli prepare-era5
python -m leona.cli download-era5
python -m leona.cli sample-lightning
python -m leona.cli lightning
python -m pytest tests -q
make BUILD=work/leona/build work/leona/build/mp_offline_driver
python scripts/check_leona_driver.py
python -m leona.cli status
```

Configuração alternativa: `python -m leona.cli --config caminho.json COMANDO`.
CDS exige credenciais locais e termos aceitos; nenhum segredo é incluído no código.
As requisições são diárias, cacheadas, verificadas e retomáveis; falhas são inventariadas.
No Wyoming, `sounding_json` é consultado por horário antes de `/wsgi/sounding`.
Ausência do inventário difere de erro do servidor e de perfil parcial.
IGRA/acervos de campanha são alternativas identificadas, sem substituição silenciosa.

## Contratos e bloqueios de execução

C0 mantém MYNN termodinâmico e não aplica suas tendências de condensado; CL aplica
apenas a parcela local; CF aplica a tendência completa. Todos usam o driver principal,
MYNN original (`mynn_edmf_type=0`) e forçantes temporais completas. Sem pulsos verticais,
saturação artificial ou microfísica adicional de pluma.

A área exige seleção meteorológica documentada, anterior à avaliação dos raios.
O caso usa uma caixa de 1 grau selecionada pelo máximo da ascendência positiva média
ERA5 em 500/700 hPa no período completo, dentro do setor visual de atividade.
`products/leona/support/selection.json` registra critério, dados, hashes e instante da
fixação anterior aos erros observacionais. Isso não garante representatividade convectiva
da coluna. Anillaco é a estação óptica; o setor vermelho é seu campo de visão.
Registrar `area=[west,south,east,north]`, justificar o suporte e fixar
`area_selection.status=fixed_before_lightning_evaluation` somente depois dessa seleção.
O recorte visual não é uma escolha automática. A área não participa do ajuste.

`observation_coverage` é uma lista booleana por intervalo de 15 min. Presença de arquivo
não certifica cobertura. Lacunas permanecem NaN e interrompem o acumulado conhecido.
A convenção 1=IC/0=CG e a cobertura completa dos dois dias foram confirmadas pelo usuário;
essas confirmações ficam registradas separadamente no JSON.

Agrupamento conjunto IC+CG, 700 ms e 10 km do primeiro evento detectado, limites
inclusivos, sem encadeamento ilimitado. Desempate: primeiro grupo elegível por tempo do
âncora e event_id. CG em qualquer evento classifica o grupo como CG. Duplicatas exatas
mantêm proveniência. Nanosegundos são inteiros; dias adjacentes são processados antes do
recorte. A variante Azambuja (1,4 s / intervalo consecutivo de 700 ms) é identificada à parte.

O leitor agora entrega `tntheta_adv` ao integrador de θ; a preparação deriva advecção
horizontal de r=q/(1-q), antes da média espacial. O driver acrescenta transporte vertical.
A pressão é prescrita; não se alega equilíbrio hidrostático exato. Forçantes acima da
superfície e do topo e fora do período não são extrapoladas. A superfície usa 2 m/10 m
mapeados ao nível 0 m AGL, decisão registrada no manifest. `rainnc` já é o equivalente
total de água; nunca somar `snownc`.

`mf_w` exporta sum(rho_interface*a*w)/sum(a), nas interfaces: contém densidade e não é
velocidade convectiva resolvida. Não é usado nos proxies que exigem velocidade convectiva
colocalizada com gelo/graupel. O método I×G usa concentrações mássicas integradas, com
r em base seca e densidade de ar seco. EP e partição fixa IC/CG ficam desativados; as
saídas legadas de raios contêm NaN. Para reproduzir a interface antiga, definir
`legacy_temperature_tendency=.true., legacy_lightning=.true.` e reutilizar a entrada
histórica original (não regenerada pelo preparador corrigido).

## Execução do caso e produtos

`leona.forcing.build` e o comando `forcing PL SL OUTPUT` reutilizam Case e convert2SCM.
`leona.experiments.run_controls` executa seis controles, valida finitude/positividade/tempo,
temperatura/pressão, vapor à entrada da microfísica e precipitação. O balanço ponderado
pela densidade nativa de entrada é registrado. Thompson recalcula densidade internamente;
esse resíduo não é um contrato de fechamento exato. `configs/budget_contract.json` define
o escopo e os limites. `scripts/check_leona_driver.py` usa GOAMAZON para testes de software e
compara diagnósticos ligados/desligados; não certifica o caso LEONA.

`leona.ensemble.prepare` exige pré-teste físico e numérico válido, seis configurações e
checksum da mesma forçante. Recebe dez perfis EDA colocalizados (T, r), gera 20 amostras
sintéticas conjuntas de T/log(r), salva os sorteios e pesos realizados e usa o mesmo membro
nas seis configurações. Nenhuma perturbação de vento é contada. `stability` registra média
e dispersão por tamanho do conjunto; a dispersão é apenas de condição inicial.

Sequência do caso: `assemble_leona_era5.py`, comando `forcing`, comando `controls`,
`collocate-eda`, `prepare-ensemble` e `run-ensemble` (ver `python -m leona.cli --help`).
`process_leona_ensemble.py` processa os membros concluídos; `--evaluate` exige o conjunto
pareado completo. `check_leona_stability.py` usa critérios previamente registrados em
`configs/ensemble_stability.json`, sem observações. `figures_leona.py` produz as quatro
figuras principais. `process_leona_controls.py` usa os processos efetivamente aplicados,
normalizados pela massa seca inicial, explicitamente distintos de fluxos instantâneos.

`leona.postprocess` diagnostica saídas sem mutá-las e integra passos nativos em intervalos
completos. `evaluate PREDICTORS OBSERVATIONS OUTPUT` exige a área e cobertura fixadas.
Preditores CSV: time,method,configuration,member,predictor (integral no suporte do intervalo).
Observações CSV: time,count. Não misturar unidades de densidade, taxa e contagem.
Coeficientes são ajustados na média dos membros perturbados, ou no controle se não houver
ensemble; um único coeficiente por método/configuração, nunca por membro. O bloco reservado
não entra no ajuste. Os cortes iniciais são pré-especificados no JSON, sem minutos aleatórios.
Preditor nulo recebe estado irrecuperável, sem atividade artificial. Uma previsão cujo
preditor continua nulo é zero para toda escala finita; preditor posterior não nulo com
escala não identificável permanece indisponível. Dados ausentes permanecem NaN.

Saídas: coefficients, heldout_predictions, metrics, fixed_scale_skill_contrasts,
physical_trajectories e ensemble_statistics, todos em CSV, mais manifest de configuração e hashes.
O critério principal é RMSE de contagens; MAE, correlação, pico e acumulado são complementares.
Mudança física e melhoria de erro são campos distintos. `leona.figures` prepara as quatro
figuras principais somente quando os produtos correspondentes existem; não fabrica séries
nem gráficos de previsão quando faltam as entradas do caso.

Não há relatórios de interpretação científica no novo fluxo. Observações e seus derivados
permanecem locais, ignorados pelo Git e sujeitos às restrições do fornecedor.
