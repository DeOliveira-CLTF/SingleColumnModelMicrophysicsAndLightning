# Ensemble ERA5 amazônico — preparado, não validado

Foram preparados 12 namelists: Thompson/WSM6 × saturação em 500–1500, 3500–5500 e 7500–9500 m × pulso adicional 0/3 m/s. Mantidos os demais parâmetros escalares do exemplo: 60 níveis, topo 20 km, primeira camada 50 m, passo 20 s, duração 2 h, largura vertical do pulso 1 km, duração 1.800 s e RH alvo 1. Os resultados futuros devem ficar separados da matriz MYNN principal. Este driver não usa MYNN nem suas forçantes temporais.

Há ERA5 local em `era5_forcing/era5_pl.nc` e `era5_forcing/era5_sl.nc`: 48 horários, grade 0,25°, latitude 6–10°S, longitude 61–65°W. O ponto amazônico do exemplo (3°S, 60°W) não está nessa cobertura. Foi escolhido um ponto de grade efetivamente disponível, 8°S, 63°W, em 05/10/2014 às 15:00 UTC (índice Fortran 16). Não se utilizou extrapolação para um ponto distante, nem se incorporou caso andino. Esse perfil ERA5 não é apresentado como equivalente ao GOAMAZON principal.

Um único teste de funcionamento, Thompson/camada baixa/pulso zero, terminou com código de retorno zero e gerou 361 instantes. Arquivos brutos, log e `summary.nc` estão na pasta `smoke_20261102T214542451445Z`. O resumo é NetCDF apesar da extensão `.csv` antiga do exemplo; os namelists preparados usam `.nc`. Pulso zero desativa somente a ascensão sintética adicional, não a velocidade derivada do ERA5.

O teste não passou no critério de aceitação científica:

- A pressão apresenta valores repetidos entre camadas, contrariando a exigência de decréscimo estrito. A interpolação/uso da altura geopotencial e o tratamento de níveis próximos à superfície precisam ser revistos. Esta é uma hipótese de origem; não se corrigiu a pressão.
- `era5_bridge.F90` lê `q` como umidade específica e atribui o valor diretamente a `col%qv`. O regradeamento do `ensemble_driver.F90` entrega esse campo como razão de mistura, sem conversão `r=q/(1−q)`. No nível mais baixo, fora da camada artificialmente saturada, entrada e saída são 0,017352733761 kg/kg, enquanto a razão de mistura correspondente é 0,017659168690 kg/kg: diferença de −1,73527%. A comparação é registrada em `status.json` com hashes.

Ausência de crash não resolve essas diferenças semânticas. Como sua correção exige revisão do caminho Fortran, proibida neste trabalho, os 12 membros não foram executados como resultados científicos e não foram geradas Tabela 4 ou Figura 6. O teste não foi incluído na Tabela 2. A advertência de camada artificial presente no registro original do smoke foi herdada do validador GOAMAZON; não é evidência de que o ERA5 tenha aquela camada. A versão atual do validador restringe essa advertência ao caso principal.

## Repetição e retomada

`python scripts/avaliar_ensemble_amazonia.py` reapresenta o diagnóstico existente, sem repetir o teste. `--rerun-smoke` executa um novo teste isolado e preserva o status anterior. Os 12 namelists existentes não são sobrescritos.

Antes de executar a matriz: resolver as inconsistências com autorização em trabalho separado; validar semântica de umidade, pressão, referência de altura, campos de nuvem ausentes e condições de superfície; repetir o smoke; verificar os dois esquemas. Só então executar cada namelist em uma pasta própria, preservar `summary.nc`/perfis/logs e construir métricas pareadas usando o mesmo perfil. A preparação atual automatiza os namelists e o teste que bloqueia a publicação, sem contornar a validação.
