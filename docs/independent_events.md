# Reprodução dos dois eventos

Configurações: `configs/cases/argentina_20191028.yaml` (calibração) e
`configs/cases/argentina_20181213.yaml` (avaliação). Dados e derivados são privados,
locais e ignorados pelo Git. Não há novo ensemble de 2019.

```powershell
python scripts/download.py era5 --config configs/cases/argentina_20191028.yaml
python scripts/download.py cache --config configs/cases/argentina_20191028.yaml
python scripts/prepare_forcing.py event select --config configs/cases/argentina_20191028.yaml
python scripts/prepare_forcing.py event build --config configs/cases/argentina_20191028.yaml
python scripts/process_lightning.py observations --config configs/cases/argentina_20191028.yaml
python scripts/process_lightning.py aggregate --config configs/cases/argentina_20191028.yaml
python scripts/run_scm.py independent run-control
python scripts/evaluate.py independent prepare
python scripts/evaluate.py independent refresh-cape
python scripts/calibrate.py independent calibrate
python scripts/evaluate.py independent evaluate
python scripts/evaluate.py meteorology download
python scripts/evaluate.py meteorology build
python scripts/plot.py article
python -m pytest tests -q
```

As aquisições ENTLN usam `entln_source`, a convenção `1=IC, 0=CG`, agrupamento
causal por primeiro evento, duração máxima 700 ms e raio 10 km. São agrupados
IC e CG conjuntamente, com remoção de duplicatas antes de selecionar a área.
O dia anterior protege o início da janela. Eventos posteriores ao fim não
alteram a criação, posição nem horário de âncoras anteriores; apenas a
multiplicidade/classificação dos flashes no limite pode ficar censurada.
O alvo principal total IC+CG não depende desses eventos posteriores. A mesma
função de agrupamento e contagem é utilizada nos dois eventos.

A área de 2019 é fixada com ERA5 antes da calibração elétrica. O critério
pré-especificado é a média espacial do produto das médias temporais da
ascensão positiva a 500/700 hPa e da precipitação horária. A busca regional
e a evolução horária permanecem em `results/argentina_20191028/support/`.
O início/fim são os dois dias completos, com observação confirmada pelo usuário.
Nenhuma seleção utiliza erro de raios.

## Convenções das três parametrizações

- CAPE×P: parcela da célula inferior, MetPy pseudoadiabática líquida. A análise
  independente usa `which_lfc=which_el=most_cape`, definido termodinamicamente
  antes da calibração; a integral assinada antiga bottom/top é preservada na
  auditoria. Não há piso positivo, mudança de T/umidade ou ajuste pelo erro ENTLN.
  Retorno ainda negativo é indisponível (NaN). CAPE em J/kg multiplica a diferença
  exata de `rainnc` dos 900 segundos seguintes, em kg/m². O produto é convertido
  de J/m² para J/km² antes de ajustar densidade de flashes.
- I×G: integral vertical de `(rho_d*qi)*(rho_d*qg)`; gelo de nuvem não é trocado
  por gelo total. Sua integral temporal tem unidade kg² s m⁻⁵.
- McCaul F3 LS: F1=0,042(w⁺ qg[g/kg]) a −15°C; F2=0,20 VII[kg/m²]; pesos
  0,95/0,05. Integração temporal dividida explicitamente por 300 s e pela área
  nominal da célula publicada (0,009° na latitude de referência de Huntsville).
  A área nominal não é uma reconstrução exata de todas as células do artigo.
  O fator adicional de 2019 é razão de máximos de densidade em intervalos de
  cinco minutos. Não substitui os coeficientes publicados de F1/F2.

CAPE×P e I×G recebem regressão proporcional sem intercepto sobre densidades de
15 minutos do controle 2019. O mapeamento de I×G é uma calibração nova do SCM:
Romps (2019) não publica uma lei universal de conversão em flashes. O fator
de McCaul é uma extensão de evento único, inspirada na comparação de máximos;
não reproduz a regressão entre máximos de várias tempestades nem a calibração
separada de F1/F2 de McCaul (2009). LS e EDMF não equivalem à velocidade
convectiva resolvida original. Não se somam as velocidades.

Os três coeficientes e a média de densidade observada de 2019 são congelados em
`results/argentina_20191028/calibration/v1/`. Não é permitido sobrescrever a
versão congelada. O manifesto registra entradas, hashes e suporte. A contagem
prevista em 2018 é densidade prevista vezes a área meteorológica fixa de 2018.
Todas as configurações e membros recebem o mesmo coeficiente por método.

## Produtos e verificações

`results/argentina_20181213/article/` contém CSVs de membros, média/desvio dos
40 membros (controle separado), máscara temporal comum, cinco tabelas principais,
seis figuras e suplementos. Os erros são previsão menos observação. Skills são
1−MSE/MSE_referência; referência zero e constante transferida da média de
densidade de 2019. Denominadores nulos e correlações indefinidas são NaN com motivo.
Picos usam o primeiro máximo; séries constantes não têm horário de pico distinto.

A validação externa usa CHIRPS v2 diário no servidor oficial atual
`https://data.chc.ucsb.edu/products/CHIRPS/v2.0/`. O produto diário é desagregado
com CFSv2 a partir de totais pentadais/mensais ajustados por satélite/estações;
não constitui verdade pluviométrica diária independente nem validação de
horários de convecção a cada 15 minutos. O CSV e o manifesto preservam essa
limitação. Precipitação simulada é somente `rainnc`, que já inclui água sólida.

Os testes cobrem bloqueio de 2018 no ajuste, referência TH_C0 membro zero,
coeficientes compartilhados, hashes nativos, conversões de área/tempo, CAPE
multicamadas, precipitação acumulada, máscaras comuns, skills e vínculos CSV/figuras.
O diagnóstico de orçamento usa densidade de entrada; o Thompson recalcula
densidade internamente. O resíduo não é um certificado de fechamento exato de
água sob pressão prescrita. Não são alteradas equações físicas ou do integrador.
