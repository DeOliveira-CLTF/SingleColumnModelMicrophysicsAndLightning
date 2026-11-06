# SCM: microfísica e raios — calibração 2019 e avaliação independente 2018

O controle determinístico Thompson C0 de 28–29/10/2019 fornece três coeficientes
fixos. O evento de 13–14/12/2018 é exclusivo de avaliação; suas observações não
podem ajustar coeficientes. O fluxo anterior está preservado em
`results/argentina_20181213/legacy/article_within_event_20261106/`.

Fluxo técnico dos dois eventos: [docs/independent_events.md](docs/independent_events.md).

Thompson/WSM6 × MYNN C0/CL/CF. Ensemble existente: 40 membros perturbados e
um controle por configuração. ENTLN e seus derivados permanecem locais.

| Pasta | Conteúdo |
|---|---|
| `src/` | SCM Fortran; instrumentação em `diagnostics/` |
| `python/scm_lightning/` | Download, forcing, observações, raios, ensemble, calibração, avaliação e figuras |
| `scripts/` | Oito pontos de entrada; implementações ficam no pacote |
| `configs/` | Configurações YAML do caso e métodos |
| `external/DEPHY-SCM/` | Dependência Git com alterações locais anteriores preservadas |
| `data/` | Cache privado, intermediários, forcing e perfis |
| `runs/argentina_20181213/` | `thompson_C0`…`wsm6_CF`, com `control/` e `members/` |
| `results/argentina_20181213/` | Figuras, tabelas, avaliação e diagnósticos |
| `archive/goamazon/` | Fontes e material histórico fora do fluxo principal |
| `build/` | Compilação local, ignorada pelo Git |

## Uso

```powershell
python -m pip install -e ".[test]"
python scripts/download.py --help
python scripts/prepare_forcing.py --help
python scripts/run_scm.py --help
python scripts/process_lightning.py --help
python scripts/calibrate.py --help
python scripts/evaluate.py --help
python scripts/plot.py --help
python scripts/reproduce.py status
python -m pytest tests -q
```

No Linux/WSL com gfortran e NetCDF-Fortran: `make -j4 BUILD=build/vertical_audit`.

```powershell
python scripts/evaluate.py vertical
python scripts/process_lightning.py diagnostics
python scripts/download.py boundaries
python scripts/plot.py all
python scripts/evaluate.py migration
python scripts/evaluate.py contracts
```

`vertical` executa as seis configurações por 48 h em pasta separada e compara
todas as variáveis nativas com os controles anteriores. Não substitui o ensemble.
Mudanças de executável ou inputs exigem uma nova pasta `--output`.

## Contratos

Omega negativo em Pa/s representa ascensão: `w=−omega/(rho*g)`. A figura 01
mostra omega sem inversão, recorta a cobertura finita do ERA5 e usa SACO/WMO
87344 em 31,32°S e 64,22°W. O ponto GPS de lançamento é distinto da coordenada
de estação.

`mf_w` é densidade vezes velocidade média das plumas nas interfaces. A saída
`.vertical.nc` exporta `w_EDMF=mf_w/rhoz`, `w_LS`, `mf_a`, densidades, alturas,
espessuras e campos da coluna no instante da chamada MYNN. Plumas inativas
têm `w_EDMF=NaN`; as velocidades não são somadas. A física principal permanece
`mynn_edmf_type=0`, `mynn_uedmf_mp=0`.

As adaptações `*_LS_1D_adaptado` e `*_EDMF_1D_adaptado` são separadas das
formulações originais. Combinar velocidade EDMF com hidrometeoros da coluna
assume representatividade da coluna nas plumas; não produz fluxos covariantes
de espécies. CAPE×P usa CAPE pseudoadiabática MetPy do perfil simulado,
amostrada a cada 15 min; não é apresentada como a CAPE isentrópica de Romps.
Esses novos diagnósticos foram calculados para os seis controles, separadamente
da avaliação principal dos 40 membros.

A área meteorológica de 2018 permanece fixa. Os dois dias completos de dezembro
são exclusivos de avaliação; a divisão antiga desenvolvimento/intervalo/avaliação
fica somente no histórico. Precipitação total vem de `rainnc`, sem somar `snownc`.

Manifests anteriores foram preservados; `resolve_path` traduz suas referências
antigas. Consulte [instalação](docs/installation.md), [fontes](docs/data_sources.md),
[metodologia](docs/methodology.md), [esquemas](docs/lightning_schemes.md) e
[reprodução](docs/reproducibility.md).

## Sensibilidade determinística do movimento vertical

`configs/experiments/updraft_sensitivity.yaml` congela o desenho meteorológico
antes dos resultados: ERA5, incremento fraco e incremento mais intenso × seis
configurações. O pulso conserva a forma cosseno temporal × gaussiana vertical
do código histórico, sem saturação artificial. O integrador principal e as
equações físicas são preservados. Não se multiplicam esses casos pelos 40 membros.

```powershell
python scripts/run_scm.py sensitivity design
python scripts/run_scm.py sensitivity run
python scripts/process_lightning.py sensitivity
```

O desenho existente é imutável; `design` recusa substituí-lo. Compilar com
`make -j4 BUILD=build/updraft_sensitivity`. Os 18 casos ficam em
`runs/argentina_20181213/updraft_sensitivity/`. NetCDF derivados, tabelas,
figuras e testes ficam em `results/argentina_20181213/updraft_sensitivity/`.

`w_LS`, `w_perturbation`, `w_effective=w_LS+w_perturbation` e `w_EDMF` são
exportados separadamente. A decomposição de quatro combinações de estado e
velocidade separa resposta do estado, efeito da velocidade e interação, sem
calibração ou leitura da ENTLN. A contabilidade de água usa massas secas iniciais
fixas e conserva separadamente os resíduos de microfísica/precipitação e as
diferenças de densidade; não certifica fechamento atmosférico completo.
A demanda de divergência para equilibrar o pulso é diagnosticada, não aplicada.
O integrador existente transporta verticalmente theta e vapor; não foi adicionado
transporte resolvido de hidrometeoros. Seus produtos de fluxo elétrico são diagnósticos.
