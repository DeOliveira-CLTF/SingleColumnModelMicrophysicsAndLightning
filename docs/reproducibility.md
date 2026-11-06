# Reprodução e contratos

`evaluate.py vertical`: saídas em `runs/argentina_20181213/_diagnostic_audit_v2/`.
Prova em `results/argentina_20181213/tests/vertical_instrumentation.json`:
todas as variáveis, nos 5761 instantes nativos, para cada configuração.

`process_lightning.py diagnostics`: 192 janelas por método/configuração,
alcance de isotermas por passo e metadados da hipótese em
`results/argentina_20181213/vertical_diagnostics/`. Arquivos `.vertical.nc`
permanecem nas execuções apontadas pelo manifest. Não são saídas dos 40 membros.

`calibrate.py adaptations` é opcional, separado, apenas para controles.
`calibrate.py fit` usa a média dos membros perturbados do ensemble principal;
uma nova avaliação exige `--evaluation-output CAMINHO_NOVO`.

`evaluate.py migration` verifica hashes das 246 saídas, 41 entradas e estado
local de DEPHY. Manifests antigos conservam caminhos originais, resolvidos por
`provenance.resolve_path`. Arquivos nativos são imutáveis, incluindo hardlinks.

Credenciais ficam fora do repositório. Não publicar ENTLN ou derivados.
