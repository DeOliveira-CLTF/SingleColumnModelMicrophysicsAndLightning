# Avaliação observacional pendente

Consulta em 02/11/2026. Nenhuma Figura 5 foi produzida. Há observações candidatas de 2014, mas não há pareamento demonstrado com a coluna simulada. Nenhum arquivo observacional foi usado como resultado ou como validação.

## Impedimento prioritário: identificar o caso

SCM e DEF locais declaram 8°S, 63°W, início 05/10/2014 12:00 e fim 06/10/2014 00:00. O eixo numérico cobre 35,255859 h, não 12 h. Não se pode escolher uma tempestade em Manaus ou Manacapuru apenas pelo nome GOAMAZON. A origem do perfil, o domínio médio da forçante e a transformação temporal precisam ser confirmados com os arquivos de origem e o responsável pelo caso. O script local `driver_DEF.py` codifica essas coordenadas/datas; isso não resolve sua validade geográfica.

## Fontes candidatas verificadas

| Fonte | Cobertura e variável | Acesso e limitação para esta tarefa |
| --- | --- | --- |
| [ARM GOAmazon — forçantes e precipitação regional](https://portal.nersc.gov/cfs/capt/ARMVAP/GOAmazon.html) | Produto de 2014–2015; domínio circular de raio 110 km centrado no radar SIPAM de Manaus; resolução de forçante de 3 h e 25 hPa | Portal público. É candidato para reconstruir a procedência e precipitação média regional; coordenadas locais ainda não reconciliadas. Campos usados para construir a forçante não são validação totalmente independente. |
| [Lopes et al. (2025), ESSD](https://essd.copernicus.org/articles/17/5489/2025/) e [dados Zenodo](https://zenodo.org/records/13732692) | Sistemas precipitantes rastreados por radar na região de Manaus em 2014–2015; associação com descargas de retorno CG da GLD360, agregadas em passos de 12 min | Arquivos públicos GeoJSON em ZIP. Disponibilidade confirmada, mas arquivos volumosos não foram baixados: sem localização/cronologia reconciliadas, não haveria seleção defensável. Strokes CG não equivalem a flashes totais do EP. |
| [Vilela (2017), tese IAG/USP](https://iag.usp.br/pos-graduacao/meteorologia/publicacoes/teses-e-dissertacoes/estrutura-eletrica-das-nuvens-na) | Descreve redes LINET e STARNET na campanha de 29/08 a 07/10/2014, região de Manacapuru | Documento público; acesso aos eventos brutos e eficiência de detecção não foi confirmado. A sobreposição de datas, sozinha, é insuficiente. |
| [GLM — inventário NOAA/NCEI](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C01527) | Arquivo listado a partir de 05/07/2017 | Incompatível com 2014; excluído da validação contemporânea. |

O vínculo local à literatura é a referência de Tian e Zhang (2025), DOI [10.1029/2024GL113920](https://doi.org/10.1029/2024GL113920), incluída nos metadados do SCM. A referência bibliográfica não basta para resolver as diferenças entre coordenadas, cronologia e séries locais.

## Dados necessários para desbloquear uma comparação

1. Identificação inequívoca do perfil e da forçante: local ou domínio médio, convenção vertical, horário UTC inicial, data de cada amostra e motivo da discrepância de duração.
2. Precipitação contemporânea com controle de qualidade: acumulado e taxa em uma área representativa da coluna ou do domínio da forçante. Documentar estação/radar, cobertura, unidades, incerteza e falhas de amostragem.
3. Série elétrica georreferenciada: horário UTC, flashes ou strokes, IC/CG, algoritmo de agrupamento e eficiência de detecção. Não converter strokes em flashes sem método documentado.
4. Janela temporal e área comuns, definidas antes de comparar; explicitar o significado espacial do diagnóstico por coluna. O fim do evento simulado é censurado, portanto a duração observada completa não é diretamente equivalente.
5. Avaliar primeiro ocorrência e organização temporal. Uma comparação quantitativa absoluta exige tratar calibração e suporte espacial; não recalibrar o EP nesta tarefa. A partição IC/CG fixa não sustenta validação independente por tipo.

Depois dessas verificações, produzir uma Figura 5 somente com dados efetivamente lidos e documentados. Até lá, o artigo deve apresentar uma comparação de sensibilidade numérica, sem alegar habilidade observacional ou superioridade de um esquema.
