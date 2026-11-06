# Verificação da entrega

- Onze testes automatizados aprovados; saída preservada em `testes.log`. Incluem reprodução das referências CF, séries nulas, censura temporal, não finitos, tempo irregular, configuração divergente, identidade inicial, parser de inteiros, arquivo ausente, análise parcial e integrais independentes dos seis casos.
- A conferência independente usa netCDF4, soma por intervalos, fechamento da chuva incremental e reconstrução geométrica da grade por bisseção, sem chamar a função de interfaces da análise principal.
- As seis simulações foram reconhecidas como reutilizáveis numa segunda chamada do executor. Nenhuma nova integração principal foi realizada nessa verificação.
- Figuras 1–4 e Figura 3 suplementar examinadas visualmente em PNG e em renderizações dos PDFs. Conferidos títulos, eixos, unidades, legendas, barras de cores, números, zeros e ausência de recortes. Ajustados contraste das isotermas, posição da indicação de campo nulo, indicação do pequeno gelo Thompson e espaço para legendas das barras.
- PNGs finais em 300 dpi; PDFs de uma página cada. Mapas de calor rasterizados; demais elementos vetoriais. Miniaturas e renderizações usadas na inspeção estão em `figuras/` nesta pasta.
- Tabelas mantêm 17 algarismos significativos; métricas ausentes e percentuais indefinidos são campos vazios, distintos de zero. O relatório apresenta arredondamento apenas para leitura.
- `auditoria_final.json` registra preservação de fontes, entrada e seis NetCDF; nenhum diff em `src/`, `data/` ou namelist principal. Alterações de build e submódulo já existentes não foram revertidas.
- Figuras 5/6 e Tabela 4 não foram fabricadas. Observações sem pareamento confirmado e smoke ERA5 reprovado estão documentados separadamente.

Limite: aprovação é dos critérios numéricos implementados. Não substitui validação observacional, análise de conservação por processos ou uma auditoria retrospectiva completa dos dois CF anteriores.
