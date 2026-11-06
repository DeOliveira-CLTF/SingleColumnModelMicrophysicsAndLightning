# Fontes e coordenadas

ERA5/EDA: CDS, com requisições, inventários e SHA256 em `data/raw/era5/`.
O domínio disponível da figura 01 é 68–50°W, 40–27°S; não há extrapolação.

SACO/WMO 87344: 31,32°S e 64,22°W, conforme o
[inventário JMA](https://www.data.jma.go.jp/tcc/tcc/products/climate/extlist_weekly/extlist_20261014.html).
O ponto GPS inicial do balão não substitui a coordenada de estação.
CSV locais em `data/raw/radiosondes/`, harmonizados em `data/processed/profiles/`.
Horários reais são preservados quando fornecidos; o perfil de 12 UTC é nominal.

ENTLN: privado, `1=IC`, `0=CG`, cobertura completa confirmada pelo usuário.
Agrupamento antecede o recorte; dados ausentes não viram zeros.

Fronteiras: Natural Earth, domínio público, admin-0 e admin-1 global.
Fontes e hashes em `data/raw/boundaries/manifest.json`.
