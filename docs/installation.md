# Instalação técnica

Python ≥3.10; dependências em `pyproject.toml` e ambiente em `environment.yml`.
Os scripts adicionam `python/` ao caminho de importação. A instalação editável
permite também `python -m scm_lightning.cli`.

Compilação Linux/WSL: gfortran, make e NetCDF-Fortran com `nf-config` no PATH.
No Windows, o executor usa `Ubuntu-24.04`. Os testes Python não exigem WSL;
a comparação das seis configurações exige o executável compilado.

Não executar reset/clean no DEPHY: suas alterações locais anteriores foram
preservadas em `external/DEPHY-SCM`.
