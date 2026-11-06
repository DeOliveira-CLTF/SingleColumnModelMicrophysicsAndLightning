#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GOAMAZON/SINGLE - segunda etapa (grade unica para o SCM), versao ERA5.

Adaptado do driver_SCM.py oficial do repositorio DEPHY-SCM (F. Couvreux,
R. Roehrig): mesma logica de regrade (convert2SCM para uma grade unica de
altura 0-20000 m com 10 m de espacamento, e a grade temporal do proprio
hfss), so trocando a fonte da temperatura de pele de superficie do arquivo
de texto original (sing_sfc_flux_278_5.txt) para o era5_ts.npz que o
driver_DEF_era5.py ja gera.

Por que esta etapa e necessaria (nao e so um "nice to have"): o
mp_offline_driver.F90/goamazon_bridge.F90 deste projeto le o arquivo
DEPHY-SCM assumindo dimensoes 'lev' e 'time' UNICAS, compartilhadas por
todas as variaveis (nf90_inq_dimid(...,'lev',...) e (...,'time',...) uma
vez so). O arquivo _DEF_ (saida do driver_DEF_era5.py) NAO tem isso: cada
variavel tem sua propria dimensao (lev_ua, lev_theta, lev_tntheta_adv...,
com tamanhos DIFERENTES entre condicao inicial e forcante, e nas alturas
nativas, irregulares, dos ~37 niveis de pressao do ERA5) -- exatamente
como o arquivo oficial do repositorio (GOAMAZON_SINGLE_DEF_driver.nc, que
tem essa mesma estrutura, conferido diretamente). Sem esta segunda etapa,
o Fortran falha ja na abertura do arquivo (dimensao 'lev' nao encontrada).

Requer a biblioteca dephycf do repositorio DEPHY-SCM (github.com/GdR-DEPHY/
DEPHY-SCM) no PYTHONPATH; nao esta no PyPI.
export PYTHONPATH=/caminho/para/DEPHY-SCM:$PYTHONPATH
"""
import numpy as np
import netCDF4 as nc
from dephycf.Case import Case

# ---- opcoes ---------------------------------------------------------------
DEF_FILE = 'GOAMAZON_SINGLE_DEF_driver.nc'   # saida do driver_DEF_era5.py
OUT_FILE = 'GOAMAZON_SINGLE_SCM_driver.nc'   # nome que goamazon_bridge.F90 espera
TS_NPZ   = 'era5_ts.npz'                     # gerado pelo driver_DEF_era5.py
# --------------------------------------------------------------------------

# 1) le a versao _DEF_ (grade nativa do ERA5, dimensoes por variavel)
case = Case('GOAMAZON/SINGLE')
case.read(DEF_FILE)

fin = nc.Dataset(DEF_FILE, 'r')

# 2) grade unica: altura 0-20000 m em passos de 10 m (a mesma que
#    goamazon_bridge.F90 espera -- ver o comentario no topo do arquivo
#    Fortran), tempo = a propria grade de hfss do arquivo _DEF_
levout = np.array(range(0, 20001, 10), dtype=np.float64)
timeout = np.array(fin['time_hfss'][:])
newcase = case.convert2SCM(time=timeout, lev=levout, levtype='altitude')

# 3) temperatura de pele (Ts): do era5_ts.npz gerado pelo driver_DEF_era5.py,
#    no lugar do sing_sfc_flux_278_5.txt do caso original
ts_data = np.load(TS_NPZ)
newcase.add_surface_temp(ts_data['ts'].astype(np.float64),
                          time=ts_data['time'].astype(np.float64), timeid='time')

newcase.set_title("Forcing and initial conditions for GOAMAZON/SINGLE case "
                   "(shallow-to-deep transition) - SCM-enabled version, ERA5-based")
newcase.set_script("driver_SCM_era5.py")

# 4) grava no nome que o driver Fortran espera
newcase.write(OUT_FILE)
