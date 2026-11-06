#!/bin/bash -x
echo git clone https://github.com/GdR-DEPHY/DEPHY-SCM.git
echo pip install netCDF4 numpy scipy matplotlib
echo export PYTHONPATH=/caminho/para/DEPHY-SCM:$PYTHONPATH
export PYTHONPATH=/p/projetos/monan_atm/paulo.kubota/pgmet/pgmet_met_348_3/caio.oliveira/mp_offline_project_02/mp_offline/era5_forcing/DEPHY-SCM:$PYTHONPATH
module swap PrgEnv-cray PrgEnv-gnu/8.6.0 
module load cray-netcdf
module load cray-python
module load cdo2
rm  era5_sl.nc
cdo merge data_stream-oper_stepType-instant.nc data_stream-oper_stepType-avg.nc  era5_sl.nc
rm era5_pl.nc
cdo mergetime  lev_05.nc  lev_06.nc era5_pl.nc
python3 driver_DEF_era5.py
python3 check.py
python3 driver_SCM_era5.py
