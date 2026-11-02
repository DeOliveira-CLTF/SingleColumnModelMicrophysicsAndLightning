import netCDF4 as nc
d = nc.Dataset('GOAMAZON_SINGLE_DEF_driver.nc')
print(d['ps'][:])       # deve ser um valor perto de 95000-96000 Pa (pressao em Manaus)
print(d['theta'][0,:5]) # theta inicial nos primeiros niveis, perto de 295-300 K
