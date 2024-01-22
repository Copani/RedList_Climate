import xarray as xr
import numpy as np


x = xr.open_dataset('../Data/ERA5/ERA5_temp_prec_wind_veg_1980_2023.nc')
# save the vegetation variable as nc file
x['tp'].to_netcdf('../Data/ERA5/ERA5_veg_1980_2023.nc')