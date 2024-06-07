import numpy as np
import xarray as xr

# load bio climatic variables
ds = xr.open_dataset('ClimateVariablesOut/bioclimatic_variables_30km.nc')

# take mean between 1990 and 2020
ds = ds.sel(year=slice('1990', '2020')).mean('year')

# store as .nc file
ds.to_netcdf('ClimateVariablesOut/bioclimatic_variables_30km_1990-2020_mean.nc')