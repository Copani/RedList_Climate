import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
# from time import time
# from datetime import date, timedelta
from calc_bioclim import *


### DATA & PREPROCESSING ###

# Load monthly min, max and avg of daily mean temperatures from ERA5 [°C]
# 0.25 x 0.25° res, 1940 - 2022
# T = xr.open_dataset('/home/claussar/ClimateVariables/ERA5monthly/ERA5_2mtemp_1940-2022_mean_dailyminmax.nc')
T = xr.open_dataset('../Data/ClimateVariables/ERA5_2mtemp_1940-2022_mean_dailyminmax.nc')

## write crs to WG84 lon,lat
T = T.rio.write_crs("EPSG:4326")

## the old calculation of calc_monthlyT_ERA5.py stored the last day of each month. Recalc to get the first day!
T['time'] = T['time'].astype('datetime64[M]')


# Load and preprocesss ERA5 precipitation data [m]
# 0.25 x 0.25° res, 1980 - 2022
# P = xr.open_dataset('/home/claussar/ClimateVariables/ERA5monthly/ERA5_prec_1940_2023.nc')
P = xr.open_dataset('../Data/ClimateVariables/ERA5_prec_1940_2023.nc')

## convert to -180, 180 longiude
P['longitude'] = (P['longitude'] + 180) % 360 - 180
P = P.sortby('longitude')

## put latitude into correct order (-90 to 90) instead of other way around
P = P.sel(latitude=P['latitude'][::-1])

## write lon lat crs
P = P.rio.write_crs("EPSG:4326")
print('loaded')

P = P.sel(time=slice(str(1970), str(2000)))
  
AP = P.tp.groupby('time.year').sum('time')[:,0] * 365.25 * 1000
AP = AP.drop_vars('spatial_ref').drop_vars('expver').rio.write_crs("EPSG:4326")

# Calculate precipitation seasonality (PS)
PS = P.tp.groupby('time.year').std('time')[:,0] / (1 + AP / 12) * 100
PS = PS.drop_vars('spatial_ref').drop_vars('expver').rio.write_crs("EPSG:4326")

# calculate bioclimatic variables
# c = calc_bioclim(T, P, years=[1970,2000])

# ds = xr.Dataset(c)
# ds.to_netcdf('ClimateVariablesOut/bioclimatic_variables_30km.nc')
PS.to_netcdf('ClimateVariablesOut/bioclimatic_variables_Bio15_30km_2.nc')
