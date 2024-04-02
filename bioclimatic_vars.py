import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
from time import time
from datetime import date, timedelta


### DATA & PREPROCESSING ###

# Load monthly min, max and avg of daily mean temperatures from ERA5 [°C]
# 0.25 x 0.25° res, 1940 - 2022
T = xr.open_dataset('/home/claussar/ERA5monthly/ERA5_2mtemp_1940-2022.nc')

## write crs to WG84 lon,lat
T = T.rio.write_crs("EPSG:4326")

# Load and preprocesss ERA5 precipitation data [m]
# 0.25 x 0.25° res, 1980 - 2022
P = xr.open_dataset('/home/claussar/ERA5monthly/ERA5_prec_1940_2023.nc')

## convert to -180, 180 longiude
P['longitude'] = (P['longitude'] + 180) % 360 - 180
P = P.sortby('longitude')

## put latitude into correct order (-90 to 90) instead of other way around
P = P.sel(latitude=P['latitude'][::-1])

## write lon lat crs
P = P.rio.write_crs("EPSG:4326")
print('loaded')
  

# calculate bioclimatic variables
def calc_bioclim(T, P, years):
  """
  Calculate bioclimatic variables for a given range of years.

  Parameters:
  T (xarray.Dataset): Dataset containing temperature data with variables 'Tavg', 'Tmax', 'Tmin'.
  P (xarray.Dataset): Dataset containing precipitation data with variable 'tp'.
  years (tuple): A tuple containing the start and end years for the calculation.

  Returns:
  dict: A dictionary containing the calculated bioclimatic variables:
      'MAT' - Mean Annual Temperature
      'MTWM' - Maximum Temperature of Warmest Month
      'MTCM' - Minimum Temperature of Coldest Month
      'AP' - Total Annual Precipitation
      'PDQ' - Precipitation of Driest Quarter
  """
  start_year, end_year = years

  # Calculate mean annual temperature (MAT) averages
  MAT = T.Tavg.groupby('time.year').mean('time')
  MAT_years = MAT.sel(year=slice(start_year, end_year))
  
  # Calculate max temperature of warmest month (MTWM)
  MTWM = T.Tmax.groupby('time.year').max('time')
  MTWM_years = MTWM.sel(year=slice(start_year, end_year))

  # Calculate min temperature of coldest month (MTCM)
  MTCM = T.Tmin.groupby('time.year').min('time')
  MTCM_years = MTCM.sel(year=slice(start_year, end_year))

  # Calculate total annual precipitation (AP)
  AP = P.tp.groupby('time.year').sum('time')[:,0]
  AP_years = AP.sel(year=slice(start_year, end_year))

  # Calculate precipitation of driest quarter (PDQ)
  Prec_quaterly = P.tp.rolling(time=4).sum()
  PDQ = Prec_quaterly.groupby('time.year').min()[:,0]
  PDQ_years = PDQ.sel(year=slice(start_year, end_year))

  return {
    'MAT': MAT_years,
    'MTWM': MTWM_years,
    'MTCM': MTCM_years,
    'AP': AP_years,
    'PDQ': PDQ_years
  }

c = calc_bioclim(T, P, years=[1940,2022])

ds = xr.Dataset(c)
ds.to_netcdf('/home/claussar/climate variables/bioclimatic_variables.nc')
