import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rasterio.features import rasterize
from time import time
from datetime import date, timedelta
from utilities_spacial_vars import *

rows = 0 # of rows = 0, all rows are loaded
year = 2022

t1 = time()


### DATA & PREPROCESSING ###

# Load species data
# species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/IUCN/AMPHIBIANS/AMPHIBIANS.shp')
if rows == 0:
  species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/GAA2/GAA2_allspecies.gdb')
  rows = len(species_data)
else:
  species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/GAA2/GAA2_allspecies.gdb', rows=rows)

t2 = time()
print('{} rows of species data loaded in {} s'.format(rows, t2-t1))

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
  
t4 = time()
print('ERA5 Temperature and Precipitation data loaded and preprocessed in {} s'.format(t4-t2))

### CALCULATIONS ###

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
  


# calc base period 1940 - 1970
bioclim_1940_1970 = xr.Dataset(calc_bioclim(T, P, (1940, 1970)))
bioclim_2002_2022 = xr.Dataset(calc_bioclim(T, P, (2002, 2022)))
# bioclim_2004_2023 = xr.Dataset(calc_bioclim_vars_for_years(T, P, (2004, 2021)))

t5 = time()
print('5 climate fields calculated in {} s'.format(t5-t4))


# For each species in this dataset, get grid, calc quantities
MAT_column = [np.nan] * rows


MATnichefrac_column = [np.nan] * rows
MTWMnichefrac_column = [np.nan] * rows
MTCMnichefrac_column = [np.nan] * rows
APnichefrac_column = [np.nan] * rows
PDQnichefrac_column = [np.nan] * rows

for i in range(rows):
  species_i = species_data.iloc[[i]]

  ## check if geometry is valid
  rmap = species_i.iloc[0].explode()
  if rmap.geometry.is_valid:
       
    ## apply grid preparations leaving time dimension intact
    MAT_1940_1970_window, geo_mask, cell_coefs = prep_for_gridcalc(bioclim_1940_1970.MAT, species_i)
    MTWM_1940_1970_window = prep_for_gridcalc(bioclim_1940_1970.MTWM, species_i, window_only=True)
    MTCM_1940_1970_window = prep_for_gridcalc(bioclim_1940_1970.MTCM, species_i, window_only=True)
    AP_1940_1970_window = prep_for_gridcalc(bioclim_1940_1970.AP, species_i, window_only=True)
    PDQ_1940_1970_window = prep_for_gridcalc(bioclim_1940_1970.PDQ, species_i, window_only=True)
    
    ## calculate fraction of 1980-2004 mean MAT outside climate niche.
    # MAT
    min_hist, max_hist = calc_niche(MAT_1940_1970_window)
    MATnichefrac_column[i], _, _ = calc_fraction_in_niche(min_hist, max_hist, bioclim_2002_2022.MAT, species_i, geo_mask=geo_mask, cell_coefs=cell_coefs)
    # MTWM
    min_hist, max_hist = calc_niche(MTWM_1940_1970_window)
    MTWMnichefrac_column[i], _, _ = calc_fraction_in_niche(min_hist, max_hist, bioclim_2002_2022.MTWM, species_i, geo_mask=geo_mask, cell_coefs=cell_coefs)
    # MTCM
    min_hist, max_hist = calc_niche(MTCM_1940_1970_window)
    MTCMnichefrac_column[i], _, _ = calc_fraction_in_niche(min_hist, max_hist, bioclim_2002_2022.MTCM, species_i, geo_mask=geo_mask, cell_coefs=cell_coefs)
    # AP
    min_hist, max_hist = calc_niche(AP_1940_1970_window)
    APnichefrac_column[i], _, _ = calc_fraction_in_niche(min_hist, max_hist, bioclim_2002_2022.AP, species_i, geo_mask=geo_mask, cell_coefs=cell_coefs)
    # PDQ
    min_hist, max_hist = calc_niche(PDQ_1940_1970_window)
    PDQnichefrac_column[i], _, _ = calc_fraction_in_niche(min_hist, max_hist, bioclim_2002_2022.PDQ, species_i, geo_mask=geo_mask, cell_coefs=cell_coefs)

    


# replace the geometry column in the original dataframe 
species_data = species_data.drop(['geometry'], axis=1)
# # or
# add data to additional columns in existing dataframe
# del species_data
# species_data = pd.read_csv('/home/claussar/gaa2_climate_2022.csv')

species_data['MAT_nichefrac_{}'.format(year)] = MATnichefrac_column
species_data['MTWM_nichefrac_{}'.format(year)] = MTWMnichefrac_column
species_data['MTCM_nichefrac_{}'.format(year)] = MTCMnichefrac_column
species_data['AP_nichefrac_{}'.format(year)] = APnichefrac_column
species_data['PDQ_nichefrac_{}'.format(year)] = PDQnichefrac_column


# print('qantities for {} species calculated in {} s'.format(rows, t4-t3))

# save new dataframe
species_data.to_csv('/home/claussar/gaa2_climate_2021_5y.csv')
# species_data.to_csv('RedList_Climate\gaa2_landuse_2021_2004_1980.csv')

t4 = time()
print('finished in {} s'.format(t4-t1))
