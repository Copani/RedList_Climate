import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rasterio.features import rasterize
from time import time
from datetime import date, timedelta
from utilities_spacial_vars import *
from calc_bioclim import *

rows = 0 # of rows = 0, all rows are loaded
year = 2004

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
# calc bioclimatic variables for base period 1940 - 1970
bioclim_1940_1970 = xr.Dataset(calc_bioclim(T, P, (1940, 1970)))

# calc bioclimatic variables for 20 years before acessmensts 2004 and 2022
bioclim20y = xr.Dataset(calc_bioclim(T, P, (year-20, year)))

# calc comparison period
if year == 2022:
  bioclim_comp = xr.Dataset(calc_bioclim(T, P, (2004-20, 2004)))
elif year == 2004:
  bioclim_comp = xr.Dataset(calc_bioclim(T, P, (1980-20, 1980)))

# For each species in this dataset, get grid, calc quantities

## summarize climate variables with dictionary of lists
var_codes = ["MAT", "MTWM", "MTCM", "AP", "PDQ"]
columns = {var_code: [np.nan] * rows for var_code in var_codes}
change_columns = {var_code: [np.nan] * rows for var_code in var_codes}
nichefrac_columns = {var_code: [np.nan] * rows for var_code in var_codes}

for i in range(rows):
  species_i = species_data.iloc[[i]]

  ## check if geometry is valid
  rmap = species_i.iloc[0].explode()
  if rmap.geometry.is_valid:  
    for var_code in var_codes:
      ## apply grid preparations leaving time dimension intact
      if var_code == 'MAT':
        var_1940_1970_window, geo_mask, cell_coefs = prep_for_gridcalc(bioclim_1940_1970[var_code], species_i)
      else:
        var_1940_1970_window = prep_for_gridcalc(bioclim_1940_1970[var_code], species_i, window_only=True)

      ## mean for 20 years before assessment
      var_mean_window = prep_for_gridcalc(bioclim20y[var_code].mean('year'), species_i, window_only=True)
      columns[var_code][i] = np.sum(geo_mask * var_mean_window * cell_coefs) / np.sum(geo_mask * cell_coefs)

      ## difference to mean of comparison period
      var_comp_window = prep_for_gridcalc(bioclim_comp[var_code].mean('year'), species_i, window_only=True)
      var_difference_windoe = var_mean_window - var_comp_window
      change_columns[var_code][i] = np.sum(geo_mask * var_difference_windoe * cell_coefs) / np.sum(geo_mask * cell_coefs)
      
      ## calculate fraction of realm outside climate niche for each variable (note that there is a time dimension in bioclim20y)
      min_hist, max_hist = calc_niche(var_1940_1970_window)
      nichefrac_columns[var_code][i], _, _ = calc_fraction_in_niche(min_hist, max_hist, bioclim20y[var_code], species_i, geo_mask=geo_mask, cell_coefs=cell_coefs)

# replace the geometry column in the original dataframe 
species_data = species_data.drop(['geometry'], axis=1)

for var_code in var_codes:
    species_data[f'{var_code}_nichefrac_{year}'] = nichefrac_columns[var_code]
    species_data[f'{var_code}_column_{year}'] = columns[var_code]
    species_data[f'{var_code}_change_column_{year}'] = change_columns[var_code]

# save new dataframe
species_data.to_csv('/home/claussar/gaa2_climate_{}.csv'.format(year))
t5 = time()
print('finished in {} s'.format(t5-t1))
