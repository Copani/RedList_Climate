import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rasterio.features import rasterize
from time import time
from datetime import date, timedelta
from utilities_spacial_vars import *

rows = 20 # of rows = 0, all rows are loaded
year = 2004
margin = 0.5 # grid resolution 0.25

t1 = time()


### DATA & PREPROCESSING ###

# Load species data
# species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/IUCN/AMPHIBIANS/AMPHIBIANS.shp')
if rows == 0:
  species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/GAA2/GAA2_allspecies.gdb')
  rows = len(species_data)
else:
  # species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/GAA2/GAA2_allspecies.gdb', rows=rows)
  species_data = gpd.read_file('Data/GAA2/GAA2_allspecies.gdb', rows=rows)

t2 = time()
print('{} rows of species data loaded in {} s'.format(rows, t2-t1))

# Load bioclimatic variavles
# bioclim = xr.open_dataset('/home/claussar/climate variables/bioclimatic_variables.nc')
bioclim = xr.open_dataset('Data/ClimateVariables/bioclimatic_variables.nc')

### CALCULATIONS ###  


# calc bioclimatic variables for base period 1940 - 1970
bioclim_1940_1970 = bioclim.sel(year=slice(1940,1970))

# calc bioclimatic variables for 20 years before acessmensts 2004 and 2022
bioclim20y = bioclim.sel(year=slice(year-20, year))

# calc comparison period
if year == 2022:
  bioclim_comp = bioclim.sel(year=slice(2004-20, 2004))
elif year == 2004:
  bioclim_comp = bioclim.sel(year=slice(1980-20, 1980))

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
    
    # initialize grid with "MAT", because it has the same shape as the other variables
    grid = GridCalc(bioclim_1940_1970['MAT'], species_i, margin=margin)
    
    for var_code in var_codes:
      ## average over 20 years prior to assessment
      var_window = grid.get_window(bioclim20y[var_code])
      columns[var_code][i] = grid.calc_average(var_window.mean('year'))

      ## difference to mean of comparison period
      var_comp = grid.get_window(bioclim_comp[var_code].mean('year'))
      var_difference = var_window.mean('year') - var_comp
      change_columns[var_code][i] = grid.calc_average(var_difference)
      
      ## calculate fraction of realm outside climate niche for each variable (note that there is a time dimension in bioclim20y)
      ## get historical window for each variable (lon, lat, time)
      hist_window = grid.get_window(bioclim_1940_1970[var_code])
      min_hist, max_hist = calc_limits(hist_window)
      
      nichefrac_columns[var_code][i] = grid.calc_fraction_in_niche(min_hist, max_hist, var_window, nichemode='minmax')
      # # calc niches
      # if var_code == 'MAT':
      #   nichefrac_columns[var_code][i] = grid.calc_fraction_in_niche(min_hist, max_hist, hist_window, nichemode='minmax')
      # elif var_code == 'MTWM':
      #   nichefrac_columns[var_code][i] = grid.calc_fraction_in_niche(min_hist, max_hist, hist_window, nichemode='max')
      # elif var_code == 'MTCM':
      #   nichefrac_columns[var_code][i] = grid.calc_fraction_in_niche(min_hist, max_hist, hist_window, nichemode='min')
      # elif var_code == 'AP':
      #   nichefrac_columns[var_code][i] = grid.calc_fraction_in_niche(min_hist, max_hist, hist_window, nichemode='minmax')
      # elif var_code == 'PDQ':
      #   nichefrac_columns[var_code][i] = grid.calc_fraction_in_niche(min_hist, max_hist, hist_window, nichemode='min')
          

  if i%1000 == 0:
    print(i, 'species done')
    
# replace the geometry column in the original dataframe 
species_data = species_data.drop(['geometry'], axis=1)

for var_code in var_codes:
    species_data[f'{var_code}_nichefrac_{year}'] = nichefrac_columns[var_code]
    species_data[f'{var_code}_{year}'] = columns[var_code]
    species_data[f'{var_code}_change_{year}'] = change_columns[var_code]

# save new dataframe
# species_data.to_csv('/home/claussar/gaa2_climate_{}.csv'.format(year))
species_data.to_csv('RedList_Climate/validation/NEWgaa2_climate_{}.csv'.format(year))
t5 = time()
print('finished in {} s'.format(t5-t1))
