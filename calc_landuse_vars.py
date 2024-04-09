import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rasterio.features import rasterize
from time import time
from utilities_spacial_vars import prep_for_gridcalc


# pop = xr.open_dataset("Data\HYDE_landuse_pop_1901_2021\population_histsoc_30arcmin_annual_1901_2021.nc")
# area=111*0.5*111*0.5*np.cos(lat*np.pi/180)
# pop_dens=np.swapaxes(np.divide(np.swapaxes(pop,1,2),area),1,2)
# pop_dens=pop_dens.filled(fill_value=np.nan)

# #log scaling + capping at 1000 people per km2
# pop_pressure=3.333*np.log(pop_dens+1)
# pop_pressure[pop_dens>1000]=10

rows = 5 # of rows = 0, all rows are loaded
year = 2022
t1 = time()

# Load and format landuse dataset
## load urban areas, cropland pastures
# urban = xr.open_dataset("../Data/HYDE_landuse_pop_1901_2021/landuse-urbanareas_histsoc_annual_1901_2021.nc", decode_coords="all", decode_times=False)
# crop = xr.open_dataset("Data/HYDE_landuse_pop_1901_2021/landuse-totals_histsoc_annual_1901_2021.nc", decode_coords="all", decode_times=False)
# pasture = xr.open_dataset("../Data/HYDE_landuse_pop_1901_2021/landuse-pastures_histsoc_annual_1901_2021.nc", decode_coords="all", decode_times=False)
pop = xr.open_dataset('../Data/HYDE_landuse_pop_1901_2021/population_histsoc_30arcmin_annual_1901_2021.nc')

## normalize the years 1901 - 2021
# urban['time'] = (urban['time'] + 1901 - pasture['time'][0]).astype(int)
# crop['time'] = (crop['time'] + 1901 - pasture['time'][0]).astype(int)
# pasture['time'] = (pasture['time'] + 1901 - pasture['time'][0]).astype(int)

# for population, sort years


## write crs WSG 84 (lon, lat)
# urban = urban.rio.write_crs("EPSG:4326")
# crop = crop.rio.write_crs("EPSG:4326")
# pasture = pasture.rio.write_crs("EPSG:4326")
pop = pop.rio.write_crs("EPSG:4326")

## extract the variables, set year
# urban = urban.sel(time=year).urbanareas
# crop = crop.sel(time=year).cropland_total
# rangeland = pasture.sel(time=year).rangeland
# managed_pasture = pasture.sel(time=year).managed_pastures
pop = pop['total-population'].sel(time=pop['time.year']==year)[0,:,:]

# take the logp1 of the population
pop = np.log(pop + 1)


t2 = time()
print('landuse data loaded in {} s'.format(t2-t1))
  
# Load species data
# species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/IUCN/AMPHIBIANS/AMPHIBIANS.shp')
if rows == 0:
  species_data = gpd.read_file('../Data/GAA2/GAA2_allspecies.gdb')
  rows = len(species_data)
else:
  species_data = gpd.read_file('../Data/GAA2/GAA2_allspecies.gdb', rows=rows)

t3 = time()
print('{} rows of species data loaded in {} s'.format(rows, t3-t2))


## for each species in this dataset, get grid, calc quantities
# urban_column = [None] * rows
# crop_column = [None] * rows
# rangeland_column = [None] * rows
# managed_pasture_column = [None] * rows
pop_density_column = [None] * rows
area_column = [None] * rows

for i in range(rows):
  species_i = species_data.iloc[[i]]

  ## check if geometry is valid
  rmap = species_i.iloc[0].explode()
  if rmap.geometry.is_valid:
       
    ## rasterize for calculation using urban dataset
    # urban_window, geo_mask, cell_coefs = prep_for_gridcalc(urban, species_i)
    pop_window, geo_mask, cell_coefs = prep_for_gridcalc(pop, species_i)

    ## rasterize other quantities
    # crop_window = prep_for_gridcalc(crop, species_i, window_only=True)
    # rangeland_window = prep_for_gridcalc(rangeland, species_i, window_only=True)
    # managed_pasture_window = prep_for_gridcalc(managed_pasture, species_i, window_only=True)

    ## calculate mean quantities in species range
    # urban_column[i] = np.sum(geo_mask * urban_window * cell_coefs) / np.sum(geo_mask * cell_coefs)
    # crop_column[i] = np.sum(geo_mask * crop_window * cell_coefs) / np.sum(geo_mask * cell_coefs)
    # rangeland_column[i] = np.sum(geo_mask * rangeland_window * cell_coefs) / np.sum(geo_mask * cell_coefs)
    # managed_pasture_column[i] = np.sum(geo_mask * managed_pasture_window * cell_coefs) / np.sum(geo_mask * cell_coefs)
    pop_density_column[i] = np.sum(geo_mask * pop_window) / np.sum(geo_mask * cell_coefs)
    area_column[i] = np.sum(geo_mask * cell_coefs)

# replace the geometry column in the original dataframe 
# species_data = species_data.drop(['geometry'], axis=1)
# or
# add data to additional columns in existing dataframe
del species_data
species_data = pd.read_csv('gaa2_l2021_lpa2004_lpa1980.csv')

# species_data['urbanareas_{}'.format(year)] = urban_column
# species_data['cropland_{}'.format(year)] = crop_column
# species_data['rangeland_{}'.format(year)] = rangeland_column
# species_data['managed_pasture_{}'.format(year)] = managed_pasture_column
species_data['population_density_{}'.format(year)] = pop_density_column
species_data['area_{}'.format(year)] = area_column
t4 = time()

print('qantities for {} species calculated in {} s'.format(rows, t4-t3))

# save new dataframe
# species_data.to_csv('/home/claussar/IUCN_range_coords/IUCN_AMPHIBIA_coords.csv')
# species_data.to_csv('gaa2_lpa2021_lpa2004_lpa1980.csv')

print('finished in {} s'.format(t4-t1))
