import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rasterio.features import rasterize
from time import time

# rows = 1
year = 2004

t1 = time()

# Load and format landuse dataset
## load urban areas, cropland pastures
urban = xr.open_dataset("Data/HYDE_landuse_pop_1901_2021/landuse-urbanareas_histsoc_annual_1901_2021.nc", decode_coords="all", decode_times=False)
crop = xr.open_dataset("Data/HYDE_landuse_pop_1901_2021/landuse-totals_histsoc_annual_1901_2021.nc", decode_coords="all", decode_times=False)
pasture = xr.open_dataset("Data/HYDE_landuse_pop_1901_2021/landuse-pastures_histsoc_annual_1901_2021.nc", decode_coords="all", decode_times=False)

## normalize the years 1901 - 2021
urban['time'] = (urban['time'] + 1901 - pasture['time'][0]).astype(int)
crop['time'] = (crop['time'] + 1901 - pasture['time'][0]).astype(int)
pasture['time'] = (pasture['time'] + 1901 - pasture['time'][0]).astype(int)

## write crs WSG 84 (lon, lat)
urban = urban.rio.write_crs("EPSG:4326")
crop = crop.rio.write_crs("EPSG:4326")
pasture = pasture.rio.write_crs("EPSG:4326")

## extract the variables, set year
urban = urban.sel(time=year).urbanareas
crop = crop.sel(time=year).cropland_total
rangeland = pasture.sel(time=year).rangeland
managed_pasture = pasture.sel(time=year).managed_pastures


t2 = time()
print('landuse data loaded in {} s'.format(t2-t1))

  
# Load species data
# species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/IUCN/AMPHIBIANS/AMPHIBIANS.shp')
species_data = gpd.read_file('Data/GAA2/GAA2_allspecies.gdb')
rows = len(species_data)

t3 = time()
print('{} rows of species data loaded in {} s'.format(rows, t3-t2))


# Calculate quantities
## for each species in this dataset, get grid, calc quantities

urban_column = [None] * rows
crop_column = [None] * rows
rangeland_column = [None] * rows
managed_pasture_column = [None] * rows

for i in range(rows):
  species_i = species_data.iloc[[i]]

  ## check if geometry is valid
  rmap = species_i.iloc[0].explode()
  if rmap.geometry.is_valid:
  
    ## get data slices for poly window
    margin = 1
    xmin, ymin, xmax, ymax = species_i.bounds.values[0]
    xmin -= margin
    ymin -= margin
    xmax += margin
    ymax += margin

    urban_window           = urban.sel(lon = slice(xmin, xmax), lat = slice(ymax, ymin))
    crop_window            = crop.sel(lon = slice(xmin, xmax), lat = slice(ymax, ymin))
    rangeland_window       = rangeland.sel(lon = slice(xmin, xmax), lat = slice(ymax, ymin))
    managed_pasture_window = managed_pasture.sel(lon = slice(xmin, xmax), lat = slice(ymax, ymin))
    
    ## rasterize the polygon according to window. 
    ## notice that all landuse data have the same grid and therefore it's enough to use on grid.
    geo_mask = rasterize([(species_i['geometry'].iloc[0])],
                    transform=urban_window.rio.transform(),
                    out_shape=urban_window.shape, 
                    all_touched=True)
    
    ## weight data cells by grid cell size
    lats = np.expand_dims(urban_window.lat,axis=1) / 360 * 2 * np.pi
    cell_coefs = np.cos(lats)

    ## for calculations set nan to 0
    urban_window0 = urban_window.fillna(0).values
    crop_window0 = crop_window.fillna(0).values
    rangeland_window0 = rangeland_window.fillna(0).values
    managed_pasture_window0 = managed_pasture_window.fillna(0).values

    ## calculate mean quantities in species range
    urban_column[i] = np.sum(geo_mask * urban_window0 * cell_coefs) / np.sum(geo_mask * cell_coefs)
    crop_column[i] = np.sum(geo_mask * crop_window0 * cell_coefs) / np.sum(geo_mask * cell_coefs)
    rangeland_column[i] = np.sum(geo_mask * rangeland_window0 * cell_coefs) / np.sum(geo_mask * cell_coefs)
    managed_pasture_column[i] = np.sum(geo_mask * managed_pasture_window0 * cell_coefs) / np.sum(geo_mask * cell_coefs)


# replace the geometry column in the original dataframe 
species_data = species_data.drop(['geometry'], axis=1)
# or
# add data to additional columns in existing dataframe
# del species_data
# species_data = pd.read_csv('RedList_Climate\gaa2_landuse_2021.csv')

species_data['urbanareas_{}'.format(year)] = urban_column
species_data['cropland_{}'.format(year)] = crop_column
species_data['rangeland_{}'.format(year)] = rangeland_column
species_data['managed_pasture_{}'.format(year)] = managed_pasture_column
t4 = time()
print('qantities for {} species calculated in {} s'.format(rows, t4-t3))

# save new dataframe
# species_data.to_csv('/home/claussar/IUCN_range_coords/IUCN_AMPHIBIA_coords.csv')
species_data.to_csv('RedList_Climate\gaa2_landuse_2021.csv')

print('finished in {} s'.format(t4-t1))