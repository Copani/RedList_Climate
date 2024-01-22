import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from netCDF4 import Dataset
import shapely
import warnings
from datetime import date, timedelta
import pickle
from time import time

t1 = time()

def range_coords(gdf, lon, lat):
  ''' # Description range coordinates
  Takes geodataframe and converts polygon geometry to df of lon, lat coordinates

  INPUTS
  - gdf: geoDataframe with "geometry" attribute
  - lon: 1d array with longitudes
  - lat: 1d array with latitudes

  OUTPUTS
  - grid_shape: np.array with (longitude, latitude) for each grid point withing the geometry
  '''
  # TODO: Visualize examples to check if everything works correctly. 
  # Check if some areas are small and not connected
  
  # get grid points from geometry
  minx, miny, maxx, maxy = gdf['geometry'].bounds.values[0]

  lon_bounded = lon[(lon >= minx) & (lon <= maxx)]
  lat_bounded = lat[(lat >= miny) & (lat <= maxy)]

  x, y = np.meshgrid(lon_bounded, lat_bounded)
  X = x.flatten(); Y = y.flatten()
  df = pd.DataFrame({'longitude': X, 'latitude': Y})
  df_geometry = gpd.points_from_xy(df.longitude, df.latitude)
  lat_lon_bound = gpd.GeoDataFrame(df, geometry=df_geometry, crs= '{}'.format(gdf.crs))
  
  gdf_grid_bool = lat_lon_bound.within(gdf.geometry.unary_union)
  # look if speedup possible
  # if there are coordinates, store as np. array (lon, lat) - otherwise take centroid
  if np.any(gdf_grid_bool):
    coords = lat_lon_bound[gdf_grid_bool].drop(['geometry'],axis=1).values
  else: 
    warnings.simplefilter('ignore') # ignore warning that "Results from 'centroid' are likely incorrect. Use 'GeoSeries.to_crs()' to re-project geometries"
    lon_centroid = lon.data[abs(lon - np.array(gdf.centroid.x)) < 0.5] # has to be adaptive & chech if min max always in grid
    lat_centroid = lat.data[abs(lat - np.array(gdf.centroid.y)) < 0.5]
    coords = np.array([[lon_centroid[0], lat_centroid[0]]])
  
  # return array of coordinates
  return coords
  


###### CODE ######

# load climate data
x = Dataset('/p/projects/climate_data_central/observation/Berkeley_Earth/gridded/monthly_land/Complete_TAVG_LatLong1.nc')
# x = Dataset('C:/Users/Claus/Master Thesis/Data/Berkely_monthly_land/Complete_TAVG_LatLong1.nc')
lon=x.variables['longitude'][:]
lat=x.variables['latitude'][:]

t2 = time()
print('cliamte data loaded in {} s'.format(t2-t1))

# load species data
species_data = gpd.read_file('/p/projects/impactee/Data/biodiversity/IUCN/AMPHIBIANS/AMPHIBIANS.shp')
# species_data = gpd.read_file('C:/Users/Claus/Master Thesis/Data/AMPHIBIANS/AMPHIBIANS.shp')

t3 = time()
print('species data loaded in {} s'.format(t3-t2))

# for each species in this dataset, calculate range grid

count = 100 #species_data['geometry'].count()
all_coords = [None] * species_data['geometry'].count()
for i in range(count):
  species_i = species_data.iloc[[i]]

  # check if geometry is valid
  rmap = species_i.iloc[0].explode()
  if rmap.geometry.is_valid:
  
    # calc range grid and exchange geometry to lon lat range grid
    coords = range_coords(species_i, lon, lat)
    all_coords[i] = coords  

# replace the geometry column in the original dataframe 
species_data = species_data.drop(['geometry'], axis=1)
species_data['lon_lat_range'] = all_coords

t4 = time()
print('coords for {} species calculated in {} s'.format(count, t4-t3))

# save new dataframe
species_data.to_csv('/home/claussar/IUCN_range_coords/IUCN_AMPHIBIA_coords.csv')

print('finished in {} s'.format(t4-t1))