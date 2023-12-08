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
  
def berkeley_data(x):
  '''# Convert berkely temperature data to dataArray

  INPUT
  x: netCDF Dataset of berkely climate data

  OUTPUT
  - xarray dataArray with maximum monthly absolute temperature [°C] in 1deg x 1deg grid (axes: longitude, latitude)
  '''

  # get imortant variables
  warnings.simplefilter("ignore")
  lon=x.variables['longitude'][:]
  lat=x.variables['latitude'][:]
  temp=x.variables['temperature'][:]
  climatology=x.variables['climatology'][:]
  time=x.variables['time'][:]

  # create useful data array of temperature deviations from month mean (climatology)
  TAVG_anomaly_data = xr.DataArray(temp.data, dims=('time', 'latitude', 'longitude'), coords=[time, lat, lon])

  # change decimal time (1891.125 := Feb. 1891) to date format
  def decimal_year_to_date(decimal_year):
    '''INPUT:    decimal_year must be a single number (no array)\\
    OUTPUT:   array of dates 
    '''

    year = int(decimal_year)
    days_in_year = (decimal_year - year) * 365.25
    time_as_date = date(year, 1, 1) + timedelta(days_in_year)

    return time_as_date

  ## apply time to array
  time_dates = [decimal_year_to_date(t) for t in time.data]
  TAVG_anomaly_data['time'] = pd.DatetimeIndex(time_dates).floor('D')

  # for absolute temperature, add climatology to each month
  months = TAVG_anomaly_data['time.month']
  TAVG_data = TAVG_anomaly_data + climatology[months - 1]

  # take max along time axis
  max_TAVG = TAVG_data.max('time')

  return max_TAVG

def max_temp_in_geo(range_grid, Tarr):
  '''# Description\\Calculate maximum temperature within range_grid 
  
  INPUTS
  - range_grid: pd df with 'latitude', 'longitude' keys for one spcies's range map
  - Tarr: xarray DataArray with temperature values for each grid cell (no temporal resolution - e.g. max. realized)
  
  OUTPUT
  - maximum temperature across gridded geometry (float) 
  '''

  # create dataframe with longitude, latitude, max temperature
  lon_mesh, lat_mesh = np.meshgrid(Tarr.longitude.data, Tarr.latitude.data)
  T_for_merging = pd.DataFrame({'longitude': lon_mesh.flatten(), 'latitude': lat_mesh.flatten(), 'temperature': Tarr.data.flatten()})

  # take only the gridcells from "grid_shape" out of the temperature array.
  T_in_geo = pd.merge(range_grid, T_for_merging, left_on=['longitude', 'latitude'], right_on=['longitude', 'latitude'])

  # maximal realized temperature (within geometry)
  mrt = T_in_geo['temperature'].max() 

  return mrt

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