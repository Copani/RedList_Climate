import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rasterio.features import rasterize

def prep_for_gridcalc(dataarray, gdf_i, window_only=False, margin=1):
  '''# Prepare for grid calculation
  This function cuts the data down to a window of rangemap bounds and
  rasterizes the range map to this window, it also gives cell coefficients important
  consistent size measutres (see [wiki](https://github.com/Copani/RedList_Climate/wiki))

  INPUTS
  - gdf_i: a single line of a geodataframe with a geometry
  - dataarray: a rioxr georeferenced dataarray with a lon, lat grid
  - window_only (opt): bool, if True, only the data_window is calculated and gdf_i rasterization aswell as coefficient calculation are skipped. Only one output if activated. 
  If you want to calculate several variables on the same grid, you may only want to do the most timeconsuming opertation in this function once.
  - margin (opt): margin around gdf_i bounds of the raster window. Should be at least twice the gridsize to avoid errors!
  
  OUTPUTS
  - data_window: np.array of data in the window with same shape as geo_mask, nan replaced by 0 (not georeferenced anymore, you may handle this afterwards, but it is usually not needed for calculations)
  - geo_mask: 1,0 np.array of shape burned into the grid (not georef., see above.)
  - cell_coefs: vector of cos(lat) coefficient for consistent area calculation
  '''

  ### SOME UNIT TESTS ###
  # Assert that gdf_i is a GeoDataFrame
  assert isinstance(gdf_i, gpd.GeoDataFrame), "gdf_i is not a GeoDataFrame"

  # Assert that gdf_i has a length of 1
  assert len(gdf_i) == 1, "gdf_i does not have a length of 1"

  # Assert that gdf_i and the dataarray are in the same crs
  assert dataarray.rio.crs == gdf_i.crs, "dataarray and gdf_i don't have the same crs" 

  # Assert if dataarray has the axes lat and lon
  assert ("lat" in dataarray.coords) & ("lon" in dataarray.coords), "The 'lat' or 'lon' axis is missing in dataarray. Make sure axes are named correctly."

  ### FUNCTION ###
  # get data slices for poly window
  xmin, ymin, xmax, ymax = gdf_i.bounds.values[0]
  xmin -= margin
  ymin -= margin
  xmax += margin
  ymax += margin

  # pick a slice from the data within geometry bounds
  data_window = dataarray.sel(lon = slice(xmin, xmax), lat = slice(ymax, ymin))
  
  # for calculations set nan to 0
  data_window0 = data_window.fillna(0).values

  if not window_only:
    # rasterize the polygon according to window.
    geo_mask = rasterize([(gdf_i['geometry'].iloc[0])],
                      transform=data_window.rio.transform(),
                      out_shape=data_window.shape, 
                      all_touched=True)

    # weight data cells by grid cell size
    lats = np.expand_dims(data_window.lat,axis=1) / 360 * 2 * np.pi
    cell_coefs = np.cos(lats)

    return data_window0, geo_mask, cell_coefs
  
  else: 
    # skip other calculations and only output data window
    return data_window0
  
