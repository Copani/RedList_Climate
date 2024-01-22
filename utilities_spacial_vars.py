import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rasterio.features import rasterize

def prep_for_gridcalc(dataarray, gdf_i, window_only=False, margin=1):
  '''This function cuts the data down to a window of rangemap bounds and
  rasterizes the range map to this window, it also gives cell coefficients important for
  consistent size measutres (see [wiki](https://github.com/Copani/RedList_Climate/wiki))

  Parameters
  ---
  - dataarray: a rioxr georeferenced dataarray with a lon, lat grid
  - gdf_i: a single line of a geodataframe with a geometry
  - window_only (opt): bool, if True, only the data_window is calculated and gdf_i rasterization aswell as cell coefs are skipped. Only one output if activated. 
  If you want to calculate several variables on the same grid, you may only want to do the most timeconsuming opertation in this function once.
  - margin (opt): margin around gdf_i bounds of the raster window. Should be at least twice the gridsize to avoid errors!
  
  Returns
  ---
  - data_window: xr.dataarray of data in the window with same shape as geo_mask, nan replaced by 0 (for values, apply data_windows.values)
  - geo_mask: 1,0 np.array of shape burned into the grid 
  - cell_coefs: vector of cos(lat) coefficient for consistent area calculation
  '''
  if 'latitude' in dataarray.coords:
    lat = 'latitude'
  elif 'lat' in dataarray.coords:
    lat = 'lat'
  if 'longitude' in dataarray.coords:
    lon = 'longitude'
  elif 'lon' in dataarray.coords:
    lon = 'lon' 

  ### SOME UNIT TESTS ###
  # Assert that gdf_i is a GeoDataFrame
  assert isinstance(gdf_i, gpd.GeoDataFrame), "gdf_i is not a GeoDataFrame"

  # Assert that gdf_i has a length of 1
  assert len(gdf_i) == 1, "gdf_i does not have a length of 1"

  # Assert that gdf_i and the dataarray are in the same crs
  assert dataarray.rio.crs == gdf_i.crs, "dataarray and gdf_i don't have the same crs" 

  # Assert if lat lon expressions exist
  assert (lat in dataarray.coords) & (lon in dataarray.coords), "there are no correct latitude/longitude coordinates in dataarray"

  # Correct if latitude and longitude are in wrong order
  if dataarray[lat].values[0] > dataarray[lat].values[-1]:
    dataarray = dataarray.sortby(lat) 
  if dataarray[lon].values[0] > dataarray[lon].values[-1]:
    dataarray = dataarray.sortby(lon)
  

  ### FUNCTION ###
  # get radian grid resolution for later
  latres = (dataarray[lat].values[1] - dataarray[lat].values[0]) / 360 * 2 * np.pi
  lonres = (dataarray[lon].values[1] - dataarray[lon].values[0]) / 360 * 2 * np.pi
  
  # get data slices for poly window
  xmin, ymin, xmax, ymax = gdf_i.bounds.values[0]
  xmin -= margin
  ymin -= margin
  xmax += margin
  ymax += margin

  # pick a slice from the data within geometry bounds
  try:
    data_window = dataarray.sel(lon = slice(xmin, xmax), lat = slice(ymin, ymax))
  except:
    data_window = dataarray.sel(longitude = slice(xmin, xmax), latitude = slice(ymin, ymax))

  
  # for calculations set nan to 0
  data_window0 = data_window.fillna(0)


  if window_only: 
    # skip other gridding of geometry and cell_coefs, only output data window
    return data_window0
  else:
    # rasterize the polygon according to window.
    geo_mask = rasterize([(gdf_i['geometry'].iloc[0])],
                      transform=data_window.rio.transform(),
                      out_shape=(data_window[lat].shape[0], data_window[lon].shape[0]),  
                      all_touched=True)

    # weight data cells by grid cell size
    try:
      lats = np.expand_dims(data_window.lat,axis=1) / 360 * 2 * np.pi
    except:
      lats = np.expand_dims(data_window.latitude,axis=1) / 360 * 2 * np.pi
      
    R = 6371 # Earth radius [km]  
    cell_coefs = R**2 * lonres * latres * np.cos(lats)

    return data_window0, geo_mask, cell_coefs

def weighted_percentile(a, q, weights=None, interpolation='step'):
    """
    Compute the qth percentile of the data a, optionally weight can be provided.
    Returns the qth percentile(s) of the array elements.
    Credit: https://gist.github.com/righthandabacus/b125c666a26a4e6e1e9ef1a19b5da4a1

    Methodology
    -----------
    If weights are not provided, we set all `a` of equal weight of 1. Then we
    normalize the weight by equal factor so that their sum is 1. Then, in sorted
    ascending order of `a`, we plot the values as a curve from 0 to 1 and lookup
    the values corresponding to `q` from the curve.
    Shape of the curve is determined by the parameter `interpolation`. If it is
    'step', the curve is cadlag steps; if 'lower', we set the leftmost edge of
    each step as the corresponding value in `a` and interpolate the adjacent
    values except the last one, which we carry the horizontal step forward to
    1.0; if 'higher', it is similar to the case of 'lower' but we set the value
    at the rightmost edge of each step instead and the horizontal step is
    preserved at the minimum value; if 'midpoint', we set the value at the
    middle of each step and the half steps at the minimum and maximum is
    preserved as horizontal.

    Parameters
    ----------
    a : array_like of float
        Input array or object that can be converted to an array.
    q : array_like of float
        Percentile or sequence of percentiles to compute, which must be between
        0 and 100 inclusive.
    weights : array_like of float, optional
        if provided, must be the same dimension as `a` and all elements are
        non-negative. This is the weights to be used
    interpolation : {'step', 'lower', 'higher', 'midpoint'}
    
    Returns
    -------
    percentile : scalar or ndarray
        If `q` is a single percentile and `axis=None`, then the result
        is a scalar. If multiple percentiles are given, first axis of
        the result corresponds to the percentiles. The other axes are
        the axes that remain after the reduction of `a`. If the input
        contains integers or floats smaller than ``float64``, the output
        data-type is ``float64``. Otherwise, the output data-type is the
        same as that of the input. If `out` is specified, that array is
        returned instead.
    """
    # sanitation check on a, q, weights
    a = np.asarray(a).flatten()
    q = np.true_divide(q, 100.0)  # handles the asarray for us too
    if q.max() > 100 or q.min() < 0:
        raise ValueError("Percentiles must be in the range [0, 100]")
    if weights is None:
        weights = np.repeat(1, a.shape)
    weights = np.asarray(weights).flatten()
    if weights.min() < 0:
        raise ValueError("Weights must be non-negative")
    if weights.max() <= 0:
        print(a)
        print(q)
        print(weights)
        raise ValueError("Total weight must be positive")
    weights = np.true_divide(weights, weights.sum())
    if weights.shape != a.shape:
        raise ValueError("Weights and input are not in the same shape")
    # sort a and weights, remove zero weights, then convert weights into cumsum
    a, weights = zip(*sorted([(a_, w_) for a_, w_ in zip(a, weights) if w_ > 0]))
    weights = np.cumsum(weights)
    # depends on the interpolation parameter, modify the vectors
    if interpolation == 'step':
        x = np.ravel(np.column_stack((a,a)))
        w = np.insert(np.ravel(np.column_stack((weights,weights)))[:-1], 0, 0)
    elif interpolation == 'lower':
        x = np.insert(a, len(a), a[-1])
        w = np.insert(weights, 0, 0)
    elif interpolation == 'higher':
        x = np.insert(a, 0, a[0])
        w = np.insert(weights, 0, 0)
    elif interpolation == 'midpoint':
        x = np.insert(np.insert(a, len(a), a[-1]), 0, a[0])
        w = [(p+q)/2 for p,q in zip(weights, weights[1:])]
        w = np.insert(np.insert(w, len(w), 1.0), 0, 0.0)
    else:
        raise NotImplementedError("Unknown interpolation method")
    # linear search of weights by each element of q
    # TODO we can do binary search instead
    output = []
    for percentile in ([q] if isinstance(q, (int, float)) else q):
        if percentile <= 0:
            output.append(x[0])
        elif percentile >= 1.0:
            output.append(x[-1])
        else:
            for i, w2 in enumerate(w):
                if w2 == percentile:
                    output.append(x[i])
                    break
                elif w2 > percentile:
                    w1 = w[i-1]
                    x1, x2 = x[i-1], x[i]
                    output.append((x2-x1)*(percentile-w1)/(w2-w1) + x1)
                    break
    return output[0] if isinstance(q, (int, float)) else np.array(output)

def calc_percentiles(percents, dataarray, gdf_i, geo_mask=None, cell_coefs=None, margin=1):
  '''Calculate quantiles from xr.dataarray and pd.GeoDataFrame. 
  
  Parameters
  ---
  percents : array_like of float
    Percentile or sequence of percentiles to compute, which must be between
    0 and 100 inclusive.
  dataarray : xarray.core.dataarray.DataArray
    With coordinates "lon", "lat" or "longitude", "latitude", with crs
  gdf_i: single lined GeoDataFrame
    if you have multiple lines, use gdf.iloc[[i]]
  geo_mask (opt) : 2D np.array
    output of prep_for_gridcalc function
  cell_coefs (opt): 2D np.array
    output of prep_for_gridcalc function
  margin (opt): float
    margin of window calculation (min. 2 * resolution recommended)

  Returns
  ---
  percentiles : np.array
    calculated values for percentiles given in percents

  '''

  # Assert that geo_mask and cell_coefs are either both values or both none
  assert ((geo_mask is None) and (cell_coefs is None)) or ((geo_mask is not None) and (cell_coefs is not None)), 'both geo_mask and cell_coefs needed as input'

  # depending on input, calc data window and perhaps geo mask and cell coefficients
  if geo_mask is None:
    data_window, geo_mask, cell_coefs = prep_for_gridcalc(dataarray, gdf_i, margin=margin)
  else:
    data_window = prep_for_gridcalc(dataarray, gdf_i, window_only=True, margin=margin)

  # reformulate in temrs of values and coefficients (prop. to grid cell size of value)
  values = data_window.values[geo_mask.astype(bool)]
  coefficients = (cell_coefs * np.ones(data_window.shape[1]))[geo_mask.astype(bool)]

  percentiles = weighted_percentile(values, percents, weights=coefficients, interpolation='step')

  return percentiles

def calc_niche(data, nsig=3):
  '''Calculate min and max excluding values above and below nsig standard seviations away from mean
  
  Input
  ---
  data: array like
    should represent data window with time and space dim. otherwise there is danger of a too narrow niche
  nsig (opt): float
    values outside nsig standard deviations from the mean are excluded
    
  Output
  ---
  '''
  # print warning if data array is not 2dimensional (later)

  # apply outlyer exclusion mask
  mask = ((data.values < data.values.mean() + nsig * data.values.std()) & 
          (data.values > data.values.mean() - nsig * data.values.std()) )
  
  data_without_outlyers = data.values[mask]
  
  data_max = np.max(data_without_outlyers)
  data_min = np.min(data_without_outlyers)
  
  return data_min, data_max
  

def calc_fraction_in_niche(min_hist, max_hist, dataarray, gdf_i, geo_mask=None, cell_coefs=None, nyears=5, margin=1):
  '''Calculate fraction of grid cells exceeding the niche.
  
  Parameters
  ---
  max_hist: float
    maximum value of niche
  min_hist: float
    minimum value of niche
  dataarray : xarray.core.dataarray.DataArray
    With coordinates "lon", "lat" or "longitude", "latitude", with crs
  gdf_i: single lined GeoDataFrame
    if you have multiple lines, use gdf.iloc[[i]]
  geo_mask (opt) : 2D np.array
    output of prec_for_gridcalc function
  cell_coefs (opt): 2D np.array
    output of prec_for_gridcalc function
  nyears (opt): int
    number of years niche must be exeeded within time slice to "count"
  margin (opt): float
    margin of window calculation (min. 2 * resolution recommended)
  
  Returns
  ---
  frac_in_niche: float
    fraction of grid cells within the niche
  frac_above_max: float
    fraction of grid cells exceeding the niche
  frac_below_min: float
    fraction of grid cells below the minimum of the niche
  '''

  # Assert that geo_mask and cell_coefs are either both values or both none
  assert ((geo_mask is None) and (cell_coefs is None)) or ((geo_mask is not None) and (cell_coefs is not None)), 'both geo_mask and cell_coefs needed as input'

  # depending on input, calc data window and perhaps geo mask and cell coefficients
  if geo_mask is None:
    data_window, geo_mask, cell_coefs = prep_for_gridcalc(dataarray, gdf_i, margin=margin)
  else:
    data_window = prep_for_gridcalc(dataarray, gdf_i, window_only=True, margin=margin)

  # calc which grid cells have exceeded threshold (at which times)
  exceed_min_bool = (data_window < min_hist)
  exceed_max_bool = (data_window > max_hist)
  
  # calc which grid cells have exceeded threshold for more then 5 years 
  exceed_min_5y_bool = exceed_min_bool.sum('year') >= nyears
  exceed_max_5y_bool = exceed_max_bool.sum('year') >= nyears
  
  # reformulate in temrs of values and coefficients (prop. to grid cell size of value)
  exceed_min_5y_bool_range = exceed_min_5y_bool.values[geo_mask.astype(bool)]
  exceed_max_5y_bool_range = exceed_max_5y_bool.values[geo_mask.astype(bool)]
  coefficients = (cell_coefs * np.ones(data_window.shape[-1]))[geo_mask.astype(bool)]

  # calc area fractions
  frac_above_max = np.sum(coefficients * exceed_max_5y_bool_range) / np.sum(coefficients)
  frac_below_min = np.sum(coefficients * exceed_min_5y_bool_range) / np.sum(coefficients)
  
  frac_in_niche = 1 - frac_above_max - frac_below_min
  
  return frac_in_niche, frac_above_max, frac_below_min

