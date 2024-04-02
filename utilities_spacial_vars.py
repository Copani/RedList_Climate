import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rasterio.features import rasterize

class GridCalc:
  def __init__(self, dataarray, gdf_i, margin=1):
    self.dataarray = dataarray
    self.gdf_i = gdf_i
    self.margin = margin
    
    # calc species range map bounds
    self.xmin, self.ymin, self.xmax, self.ymax = self.gdf_i.bounds.values[0]
    self.xmin -= self.margin
    self.ymin -= self.margin
    self.xmax += self.margin
    self.ymax += self.margin
    
    # get lat lon expressions
    if 'latitude' in self.dataarray.coords:
      self.lat = 'latitude'
    elif 'lat' in self.dataarray.coords:
      self.lat = 'lat'
    elif 'y' in self.dataarray.coords:
      self.lat = 'y'
    if 'longitude' in self.dataarray.coords:
      self.lon = 'longitude'
    elif 'lon' in self.dataarray.coords:
      self.lon = 'lon'
    elif 'x' in self.dataarray.coords:
      self.lon = 'x'
    
      
    # get window, mask and coefficients
    self.data_window = self.get_window(self.dataarray)
    self.mask = self.get_mask()
    self.cell_coefs = self.get_coefs()
    
    
    ### SOME UNIT TESTS ###
    # Assert that gdf_i is a GeoDataFrame
    assert isinstance(self.gdf_i, gpd.GeoDataFrame), "gdf_i is not a GeoDataFrame"

    # Assert that gdf_i has a length of 1
    assert len(self.gdf_i) == 1, "gdf_i does not have a length of 1"

    # Assert that gdf_i and the dataarray are in the same crs
    assert self.dataarray.rio.crs == self.gdf_i.crs, "dataarray and gdf_i don't have the same crs" 

    # Assert if lat lon expressions exist
    assert (self.lat in self.dataarray.coords) & (self.lon in dataarray.coords), "there are no correct latitude/longitude coordinates in self.dataarray"

    # Correct if latitude and longitude are in wrong order
    if self.dataarray[self.lat].values[0] > self.dataarray[self.lat].values[-1]:
      self.dataarray = self.dataarray.sortby(self.lat) 
    if self.dataarray[self.lon].values[0] > self.dataarray[self.lon].values[-1]:
      self.dataarray = self.dataarray.sortby(self.lon)
  
     

  def get_window(self, dataarray):
    '''This function cuts the data down to a window of rangemap bounds. 
    Note that any dataarray can be used as an input.

    Returns
    - data_window: xr.dataarray of data, nan replaced by 0 (for values, apply data_windows.values)
    '''

    # pick a slice from the data within geometry bounds
    coordinate_names = [('lon', 'lat'), ('longitude', 'latitude'), ('x', 'y')]

    for coord in coordinate_names:
      try:
        data_window = dataarray.sel(**{coord[0]: slice(self.xmin, self.xmax), coord[1]: slice(self.ymin, self.ymax)})
        break
      except:
        continue
    # for calculations set nan to 0
    data_window = data_window.fillna(0)
    
    return data_window
  

  def get_mask(self):
    '''This function rasterizes the range map to the window by creating a 1,0 mask of the window shape.
    Every cell that touches the geometry is set to 1.
    
    Returns
    ---
    geo_mask: 1,0 np.array of shape burned into the grid 
    '''
    # rasterize the polygon according to window.
    if self.lon == 'x':
      geo_mask = rasterize([(self.gdf_i['geometry'].iloc[0])],
                      transform=self.data_window.rio.transform('xy'),
                      out_shape=(self.data_window[self.lat].shape[0], self.data_window[self.lon].shape[0]),  
                      all_touched=True)
      return geo_mask
      
    else:
      geo_mask = rasterize([(self.gdf_i['geometry'].iloc[0])],
                        transform=self.data_window.rio.transform(),
                        out_shape=(self.data_window[self.lat].shape[0], self.data_window[self.lon].shape[0]),  
                        all_touched=True)
      return geo_mask
    
  
  def get_coefs(self):
    '''This function gives cell coefficients important for
    consistent size measutres (see [wiki](https://github.com/Copani/RedList_Climate/wiki))

    Returns
    ---
    cell_coefs: 2D np.array of shape burned into the grid, with coefficients for each cell
    '''
    # get radian grid resolution
    latres = (self.dataarray[self.lat].values[1] - self.dataarray[self.lat].values[0]) / 360 * 2 * np.pi
    lonres = (self.dataarray[self.lon].values[1] - self.dataarray[self.lon].values[0]) / 360 * 2 * np.pi
    
    # weight data cells by grid cell size
    attribute_names = ['lat', 'latitude', 'y']

    for attr in attribute_names:
        try:
            lats = np.expand_dims(getattr(self.data_window, attr), axis=1) / 360 * 2 * np.pi
            break
        except:
            continue
      
    R = 6371 # Earth radius [km]  
    cell_coefs = R**2 * lonres * latres * np.cos(lats)

    return cell_coefs

  def calc_average(self, dataarray):
    return np.sum(self.mask * dataarray.values * self.cell_coefs) / np.sum(self.mask * self.cell_coefs)
    

  def calc_fraction_in_niche(self, min_hist, max_hist, data_window, nyears=5, nichemode="minmax"):
    '''Calculate fraction of grid cells exceeding the niche.
    
    Parameters
    ---
    max_hist: float
      maximum value of niche
    min_hist: float
      minimum value of niche
    dataarray : xarray.core.dataarray.DataArray
      With coordinates "lon", "lat" or "longitude", "latitude", with crs
    nyears (opt): int
      number of years niche must be exeeded within time slice to "count"
    nichemode (opt): str
      "minmax" (default), "min" or "max" for returning fraction of grid cells exceeding the niche, below the minimum or above the maximum
    
    Returns
    ---
    frac_in_niche: float
      fraction of grid cells within the niche
    frac_above_max: float
      fraction of grid cells exceeding the niche
    frac_below_min: float
      fraction of grid cells below the minimum of the niche
    '''

    # calc which grid cells have exceeded threshold (at which times)
    exceed_min_bool = (data_window < min_hist)
    exceed_max_bool = (data_window > max_hist)
    
    # calc which grid cells have exceeded threshold for more then 5 years within 20y window
    exceed_min_5y_bool = exceed_min_bool.sum('year') >= nyears
    exceed_max_5y_bool = exceed_max_bool.sum('year') >= nyears
    
    # reformulate in temrs of values and coefficients (prop. to grid cell size of value)
    exceed_min_5y_bool_range = exceed_min_5y_bool.values[self.mask.astype(bool)]
    exceed_max_5y_bool_range = exceed_max_5y_bool.values[self.mask.astype(bool)]
    coefficients = (self.cell_coefs * np.ones(data_window.shape[-1]))[self.mask.astype(bool)]

    # calc area fractions, handel exceptions when there is no geometry
    
    frac_above_max = np.sum(coefficients * exceed_max_5y_bool_range) / np.sum(coefficients)
    frac_below_min = np.sum(coefficients * exceed_min_5y_bool_range) / np.sum(coefficients)
    
    if nichemode == "minmax":
      return 1 - frac_above_max - frac_below_min
    elif nichemode == "min":
      return  1 - frac_below_min
    elif nichemode == "max":
      return 1 - frac_above_max


def calc_limits(data, nsig=3):
  '''Calculate min and max excluding values above and below nsig standard deviations away from mean
  
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

  # apply outlier exclusion mask
  mask_no_outliers = ((data.values < data.values.mean() + nsig * data.values.std()) & 
          (data.values > data.values.mean() - nsig * data.values.std()) )
  
  data_without_outliers = data.values[mask_no_outliers]
  
  data_max = np.max(data_without_outliers)
  data_min = np.min(data_without_outliers)
  
  return data_min, data_max

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


    
  

