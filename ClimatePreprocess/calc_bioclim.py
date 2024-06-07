import numpy as np
import xarray as xr

def calc_bioclim(T, P, years):
  """
  Calculate bioclimatic variables for a given range of years according to 
  O’Donnell, M.S., and Ignizio, D.A., 2012, Bioclimatic predictors for supporting ecological applications in the conterminous United States: U.S. Geological Survey Data Series 691, 10 p.
  
  !! Note that since quarters are evaluated for three consecutive months, the last year is not reliably evaluated (two months of the following year would be required). 

  Parameters:
  T (xarray.Dataset): Dataset containing temperature data with variables 'Tavg', 'Tmax', 'Tmin'.
  P (xarray.Dataset): Dataset containing precipitation data with variable 'tp'.
  years (tuple): A tuple containing the start and end years for the calculation.

  Returns:
  dict: A dictionary containing the calculated bioclimatic variables:
      'MAT' - Mean Annual Temperature
      'TS' - Temperature Seasonality
      'MTWM' - Maximum Temperature of Warmest Month
      'MTCM' - Minimum Temperature of Coldest Month
      'AP' - Total Annual Precipitation
      'PS' - Precipitation Seasonality
      'PWeQ' - Precipitation of Wettest Quarter
      'PDQ' - Precipitation of Driest Quarter
      'PWaQ' - Precipitation of Warmest Quarter
  """
  start_year, end_year = years
  
  T = T.sel(time=slice(str(start_year), str(end_year)))
  P = P.sel(time=slice(str(start_year), str(end_year)))

  # Calculate mean annual temperature (MAT) averages
  MAT = T.Tavg.groupby('time.year').mean('time')
  MAT = MAT.drop('spatial_ref').rio.write_crs("EPSG:4326")
  
  # Calculate temperature seasonality (TS)
  TS = T.Tavg.groupby('time.year').std('time')
  TS = TS.drop('spatial_ref').rio.write_crs("EPSG:4326")
  
  # Calculate max temperature of warmest month (MTWM)
  MTWM = T.Tmax.groupby('time.year').max('time')
  MTWM = MTWM.drop('spatial_ref').rio.write_crs("EPSG:4326")

  # Calculate min temperature of coldest month (MTCM)
  MTCM = T.Tmin.groupby('time.year').min('time')
  MTCM = MTCM.drop('spatial_ref').rio.write_crs("EPSG:4326")

  # Calculate total annual precipitation (AP)
  AP = P.tp.groupby('time.year').sum('time')[:,0]
  AP = AP.drop('spatial_ref').drop('expver').rio.write_crs("EPSG:4326")
  
  # Calculate precipitation seasonality (PS)
  PS = P.tp.groupby('time.year').std('time')[:,0] / (1 + AP / 12) * 100
  PS = PS.drop('spatial_ref').drop('expver').rio.write_crs("EPSG:4326")
  
  # Calculate precipitation of wettest quarter (PWeQ)
  Prec_quaterly = P.tp.rolling(time=3).sum().shift(time=-2) # shift to sum over the next (instead of previous) 3 months  
  PWeQ = Prec_quaterly.groupby('time.year').max()[:,0]
  PWeQ = PWeQ.drop('spatial_ref').drop('expver').rio.write_crs("EPSG:4326")

  # Calculate precipitation of driest quarter (PDrQ)
  PDQ = Prec_quaterly.groupby('time.year').min()[:,0]
  PDQ = PDQ.drop('spatial_ref').drop('expver').rio.write_crs("EPSG:4326")
  
  # Calculate precipitation of warmest quarter (PWaQ)
  QTmax = T.Tavg.rolling(time=3).sum().shift(time=-2) # get T quarters
  QTmax_idxfirstmonth = QTmax.groupby('time.year').apply(xr.DataArray.idxmax, dim='time')
  PWaQ = Prec_quaterly[:,0].sel(time=QTmax_idxfirstmonth)
  PWaQ = PWaQ.drop('time').drop('expver').drop('spatial_ref').rio.write_crs("EPSG:4326")
  

  return {
    'MAT': MAT,
    'TS': TS,
    'MTWM': MTWM,
    'MTCM': MTCM,
    'AP': AP,
    'PS': PS,
    'PWeQ': PWeQ,
    'PDQ': PDQ,
    'PWaQ': PWaQ
  }