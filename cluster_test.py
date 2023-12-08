from netCDF4 import Dataset
import numpy as np


x = Dataset('/p/projects/climate_data_central/observation/Berkeley_Earth/gridded/monthly_land/Complete_TAVG_LatLong1.nc')
lon=x.variables['longitude'][:].data
np.save('/home/claussar/test_lon', lon)
