import numpy as np
import xarray as xr
import glob
from datetime import date
from time import time 
t1 = time()

# get all temperature files 1940-2022
varn='tas'
folder='/p/projects/compacts/data/climate/ERA5/' + varn + '/'
stem=varn + '_day_ECMWF-ERA5_observation_'
files=np.sort(glob.glob(folder + stem + '*'))

def bring(filename):
    '''Opens and preprocesses file'''
    T_year = xr.open_dataset(filename)

    # convert to -180, 180 longiude
    T_year['longitude'] = (T_year['longitude'] + 180) % 360 - 180
    T_year = T_year.sortby('longitude')

    # put latitude in correct order (-90 to 90) instead of other way around
    T_year = T_year.sel(latitude=T_year['latitude'][::-1])

    # from K to °C
    T_year['tas'] = T_year['tas'] - 273.15
    return T_year

    # write crs

def monthly_stat(T_year, year, stat):
    '''Calculates monthly statistics'''
    # Calc monthly stat
    if stat == 'min':
        Tstat = T_year.groupby('time.month').min('time') #(F) USGS2012 suggests mean of daily minima!
    elif stat == 'avg':
        Tstat = T_year.groupby('time.month').mean('time')
    elif stat == 'max':
        Tstat = T_year.groupby('time.month').max('time')
        
    # combine month and year to "time" dimension
    Tstat['time'] = np.array([np.datetime64(date(year, m, 1)) for m in Tstat.month.values])
    
    # drop unnecessary month dimension
    Tstat = Tstat.drop('month')
    Tstat = Tstat.isel(month=0)
    
    return Tstat

    
# loop through years and compile one dataset of all monthly statistics
for f, filename in enumerate(files):
    
    T_year = bring(filename)
    year = int(filename.split('.')[0].split('_')[4])

    Tmin = monthly_stat(T_year, year, 'min')
    Tavg = monthly_stat(T_year, year, 'avg')
    Tmax = monthly_stat(T_year, year, 'max')
    
    # if f=0: merge into one dataset, if f>0: append to existing dataset
    if f == 0:
        Tmin_all = Tmin
        Tavg_all = Tavg
        Tmax_all = Tmax
    else:
        Tmin_all = xr.concat([Tmin_all, Tmin], dim='time')
        Tavg_all = xr.concat([Tavg_all, Tavg], dim='time')
        Tmax_all = xr.concat([Tmax_all, Tmax], dim='time')
        
        
# rename dimensions from tas to min, avg, max
Tmin_all = Tmin_all.rename({'tas':'Tmin'})
Tavg_all = Tavg_all.rename({'tas':'Tavg'})
Tmax_all = Tmax_all.rename({'tas':'Tmax'})

# merge all T's into one dataset
T_all = xr.merge([Tmin_all.Tmin, Tavg_all.Tavg, Tmax_all.Tmax])
T_all.to_netcdf('/home/claussar/ERA5monthly/ERA5_2mtemp_1940-2022.nc')

t2 = time()
print('finished in', t2-t1)