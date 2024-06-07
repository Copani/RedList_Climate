import numpy as np
import xarray as xr
import glob
from datetime import date
from time import time 
import dask
import sys

t1 = time()
# dataset = "30km"
dataset = "9km"

# process data in chunks
start_index = int(sys.argv[1])
end_index = int(sys.argv[2])


if dataset == "30km":
    # get all temperature files 1940-2022
        varn='tas'
        short_varn='tas'
        folder='/p/projects/climate_data_central/reanalysis/ERA5/' + varn + '/'
        stem=varn + '_1hr_ECMWF-ERA5_observation_'
        files=np.sort(glob.glob(folder + stem + '*'))
        
if dataset == "9km":
    # get all temperature files 1940-2022
        varn='2m_temperature'
        short_varn = 't2m'
        folder='/p/projects/climate_data_central/reanalysis/ERA5-Land/' + varn + '/'
        stem=varn + '_ERA5-Land_'
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
    T_year = T_year.chunk({'time': 40})  
    T_year[short_varn] = T_year[short_varn] - 273.15
    return T_year

    # write crs
    
# loop through years and compile one dataset of all monthly statistics
datasets = []
for f, filename in enumerate(files[start_index:end_index]):
    
    T_year = bring(filename)
    
    # calc monthly mean
    Tavg = T_year.resample(time='M', label='left').mean().compute()
    
    # calc monthly mean of daily max
    daily_max = T_year.resample(time='D', label='left').max().compute()
    Tmax = daily_max.resample(time='M', label='left').mean().compute()
    
    # calc monthly mean of daily min
    daily_min = T_year.resample(time='D', label='left').min().compute()
    Tmin = daily_min.resample(time='M', label='left').mean().compute()
    
    # rename dimensions from tas to min, avg, max
    Tmin = Tmin.rename({short_varn:'Tmin'})
    Tavg = Tavg.rename({short_varn:'Tavg'})
    Tmax = Tmax.rename({short_varn:'Tmax'})

    # merge all T's into one dataset
    T_all = xr.merge([Tmin.Tmin, Tavg.Tavg, Tmax.Tmax])
    datasets.append(T_all)
    
    t2 = time()
    print('{}/{} finished in {}s'.format(f, len(files), t2-t1))

# concatenate all datasets along the time dimension
final_dataset = xr.concat(datasets, dim='time')


# save the final dataset to a netCDF file
final_dataset.to_netcdf('/home/claussar/ClimateVariablesOut/ERA5_9km_monthly/ERA5_2mtemp_1940-2022_mean_dailyminmax_{}_{}_{}.nc'.format(dataset, start_index, end_index))