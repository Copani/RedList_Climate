import numpy as np
from netCDF4 import Dataset

#constructing a human footprint index based on HYDE lannd-use and population data from 1901-2021 using the methods from Watson et al 2016:
#https://www.nature.com/articles/sdata201667#Sec1

#urban areas
urban=Dataset('/p/projects/isimip/isimip/ISIMIP3a/InputData/socioeconomic/landuse/histsoc/landuse-urbanareas_histsoc_annual_1901_2021.nc')

#cropland
totals=Dataset('/p/projects/isimip/isimip/ISIMIP3a/InputData/socioeconomic/landuse/histsoc/landuse-totals_histsoc_annual_1901_2021.nc')

#pastures
pastures=Dataset('/p/projects/isimip/isimip/ISIMIP3a/InputData/socioeconomic/landuse/histsoc/landuse-pastures_histsoc_annual_1901_2021.nc')

#human footprint: urban areas assigned 10, cropland 7 and pastures 3
LU_HF=urban['urbanareas'][:]*10 + totals['cropland_total'][:]*7 + totals['pastures'][:]*3
LU_HF=LU_HF.filled(fill_value=np.nan)

#population
pop=Dataset('/p/projects/isimip/isimip/ISIMIP3a/InputData/socioeconomic/pop/histsoc/population_histsoc_30arcmin_annual_1901_2021.nc')
lon=pop['lon'][:]
lat=pop['lat'][:]
pop=pop['total-population'][:]

#estimate population density
area=111*0.5*111*0.5*np.cos(lat*np.pi/180)
pop_dens=np.swapaxes(np.divide(np.swapaxes(pop,1,2),area),1,2)
pop_dens=pop_dens.filled(fill_value=np.nan)

#log scaling + capping at 1000 people per km2
pop_pressure=3.333*np.log(pop_dens+1)
pop_pressure[pop_dens>1000]=10

#total human footprint pressure
HF=LU_HF+pop_pressure
np.save('human_footprint/HYDE_landuse+pop_footprint.npy',HF)


