import numpy as np
from netCDF4 import Dataset
import pickle
import glob
import pandas as pd
import geopandas as gpd
import sys
from scipy.ndimage import uniform_filter1d

#load daily climate data
def bring(file_name,var_name,yr0):
        x=Dataset(file_name)
        try:
                lon=x.variables['lon'][:]
        except:
                lon=x.variables['longitude'][:]
        try:
                lat=x.variables['lat'][:]
        except:
                lat=x.variables['latitude'][:]
        var=x.variables[var_name][:]

        #convert to -180>180 lon grid
        lon-=180
        var2=np.zeros(var.shape)
        var2[...,:int(var.shape[-1]/2)]=var[...,int(var.shape[-1]/2):]
        var2[...,int(var.shape[-1]/2):]=var[...,:int(var.shape[-1]/2)]

        #convert time to helpful datetime index
        dates=np.array(str(yr0) + '-01-01',dtype=np.datetime64)
        dts=np.linspace(0,len(x.variables['time'][:])-1,len(x.variables['time'][:]))
        dates=dates+dts.astype(int)
        dates=pd.DatetimeIndex(dates)

        return(lon,lat,var2,dates)

#main files
varn=sys.argv[1]
folder='/p/projects/compacts/data/climate/ERA5/' + varn + '/'
stem=varn + '_day_ECMWF-ERA5_observation_'
files=np.sort(glob.glob(folder + stem + '*'))

years=np.linspace(1940,2021,82).astype(int)
#which quart  of the globe to calculate for?
quart_i=int(sys.argv[2])
N=int(sys.argv[3])
abs0=273.15

#number of days exceeding these thresholds per year (wet days at lowest)
for f, file_name in enumerate(files[:30]):

	yr=int(file_name.split('observation_')[1][0:4])
	[lon,lat,T,date]=bring(file_name,varn,yr)
	if varn=='tas':
		T-=abs0
	if varn=='pr':
		T*=(24*1000)
	quart=int(T.shape[-1]/4)

	if f==0:
		Temp=T[...,quart_i*quart:(quart_i+1)*quart]
		dates=date
	else:
		Temp=np.concatenate((Temp,T[...,quart_i*quart:(quart_i+1)*quart]),axis=0)
		dates=np.concatenate((dates,date))
	del T

dates=pd.DatetimeIndex(dates)
threshs=[90,95,99,99.9]
#threshs=[0.1,1,5]

#estimate moving average prior to percentiles
if N>1:
	Temp=uniform_filter1d(Temp,size=N,mode='mirror',axis=0)

##below code to be used to select only hottest months to calculate percentiles
##choose 3 months around the hottest month to estimate percentiles in
#x=[]
#for m in range(12):
#	x.append(np.mean(Temp[dates.month==m+1,:,:],axis=0))
#x=np.array(x)
#hotm=np.argmax(x,axis=0)
#
##select data from 3 hottest months
##don't know how to do this better, broadcasting a condition based on a comparison of one axis with the values in another two
#prk=np.zeros((90*40,pr.shape[1],pr.shape[2]))
#for i in range(pr.shape[1]):
#	for j in range(pr.shape[2]):
#		prk[:,i,j]=Temp[np.isin(dates.month,[(hotm[i,j]-1+x)%12+1 for x in range(3)]),i,j][:90*40]
#del Temp

#estimate percentiles
threshs=np.percentile(Temp,threshs,axis=0)
print('calculated thresholds')

#save quartiled thresholds
np.save('threshs/' + varn + '_threshs_1940-1970_' + str(N) + '_quart_' + str(quart_i) + '.npy',threshs)


