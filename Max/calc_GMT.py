import numpy as np
from netCDF4 import Dataset
import pickle
import glob
import pandas as pd
import scipy.interpolate
import sys

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

def calc_annual_extreme_stat(var,dates,threshs,stat):
	years=pd.unique(dates.year)
	var2=np.zeros((threshs.shape[0],len(years),var.shape[1],var.shape[2]))

	for i in range(len(years)):
		datum=np.where(dates.year==years[i])
		if stat=='num':
			for t in range(len(threshs)):
				var2[t,i,:,:]=np.sum(np.squeeze(var[datum,:,:])>threshs[t],axis=0)
		elif stat=='am':
			for t in range(len(threshs)):
				diff=var[datum,:,:]-threshs[t]
				var2[t,i,:,:]=np.sum(np.squeeze(diff),where=np.squeeze(diff>0),axis=0)
		else:
			print('INCORRECT STATISTIC REQUESTED')
	return(var2)

#main files
varn='tas'
folder='/p/projects/compacts/data/climate/ERA5/' + varn + '/'
stem=varn + '_day_ECMWF-ERA5_observation_'
files=np.sort(glob.glob(folder + stem + '*'))
if varn=='tas':
	files=files[:-1]

years=np.linspace(1940,2021,82).astype(int)
GMTs=[]

#number of days exceeding these thresholds per year 
for f, file_name in enumerate(files):

	yr=int(file_name.split('observation_')[1][0:4])
	[lon,lat,T,dates]=bring(file_name,varn,yr)
	
	T-=273.15
	areaweights=np.cos(np.pi*lat/180)
	
	GMT=np.mean(np.tensordot(np.mean(T,0),areaweights,axes=(0,0))/np.sum(areaweights))
	GMTs.append(GMT)
	del T
	print('done year ' + str(yr))

df=pd.DataFrame()
df['Year']=years
df['GMT']=GMTs
df.to_csv('climate_merged/ERA5_GMT.csv')

