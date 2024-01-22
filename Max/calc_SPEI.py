import numpy as np
from netCDF4 import Dataset
import pickle
import glob
import pandas as pd
import geopandas as gpd
import scipy.interpolate
import sys

def bring(file_name,var_name,yr0):
	x=Dataset(file_name)
	lon=x.variables['lon'][:]
	lat=x.variables['lat'][:]
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

ts=sys.argv[1]
data=Dataset('/home/maxkotz/ECB_inflation/climate_gridded/spei' + str(ts).zfill(2) + '.nc')

#need to flip spei axis to make it the same as the ISIMIP grid
lon=data['lon'][:]
lat=data['lat'][:]
spei=data['spei'][:]
spei=np.reshape(spei,(int(spei.shape[0]/12),12,spei.shape[1],spei.shape[2]),order='C')

#flip latitude axis to make comparable to ISIMIP
lat=np.flip(lat)
spei=np.flip(spei,axis=2)

#fill nans
spei=spei.filled(np.nan)

years=np.linspace(1901,2020,120).astype(int)

#get living planet data to merge 
livingplanet=pd.read_csv('/p/projects/compacts/data/biodiverse/LivingPlanetIndex/LPD2022_public/LPD2022_formatted.csv',encoding='latin1',dtype={'M_biome': 'str'})
livingplanet['ID']=pd.to_numeric(livingplanet['ID'],errors='coerce')
LPcoords=livingplanet[['ID','Latitude','Longitude']].groupby('ID').mean()
del livingplanet

#average spei
spei_l=[]
#over thresholds
threshs=[1.5,2,2.5,3]
spei_gt_l=[]
spei_ltm_l=[]
for t, thresh in enumerate(threshs):
	spei_gt_l.append([])
	spei_ltm_l.append([])
ID_list=[]
year_list=[]

for i, ID in enumerate(LPcoords.index):

	#population coordinates and indices in the grid
	latlon=np.array(LPcoords.loc[LPcoords.index==ID,['Latitude','Longitude']])[0]	
	lati=np.where(abs(lat-latlon[0])==np.min(abs(lat-latlon[0])))
	loni=np.where(abs(lon-latlon[1])==np.min(abs(lon-latlon[1])))

	#average over lati/loni in case of times with multiple lati/loni
	speis=np.mean(spei[:,:,lati,loni],axis=(2,3))

	#select climate data from grid
	#also average over the first two axes in case more than one lati/loni are selected
	#average spei data
	spei_l+=list(np.mean(speis[:,:],axis=(1)))
	#spei months greater/less than minnus a certain extreme threshold
	for t, thresh in enumerate(threshs):
		spei_gt_l[t]+=list(np.sum(speis[:,:]>thresh,axis=(1)))
		spei_ltm_l[t]+=list(np.sum(speis[:,:]<-thresh,axis=(1)))

	ID_list+=[ID]*len(years)
	year_list+=list(years)

	if i%100==0:
		print('done ' + str(i) + ' populations')

data=pd.DataFrame()
data['ID']=ID_list
data['year']=year_list
data['spei' + str(ts)]=spei_l
for t, thresh in enumerate(threshs):
	data['spei' + str(ts) + '_gt_' + str(thresh)]=spei_gt_l[t]
	data['spei' + str(ts) + '_ltm_' + str(thresh)]=spei_ltm_l[t]

data.to_csv('climate_merged/livingplanet_CRU_SPEI' + str(ts).zfill(2) + '.csv')





