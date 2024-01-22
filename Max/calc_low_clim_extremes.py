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

#calc monthly statistic of daily climate variable
def calc_month_stat(var,dates,stat):
	years=pd.unique(dates.year)
	var2=np.zeros((var.shape[1],var.shape[2],len(years),12))
        
	for i in range(len(years)):                
		for j in range(12):
			if stat=='std':
				var2[:,:,i,j]=np.std(np.squeeze(var[np.where((dates.year==years[i]) & (dates.month==j+1)),:,:]),axis=0)
			elif stat=='VAR':
				var2[:,:,i,j]=np.square(np.std(np.squeeze(var[np.where((dates.year==years[i]) & (dates.month==j+1)),:,:]),axis=0))
			elif stat=='mean':
				var2[:,:,i,j]=np.mean(np.squeeze(var[np.where((dates.year==years[i]) & (dates.month==j+1)),:,:]),axis=0)
			elif stat=='sum':
				var2[:,:,i,j]=np.sum(np.squeeze(var[np.where((dates.year==years[i]) & (dates.month==j+1)),:,:]),axis=0)
			else:
				sys.exit('Specified incorrect statistic to calculate: use "std" or "mean" or "var"')        
	return(var2)

def calc_annual_extreme_stat(var,dates,threshs,stat):
	years=pd.unique(dates.year)
	var2=np.zeros((threshs.shape[0],len(years),var.shape[1],var.shape[2]))

	for i in range(len(years)):
		datum=np.where(dates.year==years[i])
		if stat=='num':
			for t in range(len(threshs)):
				var2[t,i,:,:]=np.sum(np.squeeze(var[datum,:,:])<threshs[t],axis=0)
		elif stat=='am':
			for t in range(len(threshs)):
				diff=-(var[datum,:,:]-threshs[t])
				var2[t,i,:,:]=np.sum(np.squeeze(diff),where=np.squeeze(diff>0),axis=0)
		else:
			print('INCORRECT STATISTIC REQUESTED')
	return(var2)

#main files
varn=sys.argv[1]
folder='/p/projects/compacts/data/climate/ERA5/' + varn + '/'
stem=varn + '_day_ECMWF-ERA5_observation_'
files=np.sort(glob.glob(folder + stem + '*'))
if varn=='tas':
	files=files[:-1]

years=np.linspace(1940,2021,82).astype(int)

threshs=[0.1,1,5]

#load percentiles
for q in range(4):
	if q==0:
		percentiles=np.load('threshs/' + varn + '_low_threshs_1940-1970_1_quart_' + str(q) + '.npy')
	else:
		percentiles=np.concatenate((percentiles,np.load('threshs/' + varn + '_low_threshs_1940-1970_1_quart_' + str(q) + '.npy')),axis=-1)

#number of days exceeding these thresholds per year (wet days at lowest)
for f, file_name in enumerate(files):

	yr=int(file_name.split('observation_')[1][0:4])
	[lon,lat,T,dates]=bring(file_name,varn,yr)
	
	#convert precipitation from kgm-2s-2 to mm, or tas from Kelvin to Centigrade
	if varn=='tas':
		T-=273.15	
	if varn=='pr':
		T*=(24*3600/1000)*1000*24
	
	interm_num=calc_annual_extreme_stat(T,dates,percentiles,'num')
	interm_am=calc_annual_extreme_stat(T,dates,percentiles,'am')

	#append these wet days measures to the previous years values			
	if f==0:
		extrm_num=interm_num
		extrm_am=interm_am
	else:
		extrm_num=np.concatenate((extrm_num,interm_num),axis=1)
		extrm_am=np.concatenate((extrm_am,interm_am),axis=1)

	del interm_num
	del interm_am
	del T
	print('done year ' + str(yr))

livingplanet=pd.read_csv('/p/projects/impactee/Data/biodiverse/LivingPlanetIndex/LPD2022_public/LPD2022_formatted.csv',encoding='latin1',dtype={'M_biome': 'str'})
IDs=livingplanet.ID.unique()
LPcoords=livingplanet[['ID','Latitude','Longitude']].groupby('ID').mean()
del livingplanet

ID_list=[]
year_list=[]
extrm_num_l=[]
extrm_am_l=[]

for t in range(len(threshs)):
	extrm_num_l.append([])
	extrm_am_l.append([])

for i, ID in enumerate(LPcoords.index):

	#population coordinates and indices in the grid
	latlon=np.array(LPcoords.loc[LPcoords.index==ID,['Latitude','Longitude']])[0]
	#select climate data from grid
	if not np.any(np.isnan(latlon)):
		lati=np.where(abs(lat-latlon[0])==np.min(abs(lat-latlon[0])))[0][0]
		loni=np.where(abs(lon-latlon[1])==np.min(abs(lon-latlon[1])))[0][0]
		for t in range(len(threshs)):
			extrm_num_l[t]+=list(extrm_num[t,:,lati,loni].flatten())
			extrm_am_l[t]+=list(extrm_am[t,:,lati,loni].flatten())
	else:
		for t in range(len(threshs)):
			extrm_num_l[t]+=[np.nan]*len(years)
			extrm_am_l[t]+=[np.nan]*len(years)

	ID_list+=[ID]*len(years)
	year_list+=list(years)

	if i%100==0:
		print('done ' + str(i) + ' populations')

data=pd.DataFrame()
data['ID']=ID_list
data['year']=year_list
for t, thresh in enumerate(threshs):
	data[varn + '_extrm_' + str(thresh) + '_num']=extrm_num_l[t]
	data[varn + '_extrm_' + str(thresh) + '_am']=extrm_am_l[t]

data.to_csv('climate_merged/livingplanet_ERA5_' + varn + '_low_extrm.csv')
 

