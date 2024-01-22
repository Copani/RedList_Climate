import numpy as np
from netCDF4 import Dataset
import pickle
import glob
import pandas as pd
import scipy.interpolate

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

#main files
varn='pr'
folder='/p/projects/compacts/data/climate/ERA5/' + varn + '/'
stem=varn + '_day_ECMWF-ERA5_observation_'
files=np.sort(glob.glob(folder + stem + '*'))

years=np.linspace(1940,2021,82).astype(int)

#number of days exceeding these thresholds per year (wet days at lowest)
for f, file_name in enumerate(files):

	if f<19:
		yr=int(file_name.split('observation_')[1][0:4])
		[lon,lat,pr,dates]=bring(file_name,varn,yr)
	else:
		yr=int(file_name.split('observation_')[1][0:4])
		[lon,lat,pr,dates]=bring(file_name,varn,yr)
	
	#convert precipitation from kgm-2s-2 to mm
	pr*=(24*1000)

	interm_m=calc_month_stat(pr,dates,'sum')

	#append these wet days measures to the previous years values			
	if f==0:
		P=interm_m
	else:
		P=np.concatenate((P,interm_m),axis=2)

	del interm_m
	del pr
	print('done year ' + str(yr))

livingplanet=pd.read_csv('/p/projects/compacts/data/biodiverse/LivingPlanetIndex/LPD2022_public/LPD2022_formatted.csv',encoding='latin1',dtype={'M_biome': 'str'})
livingplanet['ID']=pd.to_numeric(livingplanet['ID'],errors='coerce')
LPcoords=livingplanet[['ID','Latitude','Longitude']].groupby('ID').mean()
del livingplanet

ID_list=[]
year_list=[]
P_suml=[]

for i, ID in enumerate(LPcoords.index):

	#population coordinates and indices in the grid
	latlon=np.array(LPcoords.loc[LPcoords.index==ID,['Latitude','Longitude']])[0]
	lati=np.where(abs(lat-latlon[0])==np.min(abs(lat-latlon[0])))	
	loni=np.where(abs(lon-latlon[1])==np.min(abs(lon-latlon[1])))

	#select climate data from grid
	#averaging over first two axes in case more than one lati/loni were selected
	P_suml+=list(np.sum(np.mean(P[lati,loni,:,:],axis=(0,1)),axis=-1).flatten())

	ID_list+=[ID]*len(years)
	year_list+=list(years)
	
	if i%100==0:
		print('done ' + str(i) + ' populations')

data=pd.DataFrame()
data['ID']=ID_list
data['year']=year_list
data['Psum']=P_suml

data.to_csv('climate_merged/livingplanet_ERA5_Psum.csv') 

