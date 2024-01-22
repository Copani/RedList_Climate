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
varn='tas'
folder='/p/projects/compacts/data/climate/ERA5/' + varn + '/'
stem=varn + '_day_ECMWF-ERA5_observation_'
files=np.sort(glob.glob(folder + stem + '*'))

years=np.linspace(1940,2021,82).astype(int)

#number of days exceeding these thresholds per year (wet days at lowest)
for f, file_name in enumerate(files[:-1]):

	if f<19:
		yr=int(file_name.split('observation_')[1][0:4])
		[lon,lat,T,dates]=bring(file_name,varn,yr)
	else:
		yr=int(file_name.split('observation_')[1][0:4])
		[lon,lat,T,dates]=bring(file_name,varn,yr)
	
	#convert precipitation from kgm-2s-2 to mm
	T-=273.15	
	interm_m=calc_month_stat(T,dates,'mean')

	#append these wet days measures to the previous years values			
	if f==0:
		Tmean=interm_m
	else:
		Tmean=np.concatenate((Tmean,interm_m),axis=2)

	del interm_m
	del T
	print('done year ' + str(yr))

livingplanet=pd.read_csv('/p/projects/impactee/Data/biodiverse/LivingPlanetIndex/LPD2022_public/LPD2022_formatted.csv',encoding='latin1',dtype={'M_biome': 'str'})
IDs=livingplanet.ID.unique()
LPcoords=livingplanet[['ID','Latitude','Longitude']]

threshs=[1,1.5,2,3]
ID_list=[]
year_list=[]
Tanoms_sum=[]
Tanoms_abssum=[]
Tanoms_nohot=[]
Tanoms_nocold=[]
for t, thresh in enumerate(threshs):
	Tanoms_nohot.append([])
	Tanoms_nocold.append([])

for i, ID in enumerate(IDs):

	#population coordinates and indices in the grid
	latlon=LPcoords.iloc[i][['Latitude','Longitude']]
	lati=np.where(abs(lat-latlon[0])==np.min(abs(lat-latlon[0])))	
	loni=np.where(abs(lon-latlon[1])==np.min(abs(lon-latlon[1])))

	#select climate data from grid
	Tlocal=np.squeeze(np.mean(Tmean[lati,loni,:,:],axis=(0,1)))
	#estimate monthly temp anomalies from historical seasonal cycle
	T_anom=Tlocal-np.mean(Tlocal[:30,:],axis=0)
	#normalise by their historical monthly variability 
	T_anom=np.divide(T_anom,np.std(Tlocal[:30,:],axis=0))
	#aggregate in a few different ways
	#sum those anomalies, and their absolute values
	Tanoms_sum+=list(np.sum(T_anom,axis=-1))
	Tanoms_abssum+=list(np.sum(abs(T_anom),axis=-1))
	for t, thresh in enumerate(threshs):
		#count number of hot anomalies above thresholds
		Tanoms_nohot[t]+=list(np.sum(T_anom>thresh,axis=-1))
		#count number of cold anomalies above thresholds
		Tanoms_nocold[t]+=list(np.sum(T_anom<-thresh,axis=-1))

	ID_list+=[ID]*len(years)
	year_list+=list(years)
	
	if i%100==0:
		print('done ' + str(i) + ' populations')

data=pd.DataFrame()
data['ID']=ID_list
data['year']=year_list
for t, thresh in enumerate(threshs):
	data['T_anoms_nohot_' + str(thresh)]=Tanoms_nohot[t]
	data['T_anoms_nocold_' + str(thresh)]=Tanoms_nocold[t]
data['Tanoms_sum']=Tanoms_sum
data['Tanoms_abssum']=Tanoms_abssum

data.to_csv('climate_merged/livingplanet_ERA5_T_anoms.csv')
 

