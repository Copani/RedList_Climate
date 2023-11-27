import pandas as pd
import numpy as np
import geopandas as gpd
from netCDF4 import Dataset
import glob
import sys
import math

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

Class=sys.argv[1]
#which quarter of data to process
q_i=int(sys.argv[2])

if Class=='BIRDS':
	range_maps=gpd.read_file('/p/projects/impactee/Data/biodiverse/IUCN/BIRDS/a0000000a.gdbtable')
else:
	rfiles=glob.glob('/p/projects/impactee/Data/biodiverse/IUCN/' + Class + '/' + Class + '*.shp')
	for r, rfilen in enumerate(rfiles):
		if r==0:
			range_maps=gpd.read_file(rfilen)
		else:
			range_maps=gpd.GeoDataFrame(pd.concat([range_maps,gpd.read_file(rfilen)],ignore_index=True))

#only select species which we actually have in Living planet
Class_match={'MAMMALS':'Mammalia','REPTILES':'Reptilia','BIRDS':'Aves','AMPHIBIANS':'Amphibia'}
LP=pd.read_csv('climate_bio_merged/LP+climate.csv')
LPspecies=LP.loc[LP.Class==Class_match[Class]].Binomial.unique()
LPspecies=[x.split('_')[0] + ' ' + x.split('_')[1] for x in LPspecies]
del LP
#the N to keep
N=math.ceil(len(LPspecies)/20)
LPspecies=LPspecies[q_i*N:(q_i+1)*N]

#climate data
varn='tas'
folder='/p/projects/compacts/data/climate/ERA5/' + varn + '/'
stem=varn + '_day_ECMWF-ERA5_observation_'
files=np.sort(glob.glob(folder + stem + '*'))

#get climate data grid
file_n=files[0]
yr=int(file_n.split('observation_')[1][0:4])
[lon,lat,T,dates]=bring(file_n,varn,yr)
del T

def get_grid_indices(poly,lon,lat):
	#get poly bounnds
	#bounds=poly.bounds
	bounds=np.array(poly.bounds)
	#get subgrid within those bounds
	subloni=np.where(np.logical_and(lon>bounds[0],lon<bounds[2]))
	sublati=np.where(np.logical_and(lat>bounds[1],lat<bounds[3]))
	sublon=lon[subloni]
	sublat=lat[sublati]

	#find locations within subgrid which fall within the polygon
	grid=np.meshgrid(sublon,sublat)
	gridi=np.meshgrid(subloni,sublati)
	lonlat_points=gpd.GeoSeries(gpd.points_from_xy(grid[0].flatten(),grid[1].flatten()))
	within=lonlat_points.within(poly)

	#indices of lon/lat coordinates within the polygon
	loni=gridi[0].flatten()[within]
	lati=gridi[1].flatten()[within]

	#if nothing comes back, take centroid
	if not np.any(within):
		loni=np.where(abs(lon-poly.centroid.x)==min(abs(lon-poly.centroid.x)))[0][0]
		lati=np.where(abs(lat-poly.centroid.y)==min(abs(lat-poly.centroid.y)))[0][0]		

		#failed attempt to first improve resolution of search
#		within=within.astype(int)
#		N=5
#		res=0.25/N
#		grid[0]-=int(N/2)*res
#		grid[1]-=int(N/2)*res
#		lonlat_points=gpd.GeoSeries(gpd.points_from_xy(grid[0].flatten(),grid[1].flatten()))
#		for i in range(N):
#			gridf1=grid[0]-(int(N/2)+i)*res
#			for j in range(N):
#				gridf2=grid[1]-(int(N/2)+i)*res
#				lonlat_points=gpd.GeoSeries(gpd.points_from_xy(gridf1.flatten(),gridf2.flatten()))			
#				within+=lonlat_points.within(poly).astype(int)
#		within=within.astype(bool)

	return([loni,lati])

coords=[]
#get longitude/latitude indices of rannge of all species
for s, species in enumerate(LPspecies):	

	rmap=range_maps.loc[range_maps.sci_name==species]
	#remove invalid geometries
	rmap=rmap.explode()
	rmap=rmap.loc[rmap.geometry.is_valid]
	if len(rmap)==0:
		#if no range map available, nan
		coords.append([species,np.nan,np.nan])
	else:
		#take union of all available geometries and then find the indices on the grid that match
		poly=rmap.geometry.unary_union
		[loni,lati]=get_grid_indices(poly,lon,lat)
		coords.append([species,loni,lati])	
	print(s)
	if s%10==0:
		print('done ' + str(s) + ' / ' + str(len(LPspecies)))

import pickle
with open('IUCN_range_coords/IUCN_ERA5_coords_' + Class + '_' + str(q_i) + '.pkl', 'wb') as f:
	pickle.dump(coords,f)


