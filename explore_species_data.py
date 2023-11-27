import pandas as pd
import numpy as np

#Livingplanet
# livingplanet=pd.read_csv('/p/projects/impactee/Data/biodiverse/LivingPlanetIndex/LPD2022_public/LPD2022_public.csv',encoding='latin1',dtype={'M_biome': 'str'})
livingplanet=pd.read_csv('Data/LPD2022_public.csv',encoding='latin1',dtype={'M_biome': 'str'})

colnames=list(livingplanet)
years=colnames[-71:]
colnames=colnames[:len(colnames)-71]

x=pd.melt(livingplanet,id_vars=colnames,value_vars=years,var_name='year')
x=x.sort_values(['ID','year'])
print('Total observations: ' + str(len(x.loc[~x.value.isnull()])))

x.to_csv('/p/projects/impactee/Data/biodiverse/LivingPlanetIndex/LPD2022_public/LPD2022_formatted.csv')

#number of observations per population (ID)
no_obs=x.groupby('ID')['value'].apply(lambda y: np.sum(~y.isnull()))
print('Median no. observations per population: ' + str(np.median(no_obs)))
#number of populations with at least X observations
thresh=10
print('No. populations with >' + str(thresh) + ' observations: ' + str(np.sum(no_obs>thresh)))
print('Median no. observations per population for populations with >' + str(thresh) + ' observations: ' + str(np.median(no_obs[no_obs>thresh])))

#if expressing population in terms of annual growth rate, how many observations would we have?
x['pop_growth']=x['value'].diff(periods=-1)
#set last entry to zero for each group
x.loc[x.year==2020]=np.nan

no_obs=x.groupby('ID')['pop_growth'].apply(lambda y: np.sum(~y.isnull()))
print('Median no. observations per population: ' + str(np.median(no_obs)))
#number of populations with at least X observations
thresh=10
print('No. populations with >' + str(thresh) + ' observations: ' + str(np.sum(no_obs>thresh)))
print('Median no. observations per population for populations with >' + str(thresh) + ' observations: ' + str(np.median(no_obs[no_obs>thresh])))

#visual inspection
x.loc[~x.pop_growth.isnull()][['ID','Class','Common_name','year','value','pop_growth']].iloc[:60]

#By class
Classes=x.Class.unique()
x.loc[x.Class.isin(Classes[4:]),'Class']='Fish'
Classes=x.Class.unique()

for c, Class in enumerate(Classes):
	print(Class + ': ' + str(len(x.loc[(x.Class==Class)&(~x.pop_growth.isnull())])) + ' observations')

#with at least X observations
thresh=5
for c, Class in enumerate(Classes):
	print(Class)
	no_obs=x.loc[x.Class==Class].groupby('ID')['pop_growth'].apply(lambda y: np.sum(~y.isnull()))
	print('Median no. observations per population: ' + str(np.median(no_obs)))
	print('No. populations with >' + str(thresh) + ' observations: ' + str(np.sum(no_obs>thresh)))
	print('Median no. observations per population for populations with >' + str(thresh) + ' observations: ' + str(np.median(no_obs[no_obs>thresh])))
	print('####')

#By system
systems=x.System.unique()
for s, system in enumerate(systems):
	print(system + ': ' + str(len(x.loc[(x.System==system)&(~x.pop_growth.isnull())])) + ' observations')

#with at least X observations
thresh=5
for s, system in enumerate(systems):
	print(system)
	no_obs=x.loc[x.System==system].groupby('ID')['pop_growth'].apply(lambda y: np.sum(~y.isnull()))
	print('Median no. observations per population: ' + str(np.median(no_obs)))
	print('No. populations with >' + str(thresh) + ' observations: ' + str(np.sum(no_obs>thresh)))
	print('Median no. observations per population for populations with >' + str(thresh) + ' observations: ' + str(np.median(no_obs[no_obs>thresh])))
	print('####')

#BIOTIME
biotime=pd.read_csv('/p/projects/impactee/Data/biodiverse/BioTime/BioTIMEQuery_24_06_2021.csv')
print('Total observatoins: ' + str(len(biotime.loc[~biotime['sum.allrawdata.ABUNDANCE'].isnull()])))

#PREDICTS
predicts=pd.read_csv('/p/projects/impactee/Data/biodiverse/PREDICTS/predicts_data/resource.csv')
len(predicts.loc[~predicts.Measurement.isnull()])

Classes=predicts.Class.unique()

#match Classes from predicts
matcher=predicts[['Class','Genus']].drop_duplicates()
biotime=pd.merge(biotime,matcher,how='left',left_on='GENUS',right_on='Genus')
print('Observations with Class: ' + str(len(biotime.loc[~biotime.Class.isnull()])))

Classes=biotime.Class.unique()
Classes=Classes[3:]

#data available by class
for c, Class in enumerate(Classes):
        print(Class + ': ' + str(len(biotime.loc[(biotime.Class==Class)&(~biotime['sum.allrawdata.ABUNDANCE'].isnull())])) + ' observations')

#generate unique identifier
biotime['GENUS_SPECIES_LAT_LON']=biotime['GENUS_SPECIES']+'_'+biotime.LATITUDE.round(decimals=1).astype(str)+'_'+biotime.LONGITUDE.round(decimals=1).astype(str)

#with at least X observations
thresh=5
for c, Class in enumerate(Classes):
        print(Class)
        no_obs=biotime.loc[biotime.Class==Class].groupby('GENUS_SPECIES_LAT_LON')['sum.allrawdata.ABUNDANCE'].apply(lambda y: np.sum(~y.isnull()))
        print('Median no. observations per population: ' + str(np.median(no_obs)))
        print('No. populations with >' + str(thresh) + ' observations: ' + str(np.sum(no_obs>thresh)))
        print('Median no. observations per population for populations with >' + str(thresh) + ' observations: ' + str(np.median(no_obs[no_obs>thresh])))
        print('####')

biotime.to_csv('/p/projects/impactee/Data/biodiverse/BioTime/BioTIMEQuery_24_06_2021.csv')

####PREDICTSDATA
pred=pd.read_csv('/p/projects/impactee/Data/biodiverse/PREDICTS/predicts_data/resource.csv')

#generate unique identifier for species/latlon
pred['GENUS_SPECIES_LAT_LON']=pred['Genus']+pred['Species']+'_'+pred.Latitude.round(decimals=1).astype(str)+'_'+pred.Longitude.round(decimals=1).astype(str)

pred.loc[pred.Class=='Mammalia'].sort_values('GENUS_SPECIES_LAT_LON')[['Sample_midpoint','Species','Measurement']]





