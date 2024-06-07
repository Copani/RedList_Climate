import cdsapi
import numpy as np

c = cdsapi.Client()
for year in np.arange(1950, 2023):

    c.retrieve(
        'reanalysis-era5-land-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': 'total_precipitation',
            'year': str(year),
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'time': '00:00',
            'format': 'netcdf',
        },
        'C:/Users/Claus/MasterThesis/Data/ClimateVariables/ERA5_monthly/ERA5_prec_{}.nc'.format(year))
