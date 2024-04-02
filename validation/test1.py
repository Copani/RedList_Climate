import pandas as pd
import numpy as np
import geopandas as gpd

# test if refactorization of the calc_climate_vars.py script yields same results

# load data
new = pd.read_csv('RedList_Climate/validation/NEWgaa2_climate_2004.csv')
val = pd.read_csv('RedList_Climate/validation/VALgaa2_climate_2004.csv')

# update keys
new_keys = ['Unnamed: 0', 'sci_name', 'presence', 'origin', 'seasonal', 'compiler',
       'yrcompiled', 'citation', 'generalisd', 'legend', 'Shape_Length',
       'Shape_Area', 'MAT_nichefrac_2004', 'MAT_2004', 'MAT_change_2004',
       'MTWM_nichefrac_2004', 'MTWM_2004', 'MTWM_change_2004',
       'MTCM_nichefrac_2004', 'MTCM_2004', 'MTCM_change_2004',
       'AP_nichefrac_2004', 'AP_2004', 'AP_change_2004', 'PDQ_nichefrac_2004',
       'PDQ_2004', 'PDQ_change_2004']
val.columns = new_keys

# compare
if new.equals(val):
    print('Dataframes are equal, refactorization is correct')
else:
    print('Dataframes are not equal! Refactorization is incorrect!')
    
# Python
diff = new.compare(val)
print(diff.keys())