"""
Filename:    utils.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Helpful generic functions
"""

## Imports

import os, sys
import yaml
import xarray as xr

def check_mkdir(filename):
    '''Checks if directory exists and if not, makes that directory'''
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
                
def preprocess(ds, pathvar, bbox, rename_dict, lev=None):
    '''keep only selected lats and lons and rename variables'''
    lonmin, lonmax, latmin, latmax = bbox
    if pathvar == 'huvq':
        subset = ds.sel(latitude=slice(latmax, latmin), longitude=slice(lonmin, lonmax), level=lev)
        subset = subset.rename(rename_dict)
    if pathvar == 'prec':
        subset = ds.sel(latitude=slice(latmax, latmin), longitude=slice(lonmin, lonmax))
        subset = subset.rename(rename_dict)
    if pathvar == 'iwv':
        subset = ds.sel(latitude=slice(latmax, latmin), longitude=slice(lonmin, lonmax))
        subset = subset.rename(rename_dict)
        subset = subset.drop(['tcrw', 'tcsw', 'tcw'])
    if pathvar == 'ivt':
        subset = ds.sel(latitude=slice(latmax, latmin), longitude=slice(lonmin, lonmax))
        subset = subset.rename(rename_dict)
        
    return subset
    
def load_era5_data(pathvar, bbox, anom=True, lev=None):
    path_to_data = '/home/nash/DATA/data/'
    if anom == True:
        filepath_pattern = path_to_data + 'ERA5/{0}/anomalies/daily_filtered_anomalies_{0}_*.nc'.format(pathvar)
    elif anom == False: 
        filepath_pattern = path_to_data + 'ERA5/{0}/daily/out.era5_hma_05dg_daily_{0}_*.nc'.format(pathvar)
    
    yaml_doc = '/home/sbarc/students/nash/repositories/AR_types/data/load_era5.yaml'
    config = yaml.load(open(yaml_doc), Loader=yaml.SafeLoader)
    rename_dict = config[pathvar]
        
    f = xr.open_mfdataset(filepath_pattern)
    f = preprocess(f, pathvar, bbox, rename_dict, lev)
    
    return f
    