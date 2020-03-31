"""
Filename:    ar_funcs.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions used to calculate AR Climatology Statistics in HMA
"""

## Imports

import os, sys
import numpy as np
import xarray as xr


## FUNCTIONS

def ar_climatology(dataarray, threshold):
    '''
    Returns list array of dates considered AR days based on the input subregion.
    
    Parameters
    ----------
    datarray : xarray dataarray object
        subregion where percentage of area of AR is given for each day
    threshold: number, float
        singular threshold used for the percentage of area covered by an AR
        
    Returns
    -------
    day_list : 1D array, float
        list of datetime objects that an AR covered threshold*100% of the subregion's area
    '''
    mask = dataarray.where(dataarray >= threshold).dropna(dim='time')
    mask = mask.resample(time='1D').mean()
    day_list = mask.dropna(dim='time').time
                           
    return day_list