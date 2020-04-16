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

def preprocess_ar_area_subregions(df, thres):
    '''
    Returns dataframe encoded for AR Days and identifies subregion where AR is present.
    
    Parameters
    ----------
    df : pandas dataframe
        dataframe that has subregion where percentage of area of AR is given for each day
    threshold: number, float
        singular threshold used for the percentage of area covered by an AR
        
    Returns
    -------
    df : pandas dataframe
        df that indicates whether each time step is an AR day and the location of the AR
    '''
    ## drop lev and ens cols
    df = df.drop(columns=['lev', 'ens'])
    ## resample to daily
    df = df.resample('1D').mean()
    # Add column of AR days based on threshold
    # (no LLJ day eq 0; LLJ day eq 1)
    df['ar'] = 0
    idx = (df['R01'] > thres) | (df['R02'] > thres) | (df['R03'] > thres)
    df.loc[idx, 'ar'] = 1

    # Add column of AR locations 
    # ('R01', 'R02', 'R03', 'R01/R02', 'R02/R03', 'R01/R03', 'R01/R02/R03', nan)
    df['location'] = np.nan

    idx = (df['R01'] >= thres) & (df['R02'] < thres) & (df['R03'] < thres)
    df.loc[idx, 'location'] = 'R01'

    idx = (df['R01'] < thres) & (df['R02'] >= thres) & (df['R03'] < thres)
    df.loc[idx, 'location'] = 'R02'

    idx = (df['R01'] < thres) & (df['R02'] < thres) & (df['R03'] >= thres)
    df.loc[idx, 'location'] = 'R03'

    idx = (df['R01'] >= thres) & (df['R02'] >= thres) & (df['R03'] < thres)
    df.loc[idx, 'location'] = 'R01/R02'

    idx = (df['R01'] < thres) & (df['R02'] >= thres) & (df['R03'] >= thres)
    df.loc[idx, 'location'] = 'R02/R03'

    idx = (df['R01'] >= thres) & (df['R02'] < thres) & (df['R03'] >= thres)
    df.loc[idx, 'location'] = 'R01/R03'

    idx = (df['R01'] >= thres) & (df['R02'] >= thres) & (df['R03'] >= thres)
    df.loc[idx, 'location'] = 'R01/R02/R03'
    
    return df


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