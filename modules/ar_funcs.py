"""
Filename:    ar_funcs.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions used to calculate AR Climatology Statistics in HMA
"""

## Imports

import os, sys
import numpy as np
import xarray as xr
import pandas as pd


## FUNCTIONS

def resample_track_id(df):
    '''
    Returns an array that has a single AR track ID for each 24 hr
    
    '''
    # stack up 3 subregions id numbers
    d = {'r01id': df['R01_id'],
         'r02id': df['R02_id'],
         'r03id': df['R03_id']}

    df_tmp = pd.DataFrame(data=d)
    ## combine into single series
    df_tmp = df_tmp.stack()

    ## resample to 1D taking maximum of new column
    level_values = df_tmp.index.get_level_values
    df_tmp = df_tmp.groupby([level_values(i) for i in [0,1]]+[pd.Grouper(freq='1D', level=0)]).max()


    # date array with all days
    start_date = str(df.index.year[0])+'-'+ str(df.index.month[0]) + '-' + str(df.index.day[0])
    end_date = str(df.index.year[-1])+'-'+ str(df.index.month[-1]) + '-' + str(df.index.day[-1])
    dates_allDays = pd.date_range(start=start_date, end=end_date, freq='1D')
    arr_allDays = np.zeros(len(dates_allDays), dtype=np.float)
    arr_allDays[:] = np.nan

    # Loop over AR days ID and match to list of ALL days 
    for i, date in enumerate(df_tmp.index.get_level_values(2)):
        idx = np.where(dates_allDays == date)
        arr_allDays[idx] = df_tmp.values[i]

    return arr_allDays

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
    ## Get single AR ID for each day
    track_ids = resample_track_id(df)
    
    # resample to daily
    df = df.resample('1D').mean()
    ## manually add column back into resampled df with area covered
    df['track_id'] = track_ids
    df = df.drop(columns=['R01_id', 'R02_id', 'R03_id'])
    
    # Add column of AR days based on threshold
    # (no AR day eq 0; AR day eq 1)
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

def add_ar_time_series(ds, df):
    '''Add AR time series to ds; set as coordinate variables'''
    ds['ar'] = ('time', df.ar)
    ds = ds.set_coords('ar')
    ds['location'] = ('time', df.location)
    ds = ds.set_coords('location')
    
    return ds

def calc_seasonal_contribution(ds_list, df, prec_var, mon_s, mon_e):
    '''
    For a list of ds, calculate the average total seasonal contribution of ARs for the given prec_vars
    
    '''
    ds_clim_lst = []
    ds_frac_lst = []
    ds_std_lst = []

    for k, ds in enumerate(ds_list):
        # Add AR time series to ds; set as coordinate variables
        ds = add_ar_time_series(ds, df)

        # Select months
        if mon_s > mon_e:
            idx = (df.index.month >= mon_s) | (df.index.month <= mon_e)
        else:
            idx = (df.index.month >= mon_s) & (df.index.month <= mon_e)
        ds = ds.sel(time=idx)
        
        # Select AR days
        idx = (ds.ar >= 1)
        ds_ar = ds.sel(time=idx)

        # calculate seasonal totals
        ds_ssn_sum = ds.resample(time='QS-DEC').sum()
        ds_ar_ssn_sum = ds_ar.resample(time='QS-DEC').sum()                                 

        # calculate average of seasonal totals
        ds_clim = ds_ssn_sum.mean(dim='time')
        ds_ar_clim = ds_ar_ssn_sum.mean(dim='time') 

        # things to output/append to final lists
        ds_frac_lst.append((ds_ar_clim[prec_var[k]].values/ds_clim[prec_var[k]].values)*100.)
        ds_std_lst.append(ds_ar_ssn_sum[prec_var[k]].std(dim='time').values)
        ds_clim_lst.append(ds_clim[prec_var[k]].values)
    
    return ds_clim_lst, ds_frac_lst, ds_std_lst