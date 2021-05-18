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
         'r03id': df['R03_id'],
         'r04id': df['R04_id']}

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
    df = df.drop(columns=['R01_id', 'R02_id', 'R03_id', 'R04_id'])
    
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

def preprocess_ar_SASIA(df, thres):
    '''
    Returns dataframe encoded for AR Days and identifies trackID of AR.
    
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
    df = df.drop(columns=['R01_id', 'R02_id', 'R03_id', 'R04_id'])
    
    # Add column of AR days based on threshold
    # (no AR day eq 0; AR day eq 1)
    df['ar'] = 0
    idx = (df['R04'] > thres)
    df.loc[idx, 'ar'] = 1
    
    return df

def preprocess_ar_bbox(ds, start_date, end_date, bbox):
        lat1, lat2, lon1, lon2 = bbox
        # select lats, lons, and dates within start_date, end_date and months
        ds = ds.sel(time=slice(start_date, end_date), lat=slice(lat1,lat2), lon=slice(lon1,lon2))
        # convert dataset to dataframe
        df = ds.kidmap.to_dataframe(dim_order=['time', 'lat', 'lon'])
        df = df.dropna(axis='rows')
        # keep only rows that have trackID
        new_df = df.groupby('time').kidmap.unique()
        trackID = new_df

        # create list of days that have an AR
        date_lst = new_df.index.normalize().unique()
        date_lst = date_lst + pd.DateOffset(hours=9)
        
        # Create AR dataframe
        arr_ARDays = np.ones(len(date_lst), dtype=np.int)
        tmp = pd.DataFrame({'ar':arr_ARDays, 'time':date_lst})

        # Create dataframe with all days
        dates_allDays = pd.date_range(start=start_date, end=end_date, freq='1D')
        arr_allDays = np.empty(len(dates_allDays), dtype=np.int)*np.nan
        ar_dates = pd.DataFrame({'time':dates_allDays}, index=dates_allDays)

        # merge dfs together
        df = pd.merge(ar_dates, tmp, how='outer', on='time')
        df['ar'] = df['ar'].fillna(0)
        # set time as index
        df = df.set_index(df.time)
        ar_dates = df.drop(['time'], axis=1)
        
        return trackID, ar_dates
        
def get_ar_days(reanalysis, start_date, end_date, subregions=False, bbox=None, thresh=None):
    path_to_data = '/home/sbarc/students/nash/data/ar_catalog/'
    if subregions == True:
        # Open AR dataset (subregions with area covered by AR)
        if reanalysis == 'era5':
            filename = 'ar_catalog_v3_ERAI_fraction_HASIAsubregions.nc'
        ## if MERRA2
        else:
            filename = 'ar_catalog_fraction_HASIAsubregions.nc'

        f1 = path_to_data + 'CH1_generated_data/' + filename
        ds = xr.open_dataset(f1)
        # Set dates
        ds = ds.sel(time=slice(start_date, end_date))
        ## Preprocess AR subregions - get dataframe of AR days based on area threshold
        ar_dates = preprocess_ar_area_subregions(df=ds.to_dataframe(), thres=thresh)
        trackID = ar_dates.track_id
    
    elif subregions == False:
    
        if reanalysis == 'era5':
            filename =  'globalARcatalog_ERA-Interim_1979-2019_v3.0.nc'
        elif reanalysis == 'merra2':
            filename = 'globalARcatalog_MERRA2_1980-2019_v3.0.nc'

        # open ds
        ds = xr.open_dataset(path_to_data + filename, engine='netcdf4')
        ds = ds.squeeze()
        # remove lev and ens coords
        ds = ds.reset_coords(names=['lev', 'ens'], drop=True)
        trackID, ar_dates = preprocess_ar_bbox(ds, start_date, end_date, bbox)

    return trackID, ar_dates

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
    mon_ct : 1D array, float
        array of count of ARs per month in the time series
    clim_ct : 1D array, float
        array of size 12 with the total number of ARs per month
    '''
    mask = dataarray.where(dataarray >= threshold).dropna(dim='time')
    mask = mask.resample(time='1D').mean()
    day_list = mask.dropna(dim='time').time
    mon_ct = mask.resample(time='1MS').count()
    clim_ct = mask.groupby('time.month').count('time')
                           
    return day_list, mon_ct, clim_ct

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

def ar_daily_df(ssn, nk, out_path):
    
    filepath = out_path + 'AR-types_ALLDAYS.csv'
    df = pd.read_csv(filepath)

    # set up datetime index
    df = df.rename(columns={'Unnamed: 0': 'date'})
    df = df.set_index(pd.to_datetime(df.date))
    
    ## Break up columns into different AR Types
    keys = []
    for k in range(nk):
        keys.append("AR_CAT{:1d}".format(k+1,))

    values = np.zeros((len(df.index)))
    dicts = dict(zip(keys, values))

    df_cat = pd.DataFrame(dicts, index=df.index)

    for k in range(nk):
        idx = (df['AR_CAT'] == k+1)
        col = "AR_CAT{:1d}".format(k+1,)
        df_cat.loc[idx, col] = 1
        
    # get total of all AR types
    df_cat['AR_ALL'] = df_cat['AR_CAT1'] + df_cat['AR_CAT2'] + df_cat['AR_CAT3']
    df_cat['AR_CAT'] = df['AR_CAT']
    
    return df_cat