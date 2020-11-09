"""
Filename:    teleconnections.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions for importing teleconnection indices
"""

# Import Python modules

import os, sys
import numpy as np
import pandas as pd

        
def build_teleconnection_df(ref, col, start_date, end_date):
    '''
    Gets all available indices and compiles them into single dataframe.
    ref - frequency ('monthly', 'daily', 'seasonal')
    col - column variable from index ('ANOM' or 'COND')
    '''
    # Pull all indices together
    df_sh = sh_index(ref, col) # SH
    df_ao = ao_index(ref, col) # AO
    df_pdo = pdo_index(ref, col) # PDO
    df_enso = enso_index(ref, col) # ENSO
    ## put all indices into list
    df_list = [df_ao, df_pdo, df_enso, df_sh]
    ## Loop through list to clip to start date and end date
    dfs = []
    for k, df in enumerate(df_list):
        ## Trim date range
        idx = (df.index >= start_date) & (df.index <= end_date + " 23:59:59")
        df = df.loc[idx]
        dfs.append(df)

    # create single df
    data = {'date':  dfs[0].index,
            'AO':    dfs[0].values,
            'PDO':   dfs[1].values,
            'ENSO':  dfs[2].values,
            'SH':    dfs[3].values}

    df_index = pd.DataFrame(data)
    df_index = df_index.set_index('date')
    
    return df_index

def fill_df_all_dates(df, start_date, end_date, ref, col):
    if ref == 'seasonal':
        freq = 'QS-DEC'
    elif ref == 'daily':
        freq = '1D'
    elif ref == 'monthly':
        freq = '1MS'

    # date array with all days
    dates_allDays = pd.date_range(start=start_date, end=end_date, freq=freq)
    arr_allDays = np.zeros(len(dates_allDays), dtype=np.float)
    arr_allDays[:] = np.nan

    # Loop over ar days and match to ar_full 
    for i, date in enumerate(df.index):
        idx = np.where(dates_allDays == date)
        arr_allDays[idx] = df[col].values[i]

    # Create dataframe
    data = {col:arr_allDays}
    df_all = pd.DataFrame(data, index=dates_allDays)
    df_all = df_all[col].ffill(axis = 0)
    # if COND switch to int
    if col == 'ANOM':
        df_all = df_all
    elif col == 'COND':
        df_all = df_all.astype(int)

    return df_all

def sh_index(ref, col):
    '''
    Creates a df of the siberian high index
    '''
    path_to_out  = '/home/nash/DATA/repositories/AR_types/out/' 
    fname = path_to_out + 'SH_index_ERA5_1979_2019.csv'
    df = pd.read_csv(fname, engine='python')
    df = df.rename(columns={'Unnamed: 0': 'date', 'SH': 'ANOM'})
    df = df.set_index(pd.to_datetime(df.date))
    # neutral = 0, positive = 1, negative = -1
    df['COND'] = 0
    df.loc[df['ANOM']>0, 'COND'] = 1
    df.loc[df['ANOM']<0, 'COND'] = -1
    
    if ref == 'seasonal':
        df = df[col]
    elif ref == 'monthly':
        df = fill_df_all_dates(df, df.index[0], df.index[-1], ref, col)
    elif ref == 'daily':
        df = fill_df_all_dates(df, df.index[0], df.index[-1], ref, col)

    return df

def ao_index(ref, col):
    '''create a pandas df of the Arctic Oscillation index 
    on daily or monthly timescale 
    with positive or negative conditions identified as a column
    '''
    path_to_data = '/home/nash/DATA/data/teleconnection_indices/' 
    fname_daily = 'norm.daily.ao.index.b500101.current.ascii'
    fname_monthly = 'AO_CPC_NOAA_monthly_index_1950_2019.txt'

    if ref == 'monthly':
        names=['year', 'month', 'ANOM']
        df = pd.read_csv(path_to_data + fname_monthly,
                        delim_whitespace=True, engine='python', header=None, names=names)
        df['date'] = pd.date_range('1950-01-01', '2019-11-01', freq='1MS')
        df = df.set_index('date')
        
        df['3_month_running'] = df.loc[:,'ANOM'].rolling(window=3).mean()
        # neutral = 0, positive = 1, negative = -1
        df['COND'] = 0
        df.loc[df['3_month_running']>0, 'COND'] = 1
        df.loc[df['3_month_running']<0, 'COND'] = -1

    elif ref == 'daily':
        names=['YEAR', 'MON', 'DAY', 'ANOM']
        df = pd.read_csv(path_to_data + fname_daily,
                        delim_whitespace=True, engine='python', header=0, names=names)

        # neutral = 0, positive = 1, negative = -1
        df['COND'] = 0
        df.loc[df['ANOM']>0, 'COND'] = 1
        df.loc[df['ANOM']<0, 'COND'] = -1

        df['date'] = pd.date_range('1950-01-02 9:00:00', '2019-02-28 9:00:00', freq='1D')
        df = df.set_index('date')
        
    elif ref == 'seasonal':
        names=['YEAR', 'MON', 'DAY', 'ANOM']
        df = pd.read_csv(path_to_data + fname_daily,
                        delim_whitespace=True, engine='python', header=0, names=names)
        df['date'] = pd.date_range('1950-01-02 9:00:00', '2019-02-28 9:00:00', freq='1D')
        df = df.set_index('date')
        
        df = df.resample('QS-DEC').mean()
        # reset conditions neutral=0, positive=1, negative=-1
        df['COND'] = 0
        df.loc[df['ANOM']>0, 'COND'] = 1
        df.loc[df['ANOM']<0, 'COND'] = -1

    return df[col]

def pdo_index(ref, col):
    path_to_data = '/home/nash/DATA/data/teleconnection_indices/' 
    fname = 'NOAA_PDO_index.csv'
    # PDO monthly index
    df = pd.read_csv(path_to_data+fname, engine='python', skiprows=1)
    df['date'] = pd.date_range('1854-01', '2020-07', freq='MS')
    df = df.rename(columns={'Value': 'ANOM'})
    df = df.set_index('date')
    # set conditions neutral=0, nino=1, nina=2
    df['COND'] = 0
    df.loc[df['ANOM']>0, 'COND'] = 1
    df.loc[df['ANOM']<0, 'COND'] = -1

    if ref == 'monthly':
        df = df[col]
    
    elif ref == 'seasonal':
        df = df.resample('QS-DEC').mean()
        # reset conditions neutral=0, positive=1, negative=-1
        df['COND'] = 0
        df.loc[df['ANOM']>0, 'COND'] = 1
        df.loc[df['ANOM']<0, 'COND'] = -1
        df = df[col]
        
    elif ref == 'daily':
        df = fill_df_all_dates(df, df.index[0], df.index[-1], ref, col)

    return df

def enso_index(ref, col):
    path_to_data = '/home/nash/DATA/data/teleconnection_indices/'
    fname = 'CPC_NCEP_NOAA_ONI.txt'
    df = pd.read_csv(path_to_data+fname, delim_whitespace=True, engine='python')
    df['date'] = pd.date_range('1949-12', '2019-08', freq='MS')
    df = df.set_index('date')
    df['NINOCOND'] = df['ANOM']
    df.loc[df['ANOM']>=0.5, 'NINOCOND'] = 1
    df.loc[df['ANOM']< 0.5, 'NINOCOND'] = 0

    df['NINACOND'] = df['ANOM']
    df.loc[df['ANOM']<=-0.5, 'NINACOND'] = 1
    df.loc[df['ANOM']> -0.5, 'NINACOND'] = 0

    df['ninocount'] = df['NINOCOND'].rolling(min_periods=1, window=9, center=True).sum()
    df['ninacount'] = df['NINACOND'].rolling(min_periods=1, window=9, center=True).sum()

    # set conditions neutral=0, nino=1, nina=-1
    df['COND'] = 0
    df.loc[df['ninocount']>=5, 'COND'] = 1
    df.loc[df['ninacount']>=5, 'COND'] = -1

    ## Fix a couple weird ones
    seasons = ['ASO', 'JFM', 'JJA', 'JAS', 'ASO', 'SON']
    yrs = [1984, 1984, 2008, 2008, 2008, 2008]

    for i in np.arange(6):
        df.loc[(df['SEAS'] == seasons[i]) & (df['YR'] == yrs[i]), 'COND'] = 0

    if ref == 'monthly':
        df = df[col]
    elif ref == 'seasonal': 
        df = df.resample('QS-DEC').mean()
        df['COND'] = df['COND'].astype(int)
        df = df[col]
    elif ref == 'daily':
        df = fill_df_all_dates(df, df.index[0], df.index[-1], ref, col)

    return df