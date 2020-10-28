"""
Filename:    teleconnections.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions for importing teleconnection indices
"""

# Import Python modules

import os, sys
import numpy as np
import pandas as pd


def AO_index(ref, path_to_data):
    '''create a pandas df of the AO index on daily or monthly timescale with positive or negative conditions identified as a column
    '''
    fname_daily = 'teleconnection_indices/norm.daily.ao.index.b500101.current.ascii'
    fname_monthly = 'teleconnection_indices/AO_CPC_NOAA_monthly_index_1950_2019.txt'

    if ref == 'monthly':
        names=['year', 'month', 'anom']
        df = pd.read_csv(path_to_data + fname_monthly,
                        delim_whitespace=True, engine='python', header=0, names=names)

        df['3_month_running'] = df.loc[:,'anom'].rolling(window=3).mean()
        df['COND'] = 'NEUTRAL'
        df.loc[df['3_month_running']>0, 'COND'] = 'POSITIVE'
        df.loc[df['3_month_running']<0, 'COND'] = 'NEGATIVE'
        df = df.loc[(df['year'] >= 1980) & (df['year'] <= 2018)]
        AO_df = df.loc[(df['month'] == 12) | (df['month'] == 3) | (df['month'] == 6) | (df['month'] == 9)]
    
    if ref == 'daily':
        names=['YEAR', 'MON', 'DAY', 'ANOM']
        df = pd.read_csv(path_to_data + fname_daily,
                        delim_whitespace=True, engine='python', header=0, names=names)
    
        df['COND'] = 'NEUTRAL'
        df.loc[df['ANOM']>0, 'COND'] = 'POSITIVE'
        df.loc[df['ANOM']<0, 'COND'] = 'NEGATIVE'
        
        df['AO'] = 0
        df.loc[df['ANOM']>0, 'AO'] = 1
        df.loc[df['ANOM']<0, 'AO'] = -1

        df['date'] = pd.date_range('1950-01-02 9:00:00', '2019-02-28 9:00:00', freq='1D')
        df = df.set_index('date')
    
    return df