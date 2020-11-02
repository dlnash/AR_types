"""
Collection of functions used to analyze time series (1-dimensional) data

"""

# Imports

import os, sys
import numpy as np
import xarray as xr
import pandas as pd
from scipy.stats import ttest_1samp, t, pearsonr


# Define functions

def persistence(x):
    """A function that tags persistent events
    
    Given a binary time series `x`, this function tags persistent events, 
    (x eq 1). A sequence of consecutive x=1 values is tagged as a single event.
    
    Parameters
    ----------
    x : array_like
        Binary time series of x=0 or x=1
        
    Returns
    -------
    event_id : array_like, int
        Array of the same length as x that contains tags for each event
    duration = array_like, int
        Values represent the number of observations for for each event. 1D array with a length 
        
        
    Example
    -------
    Given:
        x          = [0,1,1,0,0,1,0,0,1,0,0,1,1,0,1,1,1,0]
    Returns:
        event_id   = [0,1,1,0,0,2,0,0,3,0,0,4,4,0,5,5,5,0]
        duration   = [  2,      1,    1,    2,    3,     ]
        
    References ******
    ----------
    Adapted from Charles Jones' IDL subroutine

    """
    
    # number of observations in x
    ntot = len(x)

    # Loop to tag persistent events
    event_id = np.zeros(ntot, dtype=int)
    tag = 1
    test = 0    
    for k in range(ntot):       
        # test for active event
        if (x[k] == 1):
            test = 1
            event_id[k] = tag
        # test for end of event    
        elif (x[k] == 0) and (test == 1):
            tag += 1
            test = 0

    # Loop to find duration of each event
    nevents = event_id.max()       # Total number of tagged events
    duration = np.empty(nevents)
    for k in range(nevents):
        # eid = event id
        eid = k+1
        # find the event indices in x
        idx = np.where(x == eid) 
        # find the length of the event and store in duration arr
        duration = len(idx[0])   
    
    return event_id, nevents

def ttest_1samp_lag(sample_obs, popmean):
    '''Wrapped stats.ttest_1samp to iterate ttest over sample in xr.dataset form
       The dataset should have coords ['time', 'lat', 'lon', 'lag']
       This variation is used specifically to run ttest for multiple vars and multiple lags
       
    Parameters
    ----------
    sample_obs : xarray ds
        grouped xarray ds object
    popmean: array_like, zeros
        array of repeating zeros that matches the dimensions of the sample_obs
        
    Returns
    -------
    tvalue : array_like, float
        xarray ds object with same vars as sample_obs that has the tvalue for the one-sample t-test
    pvalue : array_like, float
        xarray ds object with same vars as sample_obs that has the pvalue for the one-sample t-test 
    
    Example
    -------
    Given:
        sample_obs          = ds.groupby('time.season')
    Returns:
        tvalue   = ds_tvalue
        pvalue   = ds_pvalue
        
    '''
    return xr.apply_ufunc(ttest_1samp,
                          sample_obs,
                          input_core_dims=[['lag', 'time']],
                          dask='allowed', output_core_dims=[['lag'], ['lag']],
                          kwargs={'popmean': popmean, 'axis': -1, 'nan_policy': 'omit'})

def independent_ttest(ds, group, alpha, df):
    '''Calculate statistical significance using 1-sample t-test of ds with lag coord
    and create mask of ds where values are considered significant'''
    nlat = len(ds.lat)
    nlon = len(ds.lon)
    nlag = len(ds.lag)
    # calculate t-statistic and p-value
    tval, pval = ttest_1samp_lag(ds.groupby(group), popmean=np.zeros([nlat, nlon, nlag]))
    # calculate the critical value
    cv = t.ppf(1.0 - alpha, df)
    print('Critical t-value: ', cv)
    # interpret via critical value
    # if abs(tvalue) >= cv, reject null hypothesis that the means are equal
    maskt_idx = (abs(tval) >= cv)
    # interpret via p-value
    # if p < alpha, reject null hypothesis that the means are equal
    maskp_idx = (pval < alpha)
    return maskt_idx, maskp_idx


#### CLEAN UP AND ADD DOC STRINGS ###
def transition_matrix(x, states):

    morder = 1                 # model order (default=1) ** make optional keyword param
    nt = len(x)-1              # number of transitions
    n = len(states)            # number of states
    transc = np.zeros((n,n), dtype=int)  # Matrix of transitions counts (initialize counts as 0)

    # Loop to count transitions
    for t in range(nt):
        i = np.where(states == x[t])[0]      # i = indice of current state
        j = np.where(states == x[t+1])[0]    # j = indice of next state    
        transc[i,j] += 1                    # add a transition count for s[i] to s[j]    

    # Compute marginal totals (sum across rows in M)
    margin = np.sum(transc, axis=1)

    # Calculate probabilities (divide each row in M by its marginal total)
    probs = transc / margin[:,np.newaxis]
    
    return transc, probs


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


def select_months(df, mon_s, mon_e):
    # Select months
    if mon_s > mon_e:
        idx = (df.index.month >= mon_s) | (df.index.month <= mon_e)
    else:
        idx = (df.index.month >= mon_s) & (df.index.month <= mon_e)

    df = df.loc[idx]
    
    return df 

def autocorr_diff(ts1, ts2):
    '''Calculate autocorrelation of the difference between two paired samples'''
    # put into pandas df
    data = {'ts1':   ts1,
            'ts2':   ts2,
            'diff_ts':   ts1-ts2}

    df = pd.DataFrame(data)
    rho1 = df.diff_ts.autocorr(lag=1)
    
    return rho1
    
def n_prime(ts1, ts2):
    '''
    calculate the effective sample size/equivalent number of independent samples 
    from two paired samples of time series data that have high autocorrelation
    
    $n'  \approx n \frac{1-\rho_1}{1+\rho_1} $
    '''
    n = len(ts1)
    rho1 = autocorr_diff(ts1, ts2)
    nprime = n * ((1-rho1)/(1+rho1))
    
    return rho1, nprime

def ttest_autocorrelation(ts1, ts2, alpha):
    '''find pvalue and tvalue based on more stringent nprime (lag1 autocorrelation of 2 time series)
       for 2 sample 2-sided t-test
    '''
    # get degrees of freedom
    if autocorr == True:
        # calculate lag-1 autocorrelation difference
        rho1, n = n_prime(ts1, ts2)
    if autocorr == False:
        n = len(ts1)
    df = n-2
    
    # calculate t-statistic and p-value
    tval, pval = ttest_ind(ts1, ts2)
    # calculate the critical value
    cv = t.ppf(1.0 - alpha, df)
#     print('Critical t-value: ', cv)
    
    # find p-value based on tval and nprime (rather than n)
    pval = t.sf(np.abs(tval), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
    # interpret via critical value
    # if abs(tvalue) >= cv, reject null hypothesis that the means are equal
#     if abs(tval) >= cv:
#         print('t-value is greater than critical value, we reject the null hypothesis that the means are equal')
#     elif abs(tval) <= cv:
#         print('t-value is less than critical value, we accept the null hypothesis that the means are equal')
#     # interpret via p-value
#     # if p < alpha, reject null hypothesis that the means are equal
#     if pval < alpha:
#         print('p-value is less than alpha, we reject the null hypothesis that the means are equal')
#     elif pval > alpha:
#         print('p-value is greater than alpha, we accept the null hypothesis that the means are equal')
    
    return tval, pval

def pearsonr_autocorrelation(ts1, ts2, alpha, autocorr=True):
    '''find pvalue and tvalue based on more stringent nprime (lag1 autocorrelation of 2 time series)
       for pearson-r correlation
    '''
    # correlation between time series
    r, p = pearsonr(ts1, ts2)
#     print('Original pearson r results: ', r, p)
    # get degrees of freedom
    if autocorr == True:
        # calculate lag-1 autocorrelation difference
        rho1, n = n_prime(ts1, ts2)
    if autocorr == False:
        n = len(ts1)
    df = n-2
    # calculate significance of correlation coefficient
    s_T = np.sqrt((1-r**2)/df)
    # get t-value
    tval = (r-0)/s_T
    # get p-value
    pval = t.sf(np.abs(tval), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
    # calculate the critical value
    cv = t.ppf(1.0 - alpha, df)
#     # interpret via critical value
#     # if abs(tvalue) >= cv, reject null hypothesis that the means are equal
#     if abs(tval) >= cv:
#         print('t-value is greater than critical value, considered statistically significant')
#     elif abs(tval) <= cv:
#         print('t-value is less than critical value, NOT considered statistically significant')
#     # interpret via p-value
#     # if p < alpha, reject null hypothesis that the means are equal
#     if pval < alpha:
#         print('p-value is less than alpha, considered statistically significant')
#     elif pval > alpha:
#         print('p-value is greater than alpha, NOT considered statistically significant')
        
    return tval, pval