"""
Collection of functions used to analyze time series (1-dimensional) data

"""

# Imports

import os, sys
import numpy as np
import xarray as xr
from scipy.stats import ttest_1samp, t


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
