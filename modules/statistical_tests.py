"""
Filename:    statistical_tests.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions for running different significance tests on tseries or ds objects
"""

# Import Python modules

import os, sys
import numpy as np
import pandas as pd
import xarray as xr
import xarray.ufuncs as xrf
from scipy import stats
from scipy.stats import ttest_1samp, t, pearsonr, norm, linregress
import scipy.stats.distributions as dist

def ttest_1samp_new(a, popmean, dim, n):
    """
    This is a two-sided test for the null hypothesis that the expected value
    (mean) of a sample of independent observations `a` is equal to the given
    population mean, `popmean`
    
    Inspired here: https://github.com/scipy/scipy/blob/v0.19.0/scipy/stats/stats.py#L3769-L3846
    
    Parameters
    ----------
    a : xarray
        sample observation
    popmean : float or array_like
        expected value in null hypothesis, if array_like than it must have the
        same shape as `a` excluding the axis dimension
    dim : string
        dimension along which to compute test
    
    Returns
    -------
    mean : xarray
        averaged sample along which dimension t-test was computed
    maskt_idx : array, bool
        Boolean array of where the tvalue is greater than the critical value
    """
    df = n - 1
    a_mean = a.mean(dim)
    d = a_mean - popmean
    v = a.var(dim, ddof=1)
    denom = xrf.sqrt(v / float(n))

    tval = d /denom
    # calculate the critical value
    cv = stats.distributions.t.ppf(1.0 - 0.05, df)
    maskt_idx = (abs(tval) >= cv)
#     prob = stats.distributions.t.sf(xrf.fabs(tval), df) * 2
#     prob_xa = xr.DataArray(prob, coords=a_mean.coords)
    return a_mean, maskt_idx

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

def independent_ttest(ds, group, alpha, n):
    '''Calculate statistical significance using 1-sample t-test of ds with lag coord
    and create mask of ds where values are considered significant'''
    nlat = len(ds.lat)
    nlon = len(ds.lon)
    nlag = len(ds.lag)
    # calculate t-statistic and p-value
    tval, pval = ttest_1samp_lag(ds.groupby(group), popmean=np.zeros([nlat, nlon, nlag]))
    # calculate the critical value
    df = n-2
    cv = t.ppf(1.0 - alpha, df)
    print('Critical t-value: ', cv)
    # interpret via critical value
    # if abs(tvalue) >= cv, reject null hypothesis that the means are equal
    maskt_idx = (abs(tval) >= cv)
    
#     # find p-value based on tval and nprime (rather than n)
#     pval = t.sf(np.abs(tval), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
#     # interpret via p-value
#     # if p < alpha, reject null hypothesis that the means are equal
#     maskp_idx = (pval < alpha)
    return maskt_idx

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

def _zstat_generic(value1, value2, std_diff, alternative, diff=0):
    '''generic (normal) z-test to save typing

    can be used as ztest based on summary statistics

    '''
    zstat = (value1 - value2 - diff) / std_diff
    if alternative in ['two-sided', '2-sided', '2s']:
        pvalue = norm.sf(np.abs(zstat))*2
    elif alternative in ['larger', 'l']:
        pvalue = norm.sf(zstat)
    elif alternative in ['smaller', 's']:
        pvalue = norm.cdf(zstat)
    else:
        raise ValueError('invalid alternative')
    return zstat, pvalue

def ztest(x1, x2=None, value=0, alternative='two-sided', usevar='pooled',
          ddof=1.):
    '''test for mean based on normal distribution, one or two samples

    In the case of two samples, the samples are assumed to be independent.

    Parameters
    ----------
    x1 : array_like, 1-D or 2-D
        first of the two independent samples
    x2 : array_like, 1-D or 2-D
        second of the two independent samples
    value : float
        In the one sample case, value is the mean of x1 under the Null
        hypothesis.
        In the two sample case, value is the difference between mean of x1 and
        mean of x2 under the Null hypothesis. The test statistic is
        `x1_mean - x2_mean - value`.
    alternative : string
        The alternative hypothesis, H1, has to be one of the following

           'two-sided': H1: difference in means not equal to value (default)
           'larger' :   H1: difference in means larger than value
           'smaller' :  H1: difference in means smaller than value

    usevar : string, 'pooled'
        Currently, only 'pooled' is implemented.
        If ``pooled``, then the standard deviation of the samples is assumed to be
        the same. see CompareMeans.ztest_ind for different options.
    ddof : int
        Degrees of freedom use in the calculation of the variance of the mean
        estimate. In the case of comparing means this is one, however it can
        be adjusted for testing other statistics (proportion, correlation)

    Returns
    -------
    tstat : float
        test statisic
    pvalue : float
        pvalue of the t-test

    Notes
    -----
    usevar not implemented, is always pooled in two sample case
    use CompareMeans instead.

    '''
    #usevar is not used, always pooled

    if usevar != 'pooled':
        raise NotImplementedError('only usevar="pooled" is implemented')

    x1 = np.asarray(x1)
    nobs1 = x1.shape[0]
    x1_mean = x1.mean(0)
    x1_var = x1.var(0)
    if x2 is not None:
        x2 = np.asarray(x2)
        nobs2 = x2.shape[0]
        x2_mean = x2.mean(0)
        x2_var = x2.var(0)
        var_pooled = (nobs1 * x1_var + nobs2 * x2_var)
        var_pooled /= (nobs1 + nobs2 - 2 * ddof)
        var_pooled *= (1. / nobs1 + 1. / nobs2)
    else:
        var_pooled = x1_var / (nobs1 - ddof)
        x2_mean = 0

    std_diff = np.sqrt(var_pooled)
    #stat = x1_mean - x2_mean - value
    return _zstat_generic(x1_mean, x2_mean, std_diff, alternative, diff=value)

def test_diff_proportion(ts1, ts2):
    '''
    Calculate the test statistic for testing 
    the difference in two population proportions
    
    ts1 : sample 1
    ts2 : sample 2
    
    return 
    Z : the test statistic
    p : the p-value
    
    '''
    n1 = len(ts1)
    n2 = len(ts2)
    y1 = ts1.sum()
    y2 = ts2.sum()
    
    p1 = y1/n1 # proportion of sample 1 
    p2 = y2/n2 # proportion of sample 2
    phat = (y1+y2)/(n1 + n2)
    
    std_err = np.sqrt(phat*(1-phat)*(1/n1 + 1/n2))
    Z = ((p1 - p2) - 0)/(std_err)
    
    # Calculate the  p-value
    # based on the standard normal distribution z-test
    pvalue = 2*dist.norm.cdf(-np.abs(Z)) # Multiplied by two indicates a two tailed testing.
    return Z, pvalue

def diff_proportion_zstat(df):
    '''
    enumerates through AR Types and Teleconnections to run z test on difference of proportion
    '''
    AR_CATS = ('AR_CAT1', 'AR_CAT2', 'AR_CAT3', 'AR_ALL')
    TELE = ['AO', 'PDO', 'ENSO', 'SH']
    zstat_array1 = []
    pval_array1 = []
    
    zstat_array2 = []
    pval_array2 = []
    
    zstat_array3 = []
    pval_array3 = []
    
    
    for i, tele in enumerate(TELE):
        
        s_positive = df.loc[(df[tele] > 0)]
        s_negative = df.loc[(df[tele] < 0)]
        s_neutral = df.loc[(df[tele] == 0)]
        
        for j, ar_type in enumerate(AR_CATS):
            zstat, pval = test_diff_proportion(s_positive[ar_type].values, s_negative[ar_type].values)
            zstat_array1.append((zstat))
            pval_array1.append((pval))
            
            zstat2, pval2 = test_diff_proportion(s_positive[ar_type].values, s_neutral[ar_type].values)
            zstat_array2.append((zstat2))
            pval_array2.append((pval2))
            
            zstat3, pval3 = test_diff_proportion(s_negative[ar_type].values, s_neutral[ar_type].values)
            zstat_array3.append((zstat3))
            pval_array3.append((pval3))
            
    
    zstat_final = [zstat_array1, zstat_array2, zstat_array3]
    pval_final = [pval_array1, pval_array2, pval_array3]
    return zstat_final, pval_final

def build_zscore_df(df):
    '''
    Creates a single df with zscore difference of proportion results
    '''
    z, p = diff_proportion_zstat(df)
    
    arrays = [['AO']*4 + ['PDO']*4 + ['ENSO']*4 + ['SH']*4,
               ['AR_CAT1', 'AR_CAT2', 'AR_CAT3', 'AR_ALL'] * 4]
    
    index = pd.MultiIndex.from_arrays(arrays, names=('Teleconnection', 'AR Type'))
    df_z1 = pd.DataFrame({'zstat': z[0],
                        'pval': p[0]},
                       index=index)
    df_z2 = pd.DataFrame({'zstat': z[1],
                        'pval': p[1]},
                       index=index)
    df_z3 = pd.DataFrame({'zstat': z[2],
                        'pval': p[2]},
                       index=index)
    
    df_z = [df_z1, df_z2, df_z3]
    return df_z