"""
Filename:    eofs.py
Author:      Tessa Montini, tmontini@ucsb.edu
Description: Functions used to calculate and anaylze EOFs
"""

## Imports

import os, sys
import numpy as np


## FUNCTIONS
def center_data(arr_list):
    ''' Remove the mean of an array along the first dimension.
    
    If *True*, the mean along the first axis of *dataset* (the
    time-mean) will be removed prior to analysis. If *False*,
    the mean along the first axis will not be removed. Defaults
    to *True* (mean is removed).
    The covariance interpretation relies on the input data being
    anomaly data with a time-mean of 0. Therefore this option
    should usually be set to *True*. Setting this option to
    *True* has the useful side effect of propagating missing
    values along the time dimension, ensuring that a solution
    can be found even if missing values occur in different
    locations at different times.
    '''
    for i, in_array in enumerate(arr_list):
        # Compute the mean along the first dimension.
        mean = in_array.mean(axis=0, skipna=False)
        # Return the input array with its mean along the first dimension
        # removed.
        arr_list[i] = in_array - mean
    return arr_list

def spatial_weights(arr_list):
    """Spatial weights
    
    Returns a 1D array of weights equal to the sqrt of the cos of latitude.
    
    Parameters
    ----------
    arr_list : list
        list of variable arrays
  
    Returns
    -------
    weighted_arrays : list of arrays
        weights equal to the sqrt of the cosine of latitude
    
    Example
    -------
    # Apply spatial weights using xarray
    wgts = spatial_weights(lats)            # compute weights
    era['wgts'] = ('latitude', wgts)        # add `wgts` to dataset
    era['uwnd_wgt'] = era.uwnd * era.wgts   # apply wgts to data variable
    
    """
    for i, in_array in enumerate(arr_list):
        latitude = in_array.lat
        # convert lats from degrees to radians
        lat_rad = np.deg2rad(latitude)
        # compute weights
        weights = np.sqrt(np.cos(lat_rad))
        # apply spatial weights to array
        arr_list[i] = in_array*weights
    return arr_list

def remove_missing_values(X):
    '''remove missing values from centered, flattened array '''
    # Find the indices of values that are missing in the whole row
    nonMissingIndex = np.where(np.logical_not(np.isnan(X[0])))[0]
    # Remove missing values from the design matrix.
    dataNoMissing = X[:, nonMissingIndex]
    return nonMissingIndex, dataNoMissing

def valid_nan(in_array):
    ''' check to see if missing values were removed properly'''
    inan = np.isnan(in_array)
    return (inan.any(axis=0) == inan.all(axis=0)).all()

def standardize_and_flatten_arrays(arr_list, mode='t'):
    '''standardize variables for EOF and put into single flattened array ''' 
    nvar = len(arr_list)
    # Data dimensions
    ntim, nlat, nlon = arr_list[0].shape
    npts = nlat*nlon
    # create empty flat array to put standarized vars in
    if mode == 't':
        Xs = np.empty((nvar*npts,ntim))
    else:
        Xs = np.empty((ntim, nvar*npts))
    
    for i, in_array in enumerate(arr_list):
        # Extract variable as numpy array
        var1 = in_array.values

        # Reshape into 2D arrays by flattening the spatial dimension
        tmp1 = np.reshape(var1, (ntim, npts))

        # Remove missing data
        tmp1_idx, tmp1_miss = remove_missing_values(tmp1)

        ## Test if the removal of nans was successful
        print(valid_nan(tmp1_miss))
        
        # Data dimensions with missing values removed
        ntim, npts = tmp1_miss.shape
        
        # if t-mode
        if mode == 't':
            # transpose to [space x time]
            X1 = tmp1_miss.T
            # Standardize by columns
            x1std = np.std(X1, axis=0)
            X1s = X1 / x1std
            # Combine variables into single data matrix Xs
            Xs[i*npts:(i+1)*npts,:] = X1s

        # if s-mode
        else:
            # keep array as [time x space]
            X1 = tmp1_miss
           # Standardize by columns
            x1std = np.std(X1, axis=0)
            X1s = X1 / x1std
            # Combine variables into single data matrix Xs
            Xs[:, i*npts:(i+1)*npts] = X1s
    
    print(Xs.shape)

    # Check that column means=0 and std dev=1
    test = np.mean(np.mean(Xs, axis=0))
    print("Column means: ", np.round(test,2))
    test = np.mean(np.std(Xs, axis=0))
    print("Column std: ", np.round(test,2))
    
    return Xs

def calc_eofs(z, mode='t'):
    """Eigenvector decomposition of covariance/correlation matrix
    
    Parameters
    ----------
    z : array_like, float
        matrix of data values [nxp];
        n = number of observations (rows);
        p = number of variables (columns)
    mode : str
        mode which you are running EOF
        't' for t-mode
        's' for s-mode
    Returns
    -------
    evals : array_like, float
        vector of eigenvalues of size [p]
    evecs : array_like, float
        pxp matrix of eigenvectors (columns)
    
    """
    if mode == 't':
        # Compute covariance/correlation matix [R]
        ntot = z.shape[0]
        R = np.matmul(z.T,z) / (ntot - 1.)

        # Eigenvector decomposition of R
        evals, evecs = np.linalg.eig(R)
    else:
        # Compute covariance/correlation matix [R]
        ntot = z.shape[1]
        R = np.matmul(z,z.T) / (ntot - 1.)

        # Eigenvector decomposition of R
        evals, evecs = np.linalg.eig(R)

    return evals, evecs


def calc_pcs(z, evecs, npcs, mode='t'):
    """Calculate principal components from eigenvalues and eigenvectors
    
    Parameters
    ----------
    z : array_like, float
        standardized data matrix
    evecs : array_like, float
        pxk matrix of eigenvectors (columns), where k<=p (may be truncated)
    npcs : scalar, int
        number of pcs to return
    mode: str
        mode for pcs - s or t
    Returns
    -------
    pcs : array_like, float
        .........
    
    """   
    if mode == 't':
        tmp = np.matmul(z, evecs[:,0:npcs])
        pcs = tmp.T
    else:
        tmp = np.matmul(evecs[0:npcs, :], z)
        pcs = tmp
    
    return pcs


def loadings(evals, evecs, neofs):
    """Calculate loadings matrix
    
    Parameters
    ----------
    evals : array_like, float
        Vector of eignvalues
    evecs : array_like, float
        Matrix of eigenvectors
    neofs : scalar, int
        number of eofs to calculate loadings for
             
    Returns
    -------
    loadings : array_like, float
        Loadings matrix
    
    """
    evals_tr = evals[0:neofs]
    evecs_tr = evecs[:, 0:neofs]
    loadings = evecs_tr * np.sqrt(evals_tr)

    return loadings


def calc_eofs_svd(z, neofs):
    """Singular value decomposition of data matrix
    
    Parameters
    ----------
    z : array_like, float
        2d array of standardized data values of size [n x p];
        n = number of observations (rows);
        p = number of variables (columns)
    
    neofs : scalar, int
        number of eofs to return (used in loadings and pcs calculation)
        
    Returns
    -------
    evals : array_like, float
        vector of eigenvalues of size [p]
    evecs : array_like, float
        array of eigenvectors (columns); size [p x p]
    loadings : array_like, float
        loadings matrix
    pcs : array_like, float
        principal components
    
    """    
    # Singular Value Decomposition of z
    U, S, Vt = np.linalg.svd(z)
    ntot = z.shape[0]

    # Compute eigenvalues
    evals = S**2.0 / (ntot-1)
    
    # Compute eigenvectors
    evecs = Vt.T
    
    # Compute loadings *** add neof
    loadings = evecs * S / np.sqrt(ntot-1.0)
    
    # Compute principal components
    tmp = U[:,0:neofs] * S[0:neofs]
    pcs = tmp.T

    return evals, evecs, loadings, pcs


def pct_variance(evals, neofs=None):
    """Explained variance of EOFs
    
    Calculates the percent of the total variance explained by each EOF.
    
    Parameters
    ----------
    evals : array_like, float
        Array of eigenvalues
    neofs : scalar, int
        Number of eigenvectors to return the percent variance for.
        Defaults to all eigenvalues.
    
    Returns
    -------
    pct_var : array_like, float
        Percent variances for each EOF
        
    """    
    slicer = slice(0, neofs)
    pct_var = evals[slicer] / np.sum(evals) * 100.
    
    return pct_var


def north_test(evals, n):
    """North Test for separation of eigenvalues
    
    Parameters
    ----------
    evals : array_like, float
        Array of eigenvalues
    n : scalar, float
        number of independent samples
        
    Returns
    -------
    error : array_like, float
        Array of errors scaled by the variance fraction (%)
    
    """   
    tmp = evals * np.sqrt(2.0/n)
    error = tmp / np.sum(evals) * 100
    
    return error


# def get_pcs_std(z, npcs, npts):
#     """Calculate standardized PC scores using SVD
    
#     **** need to verify equations... 
#     different from pc calculatin in calc_eofs_svd **** 
    
#     """
#     U, S, Vt = np.linalg.svd(z)
#     tmp = np.sqrt(npts-1) * U[:,0:npcs]
#     pcs_std = tmp.T
    
#     return pcs_std
