"""
Filename:    eofs.py
Author:      Tessa Montini, tmontini@ucsb.edu
Description: Functions used to calculate and anaylze EOFs
"""

## Imports

import os, sys
import numpy as np
import xarray as xr


## FUNCTIONS
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

def flatten_array(arr_list):
    '''Flatten data arrays
    
    Parameters
    ----------
    arr_list : list
        list of variable arrays
  
    Returns
    -------
    arr_list : list of flattened arrays
        arrays are flattened
    
    '''
    for i, in_array in enumerate(arr_list):
        # Extract variable as numpy array
        var1 = in_array.values
        ntim, nlat, nlon = var1.shape
        # flatten spatial locations to 1D
        arr_list[i] = np.reshape(var1, [ntim, nlat*nlon])
    return arr_list

def valid_nan(in_array):
    ''' check to see if missing values were removed properly'''
    inan = np.isnan(in_array)
    return (inan.any(axis=0) == inan.all(axis=0)).all()

def remove_nans(arr_list):
    '''Flatten data arrays and remove missing values
    
    Parameters
    ----------
    arr_list : list
        list of variable arrays
  
    Returns
    -------
    arr_list : list of arrays with missing data removed
        arrays are flattened and any spatial points that have missing data at any time step are removed
    
    '''
    arr_list_nan = []
    for i, in_array in enumerate(arr_list):
        # Extract variable as numpy array
        var1 = in_array
        ntim, nlocs = var1.shape
        locs = np.arange(nlocs)
        times = np.arange(ntim)
        # put flattened array into xarray dataarray object
        X = xr.DataArray(var1, coords=[times, locs], dims=['time', 'space'])
        # remove nans from the space dimension
        X_nan = X.dropna(dim='space', how='any')
        ## Test if the removal of nans was successful
        print('Nans removed success is ', valid_nan(X_nan.values))

        arr_list_nan.append(X_nan.values)
        arr_list[i] = X.values
        
    return arr_list, arr_list_nan

def standardize_arrays(arr_list, mode='t', dispersion_matrix='cor'):
    '''put variables into single flattened array and then standardize
     TO DO: change so have remove time-mean function again
     s-mode does not need a second removal of time-mean 
     Parameters
     ----------
     arr_list : list
        list of variable arrays
     
     mode : str
         mode of EOF - t or s
     
     dispersion_matrix : str
         type of dispersion matrix - cor or cov
         
     Returns
     -------
     X : single data matrix with all variables stacked
        arrays are standardized by the mode and dispersion matrix types
     
     ''' 
    print('EOF mode: ', mode)
    print('Dispersion Matrix: ', dispersion_matrix)
    nvar = len(arr_list)
    
    for i, var1 in enumerate(arr_list):
        # Data dimensions with missing values removed
        ntim, npts = var1.shape
        
        # if t-mode
        if mode == 't':
            # empty flat array to put variables in
            Xs = np.empty((nvar*npts,ntim))
            # transpose to [space x time]
            X1 = var1.T
            # Combine variables into single data matrix Xs
            Xs[i*npts:(i+1)*npts,:] = X1
        # if s-mode
        else:
            # empty flat array to put variables in
            Xs = np.empty((ntim, nvar*npts))
            # keep array as [time x space]
            X1 = var1    
            # Combine variables into single data matrix Xs
            Xs[:, i*npts:(i+1)*npts] = X1
            
    # Standardize by columns and remove column mean for ALL variables
    x1mean = np.mean(Xs, axis=0)
    x1std = np.std(Xs, axis=0)
    if dispersion_matrix == 'cor':
        # Standardize by columns (if correlation)
        # remove column mean
        X = (Xs - x1mean) / x1std
    else: ## dispersion matrix == cov (covariance)
        # remove column mean
        X = (Xs - x1mean)
                       
    print(X.shape)

    # Check that column means=0 and std dev=1
    test = np.mean(np.mean(X, axis=0))
    print("Column means: ", np.round(test,2))
    test = np.mean(np.std(X, axis=0))
    print("Column std: ", np.round(test,2))
    
    return X

def standardize_arrays_new(arr_list, mode='t', dispersion_matrix='cor'):
    '''standardize variables then put in single flattened array
     
     Parameters
     ----------
     arr_list : list
        list of variable arrays
     
     mode : str
         mode of EOF - t or s
     
     dispersion_matrix : str
         type of dispersion matrix - cor or cov
         
     Returns
     -------
     X : single data matrix with all variables stacked
        arrays are standardized by the mode and dispersion matrix types
     
     ''' 
    print('EOF mode: ', mode)
    print('Dispersion Matrix: ', dispersion_matrix)
    nvar = len(arr_list)
    ntim, npts = arr_list[0].shape
    
    # empty flat array to put variables in
    if mode == 't':
        Xs = np.empty((nvar*npts,ntim))
    else: # mode is s
        Xs = np.empty((ntim, nvar*npts))
    
    for i, var1 in enumerate(arr_list):
        # if t-mode
        if mode == 't':
            # transpose to [space x time]
            X1 = var1.T
            # Standardize by columns and remove column mean for ALL variables
            x1mean = np.mean(X1, axis=0)
            x1std = np.std(X1, axis=0)
            if dispersion_matrix == 'cor':
                # Standardize by columns (if correlation)
                # remove column mean
                X = (X1 - x1mean) / x1std
            else: ## dispersion matrix == cov (covariance)
                # remove column mean
                X = (X1 - x1mean)
            
            # Combine variables into single data matrix Xs
            Xs[i*npts:(i+1)*npts,:] = X
        
        # if s-mode
        else:
            # keep array as [time x space]
            X1 = var1
            # Standardize by columns and remove column mean for ALL variables
            x1mean = np.mean(X1, axis=0)
            x1std = np.std(X1, axis=0)
            if dispersion_matrix == 'cor':
                # Standardize by columns (if correlation)
                # remove column mean
                X = (X1 - x1mean) / x1std
            else: ## dispersion matrix == cov (covariance)
                # remove column mean
                X = (X1 - x1mean)
        
            # Combine variables into single data matrix Xs
            Xs[:, i*npts:(i+1)*npts] = X
                       
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
        # results in [time x time]
        # Eigenvector decomposition of R
        evals, evecs = np.linalg.eig(R)
    else:
        # Compute covariance/correlation matix [R]
        ## Should this be??
        # R = np.matmul(z.T,z) / (ntot - 1.) [space x space]
        ntot = z.shape[1]
        R = np.matmul(z,z.T) / (ntot - 1.)

        # Eigenvector decomposition of R
        evals, evecs = np.linalg.eig(R)
    
    return R, evals, evecs


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
        tmp = np.matmul(evecs[:, 0:npcs].T, z)
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
