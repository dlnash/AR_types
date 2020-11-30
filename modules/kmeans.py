""" Functions used in K-means Clustering """

# IMPORTS

import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns


# FUNCTIONS

def plot_optimal_k(data, kmax, create_plot=False, filename=None):
    """ Elbow plot to determine optimal number of clusters (k) 
    
    Parameters
    ----------
    data : array_like
        data to perform clustering on where rows= obs and cols = variables
    kmax : scalar, int
        maximum number of clusters; iteration over the interval k=[1,kmax]
    filename : string (optional)
        filename or path to save figure. Default is None (fig not saved).
    """
    
    # Cohesion = sum of sq dist of samples to their cluster center
    kclusters = np.arange(kmax, dtype=int) + 1
    cohesion = np.empty(kmax)
    
    # iterate over kclusters (k=1,kmax)
    for i in range(kmax):
        km = KMeans(n_clusters=kclusters[i])
        km = km.fit(data)
        cohesion[i] = (km.inertia_)
    if create_plot == True:
        # Elbow plot
        fig, ax = plt.subplots()
        sns.set_style("whitegrid")
        ax.plot(kclusters, cohesion, marker='o', linewidth=2.0, markersize=7.0)
        ax.set_title('Elbow Plot for Optimal K')
        ax.set_ylabel('Sum of Sq Dist (cohesion)')
        ax.set_xlabel('k (# of clusters)')
        ax.set_xticks(kclusters)

        if filename:
            plt.savefig(filename)
    else:
        return kclusters, cohesion

    plt.show()
