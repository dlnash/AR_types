""" Functions used in K-means Clustering """

# IMPORTS

import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns


# FUNCTIONS

def plot_optimal_k(data, neofs, kmax, create_plot=False, filename=None):
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
    kclusters = np.arange(1, kmax, dtype=int) + 1
    cohesion = np.empty((neofs-1, kmax-1))
    silhouette = np.empty((neofs-1, kmax-1))
    upper = np.empty((neofs-1, kmax-1))
    lower = np.empty((neofs-1, kmax-1))
    
    # iterate over EOFs
    for i in range(neofs-1):
        X = data[:,0:i+2]
        
        # iterate over kclusters
        for j in range(len(kclusters)):
            ki = kclusters[j]
            km = KMeans(n_clusters=kclusters[j])
            cluster_labels = km.fit_predict(X)
            km = km.fit(X)
            cohesion[i, j] = (km.inertia_)

            # The silhouette_score gives the average value for all the samples.
            # This gives a perspective into the density and separation of the formed
            # clusters
            silhouette[i, j] = silhouette_score(X, cluster_labels)
            # Compute the silhouette scores for each sample
            sample_silhouette_values = silhouette_samples(X, cluster_labels)
            
            spread = []
            for m in range(ki):
                # Aggregate the silhouette scores for samples belonging to
                # cluster i, and sort them
                ith_cluster_silhouette_values = \
                    sample_silhouette_values[cluster_labels == m]

                ith_cluster_silhouette_values.sort()
                spread.append(ith_cluster_silhouette_values.max())
            
            lower[i, j] = min(spread)
            upper[i, j] = max(spread)

            
    if create_plot == True:
        fig = plt.figure()
        fig.set_size_inches((6.0,10.0))
        fig.dpi = 300
        sns.set_style("whitegrid")
        colors = ['tab:blue', 'tab:red', 'tab:green', 'k']
        
        # Elbow plot
        ax = plt.subplot(2, 1, 1)
        for k in range(neofs-1):
            coh = cohesion[k, :]
            ax.plot(kclusters, coh, c=colors[k], marker='o', linewidth=2.0, markersize=7.0)
#             ax.set_title('Elbow Plot for Optimal K')
            ax.set_ylabel('Sum of Sq Dist (cohesion)')
            ax.set_xlabel('k (# of clusters)')
            ax.set_xlim(1., kmax+1)
            ax.set_ylim(0, 210)
            ax.set_xticks(kclusters)
        
        # Silhouette plot
        ax2 = plt.subplot(2, 1, 2)
        for k in range(neofs-1):
            sil = silhouette[k, :]
            ax2.plot(kclusters, sil, c=colors[k], marker='o', linewidth=2.0, markersize=7.0)
            ax2.fill_between(kclusters, lower[k, :], upper[k, :], fc=colors[k], ec=None, alpha=0.2)
            ax2.set_ylabel('Mean Silhouette')
            ax2.set_xlabel('k (# of clusters)')
            ax2.set_xlim(1., kmax+1)
            ax2.set_ylim(0, 1.)
            ax2.set_xticks(kclusters)
        
        if filename:
            # Save the figure
            fig.savefig(filename, bbox_inches='tight', dpi=fig.dpi)
            fig.clf()
    else:
        return kclusters, cohesion, silhouette, lower, upper
