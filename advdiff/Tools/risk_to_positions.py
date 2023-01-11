"""
This module includes riskmap treatment tools.

@author: Ketil & Guttorm
"""
import logging
import xarray              as xr
import numpy               as np
import matplotlib.pyplot   as plt
from os                    import remove
from sklearn.cluster       import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.cluster       import DBSCAN
from Tools.coord_transform import add_UTM_coords

#TODO: implement storing figures of risk maps and resulting sources.

def risk_to_sources(filenm, threshold=0.5, inpath='../Indata/', outpath='.',
    pname='location_probability', cluster='dbscan', eps=0.1, min_samples=10, n_clusters=4,
    store=False, verbose=False, plot=False):
    """
    Takes a riskmap and generates an xarray dataset containing source locations with corresponding location_probabilities.

    Parameters
    ----------
    filenm    : str
                Filename to read
    threshold : float
                Threshold of what to define as a potential source in the riskmap.
    inpath    : str
                Directory of where to look for the riskmap
    outpath   : str, optional
                Directory of where to store the result.
    pname     : str, optional
                Name of probability variable. Default is 'location_probability'
    cluster     : str
                  Type of clustering method to use. Allowable: dbscan or kmeans
    eps         : float, optional
                  Parameter to use with dbscan. Default is 0.1.
    min_samples : int, optional
                Parameter for use with dbscan. Default is 10.
    n_clusters  : int, optional
                Parameter for use with kmeans. Number of clusters. Default is 4.
    store   : bool, optional
              Whether or not to store the results. Default is False.
    verbose : bool, optional
              Whether or not to be verbose. Default is False.
    plot    : bool, optional
              Whether or not to plot. Default is False.

    Returns
    -------
    tuple
        (dataset, point_tags): dataset contains source locations. point_tags contain all point containing possible sources.
    """
    

    risk_map=xr.open_dataset(inpath + filenm).load()
    risk_map.close()

    #If riskmap is defined with lon-lat coordinates, convert to UTM
    risk_map = add_UTM_coords(risk_map, risk=True, coords={'x'   : 'coordX',
                                                           'y'   : 'coordY',
                                                           'lon' : 'coordLON',
                                                           'lat' : 'coordLAT'})
    
    # Move coordinates from data variables to coordinates
    risk_map=risk_map.assign_coords(xpos=risk_map.coordX)
    risk_map=risk_map.assign_coords(ypos=risk_map.coordY)
    risk_map=risk_map.drop(['coordX','coordY'])
    
    risk_map.attrs['max_prob']  = risk_map[pname].max().values
    risk_map.attrs['sum_prob']  = risk_map[pname].sum().values
    risk_map.attrs['risk file'] = filenm
    risk_map.attrs['threshold'] = threshold
    risk_map.attrs['cluster']   = cluster
    #Normalize
    risk_map = risk_map/risk_map[pname].max().values
    
    
    out=risk_map.where(risk_map[pname]>threshold, drop=True).stack(num=('x','y')).dropna('num').drop('num')
    out=out.rename({"xpos":"x","ypos":"y"}).copy()
    out=out.assign_attrs(risk_map.attrs)

    
    if cluster == 'dbscan':
        if verbose:
            logging.info('Clustering with DBSCAN\n')
        out, point_tags = risk_to_clusters(risk=out, cluster=cluster, eps=eps, min_samples=min_samples, verbose=verbose, plot=plot)
        
    elif cluster == 'kmeans':
        if verbose:
            logging.info('Clustering with KMeans\n')
        out, point_tags = risk_to_clusters(risk=out, cluster=cluster, n_clusters=n_clusters, verbose=verbose, plot=plot)
   
    if store:
        if verbose:
            logging.info('Storing in {}'.format(outpath+filenm.replace('.nc','_sources.nc')))
        # Hard coding deleting the file before storing. Get permision denied error if not. 
        try:
            remove(outpath+filenm.replace('.nc','_sources.nc'))
        except OSError:
            pass
        out.to_netcdf(outpath+filenm.replace('.nc','_sources.nc'), mode='w')

    return out, point_tags



def risk_to_clusters(risk, cluster='dbscan', eps=0.1, min_samples=10, n_clusters=4, verbose=False, plot=False, store=None):
    """
    Clusters together nearby riskmap points into sources.

    Parameters
    ----------
    risk        : xarray dataset
                  An xarray containing the riskmap. 
    cluster     : str
                  Type of clustering method to use. Allowable: dbscan or kmeans
    eps         : float, optional
                  Parameter to use with dbscan. Default is 0.1.
    min_samples : int, optional
                Parameter for use with dbscan. Default is 10.
    n_clusters  : int, optional
                Parameter for use with kmeans. Number of clusters. Default is 4.
    store   : bool, optional
              Whether or not to store the results. Default is False.
    verbose : bool, optional
              Whether or not to be verbose. Default is False.
    plot    : bool, optional
              Whether or not to plot. Default is False.

    Returns
    -------
    tuple
        (dataset, point_tags): dataset contains source locations. point_tags contain all point containing possible sources.
    """

    positions=np.c_[risk.x.data, risk.y.data]
    scaler = StandardScaler()
    X=scaler.fit_transform(positions)
    
    if cluster == 'dbscan':
        db     = DBSCAN(eps=eps, min_samples=min_samples).fit(X)         # Fit points into clusters
        labels = db.labels_                                              # Get what cluster each point belongs to
        
    elif cluster == 'kmeans':
        km     = KMeans(n_clusters=n_clusters, random_state=999).fit(X)  # Fit points into clusters (Also set a predictable seed to avoid confusion in the future)
        labels = km.labels_                                              # Get what cluster each point belongs to

    unique_labels = set(labels)

    risk = risk.assign({'source':('num', labels)})                              # Assign new variable with dimension num and value label
    out  = risk.location_probability.groupby(risk.source).sum().to_dataset()   # Group location probabilities by source labels and sum each group
    out  = out.assign_attrs(risk.attrs)

    x = [np.average(risk.where(risk.source==ii).dropna(dim='num').x, 
         weights=risk.where(risk.source==ii).dropna(dim='num').location_probability) 
         for ii in unique_labels]                                                      # Find the weighted average x pos for each group

    y = [np.average(risk.where(risk.source==ii).dropna(dim='num').y, 
         weights=risk.where(risk.source==ii).dropna(dim='num').location_probability) 
         for ii in unique_labels]                                                      # Find the weighted average y pos for each group

    out = out.assign_coords(x=('source', x), y=('source', y))
    
    if verbose and cluster == 'dbscan':
    # Number of clusters in labels, ignoring noise if present. 
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)  # Count number of clusters ignoring noisy points (w label = -1)
        n_noise_ = list(labels).count(-1)                            # Count number of noisy points
        logging.info('Estimated number of clusters:     %d' % n_clusters_)
        logging.info('Estimated number of noise points: %d' % n_noise_)
    
    elif verbose and cluster == 'kmeans':
        # Number of clusters in labels, ignoring noise if present.       
        logging.info('Number of clusters: %d' % n_clusters)
        
    
    point_tags = xr.Dataset(data_vars = {"source"               : (['num'], labels),
                                         "location_probability" : (['num'], risk.dropna(dim='num').location_probability.data)}, 
                              coords  = {"x"                    : (['num'], positions[:,0]), 
                                         "y"                    : (['num'], positions[:,1]),})
    point_tags = point_tags.assign_attrs(risk.attrs)
        

    if plot: 
        plt.figure(figsize=(9,9))
        out.plot.scatter(x='x',y='y',markersize='location_probability',hue='location_probability')

    return out, point_tags

