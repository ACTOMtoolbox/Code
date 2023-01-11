# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:13:52 2022

This module includes all methods which load and initialize files based on information gathered from the config file.

@author: Ketil
"""
import logging
import random
import os
import xarray as xr
import numpy  as np
import pandas as pd
from Tools.risk_to_positions     import risk_to_sources
from Tools.coord_transform       import coord_converter as cc
from Tools.coord_transform       import add_UTM_coords
from Tools.interpolator          import meshA_to_meshB
from Tools.data_processing       import get_mean_dt
from Tools.parser_tools          import strornum, strtoBool
from Tools.post_processing       import store_xarray


def concat_sources(sources_A, sources_B, point_tags_A, point_tags_B):
    """
    Concatenate two source and two point tag files together.
    The labels of set B are numbered after the labels of set A.

    Parameters
    ----------
    sources_A    : xarray dataset
                   An xarray dataset containing sources
    sources_B    : xarray dataset
                   An xarray dataset containing sources
    point_tags_A : xarray dataset
                   An xarray dataset containing point tags
    point_tags_B : xarray dataset
                   An xarray dataset containing point tags

    Returns
    -------
    Tuple of xarray datasets
        First element of the tuple contains sources, the second contains the corresponding point tags.
    """

    labels_B  = sources_B.source.data + int(np.max(sources_A.source.data + 1))
    sources_B = sources_B.assign({'source': (['source'], labels_B)}) # Wells are labled after riskmap labels

    point_tags_B_labels = point_tags_B.source.data + int(np.max(point_tags_A.source.data + 1))
    point_tags_B        = point_tags_B.assign({'source' : (['num'], point_tags_B_labels)})

    sources    = xr.concat([sources_A,    sources_B],    dim='source')
    point_tags = xr.concat([point_tags_A, point_tags_B], dim='num')
    
    return sources, point_tags


def interpolate_velocity_over_time(velocity, freq='15T', method='linear'):
    """
    Interpolate the global velocity-field over time. This is used to unify datsets to specific time-step intervals. 

    Parameters
    ----------
    velocity : xarray dataset
        An xarray dataset containing the global velocity-field stored on an unstructured grid.
    freq : str
        A string describing what frequency to interpolate at, default 15T which equals intervals of 15 minutes. Check pandas date_range() documentation for alternatives.

    Returns
    -------
    xarray dataset
        An xarray dataset containing the global velocity-field interpolated over time with a regular time-step interval.
    """
    time_range = pd.date_range(start=np.min(velocity.time.data), end=np.max(velocity.time.data), freq=freq)
    velocity   = velocity.interp(time=time_range, method=method)
    return velocity


def is_unstructured(dataset):
    """
    Check if dataset is unstructured or structured (According to the conventions I have dealt with so far.)

    Parameters
    ----------
    dataset : xarray dataset
              An xarray dataset.

    Returns
    -------
    bool
         True if dataset is unstructured, False if dataset is structured.
    """
    if 'node' in dataset.dims or 'num' in dataset.dims:
        return True   # Dataset in terms of f(t,node)
    else:
        return False  # Dataset in terms of f(t,y,x)


def select_depth(dataset, depth='average'):
    """
    Select depth from dataset.

    Parameters
    ----------
    dataset : xarray dataset
              An xarray dataset containing depth layers
    depth   : int or str, optional
              Int describing what layer to extract. If depth = 'average', take average of all layers, by default 'average'.
              If depth coordinate not found, do nothing.

    Returns
    -------
    xarray dataset
           An xarray dataset
    """
    if 'depth' in dataset.dims:
        if depth == 'average':
            dataset = dataset.mean(dim='depth')
        else:
            dataset = dataset.isel(depth=depth)
    else: # No depth found...
        dataset = dataset
    return dataset


def treat_structured_velocity(velocity, depth='average'):
    """
    Prepare raw velocity file of the form f(t,x,y) for use with AdvDiff.
    Transforms velocity file from f(t,x,y) to f(t,node). This is to unify the inputs to AdvDiff to unstructured form.
    Any NANs in a Raw velocity file will be interpreted as coastline or walls with zero velocity.

    Parameters
    ----------
    velocity : xarray dataset
               An xarray dataset object containing a unstructured velocity file of the form f(t,x,y).
    depth    : int, optional
               An int describing what depth layer to extract. If depth = 'average', then take an average of all layers, by default 'average'.
               If no depth layer found, do nothing.

    Returns
    -------
    xarray dataset
           An xarray dataset containing velocities of the form f(t,node)
    """

    velocity = add_UTM_coords(velocity)             # Incase velocity is defined in longitude latitude coordinates, redefine into xy UTM
    velocity = velocity.fillna(0.0)                 # Any NANs in a Raw velocity file will be interpreted as coastline or walls with zero velocity.
    velocity = velocity.transpose('time', 'y', 'x') # VERY IMPORTANT: If this is not done, then files of the form f(t,x,y) will give wrong output!
    
    # Select relevant depth or take average of all depths
    velocity = select_depth(velocity, depth)

    x = velocity.x.data
    y = velocity.y.data
    t = velocity.time.data

    u = velocity.u.data.reshape(t.shape[0],x.shape[0]*y.shape[0])
    v = velocity.v.data.reshape(t.shape[0],x.shape[0]*y.shape[0])

    X,Y = np.meshgrid(x,y)
    x = X.flatten()
    y = Y.flatten()

    data_vars = {"u" : (['time','node'], u), 
                 "v" : (['time','node'], v)}
            
    velocity = xr.Dataset(data_vars = data_vars, 
                        coords={"x"    : (['node'],  x), 
                                "y"    : (['node'],  y), 
                                "time" : (['time'],  t),})
    return velocity


def treat_unstructured_velocity(velocity, depth='average'):
    """
    Prepare raw unstructured velocity file of the form f(t,node) for use with AdvDiff.
    Any NANs in a Raw velocity file will be interpreted as coastline or walls with zero velocity.

    Parameters
    ----------
    velocity : xarray dataset
               An xarray dataset object containing a unstructured velocity file of the form f(t,node).
    depth    : int, optional
               An int describing what depth layer to extract. If depth = 'average', then take an average of all layers, by default 'average'.
               If no depth layer found, do nothing.

    Returns
    -------
    xarray dataset
           An xarray dataset containing velocities of the form f(t,node)
    """

    velocity = add_UTM_coords(velocity) # Incase velocity is defined in longitude latitude coordinates, redefine into xy UTM
    velocity = velocity.fillna(0.0)     # Any NANs in a Raw velocity file will be interpreted as coastline or walls with zero velocity.

    # Select relevant depth or take average of all depths
    velocity = select_depth(velocity, depth)

    u = velocity.u.data
    v = velocity.v.data
    x = velocity.x.data
    y = velocity.y.data 
    t = velocity.time.data

    data_vars = {"u" : (['time','node'], u), 
                 "v" : (['time','node'], v)}
            
    velocity = xr.Dataset(data_vars = data_vars, 
                          coords={"x"    : (['node'],  x), 
                                  "y"    : (['node'],  y), 
                                  "time" : (['time'],  t),})
    return velocity


def unstructured_to_structured(velocity, extrap_nans=True, fill_type='nearest extrapolation', Lxy=[None,None]):
    """
    Function which allows us to convert a unstructured xarray velocity file to a structured xarray file.
    It uses linear interpolation to compute the velocity for all (x,y) coordinates on a (100,100) grid. Also converts
    coordinates from longitude-latitude to x-y UTM.

    Parameters
    ----------
    velocity    : xarray dataset
                  An xarray unstructured velocity file (for example a GOM velocity file).
    extrap_nans : bool
                  Whether or not to extrapolate missing velocities with nearest neighbour interpolation. The default is True.
    fill_type   : str
                  What fill type to use if extrap_nans is False.
    Lxy         : array
                  Array containing local grid width Lx and local grid height Ly. If none, then bounding box will be the minimal size.

    Returns
    -------
    velocity : xarray dataset
               An xarray velocity file suitable for our solver.

    """
    # Construct meshes to interpolate between unstructured to structured
    meshA_points, meshB_points, X, Y, Xl, Yl = construct_meshes_from_unstructured(X=velocity.x.data, Y=velocity.y.data, outx=100, outy=100, Lxy=Lxy)  
    
    unstructured_mesh_to_structured_mesh = meshA_to_meshB(meshA_points=meshA_points, meshB_points=meshB_points, extrap_qhull=False)

    if fill_type == 'nearest extrapolation':
        extrap_nans   = True
        fill_value_u  = np.nan
        fill_value_v  = np.nan

    elif fill_type == 'zeroes':
        extrap_nans  = False
        fill_value_u = 0.0
        fill_value_v = 0.0

    elif fill_type == 'average':
        extrap_nans  = False
        fill_value_u = None
        fill_value_v = None
    
    tmp = []
    for ii, time in enumerate(velocity.time.data):
    
        u = velocity.isel(time=ii).u.data.flatten()
        v = velocity.isel(time=ii).v.data.flatten()
        
        U = unstructured_mesh_to_structured_mesh.interpolate(u, extrap_nans=extrap_nans, extrap_from_A=False, fill_value=fill_value_u)
        V = unstructured_mesh_to_structured_mesh.interpolate(v, extrap_nans=extrap_nans, extrap_from_A=False, fill_value=fill_value_v)
        
        U = np.transpose(U.reshape(X.shape[0], X.shape[1]))
        V = np.transpose(V.reshape(Y.shape[0], Y.shape[1]))
        
        out=xr.Dataset({
                "u":(['x','y'], U),
                "v":(['x','y'], V)
                }, 
                coords={"x":(['x'],Xl),
                        "y":(['y'],Yl),
                        "time":(['time'],[time]),
                })
    
        tmp.append(out)
        
    velocity = xr.concat(tmp, dim='time')
    velocity = velocity.transpose('time', 'y', 'x') # VERY IMPORTANT: If this is not done, then files of the form f(t,x,y) will give wrong output!
    
    return velocity


def construct_meshes_from_unstructured(X, Y, outx=100, outy=100, Lxy=[None,None]):
    """
    Generates meshpoints from the dataset to a [outx,outy] grid

    Parameters
    ----------
    X    : array
           Array containing x coordinates for each point.
    Y    : array
           Array containing y coordinates for each point.       
    outx : int
           Number of points along x. The default is 100.
    outy : int
           Number of points along y. The default is 100.
    Lx   : float
           Local grid width Lx
    Ly   : float
           Local grid width Ly
             

    Returns
    -------
    meshA_points : array ([n_points,2])
                   Array of points from dataset.
    meshB_points : array ([n_points,2])
                   Array of points from [outx,outy] grid.
    X  : array ([x_points,y_points])
         Meshgrid X
    Y  : array ([x_points,y_points])
         Meshgrid Y
    Xl : ([x_points])
         np.linspace(x)
    Yl : ([y_points])
         np.linspace(y)

    """
    
    x = X.flatten()
    y = Y.flatten()

    Lx, Ly = Lxy[0], Lxy[1]

    Lx = 0.0 if Lx is None else Lx
    Ly = 0.0 if Ly is None else Ly

    Xl = np.linspace(min(x)-Lx/2, max(x)+Lx/2, num=outx) # Generate bounding box
    Yl = np.linspace(min(y)-Ly/2, max(y)+Ly/2, num=outy)

    X, Y = np.meshgrid(Xl, Yl)  # 2D grid for interpolation

    meshA_points = (np.stack((x,y),axis=1))
    meshB_points = (np.stack((X.flatten(),Y.flatten()),axis=1))
    
    return meshA_points, meshB_points, X, Y, Xl, Yl


def get_velocity(inpath, filenm, drop_variables=None, time_setup=None, verbose=True):
    """
    Get velocity from netcdf file. 

    Parameters
    ----------
    inpath         : str
                     Directory for where to find file.
    filenm         : str
                     Filename to read.
    drop_variables : str
                     What variables to drop if any are problematic. The default is None. ('depth' causes a lot of trouble)
    time_setup     : dict
                     Dictionary describing what times to extract. The default is None. If none, the entire file is extracted.
    verbose        : bool
                     Sets function to verbose.

    Raises
    ------
    ERROR: If time start is larger than max time in velocity file, raise error and terminate program.

    Returns
    -------
    out : xarray dataset
          An xarray velocity file.
    """
        
    df  = xr.open_dataset(inpath+filenm, drop_variables=drop_variables).load()
    out = df.copy()
    df.close() 
    
    if 'time' not in out.coords:       
        out = out.rename({'ocean_time':'time'}) # Rename time coordinate if it is incorrectly named.
    
    if time_setup is not None:
        time_start      = time_setup['time_start']
        time_delta      = time_setup['time_delta']
        time_delta_unit = time_setup['time_delta_unit']
        time_seed       = time_setup['time_seed']
        
        ### SET START TIME ###
        if time_start == 'Start':
            time_start = out['time'].data.min()   

        elif time_start == 'Random':
            random.seed(time_seed)
            min_time     = out['time'].data.min()
            max_time     = out['time'].data.max() - pd.Timedelta(time_delta, time_delta_unit) if time_delta != 'End' else out['time'].data.max()
            random_start = min_time + (max_time - min_time) * random.random()
            time_start   = out['time'].sel(time=random_start, method='nearest').data

        else:
            time_start = out['time'].sel(time=time_start, method='nearest').data.min()
        
        ### SET END TIME ###
        if time_delta != 'End' and time_delta > 0.0:
            time_end = time_start + pd.Timedelta(time_delta, time_delta_unit)

        else:
            time_end = out['time'].data.max()

        # Handle exceptions where the user is sloppy with defining time_start and time_delta...
        if time_end > out['time'].data.max():
            logging.warning('Time end exceeds maximum allowed time. Will set time delta according to velocity file...\n') if verbose else None
            time_end = out['time'].data.max()

        if time_start < out['time'].data.min():
            logging.warning('Time start exceeds minimum allowed time. Will set time start according to velocity file...\n') if verbose else None
            time_start = out['time'].data.min()

        elif time_start > out['time'].data.max():
            logging.exception('Time start cannot be higher than the maximum allowed time...\n') if verbose else None
            raise Exception
        
        out = out.sel(time=slice(time_start, time_end))
    
    return out


def get_sources(inpath, filenm):
    """
    Read and return source netcdf files with point tags.

    Parameters
    ----------
    inpath : str
             Directory path of where to look for sources.
    filenm : str
             Filename of what file to read.

    Returns
    -------
    Tuple of xarray datasets
        First element of the tuple contains sources, the second contains the corresponding point tags.
    """
    sources    = xr.open_dataset(inpath+filenm).load()
    sources.close()
    sources    = add_UTM_coords(dataset=sources)
    point_tags = sources.drop('source').rename({'source':'num'}).assign({'source':(['num'],sources.source.data)})
    return sources, point_tags


def get_wells(inpath, filenm):
    """
    Get wells from either a .nc file or a .csv file. 

    Parameters
    ----------
    inpath : str
             Directory path for where to look for file.
    filenm : str
             Filename of what file to read.

    Returns
    -------
    Tuple of xarray datasets
        First element of the tuple contains sources, the second contains the corresponding point tags.
    """

    if '.nc' in filenm:
        sources, point_tags = get_sources(inpath, filenm)
    elif '.csv' in filenm or '.CSV' in filenm:
        sources, point_tags = get_wells_CSV(in_path=inpath, file_nm=filenm, lp_wells=1)
    
    return sources, point_tags


def get_wells_CSV(in_path, file_nm, lp_wells='Random', lp_min=0.1, lp_max=100):
    """
    Loads well locations from CSV file and returns an xarray dataset which contains x,y coordinates for each well, similar to the source datasets.
    Also returns a point_tags dataset.
    You can set random location probability or a constant location probability for all wells.

    Parameters
    ----------
    in_path  : str
               Directory path for where to look for csv file.
    file_nm  : str
               Filename of csv to load.
    lp_wells : str or float, optional
               If 'Random' then each well will be assigned a random location probability between the two bounds. 
               If a float then every well will be given the same location probability, by default 'Random'
    lp_min   : float, optional
               Lower bound for random location probability, by default 0.1
    lp_max   : int, optional
               Upper bound for random location probability, by default 100
    """

    def reject_outliers(data, m=3.5):
        data_norm  = np.abs(data - np.median(data))
        mdata_norm = np.median(data_norm)
        s = data_norm / mdata_norm if mdata_norm else 0.
        return data[s < m]
    
    wells    = pd.read_csv(os.path.join(os.path.dirname(__file__), '../'+in_path+file_nm))

    well_lonlat_83 = np.unique(np.stack([wells['Long83'].values,wells['Lat83'].values], axis=-1), axis=0)
    well_lonlat_27 = np.unique(np.stack([wells['Long27'].values,wells['Lat27'].values], axis=-1), axis=0)
    well_lon = np.nanmean(np.array([well_lonlat_83[:,0], well_lonlat_27[:,0]]), axis=0)
    well_lat = np.nanmean(np.array([well_lonlat_83[:,1], well_lonlat_27[:,1]]), axis=0)
    well_lon = reject_outliers(well_lon, m=5.0)
    well_lat = reject_outliers(well_lat, m=5.0)
    
    coord_converter = cc((np.min(well_lon), np.max(well_lon), np.min(well_lat), np.max(well_lat)), store_AOI=False)
    well_x, well_y  = coord_converter.lonlat_to_xy(well_lon, well_lat)
    positions_wells = np.stack([well_x, well_y], axis=-1)
    labels_wells    = np.arange(0, positions_wells.shape[0], 1)
    
    if lp_wells == 'Random' or lp_wells == 'random':
        location_probability_wells = np.random.uniform(lp_min, lp_max, positions_wells.shape[0])
    else:
        location_probability_wells = lp_wells * np.ones((positions_wells.shape[0]))
    
    wells_ds_rm = xr.Dataset(data_vars = {"location_probability" : (['source'], location_probability_wells)}, 
                             coords    = {"source"               : (['source'], labels_wells),
                                          "x"                    : (['source'], positions_wells[:,0]), 
                                          "y"                    : (['source'], positions_wells[:,1]),})

    wells_ds_pt = xr.Dataset(data_vars = {"source"               : (['num'], labels_wells),
                                          "location_probability" : (['num'], location_probability_wells)}, 
                             coords    = {"x"                    : (['num'], positions_wells[:,0]), 
                                          "y"                    : (['num'], positions_wells[:,1]),})

    return wells_ds_rm, wells_ds_pt


def load_velocity(config, return_unstructured=True, verbose=True):
    """
    Load velocity by name given in config file.

    config file requires:
        
        * (str)   : config['paths']['indata_path']
        * (str)   : config['paths']['velocity_path']
        * (str)   : config['velocity']['velocity_file']
        * (int)   : config['velocity']['depth']
        * (str)   : config['velocity']['fill_type']
        * (str)   : config['setup']['time_start']
        * (str)   : config['setup']['time_delta']
        * (str)   : config['setup']['time_delta_unit']
        * (float) : config['setup']['time_seed]

    Parameters
    ----------
    config              : setup.ini
                          Config file describing file path, filename etc to read.
    return_unstructured : bool
                          If True, returns the xarray dataset in unstructured form. This is what AdvDiff needs to run. Structured form is better
                          for when you want to visualize the loaded dataset for example during debugging.
    verbose             : bool
                          If True, print information about the loading process to the logger, else stay quiet.

    Returns
    -------
    velocities      : xarray dataset
                      An xarray dataset containing u and v components of velocity for (t,x,y) (and possibly multiple velocities)
    velocity_inpath : str
                      Pathname for velocity file.
    file            : str
                      Filename for velocity file.
    mean_dts        : list
                      List of mean time-step sizes in seconds. List has only one element for one velocity file.

    """
    velocity_inpath  = config['paths']['indata_path']+config['paths']['velocity_path']
    file_name        = config['velocity']['velocity_file']
    depth            = strornum(config['velocity']['depth'])
    use_custom_timer = strtoBool(config['setup']['use_custom_timer'])

    fill_type = str(config['velocity']['fill_type'])
    freq      = str(config['setup']['dt_size']+config['setup']['dt_unit'])

    if use_custom_timer:
        time_start      = config['setup']['time_start']        
        time_delta      = strornum(config['setup']['time_delta']) 
        time_delta_unit = str(config['setup']['time_delta_unit'])
        time_seed       = float(config['setup']['time_seed'])
        time_setup      = {'time_start'     : time_start,
                           'time_delta'     : time_delta,
                           'time_delta_unit': time_delta_unit,
                           'time_seed'      : time_seed}
    else:
        time_start      = 'Start'
        time_delta      = 'End'
        time_delta_unit = str(config['setup']['time_delta_unit'])
        time_seed       = float(config['setup']['time_seed'])
        time_setup      = None
    
    velocity = get_velocity(inpath=velocity_inpath, filenm=file_name, drop_variables='depth', time_setup=time_setup, verbose=verbose)

    # AdvDiff can utilize two different velocity formats, unstructured and structured. However they need different treatments to be unified.
    if is_unstructured(dataset=velocity):
        logging.info(f'Unstructured velocity file detected {file_name}\n') if verbose else None
        velocity = treat_unstructured_velocity(velocity, depth)   
    else:
        logging.info(f'Structured velocity file detected {file_name}\n') if verbose else None
        velocity = treat_structured_velocity(velocity, depth)   

    # Interpolate the velocity over time with a set stepsize.
    velocity = interpolate_velocity_over_time(velocity=velocity, freq=freq) 

    # Mostly for debugging. Not actually relevant for the AdvDiff module.
    if not return_unstructured: 
        velocity = unstructured_to_structured(velocity, extrap_nans=True, fill_type=fill_type, Lxy=[None,None])

    # Lower memory requirement
    velocity = velocity.astype('float32')   

    # Add attributes and delete all other attributes as GOM files are bloated
    velocity.attrs = {'Velocity file'   : file_name,
                      'Depth'           : depth,
                      'Start time'      : time_start,
                      'Time delta'      : time_delta,
                      'Time delta unit' : time_delta_unit,
                      'Time seed'       : time_seed if time_start == 'Random' else 'None',
                      'Start date'      : str(velocity.time.min().data),
                      'End date'        : str(velocity.time.max().data),
                      'mean_dt'         : get_mean_dt(dataset=velocity)}  

    return velocity, velocity_inpath, file_name


def load_sources(config):
    """
    Load sources by filename given in config file. (Either through riskmap or sources).

    config file requires:
        
        * (str)   : config['paths']['indata_path']
        * (str)   : config['paths']['sources_path']
        * (str)   : config['sources']['source_file']
        * (bool)  : config['sources']['get_source_from_file']
        * (str)   : config['paths']['riskmap_path']
        * (str)   : config['riskmaps']['risk_file']
        * (str)   : config['riskmaps']['cluster']
        * (float) : config['riskmaps']['threshold']
        * (float) : config['riskmaps']['eps']
        * (int)   : config['riskmaps']['min_samples']
        * (int)   : config['riskmaps']['n_clusters']

    Parameters
    ----------
    config : setup.ini
             Config file describing file path, filename etc to read.

    Returns
    -------
    sources     : xarray dataset 
                  An xarray dataset containing (x,y) coordinates for each source.
    source_file : str 
                  Filepath and filename for printing.
    point_tags  : xarray dataset
                  An xarray dataset containing (x,y) coordinates and corresponding source tags for each point in the riskmap above a certain threshold.

    """
    get_source_from = config['setup']['get_source_from'].lower() # Force lowercase to ensure case neutrality

    # Get the parameters for the different types of files.
    if 'riskmap' in get_source_from:
        risk_inpath         = config['paths']['indata_path']+config['paths']['riskmap_path']
        risk_file           = config['riskmaps']['risk_file']
        outpath             = config['paths']['outdata_path']
        cluster             = config['riskmaps']['cluster']
        threshold           = float(config['riskmaps']['threshold'])
        eps                 = float(config['riskmaps']['eps'])
        min_samples         = int(config['riskmaps']['min_samples'])
        n_clusters          = int(config['riskmaps']['n_clusters'])
    if 'well' in get_source_from:
        wells_inpath        = config['paths']['indata_path']+config['paths']['well_path']
        wells_file          = config['wells']['well_file']
    if 'source' in get_source_from:
        source_inpath       = config['paths']['indata_path']+config['paths']['sources_path']
        source_file         = config['sources']['source_file']

    # Load sources and point tags depending on file type chosen
    if get_source_from == 'sources': # Incase we already have a source file we want to use.
        sources, point_tags = get_sources(inpath=source_inpath, filenm=source_file)

    elif get_source_from == 'riskmap': # Get sources from a riskmap
        sources, point_tags = risk_to_sources(filenm=risk_file, inpath=risk_inpath, outpath=outpath, threshold=threshold, 
                                              cluster=cluster, eps=eps, min_samples=min_samples, n_clusters=n_clusters)

        source_inpath = risk_inpath
        source_file   = risk_file

    elif get_source_from == 'wells': # Get sources from a well file (.nc or .csv)
        sources, point_tags = get_wells(inpath=wells_inpath, filenm=wells_file)

        source_inpath = wells_inpath
        source_file   = wells_file

    elif 'riskmap' in get_source_from and 'wells' in get_source_from: # Get sources from riskmap and wells
        RM_sources, RM_point_tags = risk_to_sources(filenm=risk_file, inpath=risk_inpath, outpath=outpath, threshold=threshold, 
                                                    cluster=cluster, eps=eps, min_samples=min_samples, n_clusters=n_clusters)

        wells, wells_point_tags = get_wells(inpath=wells_inpath, filenm=wells_file)
        sources, point_tags = concat_sources(RM_sources, wells, RM_point_tags, wells_point_tags)

        source_inpath = [risk_inpath, wells_inpath]
        source_file   = [risk_file,   wells_file]
    

    elif 'riskmap' in get_source_from and 'sources' in get_source_from: # Get sources from riskmap and sources
        RM_sources, RM_point_tags = risk_to_sources(filenm=risk_file, inpath=risk_inpath, outpath=outpath, threshold=threshold,
                                                    cluster=cluster, eps=eps, min_samples=min_samples, n_clusters=n_clusters)
        
        sources, point_tags = get_sources(inpath=source_inpath, filenm=source_file)
        sources, point_tags = concat_sources(RM_sources, sources, RM_point_tags, point_tags)

        source_inpath = [risk_inpath, source_inpath]
        source_file   = [risk_file,   source_file]
                
    store_xarray(dataset=point_tags, outpath=config['paths']['outdata_path'], filenm=config['paths']['point_tags']) # Storing point tags in config['paths']['outdata_path']

    # Add attributes
    sources.attrs['Source file'] = source_file
    sources.attrs['Num sources'] = sources.dims['source']
    try:
        del sources.attrs['risk file'] # No need to specify what riskfile was used again
    except:
        pass

    return sources, source_inpath, source_file
    
        
def load_probes(config):
    """
    Load probes by name given in config file.
    
    config file requires:
        
        * (str)   : config['paths']['indata_path']
        * (str)   : config['paths']['probes_path']
        * (str)   : config['probes']['probes_file']

    Parameters
    ----------
    config : setup.ini
             Config file describing file path, filename etc to read.
       
    Returns
    -------
    probes        : xarray dataset
                    An xarray dataset containing (x,y) coordinates for probes.
    probes_inpath : str
                    Pathname for probe file.
    probes_file   : str
                    Filename for probe file.

    """
    probes_inpath = config['paths']['indata_path']+config['paths']['probes_path']
    probes_file   = config['probes']['probes_file']
    try:
        probes = xr.open_dataset(probes_inpath+probes_file).load()
        probes.close()
        probes = add_UTM_coords(dataset=probes)
        # Add attributes
        probes.attrs['Probe file']  = probes_file
        probes.attrs['Num probes']  = probes.dims['probe']
    except:
        probes = None
        logging.info('No file to define probe locations. Continue with no probes...\n')
    
    return probes, probes_inpath, probes_file