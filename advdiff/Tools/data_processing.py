# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 09:56:47 2022

This module includes some data processing methods.

@author: Ketil
"""

import numpy            as np
import pandas           as pd
import xarray           as xr
from Tools.interpolator import meshA_to_meshB


def concat_dict(attrs_list, liststostr=True):
    """
    Concatenate values from similar dictionaries with different variables.
    Equal variables will be ignored.

    Parameters
    ----------
    attrs_list : list
                 List of dictionaries with the same keys.
    liststostr : bool
                 If one value becomes a list, turn it into a string. Netcdf files cannot store multidimensional lists as attributes,
                 so this is a compromise.

    Returns
    -------
    dict
        Dictionary containing the concatenated values.
    """
    temp_dict = attrs_list[0]
    for item in attrs_list[1:]:
        for key, value in item.items():
            try: # Try to add value to key in dict
                if np.any(temp_dict[key] == value):
                    pass
                else:
                    if not isinstance(temp_dict[key], list):
                        temp_dict[key] = [temp_dict[key]]
                    temp_dict[key].append(value)
            except KeyError: # Key doesnt exist, add key-value pair
                temp_dict[key] = value
    if liststostr: # Netcdf cant store lists in attributes, make the list into strings instead...
        for key, value in temp_dict.items():
            if isinstance(value, list):
                temp_dict[key] = str(value)
    return temp_dict


def get_source_location_and_probability(source):
    """
    Extract source location and corresponding location_probability from a source xarray object.

    Parameters
    ----------
    source : xarray dataset
             An xarray dataset containing a source

    Returns
    -------
    tuple
        A tuple containing location and location_probability
    """
    xs       = source.x.values
    ys       = source.y.values
    location = np.array([xs,ys])                               # Local origin in global grid
    location_probability = source['location_probability'].data # Get location probability from current source
    return location, location_probability
    

def weighted_sum(field, weights=None, sum_dim='source'):
    """
    Weight a xarray variable with respect to a set of weights and sum over a specific dimension. If weights are None, do not weigh, only sum

    Parameters
    ----------
    field   : xarray dataset
              An xarray data array.
    weights : array, optional
              Array of weights, by default None
    sum_dim : str, optional
              dimension to sum over, by default 'source'

    Returns
    -------
    xarray dataset
        A xarray data array.
    """
    
    if weights is not None:
        field = field.weighted(weights) # Weight the variable of interest w.r.t location probabilities
        field = field.sum(dim=sum_dim)  # Sum the variable over the sources
    else:
        field = field.sum(dim=sum_dim)   
    
    return field


def get_source_attributes(sources):
    """
    Extract the attributes from the source xarray file.

    Parameters
    ----------
    sources : xarray dataset
              A xarray dataset.

    Returns
    -------
    tuple
        Tuple containing a dictionary of source attributes and number of sources.
    """
    sources_attrs = sources.attrs

    return sources_attrs


def get_velocity_attributes(velocity):
    """
    Extract the attributes from the velocity xarray file.
    Only extracts the filename and path as some velocity files contain too many unimportant attributes.

    Parameters
    ----------
    velocity : xarray dataset
               A xarray dataset.

    Returns
    -------
    dict
        Dictionary containing velocity file name.
    """
    velocity_attrs = velocity.attrs

    return velocity_attrs


def get_probe_attributes(probes):
    """
    Extract the attributes from the probe xarray file.

    Parameters
    ----------
    probes : xarray dataset
             A xarray dataset.

    Returns
    -------
    tuple
        Tuple containing a dictionary of probe attributes, probe type and number of probes
    """
    try:
        probes_attrs = probes.attrs
    except:
        probes_attrs = {}
    return probes_attrs


def get_global_attributes(config, sources_attrs, velocity_attrs, probes_attrs):
    """
    Generates the global attributes which are in common for all output. It is a combination of source attributes, probe attributes and velocity attributes.

    Parameters
    ----------
    config         : config file
                     A config file containing relevant attributes.
    sources_attrs  : dict
                     Source attributes.
    velocity_attrs : dict
                     Veocity attributes.
    probes_attrs   : dict
                     Probe attributes.

    Returns
    -------
    dict
        Dictionary containing global attributes.
    """
    attrs      = sources_attrs
    attrs.update(probes_attrs)
    attrs.update(velocity_attrs)
    attrs.update({'Sweeps'              : int(config['setup']['sweeps'])   if config['setup']['sweeps'].isdigit() else config['setup']['sweeps'],
                  'Boundary condition'  : config['setup']['BC_type'],
                  'Convection type'     : config['setup']['convection_type'],
                  'Diffusion type'      : config['setup']['diffusion_type'],
                  'D'                   : float(config['setup']['D']),
                  'Comp mesh'           : config['grid']['type'],
                  'Out mesh'            : config['grid']['type_out'],
                  'Extrapolate NANs'    : config['grid']['extrap_nans'],
                  'Lx'                  : float(config['grid']['Lx']),
                  'Ly'                  : float(config['grid']['Ly']),
                  'nx'                  : int(config['grid']['nx']),
                  'ny'                  : int(config['grid']['ny']),
                  'maxvol'              : float(config['grid']['maxvol']),
                  'minvol'              : float(config['grid']['minvol']),
                  'Fill type'           : config['velocity']['fill_type'],
                  'Downsampling'        : config['output']['downsample'],
                  'Downsampling factor' : float(config['output']['downsampling_factor'])})
    return attrs

    
def get_mean_dt(dataset):
    """
    Generates the mean time-step length (dt) of a dataset in seconds.

    Parameters
    ----------
    dataset : xarray dataset
              An xarray dataset containing a time dimension.

    Returns
    -------
    mean_dt : float 
              A float describing the mean time-step length (dt) in seconds.
    """
    mean_dt = (dataset['time'].diff(dim='time')/pd.Timedelta(1, 's')).mean().data
    return mean_dt


def get_sweeps(config, mean_dt):
    """
    Computes the number of sweeps to perform per time-step dt.
    If the config file defines the number of sweeps as 'Dynamic', then information about mean_dt and max_dt
    will be used to compute the number of required sweeps.

    config file requires:
        
        * (int/str) : config['setup']['sweeps']
        * (float)   : config['setup']['max_dt']
        * (str)     : config['setup']['max_dt_unit']

    Parameters
    ----------
    config  : setup.ini
              Config file with relevant variables and setup configurations.    
    mean_dt : float 
              Float describing the mean dt length of the velocity dataset.

    Returns
    -------
    sweeps : int
             An integer describing how many swipes to perform per dt.
    """
    # Dynamic sweeps
    if config['setup']['sweeps'] == 'Dynamic' or config['setup']['sweeps'] == 'dynamic':
        max_dt      = float(config['setup']['max_dt'])
        max_dt_unit = str(config['setup']['max_dt_unit'])
        sweeps = int(np.ceil(pd.Timedelta(float(mean_dt), 's') / pd.Timedelta(float(max_dt), max_dt_unit)))
    # Set sweeps
    else:
        sweeps = int(config['setup']['sweeps'])
    return sweeps


def from_unstructured_to_structured(dataset, structured_mesh):
    """
    Transform an unstructured dataset in the form f(...,'num') to a structured dataset of the form f(...,'x','y')

    Parameters
    ----------
    dataset : xarray dataset
        A xarray dataset in unstructured form.
    structured_mesh : fp.Grid2D
        A structured mesh

    Returns
    -------
    xarray dataset
        A xarray dataset in structured form.
    """
    
    dataset = dataset.transpose(..., 'num')

    x, y = structured_mesh.cellCenters
    xl   = np.unique(x)
    yl   = np.unique(y)
    x_shape = xl.shape[0]
    y_shape = yl.shape[0]
    
    coords = {'x' : (['x'], xl), 
              'y' : (['y'], yl),}
    
    keys = list(dataset.keys())
    for coord in list(dataset.coords):
        if coord == 'x' or coord == 'y':
            pass
        else:
            coords[coord] = ([coord], dataset[coord].data)
    
    unstructured_points = np.stack([dataset.x.data, dataset.y.data], axis=-1)
    structured_points   = np.stack([x,              y],              axis=-1)
    interpolator = meshA_to_meshB(meshA_points=unstructured_points, 
                                  meshB_points=structured_points)
    
    data_vars = {}
    for key in keys:
        dataarray = dataset[key]
        data      = dataarray.data
        
        out_dims = list(dataarray.dims) + ['y', 'x']
        out_dims.remove('num')
        
        shape = list(data.shape)
        if len(shape[:-1]) > 0:
            in_shape  = tuple([np.prod(shape[:-1]), shape[-1]])
            out_shape = tuple(shape[:-1] + [y_shape, x_shape])
        else:
            in_shape  = tuple([1] + shape)
            out_shape = tuple([y_shape, x_shape])

        data = data.reshape(in_shape)
        slices = []
        for data_slice in data[:]:
            slices.append(interpolator.interpolate(data_slice, 
                                                   extrap_nans=True, 
                                                   extrap_from_A=True))
        
        data_int = np.array(slices).reshape(out_shape)
        data_vars[key] = (out_dims, data_int)
        
    out = xr.Dataset(data_vars = data_vars, 
                        coords = coords)
    out.attrs = dataset.attrs
    out = out.transpose(...,'x','y')
            
    return out