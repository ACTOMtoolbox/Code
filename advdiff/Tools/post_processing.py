# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:38:11 2022

This module includes all post-dataprocessing methods.

@author: Ketil
"""

import logging
import numpy            as np
import xarray           as xr
import time             as TIME
from scipy.spatial         import KDTree
from statistics            import mode
from Tools.directory_tools import get_files, make_directory
from Tools.data_processing import concat_dict
from Tools.parser_tools    import strtoBool


def set_units(dataset):
    """
    Sets units to dataset variables and coordinates.

    Parameters
    ----------
    dataset : xarray dataset
              A xarray dataset.

    Returns
    -------
    xarray dataset
        A xarray dataset containing units.
    """
    dataset_dims = dataset.dims
    dataset_vars = list(dataset.keys())
    if 'x' in dataset_dims:
        dataset['x'].attrs        = {'unit':'m'}
    if 'y' in dataset_dims:
        dataset['y'].attrs        = {'unit':'m'}
    if 't' in dataset_dims:
        dataset['t'].attrs        = {'unit':'datetime'}
    if 'time' in dataset_dims:
        dataset['time'].attrs     = {'unit':'datetime'}
    if 'dt' in dataset_dims:
        dataset['dt'].attrs       = {'unit':'s'}
    if 'delta_t' in dataset_dims:
        dataset['delta_t'].attrs  = {'unit':'s'}
    if 'source' in dataset_dims:
        dataset['source'].attrs   = {'unit':'label'}
    if 'x_source' in dataset_dims:
        dataset['x_source'].attrs = {'unit':'m'}
    if 'y_source' in dataset_dims:
        dataset['y_source'].attrs = {'unit':'m'}
    if 'x_grid' in dataset_dims:
        dataset['x_grid'].attrs   = {'unit':'m'}
    if 'y_grid' in dataset_dims:
        dataset['y_grid'].attrs   = {'unit':'m'}
    if 'x_pos' in dataset_dims:
        dataset['x_pos'].attrs    = {'unit':'m'}
    if 'y_pos' in dataset_dims:
        dataset['y_pos'].attrs    = {'unit':'m'}
    if 'dist' in dataset_dims:
        dataset['dist'].attrs     = {'unit':'m'}

    if 'C' in dataset_vars:
        dataset['C'].attrs    = {'unit':'kg/m^2'}
    if 'dC' in dataset_vars:
        dataset['dC'].attrs   = {'unit':'kg/m^2'}
    if 'u' in dataset_vars:
        dataset['u'].attrs    = {'unit':'m/s'}
    if 'v' in dataset_vars:
        dataset['v'].attrs    = {'unit':'m/s'}
    if 'mean' in dataset_vars:
        dataset['mean'].attrs = {'unit':'kg/m^2'}
    if 'max' in dataset_vars:
        dataset['max'].attrs  = {'unit':'kg/m^2'}
    if 'var' in dataset_vars:
        dataset['var'].attrs  = {'unit':'kg^2/m^4'}
    if 'delta_mean' in dataset_vars:
        dataset['delta_mean'].attrs     = {'unit':'kg/m^2'}
    if 'delta_max' in dataset_vars:
        dataset['delta_max'].attrs      = {'unit':'kg/m^2'}
    if 'delta_min' in dataset_vars:
        dataset['delta_min'].attrs      = {'unit':'kg/m^2'}
    if 'delta_mean_abs' in dataset_vars:
        dataset['delta_mean_abs'].attrs = {'unit':'kg/m^2'}
    if 'delta_max_abs' in dataset_vars:
        dataset['delta_max_abs'].attrs  = {'unit':'kg/m^2'}
    if 'delta_var' in dataset_vars:
        dataset['delta_var'].attrs      = {'unit':'kg^2/m^4'}

    if 'Mass' in dataset_vars:
        dataset['Mass'].attrs             = {'unit':'kg'}
    if 'Mass Theoretical' in dataset_vars:
        dataset['Mass Theoretical'].attrs = {'unit':'kg'}

    if 'x_source' in dataset_vars:
        dataset['x_source'].attrs   = {'unit':'m'}
    if 'y_source' in dataset_vars:
        dataset['y_source'].attrs   = {'unit':'m'}
    if 'source tag' in dataset_vars:
        dataset['source tag'].attrs = {'unit':'label'}
    
    return dataset


def generate_global_fields(files=None, data_dir='Outdata/', prefix='statistics', variables='ALL', 
                           nx=100, ny=100, weight_source=True, sum_source=False, fillna=0):
    """
    Generate global-field dataset given the statistics from each local output.

    Parameters
    ----------
    files     : list
                List of filepaths of netcdf files to load
    data_dir  : str 
                If files=None. Directory path for where to look for files. The default is 'Outdata/'.
    prefix    : str
                If files=None. Prefix of files to look for. The default is 'statistics'.
    variables : list
                List of str of variables in local output to generate global fields for. The default is 'All'.
    nx       : int
                The number of discretizations along the x coordinate. The default is 100.
    ny       : int
                The number of discretizations along the y coordinate. The default is 100.
    quantile  : float
                What quantile between fmax and fmin to cut. The default is 0.0.
    weight_source : bool
                    Whether or not to use location_probability to weight the contributions from each source. The default is True.
    sum_source    :  bool
                    Whether or not to sum the contributions from each source together, or to keep them separated by a source dimension. The default is False.
    fillna        : float or np.nan
                    What to fill missing values with when local grids are placed on the global grid.

    Returns
    -------
    global_fields : xarray dataset
                    An xarray dataset containing the global field.
    attrs         : dict
                    Corresponding attributes for the dataset.

    """
    
    if files is None:
        files = get_files(data_dir=data_dir, prefix=prefix, ignore='global', filetype='nc', verbose=False)
        if files is None:
            return None, None
    
    results = []
    attrs_list = []
    for file in files:
        tmp = xr.open_dataset(file).load()
        tmp.close()
        attrs = tmp.attrs 
        attrs_list.append(attrs)
        if variables == 'All':
            tmp.attrs = attrs
            try:
                tmp = tmp.to_dataset()
            except:
                pass
            results.append(tmp)
        else:
            tmp = tmp[variables]
            try:
                tmp = tmp.to_dataset()
            except:
                pass
            tmp.attrs = attrs
            results.append(tmp)
        

    xmax = np.max([results[ii]['x'].max() for ii in range(len(results))])
    xmin = np.min([results[ii]['x'].min() for ii in range(len(results))])
    ymax = np.max([results[ii]['y'].max() for ii in range(len(results))])
    ymin = np.min([results[ii]['y'].min() for ii in range(len(results))])
    x = np.linspace(xmin, xmax, num=nx) 
    y = np.linspace(ymin, ymax, num=ny)
    
    global_tmp = xr.Dataset()
    global_tmp.coords['x'] = x
    global_tmp.coords['y'] = y
    
    interp_results       = [results[ii].interp_like(global_tmp).fillna(fillna) for ii in range(len(results))]
    x_source             = [results[ii].attrs['Origin'][0]                     for ii in range(len(results))]
    y_source             = [results[ii].attrs['Origin'][1]                     for ii in range(len(results))]
    tags                 = [results[ii].attrs['source tag']                    for ii in range(len(results))]
    location_probability = [results[ii].attrs['location_probability']          for ii in range(len(results))]


    global_fields = xr.concat(interp_results, dim='source')
    global_fields = global_fields.assign({'location_probability': (('source'), location_probability)})
    global_fields = global_fields * global_fields['location_probability']       if weight_source   else global_fields
    global_fields = global_fields.sum(dim='source')                             if sum_source      else global_fields
        
    global_fields = global_fields.assign({'location_probability': (('source'), location_probability),
                                          'x_source':             (('source'), x_source),
                                          'y_source':             (('source'), y_source),
                                          'source tag':           (('source'), tags)})    
    global_attrs = concat_dict(attrs_list)
    global_fields.attrs = global_attrs
    return global_fields, global_attrs
    

def cumulative_probes(files=None, data_dir='Outdata/', prefix='probe', variables='C', weight_source=True):
    """
    This function can be used to sum together the contributions from each source for the probes.

    Parameters
    ----------
    files         : list
                    List of pathnames to files to read (Optional).
    data_dir      : str
                    String of directiory for where to read files.
    prefix        : str
                    String of prefix of files to read.
    variables     : str
                    What variable to sum together. Allowable: {'variable name' , 'All'}
    weight_source : bool
                    Whether or not to weight the results from each source w.r.t location_probability.

    Returns
    -------
    probes : xarray dataset
             An xarray dataset containing the sum of contributions from each source for the probes.
    attrs  : dict
             Dictionary containing the attributes.
    """
    
    if files is None:
        files = get_files(data_dir=data_dir, prefix=prefix, ignore='cumulative', filetype='nc', verbose=False)
        if files is None:
            return None, None
    
    results = []
    attrs_list = []
    for file in files:
        tmp=xr.open_dataset(file).load()
        tmp.close()
        results.append(tmp)
        attrs = tmp.attrs
        attrs_list.append(attrs)
        
    location_probabilities = [results[ii].attrs['location_probability'] for ii in range(len(results))]
    
    tmp = xr.concat(results,dim='source').fillna(0.0) 
    tmp = tmp.assign({'location probabilities': (('source'), location_probabilities)})
    
    if variables == 'All':
        probes = tmp
    else:
        probes = tmp[variables]
        
    if weight_source:
        probes = probes.weighted(tmp['location probabilities'])               # Weight the variable of interest w.r.t location probabilities
        probes = probes.sum(dim='source')                                     # Sum the variable over the sources
    else:
        probes = probes.sum(dim='source') 
    
    global_attrs = concat_dict(attrs_list)
    return probes, global_attrs


def reduce_point_tags(point_tags, glob_stats):
    """
    Reduces the amount of points in the point_tags.nc file that is stored in Outdata/ on termination of a simulation.
    The new set of points now lie on the gricells of the global grid. Points which previously lie in the same gridcell are summed together.

    Parameters
    ----------
    point_tags : xarray dataset
                 An xarray dataset containing the point_tags.    
    glob_stats : xarray dataset
                 An xarray dataset containing the global_statistics. This is used to find where the gridcells are.

    Returns
    -------
    out : xarray dataset
          An xarray dataset that has a reduced number of point tags.
    """

    def get_mode(array): # If one pixel on the global grid contains points from multiple sources, select the source tag with most points in the pixel.
        """ 
        If one gridcell on the global grid contains points from multiple sources, then select the source tag with the most points in the pixel. (Gets the mode of the gridcell.)

        Parameters
        ----------
        array : xarray dataset
                An xarray containing the source tags for a given gridcell.    
        """
        data_mode = mode(array.data)
        x_mean    = np.mean(array.x)
        y_mean    = np.mean(array.y)
        out = xr.DataArray(data = data_mode,
                        dims   = ['num'],
                        coords = {"x" : (['num'], [x_mean]),
                                  "y" : (['num'], [y_mean])},
                        name = 'source')
        return out
    
    xp = point_tags.x.data
    yp = point_tags.y.data
    xg = glob_stats.x.data
    yg = glob_stats.y.data

    XG, YG = np.meshgrid(xg, yg)

    Tree = KDTree(list(zip(XG.ravel(), YG.ravel())))
    indx = Tree.query(np.transpose(np.array([xp, yp])))

    unique_idx = np.unique(indx[1])
    x = [XG.ravel()[ii] for ii in unique_idx]                                                     
    y = [YG.ravel()[ii] for ii in unique_idx] 

    point_tags = point_tags.assign({'indx':('num', indx[1])})
    source_tag = point_tags.source.groupby(point_tags.indx).map(get_mode).data.astype(int)

    out = point_tags.location_probability.groupby(point_tags.indx).sum().to_dataset()
    out = out.assign_attrs(point_tags.attrs)
    out = out.assign_coords(x=('num', x), y=('num', y))
    out = out.drop('indx')
    out = out.assign({'location_probability':('num', out.location_probability.data),
                      'source':              ('num', source_tag)})
    return out


def post_process_output(config, weight_source=True, sum_source=False):
    """
    Post processes the data generated in Outdata\ from one simulation.
    Constructs global field from multiple local fields.
    Constructs global statistics from multiple local statistics.
    Cumulates the contributions of each source to each probe.
    Reduces the number of point tags to fit on the global grids.
    Stores files in Outdata\ if relevant.

    config file requires:
        
        * (str)     : config['paths']['outdata_path']
        * (str)     : config['paths']['statistics']
        * (str)     : config['paths']['fields']
        * (str)     : config['velocity']['velocity_file']
        * (bool)    : config['output']['global_fields']
        * (bool)    : config['output']['global_fields']
        * (int)     : config['output']['glob_nx']
        * (int)     : config['output']['glob_ny']

    Parameters
    ----------
    config        : setup.ini
                    Config file with relevant variables and setup configurations.    
    weight_source : bool 
                    Whether or not to weight the results from each source w.r.t location_probability.
    sum_source    : bool
                    Whether or not to sum over the results w.r.t each source.

    """
    start_time = TIME.time()
    
    nx = int(config['output']['glob_nx'])
    ny = int(config['output']['glob_ny'])

    ### CONSTRUCT GLOBAL FIELD FROM STATISTICS ###
    glob_stats = None
    if strtoBool(config['output']['global_stats']): 
        outpath = config['paths']['outdata_path']
        filenm  = config['paths']['statistics']+'global'   
        glob_stats, *_ = generate_global_fields(data_dir=outpath, prefix='statistics', variables='All', 
                                                nx=nx, ny=ny, weight_source=weight_source, sum_source=sum_source)
        store_xarray(dataset=glob_stats, outpath=outpath, filenm=filenm, verbose=True)
        
    
    ### CONSTRUCT GLOBAL FIELD FROM FIELDS ###
    if strtoBool(config['output']['global_fields']):
        outpath = config['paths']['outdata_path']
        filenm  = config['paths']['fields']+'global'
        glob_fields, *_ = generate_global_fields(data_dir=outpath, prefix='fields', variables='All',  
                                                 nx=nx, ny=ny, weight_source=weight_source, sum_source=sum_source)
        store_xarray(dataset=glob_fields, outpath=outpath, filenm=filenm, verbose=True)
        
    
    ### CUMULATE THE PROBE RESULTS FROM EACH LOCAL FIELD ###
    if strtoBool(config['output']['cumulate_probes']):
        outpath = config['paths']['outdata_path']
        filenm  = config['paths']['time_series_out']+'cumulative'
        cumul_probes, *_ = cumulative_probes(data_dir=outpath, prefix='probe', variables='All', weight_source=weight_source)
        store_xarray(dataset=cumul_probes, outpath=outpath, filenm=filenm, verbose=True)


    ### REDUCE POINT TAGS TO LIE ON GLOBAL GRIDCELLS ###
    if strtoBool(config['output']['reduce_point_tags']) and glob_stats is not None:
        infile  = config['paths']['outdata_path']+config['paths']['point_tags']+'.nc'
        try:
            point_tags = xr.open_dataset(infile)
            point_tags_reduced = reduce_point_tags(point_tags=point_tags, glob_stats=glob_stats)
            outpath = config['paths']['outdata_path']
            filenm  = config['paths']['point_tags']+'_reduced'
            store_xarray(dataset=point_tags_reduced, outpath=outpath, filenm=filenm, verbose=True)
        except:
            logging.warning(f'No file of the form {infile} found.')
            pass
    elif strtoBool(config['output']['reduce_point_tags']):
        logging.warning(f'No global results found to generate reduced point tags...')

    end_time = TIME.time()

    logging.info('All post-processing done. The total post-processing time was: {0:.2f}s\n'.format(end_time-start_time))


def store_stats_to_netcdf(config, stats_out, tag, meshdict=None, attrs={}):  
    """
    Sets attributes, sets units and stores the statistics to a netcdf file.

    Parameters
    ----------
    config    : config
                Config file containing directory paths.
    stats_out : xarray dataset
                An xarray dataset containing the statistics to be stored.
    tag       : int
                Integer denoting the source tag.
    attrs     : dict, optional
                Dictionary containing attributes, by default {}
    """ 
    # If outdata is unstructured, structurize them now.
    # This is only meant to be used if the user wants the advdiff module to return unstructured output by default. 
    # Then the from_unstructured_to_structured() function can revert the format back to structured.
    # This does require you to set return_unstructured to true in tracer_transport.
    #if meshdict is not None:
    #    stats_out = from_unstructured_to_structured(stats_out, structured_mesh=meshdict['mesh_out']) # If outdata is unstructured, structurize them now.
    
    attrs.update(stats_out.attrs)                            # Add any existing attributes to global attributes
    stats_out.attrs = attrs                                  # Set attributes        
    
    outpath = config['paths']['outdata_path']
    filenm  = config['paths']['statistics']+str(tag).zfill(2)
    store_xarray(dataset=stats_out, outpath=outpath, filenm=filenm)


def store_fields_to_netcdf(config, fields_out, tag, meshdict=None, attrs={}):
    """
    Sets attributes, sets units and stores the fields to a netcdf file.

    Parameters
    ----------
    config     : config
                 Config file containing directory paths.
    fields_out : xarray dataset
                 An xarray dataset containing the fields to be stored.
    tag        : int
                 Integer denoting the source tag.
    attrs      : dict, optional
                 Dictionary containing attributes, by default {}
    """ 
    # If outdata is unstructured, structurize them now.
    # This is only meant to be used if the user wants the advdiff module to return unstructured output by default. 
    # Then the from_unstructured_to_structured() function can revert the format back to structured.
    # This does require you to set return_unstructured to true in tracer_transport.
    #if meshdict is not None:
    #    fields_out = from_unstructured_to_structured(fields_out, structured_mesh=meshdict['mesh_out']) 
    
    attrs.update(fields_out.attrs)                                  # Add any existing attributes to global attributes
    fields_out.attrs = attrs                                        # Set attributes          
    
    downsample       = strtoBool(config['output']['downsample'])                  # When storing the fields, we may want to downsample the data to reduce memory.
    if downsample: 
        downsampling_factor = float(config['output']['downsampling_factor'])      # To what degree to downsample                                                               # Downsample the fields before saving
        skips      = int(np.ceil(np.sqrt(downsampling_factor)))                   # Skip length
        fields_out = fields_out.coarsen(x=skips, y=skips, boundary='trim').mean() # Take mean of every skip

    outpath = config['paths']['outdata_path']
    filenm  = config['paths']['fields']+str(tag).zfill(2)
    store_xarray(dataset=fields_out, outpath=outpath, filenm=filenm)


def store_probes_to_netcdf(config, probes_out, tag, attrs={}):
    """
    Sets attributes, sets units and stores the probes to a netcdf file.

    Parameters
    ----------
    config    : config
                Config file containing directory paths.
    stats_out : xarray dataset
                An xarray dataset containing the probes to be stored.
    tag       : int
                Integer denoting the source tag.
    attrs     : dict, optional
                Dictionary containing attributes, by default {}
    """ 
    
    attrs.update(probes_out.attrs)                            # Add any existing attributes to global attributes
    probes_out.attrs = attrs                                  # Set attributes            

    outpath = config['paths']['outdata_path']
    filenm  = config['paths']['time_series_out']+str(tag).zfill(2)
    store_xarray(dataset=probes_out, outpath=outpath, filenm=filenm)


def store_mass_to_netcdf(config, mass_out, tag, attrs={}):
    """
    Sets attributes, sets units and stores the mass to a netcdf file.

    Parameters
    ----------
    config    : config
                Config file containing directory paths.
    stats_out : xarray dataset
                An xarray dataset containing the mass to be stored.
    tag       : int
                Integer denoting the source tag.
    attrs     : dict, optional
                Dictionary containing attributes, by default {}
    """ 
    
    attrs.update(mass_out.attrs)                              # Add any existing attributes to global attributes
    mass_out.attrs = attrs                                    # Set attributes             
    
    outpath = config['paths']['outdata_path']
    filenm  = config['paths']['mass']+str(tag).zfill(2)
    store_xarray(dataset=mass_out, outpath=outpath, filenm=filenm)


def store_xarray(dataset, outpath, filenm, filetype='nc', verbose=False):
    """
    Stores an xarray dataset as a netcdf file.

    Parameters
    ----------
    dataset  : xarray dataset
               An xarray dataset to store.
    outpath  : str
               Directory path of where to store file.
    filenm   : str
               Filename.
    filetype : str, optional
               Filetype.
    """
    make_directory(dir_name=outpath, verbose=False)  # Ensure path exists
    outfile = outpath+filenm+'.'+filetype            # Set directory and filename of output
    try:
        dataset = set_units(dataset=dataset)         # Set units
        dataset.to_netcdf(outfile, mode='w')
        logging.info(f'Writing file: {outfile}') if verbose else None
    except:
        logging.warning(f'Could not generate file: {outfile}') if verbose else None