"""
This module includes all visualizer functions.

@author: Guttorm & Ketil
"""

import logging
import sys
import glob
import configparser
import fipy   as fp
import xarray as xr
import numpy  as np
import time   as TIME
import matplotlib
import matplotlib.pyplot    as plt
import matplotlib.gridspec  as gridspec
from pathlib                import Path
from tqdm                   import tqdm
from Tools.mesh_generator   import initialize_meshes, get_NSR
from Tools.directory_tools  import make_directory, get_files
from Tools.post_processing  import generate_global_fields
from Tools.visualizer_tools import plot_field, plot_quiver, get_vbounds, make_movie
from Tools.visualizer_tools import plot_global_velocity_location, plot_probe_locations, plot_source_locations
from Tools.parser_tools     import strtoBool
from Tools.file_loader      import load_sources, load_velocity
from Tools.coord_transform  import coord_converter as cc


#%%
def plot_riskmap(data_dir='Indata/Riskmaps/', prefix='risk_', variable='location_probability', 
    out_dir='Figures/', dpi=200, filetype='jpg'):
    """
    Plots every riskmap contained in the given directory. 
    The size of the scatterplot markers will be given by the variable of choice.

    Parameters
    ----------
    data_dir : str, optional
               Directory of where to look for files, by default 'Indata/Riskmaps/'
    prefix   : str, optional
               Prefix of what files to look for, by default 'risk_'
    variable : str, optional
               What variable to use for the size of the scatterplot markers, by default 'location_probability'
    out_dir  : str, optional
               Where to store the figures, by default 'Figures/'
    dpi      : int, optional
               DPI of the figures, by default 200
    filetype : str, optional
               Filetype of the figures, by default 'jpg'
    """
    logging.info(f'Plotting riskmaps from {data_dir}')
    plt.close("all")
    files = get_files(data_dir=data_dir, prefix=None, filetype='nc')
    if files == None:
        return

    for file in files:
        file_out = out_dir+prefix+Path(file).stem+'.'+filetype
        linebreak = '\n' if file == files[-1] else '' # If saving last file in list, create a linebreak
        logging.info(f'Saving: {file_out}{linebreak}')
        tmp=xr.open_dataset(file).load()
        plt.contour(tmp[variable])
        plt.savefig(file_out, dpi=dpi)
        plt.close("all")
    


def plot_sources(data_dir='Indata/Sources/', prefix='source_', variable='location_probability',
    out_dir='Figures/', dpi=200, filetype='jpg'):
    """
    Plot all sources found i the given directory.
    The size of the scatterplot markers will be given by the variable of choice.

    Parameters
    ----------
    data_dir : str, optional
               Directory of where to look for files, by default 'Indata/Riskmaps/'
    prefix   : str, optional
               Prefix of what files to look for, by default 'risk_'
    variable : str, optional
               What variable to use for the size of the scatterplot markers, by default 'location_probability'
    out_dir  : str, optional
               Where to store the figures, by default 'Figures/'
    dpi      : int, optional
               DPI of the figures, by default 200
    filetype : str, optional
               Filetype of the figures, by default 'jpg'
    """
    logging.info(f'Plotting sources from {data_dir}')
    plt.close("all")
    files = get_files(data_dir=data_dir, prefix=None, filetype='nc')
    if files == None:
        return
    
    for file in files:
        file_out = out_dir+prefix+Path(file).stem+'.'+filetype
        tmp=xr.open_dataset(file).load()
        if all(var in list(tmp.keys()) + list(tmp.coords) for var in ['x','y']): # Check if x,y coordinates available
            coordH = 'x'
            coordV = 'y'
        elif all(var in list(tmp.keys()) + list(tmp.coords) for var in ['lon','lat']): # Check if lon,lat coordinates available
            coordH = 'lon'
            coordV = 'lat'
        try:
            tmp.plot.scatter(x=coordH, y=coordV, s=tmp[variable].data)
        except:
            logging.warning('The file {} does not contain the variable {}... Only coordinates will be plotted'.format(file, variable))
            tmp.plot.scatter(x=coordH, y=coordV, s=5.0)
        linebreak = '\n' if file == files[-1] else '' # If saving last file in list, create a linebreak
        logging.info(f'Saving: {file_out}{linebreak}')
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0)) 
        plt.savefig(file_out, dpi=dpi)
        plt.close("all")



def plot_statistics(data_dir='Outdata/', prefix='statistics', variables=['max'],
                    weight_source=True, plot_source_loc=True, plot_probe_loc=True, plot_velocity_loc=False,
                    out_dir='Figures/', filetype='jpg', logscale=False, robust=True, levels=100, nx=250, ny=250, dpi=200, figsize=(6.4,4.8),
                    markersize=10, fontsize=10, config=None):
    """
    Plot the global statistics variables.

    Parameters
    ----------
    data_dir      : str, optional
                    Directory of where to look for statistics files, by default 'Outdata/'
    prefix        : str, optional
                    Prefix of files to load, by default 'statistics'
    variables     : list, optional
                    What variables to plot, by default ['max']
    weight_source : bool, optional
                    If True, weight the statistics by location_probability. If False do not weigh the statistics, by default True
    plot_source_loc   : bool, optional
                        If True, plot the source locations, by default True
    plot_probe_loc    : bool, optional
                        If True plot the probe locations, by default True
    plot_velocity_loc : bool, optional
                        If True plot the convex hull boundary of the velocity location, by default False
    out_dir  : str, optional
               Where to store the figures, by default 'Figures/'
    filetype : str, optional
               What filetype to store, by default 'jpg'
    logscale : bool, optional
               If True plot figures using a logscale, by default False
    robust : bool, optional
             If True clip the colorbar such that the plots become more visible, by default True
    levels : int, optional
             How many color-levels to assing to the contourf plot, by default 100
    nx : int, optional
         Number of discretizations along x-axis, by default 250
    ny : int, optional
         Number of discretizations along y-axis, by default 250
    dpi : int, optional
          DPI of figures, by default 200
    figsize : tuple, optional
              Figsize, by default (6.4,4.8)
    markersize : int, optional
        Markersizes of sources and probes, by default 10
    fontsize : int, optional
        Fontsizes of source and probe labels, by default 10
    config : setup.ini, optional
        Config file to load parameters from, by default None
    """
    
    logging.info(f'Plotting statistics from {data_dir}')
    plt.close("all")
    files = get_files(data_dir=data_dir, prefix=prefix, ignore='global', filetype='nc')
    if files is None:
        return
    statistics, attrs = generate_global_fields(files=files, variables=variables, nx=nx, ny=ny,
                                                weight_source=weight_source, sum_source=True)

    if plot_source_loc:
        tags                 = statistics['source tag'].data
        source_locations     = np.stack([statistics['x_source'].data, statistics['y_source'].data], axis=-1)
        location_probability = statistics['location_probability'].data
        plt_src_loc          = plot_source_locations(tags=tags, source_locations=source_locations, location_probability=location_probability,
                                                     markersize=markersize, fontsize=fontsize)
            
    if plot_probe_loc:
        plt_prb_loc = plot_probe_locations(data_dir=data_dir, markersize=markersize, fontsize=fontsize)

    if plot_velocity_loc:
        plt_glb_loc = plot_global_velocity_location(config=config, attrs=attrs)

    if variables == 'All':
        variables = list(statistics.keys())
        ignore_list = ['x_source','y_source','location_probability','source tag']
        for ignore in ignore_list:
            while ignore in variables:
                variables.remove(ignore)

    for variable in variables:
        try:
            field_variable = statistics[variable]

            if 'dt' not in field_variable.dims: 
                select_dim = [{}]
                label      = ['']
            else: # If theres a delta_t coordinate, loop over them and add label to file name.
                select_dim = [{'dt':ii}      for ii in field_variable['dt'].data]
                label      = [f'_{str(ii)}s' for ii in field_variable['dt'].data]

            for ii, dim in enumerate(select_dim):
                plt.figure(figsize=figsize, dpi=dpi)
                field = field_variable.sel(dim)
                plot_field(field, robust=robust, levels=levels, logscale=logscale, method='contourf')
                plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
                
                if plot_source_loc:
                    plt_src_loc.plot()

                if plot_probe_loc:
                    plt_prb_loc.plot()

                if plot_velocity_loc:
                    plt_glb_loc.plot()

                file_out = out_dir + prefix + '_' + variable + label[ii] + '.' + filetype
                linebreak = '\n' if (variable == variables[-1] and dim == select_dim[-1]) else '' # If saving last file in list, create a linebreak
                logging.info(f'Saving: {file_out}{linebreak}')
                plt.title(variable+label[ii])
                plt.savefig(file_out, dpi=dpi)
                plt.close('all')
        except:
            pass



def plot_probes(data_dir='Outdata/', prefix='probe', variable='C', weight_source=True, plot_contrib=False,
    out_dir='Figures/', dpi=200, filetype='jpg'):
    """
    Plot probe time-series data.

    Parameters
    ----------
    data_dir      : str, optional
                    Directory of where to look for probe files, by default 'Outdata/'
    prefix        : str, optional
                    Prefix of files to load, by default 'probe'
    variable      : str, optional
                    What variable to plot, by default 'C'
    weight_source : bool, optional
                    If True weigh the results from each source contribution by location_probability, by default True
    plot_contrib  : bool, optional
                    Plot each source contribution individually, by default False
    out_dir       : str, optional
                    Where to store figures, by default 'Figures/'
    dpi : int, optional
        DPI of figures, by default 200
    filetype : str, optional
        Filetype of figures, by default 'jpg'
    """
    logging.info(f'Plotting probes from {data_dir}')
    plt.close("all")
    files = get_files(data_dir=data_dir, prefix=prefix, ignore='cumulative', filetype='nc')
    if files == None:
        return

    results = []
    for file in files:
        tmp=xr.open_dataset(file).load()
        tmp.close()
        if weight_source:
            tmp = tmp * tmp.attrs['location_probability']
        if plot_contrib:
            file_out = out_dir+Path(file).stem+'_'+variable+'.'+filetype
            linebreak = '\n' if file == files[-1] else '' # If saving last file in list, create a linebreak
            logging.info(f'Saving: {file_out}{linebreak}')
            tmp[variable].plot.line(x='time',hue='num')
            plt.title(Path(file).stem)
            plt.savefig(file_out, dpi=dpi)
            plt.close("all")
        results.append(tmp)
    

    location_probabilities = [results[ii].attrs['location_probability'] for ii in range(len(results))]
    
    tmp = xr.concat(results,dim='source').fillna(0.0) 
    tmp = tmp.assign({'location probabilities': (('source'), location_probabilities)})
    
    probes = tmp[variable]
    probes = probes.sum(dim='source')
       
    file_out = out_dir+'probe_total'+'_'+variable+'.'+filetype
    linebreak = '\n' if file == files[-1] else '' # If saving last file in list, create a linebreak
    logging.info(f'Saving: {file_out}{linebreak}') 
    probes.plot.line(x='time',hue='num')
    plt.title('Probe Total')
    plt.savefig(file_out, dpi=dpi)
    plt.close("all")


      
def plot_fields(data_dir='Outdata/', prefix='fields', variable='C',
    out_dir='Figures/', filetype='jpg', robust=True, logscale=False, movie=True, levels=25, dpi=200, figsize=(6.4,4.8)):
    """
    Plot each local field time-series.

    Parameters
    ----------
    data_dir      : str, optional
                    Directory of where to look for field files, by default 'Outdata/'
    prefix        : str, optional
                    Prefix of files to load, by default 'fields'
    variable : str, optional
        What variable to plot, by default 'C'
    out_dir : str, optional
        Where to store figures, by default 'Figures/'
    filetype : str, optional
        Filetype of figures, by default 'jpg'
    robust : bool, optional
        If True clip colorbar values such that plot is easier to see, by default True
    logscale : bool, optional
        If True plot with a logscale, by default False
    movie : bool, optional
        If True make a movie, by default False
    levels : int, optional
        How many colorbar-levels to assign to the contorf plot, by default 25
    dpi : int, optional
        DPI of figures, by default 200
    figsize : tuple, optional
        Figsize, by default (6.4,4.8)
    """
    logging.info(f'Plotting global {variable}-fields from {data_dir}')
    plt.close("all")
    files = get_files(data_dir=data_dir, prefix=prefix, ignore='global', filetype='nc')
    if files == None:
        return

    for ii, file in enumerate(files):
        my_dir=out_dir+Path(file).stem.capitalize()
        make_directory(dir_name=my_dir)

        field = xr.open_dataset(file).load()
        field.close()
        field = field[variable]
        
        linebreak = '\n' if file == files[-1] else '' # If saving last file in list, create a linebreak
        logging.info('Generating plots for {}-{} around source {:02d}...{}'.format(str(variable), str(prefix), ii, linebreak))
        
        vmin, vmax = get_vbounds(field=field, robust=robust)

        times  = field.time.values
        for num, time in enumerate(tqdm(times, leave=False, miniters=4)):
            filenm = my_dir + '/'+variable+str(num).zfill(5)+'.'+filetype
            
            plt.figure(figsize=figsize)
            plot_field(field.sel(time=time), vmin=vmin, vmax=vmax, levels=levels, logscale=False, method='contourf')
            plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
            plt.savefig(filenm, dpi=dpi)
            plt.close("all")
        
        if movie:
            make_movie(dir_name=my_dir, name=variable, duration=0.1, filetype=filetype)



def plot_global_fields(data_dir='Outdata/', prefix='fields', variable='C',
                       weight_source=True, plot_source_loc=True, plot_probe_loc=True, robust=True, movie=True,
                       out_dir='Figures/', filetype='jpg', logscale=False, levels=100, nx=250, ny=250, dpi=100, figsize=(6.4,4.8),
                       markersize=10, fontsize=10):
    """
    Plots the global field time-series.

    Parameters
    ----------
    data_dir          : str, optional
                        Directory of where to look for field files, by default 'Outdata/'
    prefix            : str, optional
                        Prefix of files to load, by default 'fields'
    variable          : str, optional
                        What variable to plot, by default 'C'
    weight_source     : bool, optional
                        If True, weight the statistics by location_probability. If False do not weigh the statistics, by default True
    plot_source_loc   : bool, optional
                        If True, plot the source locations, by default True
    plot_probe_loc    : bool, optional
                        If True plot the probe locations, by default True
    out_dir           : str, optional
                        Where to store figures, by default 'Figures/'
    filetype          : str, optional
                        Filetype of figures, by default 'jpg'
    robust            : bool, optional
                        If True clip colorbar values such that plot is easier to see, by default True
    logscale          : bool, optional
                        If True plot with a logscale, by default False
    movie             : bool, optional
                        If True make a movie, by default False
    levels            : int, optional
                        How many colorbar-levels to assign to the contorf plot, by default 25
    dpi               : int, optional
                        DPI of figures, by default 200
    figsize           : tuple, optional
                        Figsize, by default (6.4,4.8)
    nx                : int, optional
                        Number of discretizations along x-axis, by default 250
    ny                : int, optional
                        Number of discretizations along y-axis, by default 250
    markersize        : int, optional
                        Markersizes of sources and probes, by default 10
    fontsize          : int, optional
                        Fontsizes of source and probe labels, by default 10
    """
    logging.info(f'Plotting global fields from {data_dir}')
    my_dir = out_dir+Path('global_field').stem.capitalize()
    make_directory(dir_name=my_dir)
    
    plt.close("all") 
    files = get_files(data_dir=data_dir, prefix=prefix, ignore='global', filetype='nc')
    if files is None:
        return
    field, *_ = generate_global_fields(files=files, variables=variable, nx=nx, ny=ny,
                                          weight_source=weight_source, sum_source=True)
    
    if plot_source_loc:    
        tags                 = field['source tag'].data
        source_locations     = np.stack([field['x_source'].data, field['y_source'].data], axis=-1)
        location_probability = field['location_probability'].data
        plt_src_loc          = plot_source_locations(tags=tags, source_locations=source_locations, location_probability=location_probability,
                                                     markersize=markersize, fontsize=fontsize)

    if plot_probe_loc:
        plt_prb_loc = plot_probe_locations(data_dir=data_dir, markersize=markersize, fontsize=fontsize)

    vmin, vmax = get_vbounds(field=field[variable], robust=robust)

    times  = field.time.values
    for num, time in enumerate(tqdm(times, leave=False, miniters=4)):
        filenm = my_dir + '/' + variable + str(num).zfill(5) + '.' + filetype
        
        plt.figure(figsize=figsize, dpi=dpi)
        plot_field(field[variable].sel(time=time), vmin=vmin, vmax=vmax, levels=levels, logscale=False, method='contourf')
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

        if plot_source_loc:
            plt_src_loc.plot()
        
        if plot_probe_loc:
            plt_prb_loc.plot()
            
        plt.title(time)
        plt.savefig(filenm, dpi=dpi)
        plt.close("all")
        
    
    if movie: 
        make_movie(dir_name=my_dir, name=variable, duration=0.1, filetype=filetype)



def plot_stream(data_dir='Outdata/',prefix='fields', 
    out_dir='Figures/', filetype='jpg', colorvariable='velocity', skip=5, logscale=False, movie=True, dpi=200, figsize=(6.4,4.8)):
    """
    Plot the local stream fields.

    Parameters
    ----------
    data_dir          : str, optional
                        Directory of where to look for field files, by default 'Outdata/'
    prefix            : str, optional
                        Prefix of files to load, by default 'fields'
    out_dir           : str, optional
                        Where to store figures, by default 'Figures/'
    filetype          : str, optional
                        Filetype of figures, by default 'jpg'
    colorvariable     : str, optional
                        What variable to base the colourmap on, by default 'velocity'
    skip              : int, optional
                        How many points to skip in x,y for each vector, by default 5
    logscale          : bool, optional
                        If True plot with a logscale, by default False
    movie             : bool, optional
                        If True make a movie, by default False
    dpi               : int, optional
                        DPI of figures, by default 200
    figsize           : tuple, optional
                        Figsize, by default (6.4,4.8)
    """
    logging.info(f'Plotting streams from {data_dir}')
    plt.close("all")
    files = get_files(data_dir=data_dir, prefix=prefix, ignore='global', filetype='nc')
    if files == None:
        return

    for ii, file in enumerate(files):
        my_dir = out_dir + Path(file).stem.capitalize() + 'stream'
        make_directory(dir_name=my_dir)

        field = xr.open_dataset(file).load()
        field.close()
        field = field.transpose('time', 'y', 'x').coarsen(x=skip, y=skip, boundary='trim').mean()
       
        variable='stream'
        linebreak = '\n' if file == files[-1] else '' # If saving last file in list, create a linebreak
        logging.info('Generating plots for {}-{} around source {:02d}...{}'.format(variable, str(prefix), ii, linebreak))

        times   = field.time.values
        for num, time in enumerate(tqdm(times, leave=False, miniters=4)):
            filenm = my_dir + '/'+variable+str(num).zfill(5) + '.' + filetype
            
            plt.figure(figsize=figsize)
            plot_quiver(field.sel(time=time), colorvariable=colorvariable, logscale=logscale)
            plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
            plt.savefig(filenm, dpi=dpi)
            plt.close("all")
        
        if movie:
            make_movie(dir_name=my_dir, name=variable, duration=0.1, filetype=filetype)



def plot_global_stream(data_dir='Outdata/', prefix='fields', 
                       plot_source_loc=True, plot_probe_loc=True, logscale=False, movie=True,
                       out_dir='Figures/', filetype='jpg', nx=100, ny=100, colorvariable='velocity', skip=5, dpi=200, figsize=(6.4,4.8),
                       markersize=10, fontsize=10):
    """
    Plot the global stream field.

    Parameters
    ----------
    data_dir          : str, optional
                        Directory of where to look for field files, by default 'Outdata/'
    prefix            : str, optional
                        Prefix of files to load, by default 'fields'
    plot_source_loc   : bool, optional
                        If True, plot the source locations, by default True
    plot_probe_loc    : bool, optional
                        If True plot the probe locations, by default True
    logscale          : bool, optional
                        If True plot with a logscale, by default False
     movie             : bool, optional
                        If True make a movie, by default False
    out_dir           : str, optional
                        Where to store figures, by default 'Figures/'
    filetype          : str, optional
                        Filetype of figures, by default 'jpg'
    nx                : int, optional
                        Number of discretizations along x-axis, by default 100
    ny                : int, optional
                        Number of discretizations along y-axis, by default 100
    colorvariable     : str, optional
                        What variable to base the colourmap on, by default 'velocity'
    skip              : int, optional
                        How many points to skip in x,y for each vector, by default 5
    dpi               : int, optional
                        DPI of figures, by default 200
    figsize           : tuple, optional
                        Figsize, by default (6.4,4.8)
    markersize        : int, optional
                        Markersizes of sources and probes, by default 10
    fontsize          : int, optional
                        Fontsizes of source and probe labels, by default 10
    """
    logging.info(f'Plotting global streams from {data_dir} \n')
    my_dir = out_dir + Path('global_stream').stem.capitalize()
    make_directory(dir_name=my_dir)
    
    plt.close("all")
    files = get_files(data_dir=data_dir, prefix=prefix, ignore='global', filetype='nc')
    if files is None:
        return
    field, *_ = generate_global_fields(files=files, variables='All', nx=nx, ny=ny,
                                       weight_source=False, sum_source=False, fillna=np.nan)

    field_C  = field['C'].sum(dim='source', skipna=True)
    field_uv = field[{'u','v'}].mean(dim='source', skipna=True)
    
    if plot_source_loc:    
        tags                 = field['source tag'].data
        source_locations     = np.stack([field['x_source'].data, field['y_source'].data], axis=-1)
        location_probability = field['location_probability'].data
        plt_src_loc          = plot_source_locations(tags=tags, source_locations=source_locations, location_probability=location_probability,
                                                     markersize=markersize, fontsize=fontsize, color='black')

    if plot_probe_loc:
        plt_prb_loc = plot_probe_locations(data_dir=data_dir, markersize=markersize, fontsize=fontsize, color='black')
    
    field = field.sum(dim='source', skipna=True)
    field = field.assign(C = field_C)
    field = field.assign(u = field_uv.u)
    field = field.assign(v = field_uv.v)
    field = field.transpose('time', 'y', 'x').coarsen(x=skip, y=skip, boundary='trim').mean()

    times   = field.time.values
    for num, time in enumerate(tqdm(times, leave=False, miniters=4)):
        filenm=my_dir + '/'+'stream'+str(num).zfill(5)+'.'+filetype
        
        plt.figure(figsize=figsize)
        plot_quiver(field.sel(time=time), colorvariable=colorvariable, logscale=logscale)
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0)) 

        if plot_source_loc:
            plt_src_loc.plot()
        
        if plot_probe_loc:
            plt_prb_loc.plot()
            
        plt.title(time)
        plt.savefig(filenm, dpi=dpi)
        plt.close("all")
        
    if movie:
        make_movie(dir_name=my_dir, name='stream', duration=0.1, filetype=filetype)



def plot_mass(data_dir='Outdata/', prefix='mass', weight_source=True, plot_contrib=False,
    out_dir='Figures/', dpi=200, filetype='jpg'):
    """
    Plot the total mass of the system wrt time.

    Parameters
    ----------
    data_dir      : str, optional
                    Directory of where to look for mass files, by default 'Outdata/'
    prefix        : str, optional
                    Prefix of files to load, by default 'mass'
    weight_source : bool, optional
                    If True weigh the results from each source contribution by location_probability, by default True
    plot_contrib  : bool, optional
                    Plot each source contribution individually, by default False
    out_dir       : str, optional
                    Where to store figures, by default 'Figures/'
    dpi : int, optional
        DPI of figures, by default 200
    filetype : str, optional
        Filetype of figures, by default 'jpg'
    """
    logging.info(f'Plotting mass from {data_dir}')
    plt.close("all")
    files = get_files(data_dir=data_dir, prefix=prefix, ignore=None, filetype='nc')
    if files == None:
        return
    
    tmp_Tot_M  = None
    tmp_Tot_MT = None
    for ii, file in enumerate(files):
        tmp=xr.open_dataset(file).load()
        tmp.close()
        if weight_source:
            tmp = tmp * tmp.attrs['location_probability']
        if plot_contrib:
            file_out = out_dir + Path(file).stem + '.' + filetype
            linebreak = '\n' if file == files[-1] else '' # If saving last file in list, create a linebreak
            logging.info(f'Saving: {file_out}{linebreak}')
            tmp['Mass'].plot()
            tmp['Mass Theoretical'].plot(linestyle='dashed')
            plt.legend(['Mass', 'Mass Theoretical'])
            plt.ylabel('Mass')
            plt.title(Path(file).stem)
            plt.savefig(file_out, dpi=dpi)
            plt.close("all")
        
        if ii == 0:
            tmp_Tot_M  = tmp['Mass']
            tmp_Tot_MT = tmp['Mass Theoretical']
        else:    
            tmp_Tot_M  += tmp['Mass']
            tmp_Tot_MT += tmp['Mass Theoretical']
       
    tmp_Tot_M.plot()
    tmp_Tot_MT.plot(linestyle='dashed')
    file_out = out_dir+'mass_total'+'.'+filetype
    logging.info(f'Saving: {file_out}\n')
    plt.legend(['Mass', 'Mass Theoretical'])
    plt.ylabel('Mass')
    plt.title('Mass_total')
    plt.savefig(file_out, dpi=dpi)
    plt.close("all")
    


def plot_meshes(data_dir='Indata/', out_dir='Figures/', dpi=100, filetype='jpg'):
    """
    Plot meshes and mesh quality.

    Parameters
    ----------
    data_dir      : str, optional
                    Directory of where to look for config file, by default 'Indata/'
    out_dir       : str, optional
                    Where to store figures, by default 'Figures/'
    dpi : int, optional
        DPI of figures, by default 200
    filetype : str, optional
        Filetype of figures, by default 'jpg'
    """
    print()
    logging.info(f'Plotting meshes')
    plt.close("all")
    config=configparser.ConfigParser(allow_no_value=True)
    config.read(data_dir+'setup.ini')
    
    logging.info('Generating meshes...\n')
    names    = ['mesh', 'mesh_out']
    meshdicts = initialize_meshes(grid_type       = str(config['grid']['type']), 
                                 grid_type_out    = str(config['grid']['type_out']),
                                 Lx               = float(config['grid']['Lx']),
                                 Ly               = float(config['grid']['Ly']),
                                 nx               = int(config['grid']['nx']),
                                 ny               = int(config['grid']['ny']),
                                 minvol           = float(config['grid']['minvol']),
                                 maxvol           = float(config['grid']['maxvol']),
                                 power            = float(config['grid']['power']),
                                 reduction_factor = float(config['grid']['reduction_factor']),
                                 verbose          = False)
    fig = plt.figure(figsize=(16,18), dpi=dpi)
    gs  = gridspec.GridSpec(3, 2)
    for name in names:
        meshdict = meshdicts.copy()
        ax1 = fig.add_subplot(gs[0, :]) 
        ax2 = fig.add_subplot(gs[1, :])
        ax3 = fig.add_subplot(gs[2, 0]) 
        ax4 = fig.add_subplot(gs[2, 1]) 
        
        cmap_1 = matplotlib.cm.get_cmap(name='jet')
        cmap_2 = matplotlib.cm.get_cmap(name='RdYlGn')
    
        rnd_val  = np.random.random(meshdict[name].cellVolumes.shape[0])
        NSR_val  = get_NSR(meshdict[name])
        
        mesh_rnd = fp.CellVariable(name=name,  mesh=meshdict[name], value=rnd_val)
        mesh_NSR = fp.CellVariable(name='NSR', mesh=meshdict[name], value=NSR_val)
        
        ax1.hist(meshdict[name].cellVolumes.data, 25, rwidth=0.85)
        ax1.set_xlabel('Volume')
        ax1.set_ylabel('Count')
        ax1.title.set_text('Cell Volumes')
        ax2.hist(NSR_val, 25, density=False, stacked=False, rwidth=0.85) 
        ax2.set_xlabel('NSR')
        ax2.set_ylabel('Count')
        ax2.title.set_text('Normalized Shape Ratio') 
        fp.Matplotlib2DViewer(vars=(mesh_rnd), datamin=0.0, datamax=1.0, cmap=cmap_1, colorbar=False,       axes=ax3, figaspect='auto')
        fp.Matplotlib2DViewer(vars=(mesh_NSR), datamin=0.0, datamax=1.0, cmap=cmap_2, colorbar='vertical',  axes=ax4, figaspect='auto')
        
        file_out = out_dir+name+'.'+filetype
        plt.savefig(file_out, dpi=dpi)
        plt.clf()
    plt.close('all')
    


def plot_local_grids(data_dir='Indata/', out_dir='Figures/', variable='quiver', time=-1, dpi=200, figsize=(6.4,4.8), filetype='jpg', 
                     markersize=10, fontsize=10, config=None):
    """
    Plot the locations of the local grids superimposed on the velocity locations.

    Parameters
    ----------
    data_dir   : str, optional
                 Where to look for a config file, by default 'Indata/'
    out_dir    : str, optional
                 Where to store plots, by default 'Figures/'
    variable   : str, optional
                 What type of velocity variable to plot, by default 'quiver'
                 Allowable:  'u', 'v', 'quiver'.
    time       : int, optional
                 What time to plot for velocity data, by default -1
    dpi        : int, optional
                 DPI for plots, by default 200
    figsize    : tuple, optional
                 Figure size, by default (6.4,4.8)
    filetype   : str, optional
                 Filetype for output, by default 'jpg'
    markersize : int, optional
                 Markersize of source locations, by default 10
    fontsize   : int, optional
                 Fontsize of source labels, by default 10
    config     : setup.ini, optional
                 Config file to load. If None, load from data_dir, by default None
    """
    logging.info(f'Plotting local grids')
    logging.info('Reading sources and velocity')
    if config is None:
        config=configparser.ConfigParser(allow_no_value=True)
        config.read(data_dir+'setup.ini')

    Lx = float(config['grid']['Lx'])
    Ly = float(config['grid']['Ly'])

    tmpS, *_ = load_sources(config=config)
    tags = tmpS.source.data

    location_probability = tmpS.location_probability.data
    source_locations     = np.stack([tmpS.x.data, tmpS.y.data], axis=-1)
    plt_glb_loc = plot_global_velocity_location(config=config)
    plt_src_loc = plot_source_locations(tags=tags, source_locations=source_locations, location_probability=location_probability,
                                        Lxy=[Lx,Ly], markersize=markersize, fontsize=fontsize, linewidth=1)
    
    tmpv = plt_glb_loc.global_velocity

    plt.figure(figsize=figsize, dpi=dpi)
    if variable == 'quiver':
        skip = 4
        tmp  = tmpv.transpose('time', 'y', 'x').coarsen(x=skip, y=skip, boundary='trim').mean()
        spd  = np.sqrt(tmp.u.values**2+tmp.v.values**2)
        plt.quiver(tmp.x.values,tmp.y.values,tmp.u.values,tmp.v.values, spd)
    elif variable == 'u':
        plot_field(tmpv['u'].isel(time=time), robust=False, method=None)
    elif variable == 'v':
        plot_field(tmpv['v'].isel(time=time), robust=False, method=None)
    plt_glb_loc.plot()
    plt_src_loc.plot()
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

    file_out = out_dir+'sources_local_grids_velocity'+'.'+filetype
    logging.info(f'Saving: {file_out}\n')
    plt.savefig(file_out, dpi=dpi)
    plt.close('all')



def plot_map(config, out_dir='Figures/', file_nm='Map', filetype='jpg'):
    """
    Plots a map containing the source locations and velocity locations in space with castlines shown.

    REQUIRES cartopy package.

    Parameters
    ----------
    config   : setup.ini
               config file. Used to load source and velocity data.
    out_dir  : str, optional
               Where to store figures, by default 'Figures/'
    file_nm  : str, optional
               What to name the output, by default 'Map'
    filetype : str, optional
               What filetype to use, by default 'jpg'
    """
    logging.info(f'Plotting map')
    try:
        import cartopy.crs as ccrs
    except:
        logging.warning('Cannot import cartopy. Cancelling map plot...\n')
        return
    fig = plt.figure(figsize=(16,9), dpi=600)
    ax  = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    velocity, *_ = load_velocity(config, verbose=False)
    sources,  *_ = load_sources(config)
    coord_conv = cc()
    lonmin, lonmax, latmin, latmax = coord_conv.bounds
    if lonmin == -180.0 and lonmax == 180.0 and latmin == -90.0 and latmax == 90.0:
        logging.warning('Cancelling map plot...\n')
        return 
    lonmin = lonmin - 0.5 if lonmin > -180.0 else lonmin # Add margin
    lonmax = lonmax + 0.5 if lonmax <  180.0 else lonmax
    latmin = latmin - 0.5 if latmin > -90.0 else latmin
    latmax = latmax + 0.5 if latmax <  90.0 else latmax

    ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())

    velocity = coord_conv.assign_xy_to_lonlat(velocity)
    sources  = coord_conv.assign_xy_to_lonlat(sources)
    
    ax.quiver(velocity.lon.data,             velocity.lat.data,
              velocity.isel(time=-1).u.data, velocity.isel(time=-1).v.data,  
              np.sqrt(velocity.isel(time=-1).u.data**2+velocity.isel(time=-1).v.data**2))
    ax.scatter(x=sources.lon.data, y=sources.lat.data, c='red')

    ax.coastlines()
    ax.gridlines(draw_labels=True)

    ax.annotate('Sources', xy=(sources.lon.mean().data, sources.lat.mean().data),  xycoords='data',
                xytext=(-10., -90.), textcoords='offset points',
                arrowprops=dict(facecolor='black', shrink=0.08),
                horizontalalignment='right', verticalalignment='top',
                )
    ax.annotate('Velocity', xy=(velocity.lon.mean().data, velocity.lat.mean()),  xycoords='data',
                xytext=(190., -60.), textcoords='offset points',
                arrowprops=dict(facecolor='black', shrink=0.08),
                horizontalalignment='right', verticalalignment='top',
                )
    file_out = out_dir+file_nm+'.'+filetype
    logging.info(f'Saving: {file_out}\n')
    plt.savefig(file_out) 
    plt.close('all')


#%%
def visualizer(config=None):
    """
    Master visualizer function. If config file provided, use the information therein to plot.
    If no config file provided, load one from the outdata path.

    Parameters
    ----------
    config : setup.ini, optional
             Config containing relevant visualizer information, by default None.
    """
    masterconfig=configparser.ConfigParser(allow_no_value=True)
    masterconfig.optionxform = str  # preserve case for letters
    masterconfig.read('Indata/setup.ini')

    levels = {'DEBUG': logging.DEBUG, 'INFO': logging.INFO, 'WARNING': logging.WARNING, 'CRITICAL': logging.CRITICAL}
    logging.basicConfig(format   = '[%(levelname)s] %(message)s', 
                        level    = levels[masterconfig['DEFAULT']['verbose'].upper()], 
                        handlers = [logging.FileHandler(masterconfig['paths']['outdata_path']+'debug.log'), logging.StreamHandler(sys.stdout)])
    
    if config is None:
        config=configparser.ConfigParser(allow_no_value=True)
        config.optionxform = str  # preserve case for letters
        config.read(masterconfig['paths']['outdata_path']+'setup.ini')

    start_time = TIME.time()
    
    nx = int(masterconfig['output']['glob_nx'])
    ny = int(masterconfig['output']['glob_ny'])

    levels          = int(masterconfig['visualizer']['levels'])
    robust          = strtoBool(masterconfig['visualizer']['robust'])
    weight_source   = strtoBool(masterconfig['visualizer']['weight_source'])
    pltSourceLoc    = strtoBool(masterconfig['visualizer']['plot_source_loc'])
    pltProbeLoc     = strtoBool(masterconfig['visualizer']['plot_probe_loc'])
    logscale        = strtoBool(masterconfig['visualizer']['logscale'])
    fontsize        = float(masterconfig['visualizer']['fontsize'])
    markersize      = float(masterconfig['visualizer']['markersize'])
    filetype        = str(masterconfig['visualizer']['filetype'])

    pltRiskmap      = strtoBool(masterconfig['visualizer']['plot_riskmaps'])
    pltSources      = strtoBool(masterconfig['visualizer']['plot_sources'])
    pltStatistics   = strtoBool(masterconfig['visualizer']['plot_stats'])
    pltGlobalFields = strtoBool(masterconfig['visualizer']['plot_global_fields'])
    pltLocalFields  = strtoBool(masterconfig['visualizer']['plot_local_fields'])
    pltProbes       = strtoBool(masterconfig['visualizer']['plot_probes'])
    pltProbeContrib = strtoBool(masterconfig['visualizer']['plot_probe_contrib'])
    pltMass         = strtoBool(masterconfig['visualizer']['plot_mass'])
    pltMassContrib  = strtoBool(masterconfig['visualizer']['plot_mass_contrib'])
    pltLocalGrids   = strtoBool(masterconfig['visualizer']['plot_local_grids'])
    pltMeshes       = strtoBool(masterconfig['visualizer']['plot_meshes'])
    pltMap          = strtoBool(masterconfig['visualizer']['plot_map'])
    pltLocalStream  = strtoBool(masterconfig['visualizer']['plot_local_stream'])
    pltGlobalStream = strtoBool(masterconfig['visualizer']['plot_global_stream'])
    
    variables   = 'All'
    figures_dir = config['paths']['figures_path']
    data_dir    = config['paths']['outdata_path']
    make_directory(dir_name=figures_dir, verbose=True) # Ensure figures path exists.

    plot_riskmap(out_dir=figures_dir, filetype=filetype)                                                                            if pltRiskmap    else None
    plot_sources(out_dir=figures_dir, filetype=filetype)                                                                            if pltSources    else None 
    plot_probes(out_dir=figures_dir, data_dir=data_dir, plot_contrib=pltProbeContrib, filetype=filetype)                            if pltProbes     else None
    plot_mass(out_dir=figures_dir, data_dir=data_dir, plot_contrib=pltMassContrib, filetype=filetype)                               if pltMass       else None    
    plot_map(out_dir=figures_dir, config=config, filetype=filetype)                                                                 if pltMap        else None    
    plot_local_grids(out_dir=figures_dir, variable='u', markersize=markersize, fontsize=fontsize, filetype=filetype, config=config) if pltLocalGrids else None   

    if pltStatistics:
        plot_statistics(out_dir=figures_dir, data_dir=data_dir, prefix='statistics', variables=variables, logscale=logscale,
                        weight_source=weight_source, plot_source_loc=pltSourceLoc, plot_probe_loc=pltProbeLoc, plot_velocity_loc=False,
                        robust=robust, levels=levels, nx=nx, ny=ny, markersize=markersize, fontsize=fontsize, filetype=filetype, config=config)
    
    if pltGlobalFields:
        plot_global_fields(out_dir=figures_dir, data_dir=data_dir, weight_source=weight_source, plot_source_loc=pltSourceLoc, plot_probe_loc=pltProbeLoc,
                           robust=robust, levels=levels, nx=nx, ny=ny, dpi=200, markersize=markersize, fontsize=fontsize, filetype=filetype)
          
    if pltLocalFields:
        plot_fields(out_dir=figures_dir, data_dir=data_dir, prefix='fields', variable='C',
                    levels=levels, dpi=200, figsize=(9,6.7), filetype=filetype)
    
    if pltLocalStream:
        plot_stream(out_dir=figures_dir, data_dir=data_dir, prefix='fields', colorvariable='C', 
                    logscale=logscale, filetype=filetype, dpi=200, skip=3)    
    
    if pltGlobalStream:
        plot_global_stream(out_dir=figures_dir, data_dir=data_dir, prefix='fields', colorvariable='C', 
                           logscale=logscale, dpi=200, filetype=filetype, skip=3)

    plot_meshes(out_dir=figures_dir, filetype=filetype) if pltMeshes else None # Keep this last. The FiPy plotters like to interfere with the other plotters for some reason...
    
    plt.close('all')

    logging.info('All plots done. The total plotting time was: {0:.2f}s\n'.format(TIME.time()-start_time))
    



if __name__ == '__main__':

    masterconfig=configparser.ConfigParser(allow_no_value=True)
    masterconfig.optionxform = str  # preserve case for letters
    masterconfig.read('Indata/setup.ini')

    suffix_list = [string.replace(masterconfig['paths']['outdata_path'].strip('/'),'').strip('\\')+str('/') for string in glob.glob(masterconfig['paths']['outdata_path']+'*/')]
    if len(suffix_list) > 0:
        print(f'Available outdata: ')
        [print(suffix) for suffix in suffix_list]
        print('Enter "All" to visualize all available outdata.')
        print('Enter nothing to visualize default outdata folder.')
        suffix = input('What outdata to visualize?: ')
        print()
    else:
        suffix = ''
    
    if suffix != 'All':
        suffix_list = [suffix]
   
    for suffix in suffix_list:
        config=configparser.ConfigParser(allow_no_value=True)
        config.optionxform = str  # preserve case for letters
        config.read(masterconfig['paths']['outdata_path']+suffix+'setup.ini')
        visualizer(config)