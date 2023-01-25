"""
This module includes all visualizer tools used in the visualizer module.

@author: Ketil
"""
import glob
import logging
import imageio
import configparser
import xarray              as xr
import numpy               as np
import matplotlib.pyplot   as plt
import matplotlib.patches  as patches
from matplotlib.colors     import LogNorm, Normalize
from scipy.spatial         import ConvexHull, Delaunay
from Tools.directory_tools import get_files
from Tools.file_loader     import load_sources, load_velocity
from Tools.parser_tools    import strtoBool


def plot_field(field, xcoord='x', ycoord='y', robust=False, vmin=None, vmax=None, levels=25, logscale=False, method=None, ax=None):

    if logscale and field.mean().data > 0.0:
        eps    = xr.where(field > 0.0, field, np.nan).quantile(0.1).data
        field  = field + eps
        norm   = LogNorm()
        method = None
    else:
        norm = Normalize()

    if field.mean().data > 0.0 and vmin is None and vmax is None:
        vmin, vmax = get_vbounds(field, robust=robust, min_quantile=0.00, max_quantile=0.98)
    elif vmin is None and vmax is None:
        vmin, vmax = get_vbounds(field, robust=robust, min_quantile=0.02, max_quantile=1.00)
    
    if method == 'contourf':
        field.plot.contourf(x=xcoord, y=ycoord, levels=levels, vmin=vmin, vmax=vmax, norm=norm, ax=ax)
    elif method == 'contour':
        field.plot.contour(x=xcoord, y=ycoord, levels=levels, vmin=vmin, vmax=vmax, norm=norm, ax=ax)
    elif method is None:
        field.plot(x=xcoord, y=ycoord, levels=levels, vmin=vmin, vmax=vmax, norm=norm, ax=ax)


def plot_quiver(field, xcoord='x', ycoord='y', uvar='u', vvar='v', colorvariable=None, logscale=False, ax=None):

    if colorvariable is None:
        spd = None
    elif colorvariable == 'velocity':
        spd = np.sqrt(field[uvar].values**2+field[vvar].values**2)
    else:
        spd = field[colorvariable]

    if logscale and spd is not None:
        spd = np.log10(spd+1e-8)

    if ax is None:
        plt.quiver(field[xcoord].values, field[ycoord].values, field[uvar].values, field[vvar].values, spd)
    else:
        ax.quiver(field[xcoord].values, field[ycoord].values, field[uvar].values, field[vvar].values, spd)
    


def get_vbounds(field, robust=True, min_quantile=0.00, max_quantile=0.98):
    """
    Get colorbar bounds given a field of data and the minimum and maximum quantile to cap the bounds at.

    Parameters
    ----------
    field  : xarray dataset
             An xarray dataset containing some data.          
    robust : bool, optional
             Wether or not to use the given quantiles, by default True
    min_quantile : float, optional
                   Quantile of data to set the min bound at, by default 0.00
    max_quantile : float, optional
                   Quantile of data to set the max bound at, by default 0.98

    Returns
    -------
    tuple
          (vmin, vmax)
    """
    if robust == True:
        vmin = field.quantile(min_quantile)
        vmax = field.quantile(max_quantile)
    else:
        vmin = field.quantile(0.00)
        vmax = field.quantile(1.00)
    if vmin > vmax:
        vmin, vmax = vmax, vmin
    if vmin == vmax: # Ensure contour bounds are always increasing. Else we get an error
        vmax = vmin + 1e-15
    return vmin, vmax


def make_movie(dir_name, name, duration, filetype):
    """
    Make a movie (.gif) out of the images in a certain directory.

    Parameters
    ----------
    dir_name : str
               Directory name. Where to look for images.
    name     : str
               Name of images to find.
    duration : float
               Duration of each frame on screen.
    filetype : str
               Filetype of images to find.
    """
    images = []
    for filename in sorted(glob.glob(dir_name + '/' + name + '*' + filetype)):
            images.append(imageio.imread(filename))
    imageio.mimsave(dir_name + '/' + name + 'movie.gif', images, duration=duration)


class plot_source_locations:
    
    def __init__(self, tags, source_locations, location_probability, Lxy=[None,None], color='white', fontsize=10, markersize=10, linewidth=2):
        """
        Plot source locations.

        Parameters
        ----------
        tags                 : array ([n_sources])
                               Array containing the tags for every source.
        source_locations     : array ([n_sources,2])
                               Array (x,y) containing the source locations for each source.
        location_probability : array ([n_sources])
                               Array containing the respective location probabilities.
        Lxy        : list, optional
                     List containing sidelengths of the local grids [Lx,Ly], by default [None,None]. If Lxy==[None,None], then do not plot local grid boxes.
        color      : str, optional
                     Color of the source marker, by default 'white'.
        fontsize   : float or int
                     Fontsize of source marker label.
        markersize : float or int
                     Markersize of source marker.
        linewidth  : float or int
                     Linewidth of source marker.
        """
        self.x_source = source_locations[:,0]
        self.y_source = source_locations[:,1]
        if np.all(np.array(Lxy) != None):
            self.Lx = Lxy[0]
            self.Ly = Lxy[1]
            self.plot_grids = True
        else:
            self.plot_grids = False
        
        self.location_probability = location_probability
        self.tags       = tags
        self.color      = color
        self.fontsize   = fontsize
        self.markersize = markersize
        self.linewidth  = linewidth
    
    def plot(self):
        """
        Plot the sources.
        """
        plt.scatter(x=self.x_source, y=self.y_source, s=self.location_probability*self.markersize/10, marker='o', edgecolors=self.color, facecolors='none')
        for i, txt in enumerate(self.tags):
            plt.annotate(txt, (self.x_source[i], self.y_source[i]), xytext=(2,2), textcoords='offset points', color=self.color, fontsize=self.fontsize)
        if self.plot_grids:
            ax = plt.gca()
            for i in range(self.tags.shape[0]):
                rect = patches.Rectangle((self.x_source[i] - self.Lx/2, self.y_source[i] - self.Ly/2), self.Lx, self.Ly, linewidth=self.linewidth, edgecolor='k', facecolor='none', linestyle=':')
                ax.add_patch(rect)   


class plot_probe_locations:
    
    def __init__(self, data_dir, color='white', fontsize=10, markersize=10):
        """
        Plot probe locations.

        Parameters
        ----------
        data_dir   : str
                     Where to look for probe files.
        color      : str, optional
                     Color of the source marker, by default 'white'.
        fontsize   : float or int
                     Fontsize of probe marker label.
        markersize : float or int
                     Markersize of probe marker.
        """
        
        files = get_files(data_dir=data_dir, prefix='probe', ignore='cumulative', verbose=False)
        try:
            tmp = xr.open_dataset(files[0]).load()
            tmp.close()
            self.tags   = tmp.probe.data
            self.x_prob = tmp.x.data
            self.y_prob = tmp.y.data
        except:
            logging.warning('Could not load probe file...')
            self.plot = self.noneplot
        self.color = color
        self.fontsize = fontsize
        self.markersize = markersize
    
    def noneplot(self):
        """
        Do not plot anything.
        """
        pass

    def plot(self):
        """
        Plot the probes.
        """
        plt.scatter(x=self.x_prob, y=self.y_prob, s=self.markersize, marker='x', facecolors=self.color)
        for i, txt in enumerate(self.tags):
            plt.annotate('${}^*$'.format(str(txt)), (self.x_prob[i], self.y_prob[i]), xytext=(-8,2), textcoords='offset points', color=self.color, fontsize=self.fontsize)
             

class plot_global_velocity_location:

    def __init__(self, config=None):
        
        if config is None:
            config=configparser.ConfigParser(allow_no_value=True)
            config.read('External-Indata/AdvDiff.ini')

         ### LOAD SOURCES ###
        sources, *_ = load_sources(config=config)

        ### LOAD VELOCITY ###
        global_velocity, *_ = load_velocity(config=config, verbose=False)

        self.config   = config
        self.sources  = sources

        self.Lx                          = int(self.config['grid']['Lx'])
        self.Ly                          = int(self.config['grid']['Ly'])
        self.fill_type                   = str(self.config['velocity']['fill_type'])
        self.allow_synthetic_translation = strtoBool(self.config['velocity']['allow_synthetic_translation'])
        self.allow_synthetic_scaling     = strtoBool(self.config['velocity']['allow_synthetic_scaling'])
        
        self.global_points      = np.stack([global_velocity.x.values, global_velocity.y.values], axis=-1) # Assuming x = X.flatten(), y = Y.flatten()
        self.global_hull        = ConvexHull(self.global_points)
        self.global_points_hull = self.global_points[np.unique(self.global_hull.simplices)]

        # Check if all sources lie in global velocity hull.
        if np.all(self.in_hull(self.global_points_hull, np.stack([self.sources.x.data, self.sources.y.data], axis=-1))): 
            self.global_velocity, *_ = load_velocity(config=config, return_unstructured=False, verbose=False)

        # If one or more sources lies within the global velocity hull, then we can directly interpolate
        # as long as we do not allow for synthetic translation. Sources outside the global velocity field will be ignored completely
        elif np.any(self.in_hull(self.global_points_hull, np.stack([self.sources.x.data, self.sources.y.data], axis=-1))) and not self.allow_synthetic_translation:
            self.global_velocity, *_ = load_velocity(config=config, return_unstructured=False, verbose=False)

        # If one or more sources lie outside the global velocity hull and synthetic translation is enabled. Then we will translate the velocities onto the source coordinates.
        elif self.allow_synthetic_translation:
            global_velocity, *_ = load_velocity(config=config, return_unstructured=False, verbose=False)
            if self.allow_synthetic_scaling:
                xmax_target = self.sources.x.max().values + self.Lx/2 # Global coords
                xmin_target = self.sources.x.min().values - self.Lx/2
                ymax_target = self.sources.y.max().values + self.Ly/2
                ymin_target = self.sources.y.min().values - self.Ly/2

                xmax = global_velocity.x.max().values
                ymax = global_velocity.y.max().values
                xmin = global_velocity.x.min().values
                ymin = global_velocity.y.min().values

                x = (xmax_target-xmin_target)*(global_velocity.x.values - xmin)/(xmax - xmin) + xmin_target # Translate and scale coorinates
                y = (ymax_target-ymin_target)*(global_velocity.y.values - ymin)/(ymax - ymin) + ymin_target 

            else:
                xmean_target = self.sources.x.mean().values # Global coords
                ymean_target = self.sources.y.mean().values

                xmean = global_velocity.x.mean().values
                ymean = global_velocity.y.mean().values

                x = (global_velocity.x.values - xmean) + xmean_target # Translate coordinates
                y = (global_velocity.y.values - ymean) + ymean_target 

            u = global_velocity.u.data
            v = global_velocity.v.data
            t = global_velocity.time.data

            data_vars = {"u" : (['time','y','x'], u), 
                         "v" : (['time','y','x'], v)}
                    
            global_velocity = xr.Dataset(data_vars = data_vars, 
                                coords={"x"    : (['x'],  x), 
                                        "y"    : (['y'],  y), 
                                        "time" : (['time'],  t),})
            
            self.global_velocity    = global_velocity
            X, Y = np.meshgrid(global_velocity.x.values, global_velocity.y.values)
            self.global_points      = np.stack([X.flatten(), Y.flatten()], axis=-1) # Assuming x = X.flatten(), y = Y.flatten()
            self.global_hull        = ConvexHull(self.global_points)
            self.global_points_hull = self.global_points[np.unique(self.global_hull.simplices)]

    def in_hull(self, hull, query_points):
        """
        Check if query points are inside the convex hull defined by a set of points.
        Thanks to the helpful people at https://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl
        for the algorithm.

        Parameters
        ----------
        hull         : array or Delaunay object
                       An MxK array of M points in K dimensions of points defining a convex hull or a scipy.spatial.Delaunay object.
        query_points : array
                       A set of points to query. Should be an NxK array coordinates of N points in K dimensions.

        Returns
        -------
        array
            Boolean array.
        """
        if not isinstance(hull, Delaunay):
            hull = Delaunay(hull)
        return hull.find_simplex(query_points)>=0

    def noneplot(self):
        """
        Do not plot anything.
        """
        pass

    def plot(self, linestyle='-', linewidth=2, color='r', ax=None):
        """
        Plot the convex hull of the velocity.

        Parameters
        ----------
        linestyle : str, optional
                    Linestyle, by default '-'
        linewidth : int, optional
                    Linewidth, by default 2
        color     : str, optional
                    Color, by default 'r'
        ax : object, optional
                    What axes object to plot to, by default None
        """
        if ax == None:
            for simplex in self.global_hull.simplices:
                plt.plot(self.global_points[simplex, 0], self.global_points[simplex, 1], linestyle=linestyle, linewidth=linewidth, color=color)
        else:
            for simplex in self.global_hull.simplices:
                ax.plot(self.global_points[simplex, 0], self.global_points[simplex, 1], linestyle=linestyle, linewidth=linewidth, color=color)
