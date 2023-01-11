# -*- coding: utf-8 -*-
"""
@author: Ketil

This module defines the global_to_local_variables class, with a number of methods. 
The global_to_local_variables class module extracts the relevant local variables from the global ones for use in AdvDiff. 

"""
import logging
import numpy  as np
import xarray as xr
import matplotlib.pyplot as plt
from Tools.parser_tools    import strtoBool
from Tools.interpolator    import meshA_to_meshB
from scipy.spatial         import ConvexHull, Delaunay


class global_to_local_variables:

    
    def velocity_viewer(data, time, tag, variable='u', name='none'):
        """
        Viewer meant for debugging.

        Parameters
        ----------
        data     : xarray dataset
                An xarray dataset containing velocity.
        time     : int
                What timestep to plot.
        tag      : int
                What source does the data belong to.
        variable : str, optional
                Which velocity component to plot, by default 'u'.
        name     : str, optional
                Name of output image, by default 'none'.
        """
        plt.close('all')
        out_dir='Figures/'
        data.isel(time=time)[variable].plot(x='x',y='y') 
        file_out = out_dir+name+str(tag)+'.'+'.jpg'    
        logging.info(f'Saving: {file_out}')
        plt.xticks(rotation=30)
        plt.savefig(file_out, dpi=200)
        plt.close("all")


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
    

    def get_local_grid_hull_points(self, tag):
        """
        Get convex hull points defining the local grid.

        Parameters
        ----------
        tag : int
              source tag

        Returns
        -------
        array
            A set of points defining the convex hull of the local grid.
        """
        local_meshdict    = self.get_local_meshes(tag=tag)
        local_points      = np.transpose(local_meshdict['mesh'].cellCenters)
        local_hull        = ConvexHull(local_points)
        local_hull_points = local_points[np.unique(local_hull.simplices)]
        return local_hull_points


        # Interpolates the velocity file from the velocity mesh onto the computational mesh.
    def extract_and_interpolate_velocity(self, xmid, ymid, extrap_from_A=True):
        """
        Interpolates the velocity dataset onto the computational mesh which will be used by the fipy solver.
        Interpolation is done using linear meshA_to_meshB interpolator. Any NAN values encountered during 
        interpolation may be extrapolated using nearest neighbor interpolator, zeroes, an average value or lead to the velocity being ignored.
        This is based on the fill type that is set by the user.

        Parameters
        ----------
        xmid : float
            x-coordinate of center of local grid.
        ymid : float
            y-coordinate of center of local grid.
        extrap_from_A : bool, optional
            Whether or not NaNs should be extrapolated using data stored on the original mesh (A) or new mesh (B), by default True

        Returns
        -------
        _type_
            _description_
        """
        
        u_arr = [] 
        v_arr = []
        
        # If our computational grid is a nonstructured triangle grid, then we have no option but to interpolate using a
        # custom linear interpolator method (meshA_to_meshB)...

        mesh = self.meshdict['mesh'] + ((xmid,),(ymid,))
        comp_mesh_points = np.transpose(mesh.cellCenters)
        
        if self.fill_type == 'nearest extrapolation':
            extrap_nans   = True
            fill_value_u  = np.nan
            fill_value_v  = np.nan

        elif self.fill_type == 'zeroes':
            extrap_nans  = False
            fill_value_u = 0.0
            fill_value_v = 0.0

        elif self.fill_type == 'average':
            extrap_nans  = False
            fill_value_u = None
            fill_value_v = None

        elif self.fill_type == 'ignore':
            extrap_nans  = False
            fill_value_u = np.nan
            fill_value_v = np.nan

        else: # Default to nearest
            extrap_nans  = True
            fill_value_u = np.nan
            fill_value_v = np.nan

        vel_mesh_to_comp_mesh = meshA_to_meshB(meshA_points=self.global_points, meshB_points=comp_mesh_points)
        for time in self.global_velocity.time.values:
                
                data_u = self.global_velocity['u'].sel(time=time).data
                data_v = self.global_velocity['v'].sel(time=time).data

                if self.fill_type == 'average':
                    fill_value_u = np.mean(data_u)
                    fill_value_v = np.mean(data_v)
                
                interp_u = vel_mesh_to_comp_mesh.interpolate(data_u, extrap_nans=extrap_nans, extrap_from_A=extrap_from_A, fill_value=fill_value_u)
                interp_v = vel_mesh_to_comp_mesh.interpolate(data_v, extrap_nans=extrap_nans, extrap_from_A=extrap_from_A, fill_value=fill_value_v)

                if self.fill_type == 'ignore' and np.any([np.isnan(interp_u), np.isnan(interp_v)]):
                    return None
                
                u_arr.append(interp_u)
                v_arr.append(interp_v)   


        data_vars = {"u" : (['time','node'], u_arr), 
                     "v" : (['time','node'], v_arr)}
            
        local_velocity = xr.Dataset(data_vars = data_vars, 
                            coords={"x"    : (['node'],  comp_mesh_points[:,0]), 
                                    "y"    : (['node'],  comp_mesh_points[:,1]), 
                                    "time" : (['time'],  self.global_velocity.time.values),})  

        return local_velocity


    # Check which sources will have missing values and need to be filled / ignored...
    def perform_fill_check(self):
        """
        Checks each local grid if any intersect with the exterior of the convex hull of the global velocity data. If so
        warn the user which sources this affects.
        """
        affected_tags = []
        for tag in self.sources.source.data:
        #for tag in self.sources.source.data:
            local_hull_points = self.get_local_grid_hull_points(tag=tag)
            if np.all(self.in_hull(self.global_points_hull, local_hull_points)):
                # Local grid is all inside global grid, no fill will be needed...
                pass
            else:
                # Local grid lies partially outside global grid, some fill treatment will be needed...
                affected_tags.append(tag)
        
        if len(affected_tags) > 0:
            if self.fill_type != 'ignore':
                logging.warning('The local grid around some of the sources may partially lie outside the global grid.')
                logging.warning('Any missing velocity values will be filled with {}...'.format(str(self.fill_type)))
                logging.warning('This affects source: {}\n'.format(str(affected_tags)))

            elif self.fill_type == 'ignore':
                logging.warning('The local grid around some of the sources may partially lie outside the global grid.')
                logging.warning('Any source with missing velocity values will be IGNORED...')
                logging.warning('This affects source: {}\n'.format(str(affected_tags)))


    def perform_source_check(self):
        """
        Check if sources are inside convex hull of the velocity data. 
        Any sources outside the convex hull will lead to a warning being printed (as long as verbose is enabled).
        """

        affected_tags = []
        for tag, x, y in zip(self.sources.source.data, self.sources.x.data, self.sources.y.data):
            if (self.in_hull(self.global_points_hull, np.array([x, y]))):
                pass # Source inside hull, check passed.
            else:
                affected_tags.append(tag)
        
        if len(affected_tags) > 0:
            logging.warning('Some sources lie outside the global grid.')
            logging.warning('Any such sources will be ignored...')
            logging.warning('This affects source: {}\n'.format(str(affected_tags)))

    
    def perform_probe_check(self):
        """
        Check if probes are inside the convex hull of the velocity data. 
        Any probes outside the convex hull will lead to a warning being printed. 
        If all probes are outside the convex hull, then probes will be set to None so that they can be ignored.
        """
        affected_tags = []
        tag = 0
        for x, y in zip(self.probes.x.data, self.probes.y.data):
            if (self.in_hull(self.global_points_hull, np.array([x, y]))):
                pass # Probe inside hull, check passed.
            else:
                affected_tags.append(tag)
            tag += 1

        if len(self.probes.x.data) > len(affected_tags) > 0:
            logging.warning('Some probes lie outside the global grid.')
            logging.warning('This affects probe: {}\n'.format(str(affected_tags)))
        elif len(self.probes.x.data) == len(affected_tags):
            logging.critical('No probes lie inside the global grid.')
            logging.critical('Probes will be ignored...\n')
            self.probes = None

    # Slice and extract a local velocity from the global velocity file...
    def get_local_velocity(self, tag, origin=None, extrap_from_A=True, view=False):
        """

        Get a local velocity xarray dataset extracted from the global velocity.

        Parameters
        ----------
        tag    : int
                 What source to extract velocities around.
        origin : array ([x,y]), optional
                 Location in space to extract velocities around, by default None.
        view   : bool, optional
                 Whether or not to plot and store the extracted velocity, by default False.

        Returns
        -------
        xarray dataset
                 An xarray dataset containing velocities from a local region around a source.
        """

        xmid = self.sources.isel(source=tag).x.data if origin==None else origin[0]
        ymid = self.sources.isel(source=tag).y.data if origin==None else origin[1]
        
        if (self.in_hull(self.global_points_hull, np.array([xmid, ymid]))):
            local_velocity = self.extract_and_interpolate_velocity(xmid, ymid, extrap_from_A=extrap_from_A)
        else:
            local_velocity = None

        if view == True and local_velocity is not None:
            self.velocity_viewer(data=local_velocity, time=0, tag=tag, variable='u', name='velocity_local_preinterp')

        return local_velocity


    def get_local_source(self, tag):
        """
        Extract single source given a tag from dataset containing multiple sources.

        Parameters
        ----------
        sources : xarray dataset
            A xarray dataset containing sources
        tag : int
            What source to extract

        Returns
        -------
        xarray dataset
            A xarray dataset containing a single source.
        """
        source = self.sources.isel(source=tag)
        return source


    def get_probe_locations(self, tag):
        """
        Get probe locations.

        Parameters
        ----------
        probes : xarray dataset
            A xarray dataset containing probe locations.
        source : xarray dataset
            A xarray dataset containing one source location.

        Returns
        -------
        xarray dataset
            A xarray dataset containing probe locations
        """
        if self.probes is not None:    
            probes = self.probes
        else: 
            probes = None
        return probes


    def get_local_meshes(self, tag, origin=None):
        """
        Get local mesh centered around source given by tag, or around coordinates given by origin.

        Parameters
        ----------
        tag    : int
                 Integer for source which source to center the local mesh around.
        origin : array ([x,y]), optional
                 Coordinate for where to center the local mesh if you do not want to use a specific source, by default None

        Returns
        -------
        dict
            dictionary containing computational and output meshes that have been translated to the relevant coordinates.
        """

        xmid     = self.sources.isel(source=tag).x.data if origin is None else origin[0]
        ymid     = self.sources.isel(source=tag).y.data if origin is None else origin[1]
        meshdict = self.meshdict.copy() # Important to copy or else self.meshdict will also change.

        meshdict['origin']   = meshdict['origin']   + np.array([xmid, ymid]) # Translate origin to source origin
        meshdict['mesh']     = meshdict['mesh']     + ((xmid,), (ymid,))     # Translate mesh to source origin
        meshdict['mesh_out'] = meshdict['mesh_out'] + ((xmid,), (ymid,))     # Translate mesh to source origin

        return meshdict


    def __init__(self, config, global_velocity, sources, probes, meshdict):
        """
        Generates the local velocities given some global velocity file and extracts other local variables such as
        source locations, probe locations and local meshes. 
        There are several contingencies which allows us to generate
        local velocities even if we do not have an overlap with source coordinates and velocity coordinates. There are 3 scenarios to consider:
            
        *   If all sources lie within the bounding box of velocity coordinates, then we directly interpolate.
            This option is what should be used in an employed model for which we are required
            that sources and velocity coordinates do indeed overlap before simulations are run.
            
        *   If one or more sources lie within the bounding box of the global velocity field, then we directly interpolate
            as long as synthetic_translation == False. Sources outside the global velocity field will be ignored completely as
            we do not have enough velocity information here to run that simulation. A warning will be raised to console if a source has been ignored.
            
        *   If one or more sources lie outside the bounding box of the global velocity field AND synthetic_translation == True,
            then we will artificially map the global velocity onto the riskmap coordinates such that all sources will be inside the 
            new global velocity field by construction. If however, no riskmap can be found, we instead artificially map the global
            velocity field onto the maximum bounding box given by source coordinates and local grid sizes (Lx,Ly).
            Note: synthetic_translation == True will lead to unrealistic simulations as we will translate and stretch the global
            velocities to fit our sources. synthetic_translate == True should only be used for testing.

        config file requires:
            
            * (str)     : config['veloctiy']['fill_type']
            * (bool)    : config['velocity']['allow_synthetic_translation']
            * (bool)    : config['velocity']['allow_synthetic_scaling']
            * (str)     : config['paths']['indata_path']
            * (str)     : config['paths']['riskmap_path']
            
        Parameters
        ----------
        config          : setup.ini
                          Config file with relevant variables and setup configurations.     
        global_velocity : xarray dataset
                          An xarray dataset containing the global velocity components u and v for (t,x,y).
        sources         : xarray dataset
                          An xarray dataset containing the (x,y) coordinates for each source.
        probes          : xarray dataset
                          An xarray dataset containing the (x,y) coordinates for each probe.
        meshdict        : dictionary
                          A dictionary containing all necessary computational- and output-mesh information.

        Returns
        -------
        global_to_local_variables : object
                         An global_to_local_variables object which can be used to generate the necessary local variables used for each instance of tracer_transport().
        """
        logging.info('Initializing global-to-local generator...\n')

        self.meshdict = meshdict
        self.config   = config
        self.sources  = sources
        self.probes   = probes

        self.Lx                          = meshdict['Lx']
        self.Ly                          = meshdict['Ly']
        self.fill_type                   = str(self.config['velocity']['fill_type'])
        self.allow_synthetic_translation = strtoBool(self.config['velocity']['allow_synthetic_translation'])
        self.allow_synthetic_scaling     = strtoBool(self.config['velocity']['allow_synthetic_scaling'])
        
        self.global_points      = np.stack([global_velocity.x.values, global_velocity.y.values], axis=-1) # Assuming x = X.flatten(), y = Y.flatten()
        self.global_hull        = ConvexHull(self.global_points)
        self.global_points_hull = self.global_points[np.unique(self.global_hull.simplices)]

        # Check if all sources lie in global velocity hull.
        if np.all(self.in_hull(self.global_points_hull, np.stack([self.sources.x.data, self.sources.y.data], axis=-1))): 
            logging.info('All sources lie within the global velocity field...')
            logging.info('No synthetic translation or scaling needed...\n')
            self.global_velocity = global_velocity

        # If one or more sources lies within the global velocity hull, then we can directly interpolate
        # as long as we do not allow for synthetic translation. Sources outside the global velocity field will be ignored completely
        elif np.any(self.in_hull(self.global_points_hull, np.stack([self.sources.x.data, self.sources.y.data], axis=-1))) and not self.allow_synthetic_translation:
            logging.warning('Only some sources lie within the global velocity field...')
            logging.warning('Synthetic translation: Disabled')
            logging.warning('Any sources outside the global velocity field will be ignored...\n')
            self.global_velocity = global_velocity

        # If one or more sources lie outside the global velocity hull and synthetic translation is enabled. Then we will translate the velocities onto the source coordinates.
        elif self.allow_synthetic_translation:
            if self.allow_synthetic_scaling: # Allow synthetic scaling of velocity to fit
                logging.warning('The sources do NOT lie within the global velocity coordinates.')
                logging.warning('Synthetic translation: Enabled')
                logging.warning('Synthetic scaling:     Enabled')
                logging.warning('Velocity will be translated and scaled to fit over sources...\n')
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

            else: # No synthetic scaling 
                logging.warning('The sources do NOT lie within the global velocity coordinates.')
                logging.warning('Synthetic translation: Enabled')
                logging.warning('Synthetic scaling:     Disabled')
                logging.warning('Velocity will be translated to sources...\n')
                xmean_target = self.sources.x.mean().values # Global coords
                ymean_target = self.sources.y.mean().values

                xmean = global_velocity.x.mean().values
                ymean = global_velocity.y.mean().values

                x = (global_velocity.x.values - xmean) + xmean_target # Translate coordinates
                y = (global_velocity.y.values - ymean) + ymean_target 

            u = global_velocity.u.data
            v = global_velocity.v.data
            t = global_velocity.time.data

            data_vars = {"u" : (['time','node'], u), 
                         "v" : (['time','node'], v)}
                    
            global_velocity = xr.Dataset(data_vars = data_vars, 
                                coords={"x"    : (['node'],  x), 
                                        "y"    : (['node'],  y), 
                                        "time" : (['time'],  t),})
            
            self.global_velocity    = global_velocity
            self.global_points      = np.stack([global_velocity.x.values, global_velocity.y.values], axis=-1) # Assuming x = X.flatten(), y = Y.flatten()
            self.global_hull        = ConvexHull(self.global_points)
            self.global_points_hull = self.global_points[np.unique(self.global_hull.simplices)]

        self.perform_source_check()
        self.perform_probe_check()
        self.perform_fill_check()