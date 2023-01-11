#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:49:28 2020

@author: Guttorm & Ketil

This module defines the tracer_transport class, with a couple of methods. 

The module is built on the FiPy package, https://www.ctcms.nist.gov/fipy/index.html
It also uses numpy, scipy, pandas, math and random. 

"""
import logging
import fipy   as fp
import numpy  as np
import xarray as xr
import time   as TIME
from scipy.spatial        import KDTree, Delaunay, distance
from Tools.mesh_generator import initialize_meshes
from Tools.interpolator   import meshA_to_meshB
from tqdm.auto            import trange

#%%   
class tracer_transport:
    """
    A class that contains all information needed and methods for 
    simulating tracer transport. 
    
    The class uses the fipy package for solving the PDE. 
    
    Velocity fields and source location and fluxes are input to the class. 
    
    """
    
    
    def _viewer(self, variable):
        """
        Viewer function meant for debugging.

        Parameters
        ----------
        variable : array ([n_cells])
            Data variable to plot using the built-in FiPy plotter.
        """
        viewer = fp.Viewer(vars=variable,
                  limits={'ymin': self.y.min(), 'ymax': self.y.max(),
                          'xmin': self.x.min(), 'xmax': self.x.max()})
        viewer.plotMesh()
             

    def closest_node(self, mesh_points, query_point):
        """
        Finds the index of the nearest node from a set of nodes ([n_nodes,2]) to the point ([x,y]).
        This is used to find the nearest cell to the source.

        Parameters
        ----------
        query_point : array ([x,y])
                      Array containing a point (x,y).
        mesh_points : array ([n_nodes,2])
                      Set of points to query.

        Returns
        -------
        closest_index : int
                        index of the nearest node in nodes to the point (x,y).
        """
        closest_index = distance.cdist([query_point], mesh_points).argmin()
        return closest_index


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


    def mean_as_xr(self):
        """
        Construct xarray dataset containing the mean value of C

        Returns
        -------
        xarray dataset
            xarray dataset
        """
        name = 'mean'
        mean = self.accu/self.N
        return self.stats_to_xarray(mean, name=name)


    def max_as_xr(self):
        """
        Construct xarray dataset containing the max value of C

        Returns
        -------
        xarray dataset
            xarray dataset
        """
        name = 'max'
        max  = self.maxval
        return self.stats_to_xarray(max, name=name)


    def var_as_xr(self):
        """
        Construct xarray dataset containing the variance of C

        Returns
        -------
        xarray dataset
            xarray dataset
        """
        name  = 'var'
        mean  = self.accu/self.N
        mean2 = self.accu2/self.N
        var   = mean2-mean*mean
        return self.stats_to_xarray(var, name=name)

    
    def delta_mean_as_xr(self):
        """
        Construct xarray dataset containing the mean value of dC with different time-step sizes

        Returns
        -------
        xarray dataset
            xarray dataset
        """
        name = 'delta_mean'
        dmean_list = []
        for ii in range(self.step_N_stats.shape[0]):
            dmean = self.delta_accu[ii]/self.N_arr[ii]
            dmean_list.append(self.stats_to_xarray(dmean, name=name))
        ds = xr.concat(dmean_list, dim=xr.DataArray(data=self.dt_size*self.step_N_stats, dims='dt'))
        return ds


    def delta_mean_abs_as_xr(self):
        """
        Construct xarray dataset containing the mean value of abs(dC) with different time-step sizes

        Returns
        -------
        xarray dataset
            xarray dataset
        """
        name = 'delta_mean_abs'
        dmean_abs_list = []
        for ii in range(self.step_N_stats.shape[0]):
            dmean_abs = self.delta_accu_abs[ii]/self.N_arr[ii]
            dmean_abs_list.append(self.stats_to_xarray(dmean_abs, name=name))
        ds = xr.concat(dmean_abs_list, dim=xr.DataArray(data=self.dt_size*self.step_N_stats, dims='dt'))
        return ds


    def delta_max_as_xr(self):
        """
        Construct xarray dataset containing the max value of dC with different time-step sizes

        Returns
        -------
        xarray dataset
            xarray dataset
        """
        name = 'delta_max'
        dmax_list = []
        for ii in range(self.step_N_stats.shape[0]):
            dmax = self.delta_max[ii]
            dmax_list.append(self.stats_to_xarray(dmax, name=name))
        ds = xr.concat(dmax_list, dim=xr.DataArray(data=self.dt_size*self.step_N_stats, dims='dt'))
        return ds


    def delta_max_abs_as_xr(self):
            """
            Construct xarray dataset containing the max value of abs(dC) with different time-step sizes

            Returns
            -------
            xarray dataset
                xarray dataset
            """
            name = 'delta_max_abs'
            dmax_abs_list = []
            for ii in range(self.step_N_stats.shape[0]):
                dmax_abs = self.delta_max_abs[ii]
                dmax_abs_list.append(self.stats_to_xarray(dmax_abs, name=name))
            ds = xr.concat(dmax_abs_list, dim=xr.DataArray(data=self.dt_size*self.step_N_stats, dims='dt'))
            return ds


    def delta_min_as_xr(self):
        """
        Construct xarray dataset containing the min value of dC with different time-step sizes

        Returns
        -------
        xarray dataset
            xarray dataset
        """
        name = 'delta_min'
        dmin_list = []
        for ii in range(self.step_N_stats.shape[0]):
            dmin = self.delta_min[ii]
            dmin_list.append(self.stats_to_xarray(dmin, name=name))
        ds = xr.concat(dmin_list, dim=xr.DataArray(data=self.dt_size*self.step_N_stats, dims='dt'))
        return ds


    def delta_var_as_xr(self):
        """
        Construct xarray dataset containing the variance of dC with different time-step sizes

        Returns
        -------
        xarray dataset
            xarray dataset
        """
        name = 'delta_var'
        dvar_list = []
        for ii in range(self.step_N_stats.shape[0]):
            dvar = self.delta_accu2[ii]/self.N_arr[ii]-(self.delta_accu[ii]/self.N_arr[ii])**2
            dvar_list.append(self.stats_to_xarray(dvar, name=name))
        ds = xr.concat(dvar_list, dim=xr.DataArray(data=self.dt_size*self.step_N_stats, dims='dt'))
        return ds
    
    
    def get_stat(self):
        """
        Retrieves the statistics of the solution and merge all statistics into one xarray dataset object.

        Returns
        -------
        stat : xarray dataset
            An xarray dataset containing mean, variance, max, dmean, dvariance, dmax and dmin as variables

        """
        # Call on the stat functions to retrieve the stats the user wants to store
        stat = xr.merge([self.stat_functions[variable]() for variable in self.stat_functions if variable in self.stat_variables])
        stat.attrs = {'N_steps' : self.N} # Number of time-steps 
        return stat
    
    
    def get_transient_term(self):
        """
        Get a fipy.TransientTerm() object.

        Returns
        -------
        fipy.TransientTerm()
            Transient term for use in defining the differential equation.
        """
        return fp.TransientTerm(var=self.tracer)

    
    def get_convection_term(self):
        """
        Get a fipy.ConvectionTerm() object.

        Returns
        -------
        fipy.ConvectionTerm()
            Convection term for use in defining the differential equation.
        """
        if self.conv_type == 'PowerLaw':
            return fp.PowerLawConvectionTerm(var=self.tracer,          coeff=self.velocity) # Powerlaw convection term
        elif self.conv_type == 'CentralDiff':
            return fp.CentralDifferenceConvectionTerm(var=self.tracer, coeff=self.velocity) # Central difference convection term
        elif self.conv_type == 'Exponential':
            return fp.ExponentialConvectionTerm(var=self.tracer,       coeff=self.velocity) # Exponential convection term
        elif self.conv_type == 'Hybrid':
            return fp.HybridConvectionTerm(var=self.tracer,            coeff=self.velocity) # Hybrid convection term
        elif self.conv_type == 'Upwind':
            return fp.UpwindConvectionTerm(var=self.tracer,            coeff=self.velocity) # Upwind convection term
        elif self.conv_type == 'ExplicitUpwind':
            return fp.ExplicitUpwindConvectionTerm(var=self.tracer,    coeff=self.velocity) # Explicit convection term
        elif self.conv_type == 'VanLeer':
            return fp.VanLeerConvectionTerm(var=self.tracer,           coeff=self.velocity) # VanLeer convection term
        else:
            logging.critical(f'{self.conv_type} is not a valid convection term type...')
            raise ValueError(f'{self.conv_type} is not a valid convection term type...')
            
    
    def get_diffusion_term(self):
        """
        Get a fipy.DiffusionTerm() object or scalar.
        If the diffusion rate (self.D) == 0.0, then return a scalar == 0.0.
        Then no diffusion term is used in the differential equation. This may not be preferable. 

        Returns
        -------
        fipy.DiffusionTerm()
            Diffusion term for use in defining the differential equation.
        """
        if (self.D == 0.0).all():
            return 0.0                                                     # No diffusion term
        elif self.diff_type == 'Implicit':
            return fp.DiffusionTerm(var=self.tracer,         coeff=self.D) # Implicit diffusion term
        elif self.diff_type == 'Explicit':
            return fp.ExplicitDiffusionTerm(var=self.tracer, coeff=self.D) # Explicit diffusion term
        else:
            logging.critical(f'{self.diff_type} is not a valid diffusion term type...')
            raise ValueError(f'{self.diff_type} is not a valid diffusion term type...')
        

    def set_boundary_conditions(self):
        """
        Sets the boundary conditions for the tracer variable. Also sets an artifical boundary source term self._BC_SourceTerm().
        The artifical source is == 0.0 for every BC type except Robin.

        The Robin condition seems to be the best overall boundary condition, but I do not yet fully understand how it works.
        I have been copying ideas and methods from:
        https://www.ctcms.nist.gov/fipy/documentation/USAGE.html,
        https://www.ctcms.nist.gov/fipy/examples/convection/generated/examples.convection.robin.html,
        https://www.ctcms.nist.gov/fipy/documentation/USAGE.html#applying-robin-boundary-conditions,
        https://www.mail-archive.com/fipy@nist.gov/msg04328.html,
        https://stackoverflow.com/questions/61632193/1d-coupled-transient-diffusion-in-fipy-with-reactive-boundary-condition,
        https://github.com/usnistgov/fipy/issues/426,
        https://stackoverflow.com/questions/42769749/general-boundary-conditions,
        https://github.com/usnistgov/fipy/issues/788,
        https://github.com/usnistgov/fipy/blob/8efb4378648c339a710ce021fb5e0c25a1899b4e/examples/convection/robin.py,
        https://charlesreid1.com/wiki/Fipy_Verification,
        https://stackoverflow.com/questions/57754060/diffusion-reaction-pde-with-nonlinear-source-term-potential-term,

        Note: The Dirichlet condition behaves quite strangely and I do not know why. Might be a bug!
        """
        if self.BC_type == 'No_flux':
            pass                                                               # No_flux BC is the default BC applied to FiPY variables

        elif self.BC_type == 'Dirichlet':
            self.tracer.constrain(0.0,          where=self.mesh.exteriorFaces) # Dirichlet BC
        
        elif self.BC_type == 'Neumann':
            self.tracer.faceGrad.constrain(0.0, where=self.mesh.exteriorFaces) # Neumann   BC
        
        elif self.BC_type == 'Robin':                                          # Roundabout, but stable way of implementing Robin BCs
            np.seterr(divide='ignore', invalid='ignore')                       # Since D is constrained to 0 on the boundaries, numpy likes to print warnings about dividing by 0. I want to ignore this...
            DC = 1.0                                                           # Arbitrary dampening constant for how quickly tracer is dampened at the boundaries                                      
            self.velocity.constrain(0., where=self.mesh.exteriorFaces)         # To get Robin to work, I had to use ImplicitSourceTerms
            self.D.constrain(0.,        where=self.mesh.exteriorFaces)         # This means constraining D and v along the boundaries and letting sources along boundaries drive the behaviour along the boundary.
            self._BC_SourceTerm = fp.ImplicitSourceTerm(var=self.tracer, coeff=-DC*(self.mesh.exteriorFaces * self.mesh.faceNormals).divergence)
        
        elif self.BC_type == 'Robin experimental':                             # Robin BCs as shown in FiPy documentation
            DC = 0.0005
            self.tracer.faceGrad.constrain(-DC*self.tracer.faceValue, where=self.mesh.exteriorFaces)        
        
        else:
            logging.critical(f'{self.BC_type} is not a valid boundary condition type...')
            raise ValueError(f'{self.BC_type} is not a valid boundary conditiæon type...')
        
        if self.BC_type != 'Robin':
            self._BC_SourceTerm = 0.0 # If the main Robin BC is not used, set the BC sources to 0


    def set_velocity(self, u=1., v=1.):
        """
        Set the velocity for the solver.
        
        Parameters
        ----------
        u : array ([n_cells])
            Array of East-West velocity for every cell. The default is 1., i.e homogeneous.
        v : array([n_cells])
            Array of North-South velocity for every cell. The default is 1., i.e. homogeneous.

        """
        self.xVelocity.setValue(u) # Set cellVariable
        self.zVelocity.setValue(v) # Set cellVariable
        self.velocity.setValue(value=(self.xVelocity.arithmeticFaceValue, 
                                      self.zVelocity.arithmeticFaceValue)) # Set faceVariable
    

    def set_sources(self, source, view=None):
        """
        Set sources in relation to the computational mesh. 
        This is done by setting self._SourceTerm to non-zero at the gridcell nearest to the source coordinate. 
        The flux assigned to the cell is 1/Cell_volume.

        Parameters
        ----------
        source : array ([x,y])
                 Array containing the (x,y) coordinate for the source.
        view   : bool
                 Optional viewer, lets the user see where the source is placed for debugging. The default is None.

        """

        self.source = source

        self.flux = 1 
        CellVols  = self.mesh.cellVolumes
    
        nearest_index = self.closest_node(mesh_points=self.mesh_points, query_point=np.array([source.x.data, source.y.data]))
    
        self._SourceTerm[nearest_index] = self.flux / CellVols[nearest_index]
        
        # Incase something goes wrong, make sure the total flux == 1
        # tmp = np.sum(np.multiply(self._SourceTerm.value, CellVol)) # Compute total flux F = f_i * v_i
        # self._SourceTerm.setValue(self._SourceTerm.value / tmp)    # Make sure whatever source we have cumulatively contains a total flux of 1
        
        if view is not None: 
            self._viewer(variable=self._SourceTerm)
            
    
    def set_probes(self, probes):
        """
        Assigns a probe to the nearest meshpoint to the position given in probes. 
        Each probe retrieves the solution at the gridcell nearest the probe coordinate. 
        Any probe outside the global grid will only probe NaNs. Variables stored in self.probes.

        Parameters
        ----------
        probes        : xarray dataset
                        xarray containing (x,y) coordinates for probes.
        """
        
        #x, y = self.mesh.cellCenters
        #hull_points = self.mesh_points

        x, y = self.mesh_out.cellCenters
        hull_points = self.mesh_out_points

        # Use a K-D Tree to find the nearest neighbor meshpoints to the positions (a,b)
        tree = KDTree(list(zip(x.ravel(), y.ravel())))
        a    = probes.x.values
        b    = probes.y.values
        indx = tree.query(np.transpose(np.array([a,b]))) 
            
        time   = self.time.value   

        self.probes = probes.assign_coords(
            indx   = (['num'],  indx[1]),
            x_grid = (['num'],  x[indx[1]]),
            y_grid = (['num'],  y[indx[1]]),
            x_pos  = (['num'],  a),
            y_pos  = (['num'],  b),
            dist   = (['num'],  indx[0]),
            time   = (['time'], [time])
        )

        # If probe lands outside the local grid... Set any result to nan (Any shape works)
        # Sets condition that is used in self.probes_to_xarray()
        self.probe_condition = self.in_hull(hull=hull_points, query_points=np.stack([self.probes.x_pos.data, self.probes.y_pos.data], axis=-1))

    
    def set_stats(self):
        """
        Set and update the statistics variables. This function is called for every timestep.
        """
        # Collect and set tracer statistics
        self.accu  += self.tracer_output                  
        self.accu2 += self.tracer_output*self.tracer_output
        self.maxval = np.maximum(self.maxval, self.tracer_output)
        
        # Set delta statistics
        for ii in range(self.step_N_stats.shape[0]):
            if np.mod(self.N, self.step_N_stats[ii]) == 0:
                self.delta_tracer_arr    = (self.tracer_output-self.tracer_old[ii])
                self.delta_accu[ii]     += self.delta_tracer_arr
                self.delta_accu_abs[ii] += np.abs(self.delta_tracer_arr)
                self.delta_accu2[ii]    += self.delta_tracer_arr*self.delta_tracer_arr
                self.delta_max[ii]       = np.maximum(self.delta_max[ii],            self.delta_tracer_arr)
                self.delta_min[ii]       = np.minimum(self.delta_min[ii],            self.delta_tracer_arr)
                self.delta_max_abs[ii]   = np.maximum(self.delta_max_abs[ii], np.abs(self.delta_tracer_arr))


    def set_stat_variables(self, variables):
        """
        Used to modify what statistics to store during the simulation.
        By default: ['mean', 'max', 'var', 'delta_mean', 'delta_max', 'delta_var', 'delta_mean_abs','delta_max_abs', 'delta_min']

        Parameters
        ----------
        variables : list of str
                    List of variablenames to store. 
        """
        self.stat_variables = variables
        

    def set_time(self, time):
        """
        Set current time value.

        Parameters
        ----------
        time : time
        """
        self.time.setValue(time)


    def set_tracer_output(self):
        """
        Update self.tracer_output with current self.tracer value.
        If no interpolation is needed between the computational and output meshes, then self.tracer_output == self.tracer.value.
        """
        self.tracer_output = self.format_data(self.tracer.value)


    def set_tracer_old(self):
        """
        Update self.tracer_old for several steps back in time.
        The shape of each row of self.tracer_old is the same as self.tracer_output.
        """
        for ii in range(self.step_N_stats.shape[0]):         # Store old value for several steps back in time
            if np.mod(self.N, self.step_N_stats[ii]) == 0:   
                self.tracer_old[ii] = self.tracer_output

        
    def step(self, dt=1, sweeps=1, u=None, v=None, view=None):
        """
        Step the solver forward in time with step size dt and a given number of sweeps per dt.
        Updates self.tracer, self.tracer_old, self.tracer_output, self.velocity, self.xVelocity, self.zVelocity and the statistics.

        Parameters
        ----------
        dt     : float
                 Size of timestep.
        sweeps : int
                 Number of sweeps per timestep.
        u      : array ([n_cells])
                 Array of East-West velocity for every cell. 
        v      : array ([n_cells])
                 Array of North-South velocity for every cell.
        view   : bool
                 Optional viewer. The default is None.

        Returns
        -------
        None.
        """
        
        loc_dt = dt/sweeps
        u_old  = self.xVelocity.value
        v_old  = self.zVelocity.value
        
        if u is None:
            u = u_old
        if v is None:
            v = v_old 

        self.set_tracer_old()
        
        for ii in range(sweeps):
            # One sweep per dt   => loc_dt == dt
            # Many sweeps per dt => Need to interpolate velocities 
            tmp_fac = 1 if sweeps == 1 else ii/(sweeps-1)  
                
            u_in = (1-tmp_fac)*u_old+tmp_fac*u      # Linearly scale from u_old to u for each sweep
            v_in = (1-tmp_fac)*v_old+tmp_fac*v      # Linearly scale from v_old to v for each sweep
            
            self.sweep(dt=loc_dt, u=u_in, v=v_in, clip=False, view=None)  # Solve for one substep in the current dt

        self.N     +=  1
        self.N_arr += [1 if np.mod(self.N, self.step_N_stats[ii]) == 0 else 0 for ii in range(self.step_N_stats.shape[0])]
        
        # Update the output tracer with new self.tracer value. Interpolate if necessary.
        self.set_tracer_output()

        # Update stats
        self.set_stats()

        # Update velocity
        self.set_velocity(u=u, v=v)
        
        if view is not None:
            self._viewer(variable=self.tracer)
            

    def sweep(self, dt=1, u=None, v=None, clip=False, view=None):
        """
        Sweep the solution forward in time by some dt. This function is called once for each sweep by the self.step() function.
        Only updates self.tracer, self.velocity, self.xVelocity and self.yVelocity. 
        
        Parameters
        ----------
        dt   : float
               Size of timestep. The default is 1.
        u    : array ([n_cells])
               East-West velocity. If None, use the stored velocity. 
        v    : array ([n_cells]) 
               North-South velocity. If None, use the stored velocity.
        clip : bool, optional
               If True, clip self.tracer to contain only positive values. The default is False.
        view : bool, optional
               If not None, plot the tracer field. The default is None.
        """
    
        if u is not None or v is not None:
            self.set_velocity(u=u, v=v)

        self.eq.sweep(dt=dt)                   # Integrates the AdvDiff equation forward in time and updates self.tracer
        self.time.setValue(self.time.value+dt) # Updates time value
        if clip:
            self.tracer.setValue(np.clip(self.tracer.value, a_min=0.0, a_max=None)) # Force tracer to be positive by clipping negatives...

        if view is not None:
            self._viewer(variable=self.tracer)     


    def stats_to_xarray(self, data, name='Anon'):
        """
        Takes a data-array and constructs an xarray dataset from it.
        If the computational mesh and the output meshes are dissimilar, an interpolation is performed first.
        The output is indexed by x,y coordinates.

        Parameters
        ----------
        data  : array ([n_cells])
            Array containing the cell-variables for each cell in the computational mesh.
        name  : str
            Name to give to dataset.

        Returns
        -------
        dataset : xarray dataset 
            An xarray dataset containing the data from the array, indexed by x,y coordinates.
        """
        
        data_vars = {name : (['num'], data),}
        out = xr.Dataset(data_vars = data_vars, 
                         coords    = {"x" : (['num'], self.x_out), 
                                      "y" : (['num'], self.y_out),})
        
        if self.return_structured:
            out = out.set_index(num=("x","y")).unstack('num')
        return out   


    def probes_to_xarray(self, timestamp=None):
        """
        Probe the solution C at the nearest meshpoint to the measurement probes.
        If a probe lies outside the local mesh, then the probe will only output NANs.

        Parameters
        ----------
        timestamp : np.datetime64 
                    Timestamp of current timestep. If none, get stored timestamp. The default is None.

        Returns
        -------
        out : xarray dataset
              An xarray dataset containing C and dC at the measurement probes at the current timestep.

        """
        
        time = self.time.value if timestamp is None else timestamp    
        indx = self.probes.indx.values                                   # Gridcell index of probe
        #C    = np.where(self.probe_condition, self.tracer.value[indx], np.nan) # Ignore probes lying outside grid
        C    = np.where(self.probe_condition, self.tracer_output[indx], np.nan) # Ignore probes lying outside grid
        
        out         = self.probes.copy()
        out['time'] = ('time', [time])
        out['C']    = ('num',  C)
        return out


    def mass_to_xarray(self, timestamp=None):
        """
        Get the total mass and theoretical mass at the current timestep.

        Parameters
        ----------
        timestamp : np.datetime64
                    Timestamp of current timestep. If none, get stored timestamp. The default is None.
        
        Returns
        -------
        out : xarray dataset
              xarray containing the mass and theoretical mass of the system at the current timestep.

        """
        time = self.time.value if timestamp is None else timestamp
            
        if self.start_time is None:
            self.start_time = time
        
        C       = self.tracer.value
        CellVol = self.mesh.cellVolumes
        Mass    = np.sum(np.multiply(C, CellVol))
        
        # Compute total theoretical mass (M = F*t no volume included as F is total flux)
        Mass_t  = self.flux * (time-self.start_time)/np.timedelta64(1, 's')
        
        out = xr.Dataset({
            "Mass"             : (['time'], [Mass]),
            "Mass Theoretical" : (['time'], [Mass_t])},
            coords = {"time"   : (['time'], [time])})
        return out
    
    
    def fields_to_xarray(self, timestamp=None, store_uv=False):
        """
        Get the solution variables C (and u,v if store_uv == True) at the current time-step and store 
        the solution into a xarray dataset.

        Parameters
        ----------
        timestamp : np.datetime64
                    Timestamp of current timestep. If None, get stored timestamp. The default is None.
        store_uv  : bool
                    If True, store the local velocities u and v in addition to C. The default is False.

        Returns
        -------
        out : xarray dataset
              An xarray dataset containing C (and u,v if store_uv == True).

        """
        time = self.time.value if timestamp is None else timestamp 
        
        C  = self.format_data(self.tracer.value)
        data_vars = {"C" : (['num'], C)}

        if store_uv:
            u  = self.format_data(self.xVelocity.value)
            v  = self.format_data(self.zVelocity.value)
            data_vars.update({"u" : (['num'], u), 
                              "v" : (['num'], v)})
        
        out = xr.Dataset(data_vars = data_vars, 
                          coords={"x"    : (['num'],  self.x_out), 
                                  "y"    : (['num'],  self.y_out), 
                                  "time" : (['time'], [time]),})
        
        if self.return_structured:
            out = out.set_index(num=("x","y")).unstack('num')
        return out


    def interpolate(self, data):
        """
        Interpolate data defined on computational mesh onto the output mesh.

        Parameters
        ----------
        data : array 
               A numpy 1xN_cells data array.

        Returns
        -------
        array
              A numpy data array.
        """
        return self.interpolator.interpolate(data, extrap_nans=self.extrap_nans)

    
    def format_data(self, data):
        """
        Format data defined on computational mesh onto output mesh.
        Return the raw input data if no interpolation is needed between the computational and output meshes.
        Return interpolated input data if interpolation IS needed between the computational and output meshes.

        Parameters
        ----------
        data : array
               A numpy 1xN_cells data array.

        Returns
        -------
        array
              A numpy data array.
        """
        return data if self.is_onemesh else self.interpolate(data)


    def run(self, velocity, sweeps=1, 
            store_stats=True, store_fields=False, store_uv=False, store_probes=False, store_mass=False, 
            enable_pbar=True, desc_pbar='Progress: ', position_pbar=0):
        """
        Run a AdvDiff simulation.
        It steps through every time in the velocity data and generates a solution for the tracer variable.

        Parameters
        ----------
        velocity      : xarray dataset
                        An xarray dataset containing u and v velocity data for every gridcell at every time. 
        sweeps        : int, optional
                        Number of sweeps to perform per time-step, by default 1.
        store_stats   : bool, optional
                        Whether or not to solve for and store the statistics, by default True.
        store_fields  : bool, optional
                        Whether or not to solve and store for the field variables for every time, by default False.
        store_uv      : bool, optional
                        Whether or not to include the velocities in the field variables (if relevant), by default False.
        store_probes  : bool, optional
                        Whether or not to solve for and store the probe data for every time, by default False.
        store_mass    : bool, optional
                        Whether or not to solve for and store the total mass in the system for every time, by default False.
        enable_pbar   : bool, optional
                        Whether or not to enable the progress bar, by default True.
        desc_pbar     : str, optional
                        What prefix to assign progress bar, by default Progress:
        position_pbar : int, optional
                        Position of the progress bar (Allowable: Int >= 0), by default 0.
        

        Returns
        -------
        tuple
            A tuple containing an outdata dictionary which stores the statistics, fields, probes and mass and the average process time for each time-step dt.
        """
        if self.probes is None: # Ensure that probes are not run if user has not set probes.
            store_probes = False
        
        ### INITIALIZE OUTPUT ###
        self.dt_size = (velocity.time.diff(dim='time')/np.timedelta64(1, 's')).mean().data # dt_size for use in delta statistics

        time = velocity.time.min()

        self.set_velocity(u=velocity.sel(time=time).u.data, v=velocity.sel(time=time).v.data) # Set initial velocity.
        
        fields_arr  = [self.fields_to_xarray(timestamp=time.values, store_uv=store_uv)] if store_fields else None # Initialize fields outdata 
        probes_arr  = [self.probes_to_xarray(timestamp=time.values)]                    if store_probes else None # Initialize probes outdata
        masses_arr  = [self.mass_to_xarray(timestamp=time.values)]                      if store_mass   else None # Initialize masses outdata
        
        ### ITERATE THROUGH TIME ###
        start_time = TIME.process_time()
        for ii in trange(np.shape(velocity.time)[0]-1, desc=desc_pbar, disable=not enable_pbar, position=position_pbar, 
                         miniters=16, lock_args=None, leave=False): # For-loop with multiprocessing progressbars.
                
            dt   = (velocity.time[ii+1]-velocity.time[ii])/np.timedelta64(1, 's')
            time =  velocity.time[ii+1].values 
            self.step(dt=dt, sweeps=sweeps, u=velocity.isel(time=ii+1).u.data, v=velocity.isel(time=ii+1).v.data)
            
            fields_arr.append(self.fields_to_xarray(timestamp=time, store_uv=store_uv)) if store_fields and np.mod(self.N, self.step_N_fields) == 0 else None
            probes_arr.append(self.probes_to_xarray(timestamp=time))                    if store_probes else None
            masses_arr.append(self.mass_to_xarray(timestamp=time))                      if store_mass   else None
                
        end_time = TIME.process_time()
        avg_process_time = (end_time-start_time)/self.N # Average time spent per time-step

        stats    = self.get_stat()                   if store_stats  else None # Get stats
        fields   = xr.concat(fields_arr, dim='time') if store_fields else None # Concatenate fields over all times  
        probes   = xr.concat(probes_arr, dim='time') if store_probes else None # Concatenate probes over all times
        mass     = xr.concat(masses_arr, dim='time') if store_mass   else None # Concatenate masses over all times

        outdata = {'stats'        : stats,
                   'fields'       : fields,
                   'probes'       : probes,
                   'mass'         : mass,
                   'process_time' : avg_process_time,
                   }

        return outdata


    def __init__(self, grid_type="equi", grid_type_out='equi', conv_type='PowerLaw', diff_type='Implicit', BC_type='Robin',
                 D=0.1, Lx=2000, Ly=2000, nx=100, ny=100, reduction_factor=0.96, maxvol=25000, minvol=250, step_N_stats=[1], step_N_fields=1,
                 source=None, probes=None, meshdict=None, extrap_nans=True, return_structured=True):
        """
        A tracer-transport (Advection-Diffusion) solver which solves for the tracer transport problem around one source located
        in the center of the domain with total flux == 1. 
        The solver utilizes the FiPY package to step forwards in time.

        Note: You do not need to specify grid and domain variables if meshdict is given.

        Parameters
        ----------
        grid_type         : str 
                            Type of computational grid. 
        grid_type_out     : str
                            Type of output grid.
        BC_type           : str
                            Type of boundary condition to apply.
        conv_type         : str
                            Type of convection.
        diff_type         : str
                            Type of diffusion.
        D                 : float
                            The diffusion coefficient. 
        Lx/Ly             : float
                            Length of the domain in the x/y direction. 
        nx/ny             : int
            	            Number of discretizations in x/y direction in a rectilinear mesh. 
        reduction_factor  : float
                            How much strong the exponential reduction of gridcell size is when using expnential grid.
        maxvol            : float
                            The maximum volume for the largest gridcell in a triangular mesh.
        minvol            : float
                            The minimum volume for the smallest gridcell in a triangular mesh.
        step_N_stats      : list
                            List of ascending integers greater or equal to 1. Each integer sets for what intervals of time-steps to store delta variables.
        step_N_fields     : int
                            Int greater or equal to 1. The integer sets for what intervals of time-steps to store the fields.
        source            : xarray dataset
                            An xarray containing location for origin of local grid. (The middle of the grid).
        probes            : xarray dataset
                            An xarray containing x,y locations for probes. Not needed if probes are set after init with .set_probes().
        meshdict          : dict
                            Dictionary containing mesh parameters. Not needed if meshes are to be generated for each instance of tracer_transport().
                            If meshdict is given, then it will override all other mesh related parameters!!!
        extrap_nans       : bool
                            To extrapolate NANs when interpolating between meshes. RECOMMENDED TRUE.
        return_structured : bool
                            To return the statistics and fields in structured form aka f(t,x,y) instead of f(t,num).
                            If false, then global stats and fields will not work correctly and you will not be able to plot the output easily.

        Returns
        -------
        tracer_transport : object
                           Tracer_transport object solving the tracer-transport problem. 

        """
        # Set variables
        self.conv_type     = conv_type     # Convection type
        self.diff_type     = diff_type     # Diffusion type
        self.BC_type       = BC_type       # Boundary condition type

        # Initialize meshes
        if meshdict is not None:           # If meshdict is given, use information contained within (FAST as you do not have to regenerate the meshes).
            self.meshdict = meshdict       # Mesh information
        else:                              # Else initialize meshes now using the input patameters (SLOW as you have to generate a new set of meshes)                                                   
            self.meshdict = initialize_meshes(grid_type        = grid_type,        # Computational mesh type
                                              grid_type_out    = grid_type_out,    # Output mesh type    
                                              Lx               = Lx,               # Grid width
                                              Ly               = Ly,               # Grid height
                                              nx               = nx,               # Number of discretizations along width  (Only relevant for rectilinear mesh)
                                              ny               = ny,               # Number of discretizations along height (Only relevant for rectilinear mesh)
                                              maxvol           = maxvol,           # Maximum volume of largest triangle     (Only relevant for triangular mesh)
                                              minvol           = minvol,           # Minimum volume of smallest triangle    (Only relevant for triangular mesh)
                                              reduction_factor = reduction_factor, # Reduction factor for exp grids
                                              origin           = np.array([source.x.data, source.y.data]), # Source location
                                              verbose          = False) 

        self.mesh              = self.meshdict['mesh']            # Computational mesh
        self.mesh_out          = self.meshdict['mesh_out']        # Output mesh
        self.is_onemesh        = self.meshdict['is_onemesh']      # Is computational mesh == output mesh
        self.is_unstructured   = self.meshdict['is_unstructured'] # Is the computational mesh unstructured
        self.return_structured = return_structured                # Return output in structured form? (Bool)
        self.x,     self.y     = self.mesh.cellCenters            # The x, y coords on computational mesh
        self.x_out, self.y_out = self.mesh_out.cellCenters        # The x, y coords on output mesh
        

        # Initialize interpolator
        self.mesh_points       = np.transpose(self.mesh.cellCenters)     # Computational mesh points
        self.mesh_out_points   = np.transpose(self.mesh_out.cellCenters) # Output mesh points
        self.extrap_nans       = extrap_nans                             # Extrapolate possible NaNs when interpolating? (Bool)
        self.interpolator      = meshA_to_meshB(meshA_points=self.mesh_points, meshB_points=self.mesh_out_points) # Set interpolator

        # Define diffusion rate as a field
        self.D         = fp.FaceVariable(mesh=self.mesh, name='diffusion', rank=0, value=D) 

        # Define the solution variable
        self.tracer    = fp.CellVariable(mesh=self.mesh, name='tracer',    rank=0, value=0.)
        
        # Set velocity 
        self.xVelocity = fp.CellVariable(mesh=self.mesh, name='Xvelocity', rank=0, value=0.)      # cellVariable
        self.zVelocity = fp.CellVariable(mesh=self.mesh, name='Zvelocity', rank=0, value=0.)      # cellVariable
        self.velocity  = fp.FaceVariable(mesh=self.mesh, name='velocity',  rank=1, value=(0.,0.)) # faceVariable

        # Set boundary conditions
        self.set_boundary_conditions()

        # Define differential terms and source term
        self._TransientTerm  = self.get_transient_term()
        self._ConvectionTerm = self.get_convection_term()
        self._DiffusionTerm  = self.get_diffusion_term()
        self._SourceTerm     = fp.CellVariable(mesh=self.mesh, name='source', rank=0, value=0.) # Value of source is given by set_source()
        
        # Define the PDE for tracer transport. The _BC_SourceTerm is only relevant for Robin BC.
        self.eq = (self._TransientTerm + self._ConvectionTerm == self._DiffusionTerm + self._SourceTerm + self._BC_SourceTerm) 

        # Initialize time
        self.start_time = None
        self.time       = fp.Variable(0.)
        
        # Initialize time-steps (N). 
        self.N             = 0
        self.step_N_stats  = np.array(step_N_stats)
        self.step_N_fields = step_N_fields
        self.N_arr         = self.N*np.ones(shape=self.step_N_stats.shape)

        # Set source field
        self.set_sources(source=source)

        # Set probes
        if probes is not None: # If probes are given, set them now
            self.set_probes(probes=probes) 
        else:
            self.probes = None

        # Set the statistics functions
        self.stat_functions = {'mean'           : self.mean_as_xr,
                               'max'            : self.max_as_xr,
                               'var'            : self.var_as_xr,
                               'delta_mean'     : self.delta_mean_as_xr,
                               'delta_max'      : self.delta_max_as_xr,
                               'delta_var'      : self.delta_var_as_xr,
                               'delta_mean_abs' : self.delta_mean_abs_as_xr,
                               'delta_max_abs'  : self.delta_max_abs_as_xr,
                               'delta_min'      : self.delta_min_as_xr}
        
        # Set the statistics variable names to store
        self.stat_variables = [variable for variable in self.stat_functions]                                   

        # Set output variable
        # This variable is the same as self.tracer, but interpolated onto the output mesh.
        # It is used for storing statistics as it is better to do interpolation first then compute statistics, instead of the other way around. 
        # Notice: self.tacer_old also gets the same shape as self.tracer_output, not self.tracer. 
        # This is because self.tracer_old is only used to compute the delta statistics.
        self.tracer_output  = self.format_data(self.tracer.value)

        # Set statistics 
        self.accu           = np.zeros(shape=self.tracer_output.shape) 
        self.accu2          = np.zeros(shape=self.tracer_output.shape) 
        self.maxval         = np.zeros(shape=self.tracer_output.shape) 

        # Set delta statistics 
        self.tracer_old     = np.zeros(shape=(self.step_N_stats.shape[0], self.tracer_output.shape[0])) 
        self.delta_accu     = np.zeros(shape=(self.step_N_stats.shape[0], self.tracer_output.shape[0]))
        self.delta_accu_abs = np.zeros(shape=(self.step_N_stats.shape[0], self.tracer_output.shape[0]))
        self.delta_accu2    = np.zeros(shape=(self.step_N_stats.shape[0], self.tracer_output.shape[0]))
        self.delta_max      = np.zeros(shape=(self.step_N_stats.shape[0], self.tracer_output.shape[0])) 
        self.delta_max_abs  = np.zeros(shape=(self.step_N_stats.shape[0], self.tracer_output.shape[0]))
        self.delta_min      = np.zeros(shape=(self.step_N_stats.shape[0], self.tracer_output.shape[0])) 
