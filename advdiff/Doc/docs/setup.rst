.. _setup:

Setup.ini
==========

The config file contains all the options and parameters describing how to run the module on the given data. 
The config file is divided into sections which specify what parts of the module are affected.

.. include:: .//setup.ini
    :literal:


[DEFAULT] 
----------
Here we can specify how many CPU cores to use in parallel. When running the module, each CPU core works on one source at a time. 
We can also specify verbosity of the module and how much information and feedback to write to stdout and the log file as the module is running.
In addition, we can enable or disable the progress bar for each CPU core.

* **max_cpus**: Number of CPU threads to run in parallel. Set to -1 to run all.

* **verbose**: Level of verbosity of module. Allowable levels: DEBUG, INFO, WARNING, CRITICAL. 

* **progress_bar**: If True, enable progression bar for each individual source. If False, no progression bars will be shown during runtime.


[setup]
--------
Here we can set where we get the sources from, the start time for the simulation, and for how long to run the simulation for. 
In addition we can set the numerical schemes to be used by our solver, and the diffusion rate.

* **get_source_from**: Where to get sources from. You can set the module to generate sources from sources, riskmaps, wells, riskmaps and wells or riskmap and sources. 
  The file formats for each type of file is shown in the :ref:`inputs section <input>`.


* **use_vel_timer / time_delta / time_delta_unit / time_seed**: If **use_vel_timer** is set to **True**, then the module uses the timer in the given velocity data. 
  In other words, the start time begins at the first time in the velocity data while the end time is at the last time in the velocity data.
  BEWARE setting this to true if the velocity file spans a very long time-scale.
  If **use_vel_timer** is set to **False**, then we can specify the start time of the simulation using **time_start**, 
  and the duration using **time_delta** together with **time_delta_unit**. 
  The **time_seed** parameter is only relevant when our start time is set to **Random**.
  

* **dt_size / dt_unit**: **dt_size** sets what time-step size to integrate the solution with, this can be set to any positive integer. 
  The **dt_unit** variable sets the corresponding temporal unit (allowbale:  S, min, H). 
  Effectively, this will also interpolate the provided velocity file with the same time-step increments.
  For example, if the velocity file steps in time with 1 hour increments, 
  we can interpolate it over 15 minute increments by setting **dt_size = 15** and **dt_unit = min**.
  The opposite case can also be achieved, if the velocity file steps in time with 15 minute increments, 
  you can set **dt_size = 1** and **dt_unit = H** to force time-stepping in 1 hour increments.
  The main purpose behind this is to unify the time-step lengths when using different velocity datasets. 
  To maintain order of accuracy between different datasets, it is best to keep time-steps the same length. 


* **Sweeps / max_dt / max_dt_unit**: **Sweeps** specifies how many sweeps the module performs between each time step **dt**. 
  The time step length **dt** is determined by the time step length between each instance of time in the velocity file. By default **sweeps = 1**, 
  but it can be set to any positive integer, or to **Dynamic**. 
  If it set to **Dynamic**, then the module will compute how many sweeps are required to meet a maximum time step **max_dt** with unit **max_dt_unit**. 
  So for instance, if the velocity file increments in time by dt=1h, but we require that the module maximally steps forward in time by 15m, 
  then we can either set sweeps = 4 or we can set sweeps = Dynamic, and then **max_dt = 15**, **max_dt_unit = m**. 
  The module will in either scenario sweep 4 times per step in time.


* **step_N_stats**: The number of time-step increments to store delta variables for. Can be set to a single integer or a list of integers.
  For example if **step_N_stats = [1,2,4]**, and if **dt_size = 15** and **dt_unit = min**, 
  then the delta variables will be computed for **1*dt_size = 15**, **2*dt_size = 30** and **4*dt_size = 60** minute increments.


* **step_N_fields**: The number of time-step increments to store the fields for. Can be set to a single integer.
  For example if **step_N_fields = 4**, and if **dt_size = 15** and **dt_unit = min**, 
  then the fields will be stored for **4*dt_size = 60** minute increments. This can help to reduce the memory requirement for storing large fields.


* **BC_type**: The **BC_type** variable sets the boundary condition to apply. The allowable types are Robin, Neumann, Dirichlet, No flux and Robin experimental. 
  See more about how the different boundary conditions behave :doc:`here <./boundary_conditions>`.


* **convection_type / diffusion_type**: The numerical schemes can be set by the variables **convection_type** and **diffusion_type**. 
  The allowable types are listed in the config file. 
  For more information about what each scheme does numerically, see the `FiPy documentation <https://www.ctcms.nist.gov/fipy/documentation/numerical/index.html>`_. 


* **D**: The diffusion rate can simply be set by letting **D** be a non-negative float.
  As of now, if D == 0.0, then the diffusion term will not be added to the numerical scheme. 
  This increases computation speed slightly. But I am not sure if it alters the solution.


* **return_structured**: Keep this to **True** by default. It tells the module to structurize the outdata. If false, the outdata will stay in unstructured form, and
  any plots and post-processing modules may stop working correctly. This is mostly here to future-proof the module in case you want to keep the output in unstructured form.


.. note:: 

  The **dt_size** variable differs from the **max_dt** variable in the sweeps solution!
  The **dt_size** variables tells for what time-step intervals the AdvDiff solver should solve and store the solution.
  However, AdvDiff does not store the solution for every sweep. 
  So, sweeps are mostly used to maintain the accuracy of solutions whenever storing the solution for refined time-steps is not desirable or necessary.



[output]
--------
Here we can specify what sorts of outputs we want when the simulation terminates.

* **store_stats**: Do we want to store the maximum and mean values and the statistics for each simulation?


* **store_fields**: Do we want to store the field variables for every time step of the simulation? 
  This may require a lot of data if you run high fidelity simulations with many sources and long time series. The field variables are mostly useful for creating movies.


* **store_probes**: Do we want to store the time series data from the measurement probes? 
  This does not require as much data, as we will only look at the field variables at a set of positions for each time.


* **store_mass**: Do we want to store the total mass generated around each source with respect to time? 
  This is useful if you want to see the effect the boundaries have on the simulation.


* **store_uv**: Do we want to store the x- and y-components of velocity together with our field variables. 
  This will only take effect if store_fields is also True.


* **stat_variables**: This is a list of types of statistics to store.


* **reduce_point_tags**: If set to true, a netcdf file called point_tags_reduced.nc will be computed. 
  It includes the point tags for every pixel in the global grid.


* **downsample / downsampling_factor**: Storing field variables requires a lot of storage, 
  especially when running long timescale simulations utilizing meshes with many gridcells. 
  To help with memory management, we have an option to **downsample** the data. Downsample allows us to reduce the memory usage whenever store_fields = True. 
  The fields will be downsampled by a factor given by the **downsampling_factor** parameter.


* **global_stats**: If True, then stitch each of the local statistics into one global statistics file.


* **global_fields**: If True, then stitch each of the local fields into one global fields file.


* **cumulate_probes**: If True, take the cumulative sum of every probe contribution and output a cumulative probe file. 
  The cumulative probe file will represent the actual measurements of each probe given that every source of the system is accounted for.


* **reduce_point_tags**: If True, then reduce the point_tags.nc file onto the global grid.


* **glob_nx / glob_ny**: The respective number of discretizations along the x and y axes of the global grid.


* **weight_source**: If True, weight the global variables by location probability.


[paths]
--------
This section sets the directories, path names and file names for inputs and outputs. There is no need to change these parameters unless completely necessary.


[coordinates]
--------------
**datum**: Sets the datum for coordinate convertions betwen longitude-latitude coordinates and UTM.


[sources]
----------
This section is only relevant if **get_source_from == Sources**

* **source_file**: Sets the filename for the desired sources to run.



[riskmaps]
-----------
This section is only relevant if **get_source_from == Riskmap** or **get_source_from == Riskmap and Wells**.

* **risk_file**: Sets the filename for the desired riskmap to use to generate source coordinates.


* **out_path**: Sets the directory for where to save the generated source coordinates. This is only relevant if you use the **risk_to_sources** module as standalone.


* **threshold**: Sets the threshold for what is considered a source in the riskmap.


* **cluster**: Sets the clustering algorithm to use to generate sources. You can either use **dbscan** or **kmeans** clustering.


* **eps / min_samples**: These are dbscan specific attributes. See `dbscan documentation <https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html>`_ for more information.


* **n_clusters**: This is a kmeans specific attribute. It specifies how many clusters (sources) you want to build out of the riskmap.



[wells]
---------
This section is only relevant if **get_source_from == Wells** or **get_source_from == Riskmap and Wells**.

* **well_file**: What .nc or .csv file to load well locations from.



[velocity]
----------
Here we set what velocity file to read.

* **velocity_file**: Sets the filename for the desired velocity to use to run the simulation.


* **depth**: We encountered that some velocity files had an additional depth coordinate. 
  Here you can specify which depth index to use, or set **depth = average** to take an average velocity of all depths.


* **fill_type**: When generating local fields from the global velocity field, 
  you may encounter NAN values wherever there is missing velocity data, such as whenever the local fields extend beyond the limits of the global field. 
  If **fill_type = nearest extrapolation** or **fill_type = fill extrapolation**, then any NANs are filled with extrapolated values. 
  If **fill_type = zeroes**, then NANs are filled with zeroes. 
  If **fill_type = average**, then NANs are filled with the average value of the field. 
  If **fill_type = ignore**, then the simulation around the source in question will be ignored and terminated if any NANs are encountered.


* **allow_synthetic_translation**: If set to **False**, then sources which land outside the limits of the global field will be ignored and terminated. 
  This should be the default option for any realistic simulation as it does not impose any stretching of the velocity data to fit onto the source coordinates. 
  If set to **True**, then we will allow the module to stretch the global field to fit onto the source coordinates, this will allow us to run velocities and sources with completely different coordinate data.


* **allow_synthetic_scaling**: Same as above, but also allows for scaling of the velocity data to fit over all sources. It will only work in tandem with **allow_synthetic_translation**.



[probes]
---------
Here we set the probe file to use.

* **probe_file**: This is optional, and if you do not want to run any probes, set **probe_file = None**. If the probe coordinates do not overlap with the 
  global grid, then probes will be ignored.



[grid]
-------
Here we can specify the parameters for the meshes used during computation and storage.

* **type**: What mesh type should be used for the computational grid. The available options are listed in the config file. 
  The triangle and exponential options allow for grids which become more and more refined towards the source.


* **type_out**: What mesh type should be used for the output data. 
  This parameter does not allow for triangle grids as of now as xarray offers no simple way for plotting data in unstructured grids, but I may add this in the future.


* **Lx / Ly**: These parameters specify the length of local grids along the x- and y-axis.


* **extrap_nans**: In rare instances when interpolating from a rectilinear mesh to a triangular mesh, NANs will appear along a few of the cells along the boundaries. 
  If **True**, then we will allow for nearest interpolation to patch these holes. If **False**, then any NANs will be filled with zeroes instead. 
  It is highly recommended to keep this to True.


* **nx / ny**: This specifies the number of discretizations the x- and y-axis should be divided into for the rectilinear grids.


* **maxvol / minvol**: This specifies maximum size of the largest triangle and the minimum size of the smallest triangle in the triangle mesh. 
  These values are approximate, and the final values may be off by about 10%.

* **reduction_factor**: How strongly the gridcells decrease in size towards the center when using an exp grid.

* **power**: How strongly the gridcells decrease in size towards the center when using a triangular grid.
             Power = 0.0 --> cell volumes are constant = maxvol.
             Power = 0.5 --> cell volumes increase proportionally to the square root of the distance from the center from minvol to maxvol.
             Power = 1.0 --> cell volumes increase linearily by the distance from the center from minvol to maxvol.
             Power = 2.0 --> cell volumes increase proportionally to the square of the distance from the center from minvol to maxvol.




[visualizer]
-------------
This section specifies the settings for the visualizer module in visualize.py.

* **plot**: When **True**, visualize the results on completion of the simulations in main.py. If **False**, you will have to run python visualizer.py separately.


* **robust**: Sets robust colorbar scaling. This is robust to outliers and points with abnormally large concentration. IF set to False, then plots will typically be hard to see.


* **levels**: Number of levels in contour field.


* **logscale**: Whether or not to plot colorbar on a logscale. (Is fairly buggy).


* **weight_source**: Whether or not to weight the source contributions with respect to location probability. This is independent of weight_source in output section.


* **plot_source_loc**: Whether or not to plot source locations ontop of visualizations as circles. The size of the circles is proportional to location probability.
  The sources are numbered from 0 to N.


* **plot_probe_loc**: Whether or not to plot probe locations optop of visualizations as crosses. The probes are numbered from 0 to N^*.


* **fontsize**: The fontsize of source and probe numbers.


* **markersize**: The markersizes of source and probe markers.


* **filetype**: What image filetype to store figures as.


* **plot_stats**: Plot the statistics?


* **plot_global_fiels**: Plot the global field as a movie?


* **plot_local_fields**: Plot each of the individual local fields as movies?


* **plot_probes**: Plot the probe measurement over time?


* **plot_probe_contrib**: Plot the contributions from each source to every probe measurement over time?


* **plot_mass**: Plot the total mass over time?


* **plot_mass_contrib**: Plot the mass contributions from each source over time?


* **plot_meshes**: Plot the meshes and mesh quality?


* **plot_local_grids**: Plot the local grid locations superimposed on the velocities? (Useful for debugging)


* **plot_map**: Plot a map containing the source locations and velocity locations in space with castlines shown. Requires the `Cartopy module <https://scitools.org.uk/cartopy/docs/latest/>`_ to work.


* **plot_global_stream**: Plots the stream fields for the global system. This will only work if **store_fields** and **store_uv** are both true.


* **plot_local_stream**: Plots the stream fields for each local system. This will only work if **store_fields** and **store_uv** are both true.


* **plot_riskmaps**: Plot the riskmaps for every file found in the riskmap directory on it's own figure?


* **plot_sources**: Plot the source locations for every file in the source directory on it's own figure?


