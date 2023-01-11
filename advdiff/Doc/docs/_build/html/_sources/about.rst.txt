About AdvDiff
==============

Introduction
-------------

**AdvDiff** is a Python module for simulating the tracer transport problem on an ensemble of point sources. 
To this end, the module utilizes parallel computing capabilities of multi-core CPUs. 
Each core/thread of the CPU works with a single source at any given time, and each source is solved for on individual local grids.
As the advection-diffusion equation is a linear partial differential equation (PDE), 
we are able to assemble the aggregate solution by simply summing the contributions from each local solution together.
As such, we can construct a global solution which describes the behaviour of all sources acting together. 
This has several advantages over solving all sources together on the same global grid.
For one, it allows us to utilize parallelization as each source is solved separately.
Moreover, it gives us the ability to generate multiple possible scenarios from running the module only once.
As each source is stored separately on individual local grids, we are able to change the resulting global outcome
by weighting the contributions from each source independently.
Therefore, we are able to account for many possible leakage scenarios without having to redo simulations.

The tracer transport problem in two dimensions is modelled by the following PDE:

.. math:: 
   
   \frac{\partial C}{\partial t} + \mathbf{v} \cdot \nabla C = D \Delta C + F(\mathbf{x}), \quad \mathbf{x} \in \mathbb{R}^2, t > 0

For which :math:`C [\frac{kg}{m^2}]` is the tracer variable, :math:`\mathbf{v} [\frac{m}{s}]` 
is the velocity field which may be space-time dependent, :math:`D [\frac{m^2}{s}]` is the diffusion rate
and :math:`F [\frac{kg}{m^2 s}]` is the source function of the form:

.. math:: 

   F(\mathbf{x}) = \sum_{i=1}^{N_s} \delta(\mathbf{x}-\mathbf{x}_i) / V_i

Here, :math:`\delta(\mathbf{x})` is the Kronecker-delta function and :math:`\mathbf{x}_i` are the source locations in the plane.
Due to discretizations, we assume that the source only occupies a single gridcell of the computational mesh.
For each source to have a total flux of 1, we scale each source by the volume of the source cell :math:`V_i`, where
the gridcell nearest :math:`\mathbf{x}_i` is the source cell. 

The solution of the tracer transport problem is computed using the finite-volumes python package FiPy.
See the `FiPy website <https://www.ctcms.nist.gov/fipy/>`_ and the `FiPy manual <https://www.ctcms.nist.gov/fipy/download/fipy-3.1.3.pdf>`_ for more information.

.. note::

   The AdvDiff module works best on the open ocean, away from islands or coastlines that could affect the solution. 
   This is because there is no special consideration of landmasses in the solver.
   At the current moment, if the velocity data contains any NANs (typical for locations where there is land), then these locations will be filled with 
   zero velocities. This is to simulate a no flux boundary on land. This is obviously only a rough approximation, but it was the simplest to implement.
   As such, if a plume of tracer moves rapidly onto the coastline, it may get 'stuck' there due to the zero velocity at these locations.
   If you want to change this behaviour, check the **fillna()** method in the **treat_structured_velocity()** and **treat_unstructured_velocity()** 
   functions in :doc:`file_loader.py module <./file_loader>`.

   One option for the future, is to generate a mask based on locations in space that contain land, then using this mask to set no-flux, Diriclet or Robin 
   boundary conditions on this mask in addition to the regular boundary condition. This should be possible by using the **where** keyword in 
   FiPy when setting boundary conditions or defining variables. However, this has never been tested.



.. _input:

Input
------

To run the tracer transport module, you require:

* A config file which specifies the module parameters. See :doc:`setup.ini <./setup>`. 
* A netcdf file containing a riskmap, source coordinates or well locations.
* A netcdf file containing velocity data.
* A netcdf file containing coordinates for measurement probes (Optional).

The repository comes pre-packaged with several netcdf files and a template for the config file, which the user can tweak and work with.

If you just want to run a test scenario, you can keep everything in the config file at default for now, and head over to :doc:`quick-start <./usage>`.

.. note::

   It is important that if you design your own riskmap, source, well, velocity or probe data, that it follows the conventions below **EXACTLY**. 
   Many of the submodules require that the coordinate-, dimension- and variable-names align with the expected conventions.



Sources & Riskmaps & Wells
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sources, wells and riskmaps are needed to determine where to place the point fluxes.
Any of these files can be used to specify where in space the sources are located. 

* **Sources**: Provides the module the exact coordinates for where the sources should be located and also gives a corresponding **location_probability**.

* **Riskmaps**: Provides the module potential zones for where sources could be located. 
  A clustering algorithm is then used to determine exact source coordinates and corresponding location probabilities.

* **Wells**: Provides the module the exact coordinates for where the sources should be located and also gives a corresponding **location_probability**.
  This is similar to the source file. However, you can also load wells from a CSV fileformat.


To determine what files to read from, set the **get_source_from** variable in :doc:`setup.ini <./setup>` to either 
**Sources**, **Riskmap**, **Wells**, **Riskmap and Wells** or **Riskmap and Sources**.
Regardless of which files you use, they need to contain coordinates for each relevant location 
and a corresponding **location_probability** variable which determines probabilities of each location.

.. note::

   The **location_probability** does not affect the flux of the sources during the simulation. All sources are run with a flux of 1. 
   However, **location_probability** can be used to scale the results in post-processing (which would be equivalent to scaling the flux due to linearity).
   If you do not want to scale the solutions, set **weight_source** to False in the ouptut section of :doc:`setup.ini <./setup>`.

A valid riskmap should look like one of the following:

.. figure:: .//images//Inputs//riskmap.png
   :width: 50 %
   :align: center

   Example of a Riskmap netcdf file in UTM coordinates.

.. figure:: .//images//Inputs//riskmap_lonlat.png
   :width: 50 %
   :align: center

   Example of a Riskmap netcdf file in WGS84 coordinates.

A valid source file should look like one of the following:

.. figure:: .//images//Inputs//sources.png
   :width: 50 %
   :align: center

   Example of a Source netcdf file in UTM coordinates.

.. figure:: .//images//Inputs//sources_lonlat.png
   :width: 50 %
   :align: center

   Example of a Source netcdf file in WGS84 coordinates.

A valid wells file should look like one of the following:

.. figure:: .//images//Inputs//wells_netcdf.png
   :width: 50 %
   :align: center

   Example of a Wells netcdf file in UTM coordinates.

.. figure:: .//images//Inputs//wells_netcdf_lonlat.png
   :width: 50 %
   :align: center

   Example of a Wells netcdf file in WGS84 coordinates.

.. note:: 

   Ensure that any x-y coordinates in any of the source, riskmap, or well files are in UTM coordinates (meters Easting and meters Northing)
   and that any lon-lat coordinates are in WGS84 coordinates.


In addition, a csv file can be read to determine well locations.
The csv files we have encountered set the coordinates of these wells with coordinates with names ['Long83','Lat83'] and ['Long27','Lat27']. 
At the moment, the csv reader is hard-coded to deal with this naming scheme only. 
If you want to change this behaviour, look at the **get_wells_CSV()** located in the :doc:`file_loader <./file_loader>` tool.
Moreover, the csv files we have encountered do not include a 'location_probability', as such this is artificially set to 1 by the **lp_well** variable in **get_wells_CSV()**.

A valid well csv file should look like the following:

.. figure:: .//images//Inputs//Wells.png
   :width: 50 %
   :align: center

   Example of a Wells csv file in WGS84 coordinates.





Velocities
~~~~~~~~~~~

Velocities are used to determine how the tracer flows in the domain with time. 

.. note::
   For best results try to ensure that the area of where the velocities are defined, covers all the source locations
   with some margin to give room for the local grids. If the module discovers that there are some sources which lie outside the region where the velocities are defined, 
   it will display a warning and proceed depending on what the user sets in the config file (see :doc:`synthetic translation <./synthetic_translation>` and :doc:`setup.ini <./setup>`).

.. note::

   Ensure that the velocities are defined in meters per second :math:`\frac{m}{s}`. 
   The module does not check for this and it will not do any scaling to account for any other possible velocity units.

There are currently two different velocity file formats that are valid inputs to our module:

* Structured velocity time-series data defined on a rectilinear grid (:math:`f(t,x,y))`.

* Unstructured velocity time-series data defined on a set of unstructured points in space (:math:`f(t,node)`).

Structured velocities
~~~~~~~~~~~~~~~~~~~~~~~

The structured velocity data may look something like:

.. figure:: .//images//Inputs//velocity_dense.png
   :width: 50 %
   :align: center

   Example of a structured velocity netcdf file.
 
The velocity is defined for every point :math:`(x,y)` on a rectilinear grid and for each time-step. 

Unstructured velocities
~~~~~~~~~~~~~~~~~~~~~~~~

The module also supports unstructured velocity data. The following example is taken from the Gulf of Mexico dataset:

.. figure:: .//images//Inputs//velocity_sparse.png
   :width: 50 %
   :align: center

   Example of an unstructured velocity netcdf file.
 
The velocity is defined only for a handful of locations :math:`(nodes)` in space.
It has the benefit of not requiring as much data in memory, however the module needs to interpolate this data to 'fill the empty gaps'.
And as such, we lose out on the high fidelity variations in the velocity field.

.. note::

   * Both types of velocity files accept a depth dimension. The depth slice is selected by the user in the config file see :doc:`setup.ini <./setup>`.

   * The time dimension can be named either 'time' or 'ocean_time'. 
   
   * Both types of files can be in terms of either UTM or WGS84 coordinates.

.. note::

   The time-stepping procedure that AdvDiff uses, is independent of the time-step length in the velocity data. 
   For example, if the velocity file steps forward in time with 1 hour increments, the module will be able to step forward in time with say 15 minute increments
   or vice versa. The AdvDiff module will interpolate the velocity data in time if needed.
   To set a time-step length, see **dt_size** and **dt_unit** in :doc:`setup.ini <./setup>`.





Probes 
~~~~~~~

The probes are an optional part of the module. Probes lets the user specify measurement locations for where to 'probe' the solution for every timestep. 
If you want to run probes, you should provide a file of the following form:

.. figure:: .//images//Inputs//probes.png
   :width: 50 %
   :align: center

   Example of a Probes netcdf file in UTM coordinates.

.. figure:: .//images//Inputs//probes_lonlat.png
   :width: 50 %
   :align: center

   Example of a Probes netcdf file in WGS84 coordinates.



Output 
-------
Currently, the module has several different outputs which serve different purposes, **Statistics**, **Fields**, **Probes**, **Mass** and **point tags**.
The user has the ability to enable or disable storing these outputs during the simulation. See :doc:`setup.ini <./setup>`.

statistics
~~~~~~~~~~~

The statistics store the maximum, mean, variance etc of the tracer field variable :math:`C(t,x,y)`. 
They can be used to show aggregation zones and the overall behaviour of the system over longer time-scales without having to store the solution for each time-step.
This allows us to observe long time-scale behaviour while preserving memory.
Currently, the AdvDiff module can output the following statistics variables:

* **max**:             
   .. math::
     
      \max_{t} {C(t)}

* **mean**:           
   .. math:: 
      
      \text{mean}_{t} {C(t)}

* **var**:             
   .. math::
      
      \text{var}_{t} {C(t)} = \text{mean}_{t} {C(t)^2} - [\text{mean}_{t} {C(t)}]^2

* **delta_max**:       
   .. math::

      \max_{t} {[C(t+dt)-C(t)]}

* **delta_min**:       
   .. math::

      \min_{t} {[C(t+dt)-C(t)]} 

* **delta_mean**:      
   .. math::
      
      \text{mean}_{t} {[C(t+dt)-C(t)]}

* **delta_var**:       
   .. math::
      
      \text{var}_{t} {[C(t+dt)-C(t)]} = \text{mean}_{t} {[C(t+dt)-C(t)]^2} - [\text{mean}_{t} {[C(t+dt)-C(t)]}]^2

* **delta_max_abs**:   
   .. math::
      
      \max_{t} {|C(t+dt)-C(t)|}

* **delta_mean_abs**:  
   .. math::
      
      \text{mean}_{t} {|C(t+dt)-C(t)|}

The user has the ability to pick and choose which of these variables to solve for and store in the output. See :doc:`setup.ini <./setup>`.
Typical statistics outputs are shown below:


.. image:: .//images//Output//results//results_max.jpg
   :width: 33 %
.. image:: .//images//Output//results//results_mean.jpg
   :width: 33 %
.. image:: .//images//Output//results//results_var.jpg
   :width: 33 %

.. image:: .//images//Output//delta_C//statistics_delta_max15.jpg
   :width: 33 %
.. image:: .//images//Output//delta_C//statistics_delta_max30.jpg
   :width: 33 %
.. image:: .//images//Output//delta_C//statistics_delta_max60.jpg
   :width: 33 %

.. note::

   * For the delta variables, the user has the ability to change the size of :math:`dt` to any integer multiple of the time-step length **dt_size**. 
     See :doc:`setup.ini <./setup>`.

   *  See the clear rectangular shaped boundaries created by the tracer? This is caused by the fact that the local grids are a subset of the global grid. 
      Any tracer that moves beyond the boundary of the local grid is no longer accounted for in the simulation, as such we get clearly defined borders around each source. 
      A solution to this issue is to just make the local grids large enough such that the tracer cannot travel beyond the boundaries. 
   
   *  The white circles in the plot show each of the locations of the sources, they are numbered from :math:`0` to :math:`N_{sources}`. 
      The size of each individual circle is related to the location_probability of the respecive source.
      In this specific example, we can see that source 3 has the largest circle, and as such it dominates the behaviour of the system, while source 6 has one of the smallest circles 
      and thus contributes an indistinguishable amount to the overall solution.
   
   *  The white crosses in the plot show the locations of where the probes are placed. They are numbered from :math:`0^*` to :math:`N_{probes}^*`.
   
   *  Plotting both source locations and probe locations can make the plots messy. As such, :doc:`setup.ini <./setup>` includes an option to include or exclude them from plots.



fields
~~~~~~~

The fields output store the field-variables for the tracer :math:`C` and the horizontal and vertical components for velocity :math:`(u,v)` for every time-step. 
They can be used to visualize the behaviours of the tracer or velocity for all times using animations. 
However, storing the fields is quite memory intensive, as fields store several solution variables for thousands of grid-cells 
over potentially hundreds of time-steps.
Typical results are shown below:

.. image:: .//images//Output//fields//Cmovie.gif
   :width: 50 %
   :align: center



probes
~~~~~~~

The probes output store the field-variables for the tracer :math:`C` and the horizontal and vertical components for velocity :math:`(u,v)` for every time-step, 
at a small set of points in space. 
They can be used to reduce the amount of memory needed in comparison to storing **fields**.

.. image:: .//images//Output//probes//probe_total_C.jpg
   :width: 50 %
   :align: center

.. note::

   * The probes are numbered here from :math:`0` to :math:`N_{probes}` as in the plots shown above, but dropping the * superscript. 

   * Each timeseries is what each probe would have measured if placed at their respective locations in the global grid.


mass
~~~~~

The mass output stores the total mass generated by the source in the local grid per time. 
It also stores the expected theoretical mass, based on the total flux of the source and time since simulation start.
Any deviation of the actual mass versus the theoretical mass means that boundary effects are affecting the solution.
This can give you a hint to whether or not they should try to increase the size of the local grids to improve 
the accuracy of the simulations.

.. image:: .//images//Output//mass//mass_total.jpg
   :width: 50 %
   :align: center

.. note::

   * If the theoretical mass > actual mass, then mass is escaping through the boundaries and can no longer be accounted for in the simulation.

   * If the theoretical mass < actual mass, then mass is being generated by the boundaries due to the Neumann boundary conditions, 
     and the results should be taken with a grain of salt.




Pipeline
---------
This following section gives a rough explanation for how the data is treated and how the AdvDiff module arrives to its solution.

Let us say we want to solve for a system containing 4 sources located at the Gulf of Mexico. 


Loading source locations
~~~~~~~~~~~~~~~~~~~~~~~~~~

First, the riskmap is analyzed by the :doc:`risk_to_positions.py <./risk_to_positions>` module. Here, a clustering algorithm is used to select all
all points in the riskmap which satisfy a threshold, and organizes them in distinct clusters. 

The resultant output is a set of coordinates with an accompanying **location_probability** variable. 
The location_probability determines the 'strength' of the source cluster.
The following figure displays a riskmap with each source discovered by the clustering algorithm given by the black circles. 
The sizes of each circle are corrolated to the **location_probability** variable, and it determines how we weight the results at the end.
This specific example shows the resulting clustering using kmeans with n=4. Using dbscan instead would result in a different set of sources. 
The desired clustering algorithm can be set in the config file, see :doc:`setup.ini <./setup>`.

.. figure:: .//images//riskmap_and_sources.png
   :width: 75 %
   :align: center

   Source locations superimposed on the riskmap.


Pre-processing velocity data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the module has aquired source locations, it begins pre-processing the velocity data.
The following figure shows the global u-velocity provided to us by the GOM velocity file **GOM-1-1993-01.nc**.

.. figure:: .//images//velocity_raw.png
   :width: 75 %
   :align: center

   Raw GOM velocity file.


Notice how the coordinates are in in terms of latitudes-longitudes values.
Our module solves for the tracer on x-y coordinates in terms of UTM coordinates only. 
As such the latitude-longitude coordinates are transformed into UTM coordinates. 
This is done by using the `pyproj package <https://pyproj4.github.io/pyproj/stable/#>`_.

.. note::

   The red boundary is the convex hull of all points in the dataset. 
   It shows the boundary of the region where we are able to interpolate
   data using a linear interpolator. Outside this region, we have no other option but to extrapolate using a nearest neighbor interpolator.
   This causes a degredation of accuracy outside the convex hull boundary. 
   This is not usually a problem unless the sources lie very close the convex hull boundary.
   The best solution is to either let the local grids be so small that no part of the local grid extends beyond the convex hull, 
   or let the input velocity span a large enough area to encompass all sources and all local grids such that linear interpolation 
   can be performed for all subsets covered by each of the local grids. 


Generating local velocity data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once source locations has been set, and the global velocity has been loaded, 
the module extracts the local velocity data around each individual source. They are then stored on seperate local grids.


.. figure:: .//images//sources_and_local_grids.png
   :width: 50 %
   :align: center

   Source locations and local grids superimposed on the u-velocity field.

Here we see the source locations as black circles, superimposed onto the global u-velocity.
The boxes show us the boundaries of the local grids.

.. note::

   Notice how the local grid around source 1 does not intersect the convex hull of the global velocity and as such all velocity values are interpolated linearly. 
   However, for all the other sources, a major section of the local grid needs to be interpolated using nearest interpolation...
 
As stated previously, each thread/core of the CPU will be given one local grid each to work with at any given time. 
So in essence, each thread/core will be individually working with the following four instances:

.. figure:: .//images//all_local_grids.png
   :width: 75 %
   :align: center

   Each source centered in a local grid superimposed on the u-velocity.


Generating the solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Once each of the 4 simulations are finished, we may end up with something like:

.. figure:: .//images//all_local_results.png
   :width: 75 %
   :align: center

   The results on each local grid at a certain time after simulation terminates.


Here we see the field solutions for the tracer :math:`C` at one instance in time. Notice that the amplitudes for each source are approximately the same.
This is because each source are assigned a total flux of 1 during computation.
We can assemble these solutions back together onto the global grid such that we end up with the following:

.. figure:: .//images//global_result.png
   :width: 75 %
   :align: center

   The resultant global solution.

Now notice how the amplitudes for each of the sources are different. 
This is because we construct the global solution using a weighted sum of contributions from each source.
The weights are given by **location_probability**.
