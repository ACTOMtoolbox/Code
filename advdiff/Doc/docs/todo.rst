Todo  
=====


Improving performance
----------------------

The default run-time configuration of FiPy is based on the available packages that are installed on the machine.
To override the default behaviour, command-line flags can be used. The --inline flag is causes 
many mathematical operations to be computed through C instead of Python, which can improve the performance.
However, for this flag to work, the Weaver package has to be installed. Unfortunately, Weaver is not available on Python 3.x,
and thus I have been unable to get --inline to work. 

In addition to this, FiPy can use several different solvers. The solvers include “petsc”, “scipy”, “pysparse”, and “trilinos”.
This version of AdvDiff runs with the scipy solver, however it is also the slowest solver of the bunch.
Unfortunately, I have been unable to get any of the other solvers to work with version 3.8 of Python.

If you are able to install Weaver and/or get any of the other solvers to work with Python 3.8, then you may gain some performance improvements.
Alternatively, you could look into porting AdvDiff over to another version of Python that supports the other solvers.


Uniformity of inputs and outputs
-------------------------------------

One possible issue with AdvDiff that shoule be rectified, is uniformity between input and output data.
FiPy solvers work on unstructured data, however xarray works much better with structured data.
You can store unstructured data with xarray just fine, but much of the nice xarray functionality will not work as well.
Such as easily being able to generate plots, interpolate multiple grids together, etc.
As such, we have constructed the AdvDiff module to generate solutions on an unstructured set of data, to then be interpolated onto a structured grid
which are then stored in xarray dataset. In short:

* User gives the AdvDiff solver input data in an unstructured form :math:`F(t,num)` where **num** is the gridcell id of an unstructured set of data.

* FiPy generates solution of the form :math:`C(t,num)`.

* The unstructured solution is then interpolated onto a structured grid of the form :math:`C(t,x,y)`.

* The structured solution is then stored to a netcdf file via xarray.


We see that it may be confusing to the user that input data needs to be given in terms of an unstructured set of data, while the output will always be in terms of 
structured dataset. To fix this, we have developed an alternate approach. In short:

* User gives the AdvDiff solver input data in an unstructured form :math:`F(t,num)` where **num** is the gridcell id of an unstructured set of data.

* FiPy generates solution of the form :math:`C(t,num)`.

* The unstructured solution is then stored to a netcdf file via xarray.

* A post-processing module is then used to interpolate the data :math:`C(t,num)` and reformat the unstructured xarray file to represent the solution in terms of :math:`C(t,x,y)`.


The former solution has the benefit of only needing to parse the data and interpolate once during runtime, while the latter solution
has the benefit of unifying input and output formats. Either solution works well, but it is up to the discretion of the future developer to decide what method
is best.



Dealing with coastlines, islands and other landmasses
---------------------------------------------------------

The AdvDiff module works best on the open ocean or away from islands or coastlines that could affect the solution. 
This is because there is no special consideration of landmasses in the solver.
At the current moment, if the velocity data contains any NANs (typical for locations where there is land), then these locations will be filled with 
zero velocities. This is to simulate a no flux boundary on land. This is obviously only a rough approximation, but it was the simplest to implement.
As such, if a significant amount of tracer moves rapidly onto the coastline, it may get stuck there due to the zero velocity at these locations.
If you want to change this behaviour, check the **fillna()** method in the **treat_xy_velocity()** and **treat_sparse_velocity()** functions in :doc:`file_loader.py module <./file_loader>`.

One option for the future developer is to generate a mask based on locations in space that contain land, then using this mask to set no-flux, Diriclet or Robin 
boundary conditions on this mask. This should be possible by using the **where** keyword in FiPy when setting boundary conditions or defining variables, but this has never been tested.


Investigating the Diriclet boundary condition
---------------------------------------------------

The Dirichlet boundary condition behaves strangely. This should be investigated as it could reveal some underlying problem.