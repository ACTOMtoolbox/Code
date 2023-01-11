Coordinate transformations
===========================

As seen in the :ref:`inputs section <input>`, the AdvDiff module can take two different coordinate systems as inputs.

* UTM coordinates in terms of meters Easting (x-axis) and meters Northing (y-axis)

* WGS84 coordinates in terms of longitude (x-axis) and latitudes (y-axis).

The AdvDiff module computes its solutions on flat 2D grids mapped on UTM coordinates.
And thus, all inputs are converted to UTM before the simulation can run. 
One problem with UTM, is that they do not offer a unique mapping to longitude-latitude coordinates. 
One set of UTM coordinates can point to several different locations on the Earth, depending on what UTM zone is implied.
To perform an accurate transformation, a UTM zone or a rough area of interest has to be given.
Unfortunately, the inputs we have been provided do not contain any UTM zones by default.
Therefore we have chosen to not include an option to transform the ouput data back in terms of longitude-latitude coordinates.

However, any set of longitude-latitude coordinates can map to a unique UTM zone. If we know what longitude-latitude coordinates we are expecting,
then a UTM zone can be extracted and hence we will be able to convert UTM to longitude-latitude coordinates.
The `pyproj package <https://pyproj4.github.io/pyproj/stable/#>`_, gives us the ability to automatically discover the correct UTM zone
given a bound for the longitude-latitude coordinates.
So, if the user defines sources, riskmaps, velocities, or probes in terms of longitude-latitude coordinates, then
the module will be able to convert UTM back into the correct longitude-latitude coordinates if needed. 
Therefore it may be beneficial to start defining inputs in terms of longitude-latitude coordinates instead of UTM as you will be able to convert UTM back
into the correct location.

To do this, have a look at the :doc:`coord_transform module <./coord_transform>`. It can transform longitude-latitude to UTM but also back, 
**if and only if** an area of interest is given. If any of the inputs given to AdvDiff are in terms of longitude-latitude coordinates during the simulation, 
an area of interest netcdf file will be generated. This file simply contains the longitude-latitude bounds of the area we are simulating. 
The bounds are enough information for `pyproj <https://pyproj4.github.io/pyproj/stable/#>`_, to recognize the correct UTM zone. This file is then stored in Indata/Coordinates/. 
Whenever coordinate transforms are called for, the coord_transform module will check if this file exists, if not then there will be a warning whenever UTM coordinates
are attempted to be transformed into longitude-latitude coordinates. 

