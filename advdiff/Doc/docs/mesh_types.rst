Different mesh-types
=====================

There are two main mesh-types used in AdvDiff; structured meshes and unstructured meshes.
Supported unstructured meshes include **triangular meshes** while structured meshes include **equidistant meshes**, **exponential meshes** and **segmented meshes**.

To score the mesh quality for the different available meshes, we have used the Normalized Shape Ratio (NSR).
We chose this measure as it is easy to generalize the score for both rectangular and triangular gridcells.

.. note::

  The user has the option to set a **grid_type** and a **grid_type_out** in :doc:`setup.ini <./setup>`.
  FiPy will use whatever grid set with **grid_type** as the computational grid, while the results will be interpolated onto the grid set with **grid_type_out**.
  This means you have the option to set what mesh to use compute the solution, and a seperate grid to store the results on.
  This is done as xarray is not very flexible with unstructured datasets. As such, it is best to always set **grid_type_out** to any of the rectilinear meshes only.
  Otherwise, you will not be able to plot the results as easily, and any of the post-processing modules may not give the correct output.



Normalized Shape Ratio
-----------------------

The Normalized Shape Ratio (more typically referred to as aspect ratio) is a measure to determine mesh quality, 
in other words how poor or how good a mesh is for use with numerical methods such as finite elements or finite volumes.
In short, when utilizing rectangular gridcells you typically want the gridcells to be as square as possible.
Similarly, when utilizing triangular meshes, it is preferred that each triangle is as equilateral as possible.
Keeping a consistent gridcell shape helps alleviate discretization errors. 

The NSR computes how similar a rectangle is to a perfect square and how similar a triangle is to a perfect equilateral triangle. This is done by computing the 
ratio between the inscribed circle and the excircle of the respective shapes.
For perfect squares and perfect equilateral triangles, the NSR is equal to 1.0, 
for any other rectangle or triangle, the NSR will lie somewhere in the range [0, 1], for which a value near 0 signifies a significantly squashed shape.

You can read more about the NSR and other mesh quality metrics here:
`link 1 <https://www.pygimli.org/_tutorials_auto/1_basics/plot_6-mesh-quality-inspection.html>`_,
`link 2 <https://cfdisrael.blog/2019/02/01/know-thy-mesh-mesh-quality-part-i/>`_,
`link 3 <https://help.solidworks.com/2021/English/SolidWorks/cworks/c_Mesh_Quality_Checks.htm>`_,
`link 4 <https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-c-mgnconcpartitionverify.htm>`_,
`link 5 <https://en.wikipedia.org/wiki/Types_of_mesh#Aspect_ratio>`_.



Unstructured meshes
---------------------

**Triangular mesh**:

The triangular mesh is the only unstructured mesh supported in AdvDiff. Delaunay triangularization is used to generate a mesh with high density in the center
and low density along the borders. The average sidelengths of each triangle decreases radially while approaching the center. 
This is done to keep accuracy near the source high while also being able to increase the size of local grids. 
To set the width and height of the mesh, set **Lx** and **Ly** respectively in :doc:`setup.ini <./setup>`.
To set the maximum size of the largest triangle and the minimum size of the smallest triangle, set **maxvol** and **minvol** respectively in :doc:`setup.ini <./setup>`.

The main benefits of the Triangular mesh are:

* The triangular mesh allows high refinement near the source and low refinement far away from the source, 
  for which the degree of refinement decreases uniformly and radially away from the source. 

* Almost all triangles in the mesh have a 'good' NSR score. In other words, the triangles are nearly all perfectly equilateral.
  And this lowers the amount of discretization errors. 

* To get a similar mesh refinement near the source as with an exponential grid, a lesser number of gridcells are required.

The main downsides of the Triangular mesh are:

* Currently, you do not have the option to store the results on triangular meshes as triangular meshes are inherently unstructured. 

* The triangular mesh takes much longer to initialize than the rectilinear meshes, 
  but overall the triangular mesh typically requires less gridcells for the same refinement of an equivalent rectilinear mesh.


 .. figure:: .//images//Mesh_types//Triangular_mesh.jpg
   :width: 90 %
   :align: center

The uppermost figure shows the distribution of gridcell volumes (areas). We see that most gridcells are small, while the rest of the space is filled with a smaller number
of large gridcells. The middle figure shows the distribution of NSR scores. Most gridcells score a 1.0, while a very small amount of gridcells score below a 0.9.
The lower-left figure shows the structure of the mesh. The gridcells are coloured randomly such they are more easily distinguished. The lower-right figure
shows the corresponding NSR score for each gridcell.



Structured meshes
------------------


**Equidistant mesh**:

As the name suggests, the equidistant mesh is a regular rectilinear mesh with equidistant grid spacing. 
To set the width and height of the mesh, set **Lx** and **Ly** respectively in :doc:`setup.ini <./setup>`.
To set the number or discretizations along each axis, set **nx** and **ny** respectively in :doc:`setup.ini <./setup>`.

The main benefits of the equidistant mesh are:

* Easy to interpret and understand.

* It is structured, so data can be stored in an intuitive way using xarray datasets in the form of f(t,x,y).

* Equidistant grids are usually preferred to avoid computational artefacts such as artifical diffusion or dispersion.

* As long as you properly scale such that **Lx** / **Ly** = **nx** / **ny**, then the NSR will be 1.0 everywhere.
  This is because every gridcell will be a perfect square. This is preferred to avoid artefacts.

The main downsides of the equidistant mesh are:

* If you want a highly refined mesh around the source, you are forced to have a higly refined mesh everywhere else also. 
  As such, you may waste a lot of computational energy on performing high density computation far away from the source.


 .. figure:: .//images//Mesh_types//Equidistant_mesh.jpg
   :width: 90 %
   :align: center

Here we see that every gridcell has the same size, and the NSR is equal to 1.0 everywhere.



**Exponential mesh**:

The exponential mesh is a rectilinear mesh for which the grid spacing decreases exponentially towards the center of the mesh where the source lies.
This is done to keep accuracy near the source high while also being able to increase the size of local grids. 
To set the width and height of the mesh, set **Lx** and **Ly** respectively in :doc:`setup.ini <./setup>`.
To set the number or discretizations along each axis, set **nx** and **ny** respectively in :doc:`setup.ini <./setup>`.
To set the degree of which grid spacing decreases, set **reduction_factor** in :doc:`setup.ini <./setup>`.

The main benefits of the exponantial mesh are:

* It is structured, so data can be stored in an intuitive way using xarray datasets in the form of f(t,x,y).

* You are able to have high refinement at the source and low refinement far away from the source. As such, you can increase the size of the meshes, while
  Still retaining accuracy near the source, while also not wasting computational resources on locations far away from the source.

The main downsides to the exponential mesh are:

* Due to how the fact that the exponential mesh is structured, you cannot avoid also having high refinement along the coordinate axes, as such wasting some computational resources.

* The mesh will have a poor NSR score. This is due to the problem above. The cells along the diagonals will be square, while all other cells will be squashed rectangles.


 .. figure:: .//images//Mesh_types//Exponential_mesh.jpg
   :width: 90 %
   :align: center

Here we see that the gridcells do indeed become smaller towards the center, however they become squashed along the coordinate axes. This leads to a poor NSR score
along the axes and a good NSR score only along the diagonals. This problem becomes a lot worse if **Lx / Ly =/= 1**. 


**Segmented mesh**

The segmented mesh is a rectilinear mesh for which the grid spacing is refined in the middle, then abruptly changes outside some distance away from the source.
This mesh is not recommended as the gridcells in the four corners will be very large.


 .. figure:: .//images//Mesh_types//Segmented_mesh.jpg
   :width: 90 %
   :align: center

Here we see that there are two main regions, one with high density and one with low density. Here we also run into the problem of squashed gridcells along the axes,
leading to four regions of very poor NSR score.
As such this mesh type is not recommended.
