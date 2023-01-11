Interpolation with meshA_to_meshB
=======================================

Credits
---------

Special thanks for the linear interpolation function goes to `Jaime <https://stackoverflow.com/users/110026/jaime>`_ who designed an excellent linear interpolator 
and posted it to `stackoverflow <https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids>`_.
I used his framework to design a more flexible version with more features and improved ease of use.



What is meshA_to_meshB
---------------------------

One of the core submodules of **AdvDiff** is the **meshA_to_meshB** interpolator module (see :doc:`interpolator.py <./interpolator_module>`). 
It is a general purpose linear- and nearest-neighbor interpolator that can interpolate data defined on any arbitrary set of points (meshA) onto any other arbitrary set of points (meshB). 

Interpolation is a central tool for enabling AdvDiff's flexibility. Data can come in many forms and be defined on many different types of structured and unstructured grids.
For example, interpolation is needed when transferring velocity data from the raw netcdf files onto the unstructured computational FiPy grid. 
The same goes for when data variables on the FiPy grid is transferred onto a structured output grid which can then be stored as a netcdf file.

As the user can choose what grid-type to utilize during the simulation, 
it is important that AdvDiff can accurately anc quickly interpolate data from one arbitrary grid to any other.


Why meshA_to_meshB
------------------------

Why design a custom interpolator when Scipy's `scipy.interpolate.griddata() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html#scipy-interpolate-griddata>`_ 
module does the exact same thing? 
There are two main reasons:

1. The Scipy griddata() module is VERY slow when interpolating between two meshes that both contain many points.
   For use with AdvDiff, it is common to interpolate between meshes containing 10.000 points or more. 
   Scipy griddata() is slow as the module reinitializes the Delaunay triangularization for both meshA and meshB for every call. 
   This is computationally expensive and wasteful, as we are interpolating between the same two static meshes for every time-step and thus there is no need to 
   recompute the Delaunay triangularization.

2. The Scipy griddata() module cannot extrapolate outside the convex hull of meshA. 
   This leads to several problems. Whenever points on meshB lie outside the convex hull of meshA, griddata() gives NAN values. 
   FiPY does not know what to do with NAN values and it will crash if any of it's variables contain any NANs.
   This issue can be overcome by filling the NANs with some value, say 0.0, but that has it's own set of problems. Another more subtle issue is that 
   unless the border gridcells of meshA and meshB align in a favourable way, the borders may be filled with NANs also. 
   This specific problem can be seen on the figure in the `griddata() documentation <https://docs.scipy.org/doc/scipy/_images/scipy-interpolate-griddata-1.png>`_,
   notice the white border around the data for linear and cubic interpolation? This problem in particular made using griddata() unreliable.


The **meshA_to_meshB** module solves both of these problems. The interpolator only computes the Delaunay triangularization once on initialization of the interpolator object. 
Then for each call of **meshA_to_meshB.interpolate()**, it will use the precomputed triangularizations to perform linear interpolation. 
Once linear interpolation is performed, the module checks the result for any NANs. If any are found, then it will perform 
nearest interpolation on only those points. Nearest interpolation can be used outside the convex hull of meshA, and as such we end up with a result with no 'holes'. 



meshA_to_meshB** vs. scipy.interpolate.griddata()
--------------------------------------------------------


I have created an artificial scenario to showcase the power of **meshA_to_meshB**. 
We want to interpolate data defined on a random set of points onto a 100x100 grid. We compare interpolation accuracy and efficiency.
The results are shown in the next two sections.



Accuracy comparison
~~~~~~~~~~~~~~~~~~~~

Figure on the left shows randomly placed points in the plane, each point is associated with a value [0,1] depicted by the hue. 
The red boundary line is the convex hull of the set of points.
The middle-left figure shows the resultant linear interpolation of the data using **scipy.interpolate.griddata()** ,
while the middle-right figure shows the resultant linear interpolation using **meshA_to_meshB**. 
The scipy interpolator cannot extrapolate values outside the convex hull.
I have set **extrap_nans = False** for a fair apples-to-apples comparison.


To ensure that **meshA_to_meshB** is indeed accurate, I also compute **np.isclose()** for every point on the interpolated grid, the result is shown on the rightmost figure.
Yellow means the two solutions are the same while purple means there is a discrepancy.
Here we see that the two interpolators do indeed produce the same results for every point. 
As such, we are confident that **meshA_to_meshB** produces the correct linear interpolation.

 .. figure:: .//images//meshA_to_meshB_vs_scipy//mA2mBvsScipy_no_extrapolation.png
   :width: 90 %
   :align: center

If however, I set  **extrap_nans = True** in **meshA_to_meshB** we get:

 .. figure:: .//images//meshA_to_meshB_vs_scipy//mAmBvsScipy_with_extrapolation.png
   :width: 90 %
   :align: center

Notice that **meshA_to_meshB** has now produced values outside of the convex hull. This is important for our case, because FiPY does not know what to do with NANs.
If local grids ever poke outside the global velocity field, we want to extrapolate values such that we can still run simulations if we want to.
On the rightmost figure, we see that the two interpolators now differ, but only outside the convex hull of the datapoints. 
This is because **scipy-interpolate-griddata()** only produce NANs outside the convex hull of points, while **meshA_to_meshB** utilizes 
nearest interpolation to compute values outside the convex hull.



Speed comparison
~~~~~~~~~~~~~~~~~~~

The main reason why **meshA_to_meshB** is preferred, is due to it's incredible speed and efficient scaling. Thanks to `Jamie's <https://stackoverflow.com/users/110026/jaime>`_ linear interpolator, we can interpolate
data between two grids over several time-steps very rapidly. The following figure shows the average time taken to perform one interpolation with **scipy-interpolate-griddata()**,
**meshA_to_meshB without extrapolation** and **meshA_to_meshB with extrapolation** over number of datapoints. We see that there is a clear speed benefit for using **meshA_to_meshB** in both cases.
Moreover, the time average time taken for performing one interpolation with **meshA_to_meshB** does not seem to increase as we increase the number of datapoints, whereas with **scipy-interpolate-griddata()**,
the time increases exponentially.

 .. figure:: .//images//meshA_to_meshB_vs_scipy//Times_comparison.png
   :width: 90 %
   :align: center

 .. figure:: .//images//meshA_to_meshB_vs_scipy//Speed_Comparison.png
   :width: 30 %
   :align: center

However, the slowest part of the **meshA_to_meshB** interpolator is indeed the nearest interpolator function. This explains the constant width gap between **meshA_to_meshB without extrapolation** and **meshA_to_meshB with extrapolation**.
So when performing interpolation, it is beneficial to try to set up scenario for which there is as little extrapolation needed as possible. 
But, the safest bet is still to set **extrap_nans=True** because, if no NANs are actually found, it will not waste time trying to extrapolate.
I have made an example notebook that you can play with and confirm my findings. You can get it from my GitHub `here (Work In Progress) <https://github.com/KetilFIversen/meshA_to_meshB-Interpolator>`_. 


.. note::

  * The **meshA_to_meshB** interpolator is only faster than the **scipy.interpolate.griddata()** for when interpolating between large grids. If the grids are relatively coarse (less than 500 points),
    then there is no benefit to **meshA_to_meshB**, aside from the extrapolation property

  * **meshA_to_meshB** is only beneficial for when interpolating between the same two meshes multiple times. If you are going to do an interpolation only once, then **scipy.interpolate.griddata()** is better.
    Also, if you want to do interpolation between two changing meshes, for example when implementing adaptively refined meshes, then do not use **meshA_to_meshB**.

  * The nearest interpolator function is very slow, as such extrapolation is slow. When you want to perform rapid interpolation, try to make sure the convex hull of meshB fully contained
    within the convex hull of meshA such that extrapolation is either not needed or only needed for a small set of points.

  * The **extrap_from_A** boolean determines if nearest extrapolation should be in respect to data defined on meshA or the interpolated data on meshB. The best results depends on the situation.
    Setting **extrap_from_A = True** is best for when the two meshes overlap well and contain the same order of magnitude number of points. However, if meshA is very coarse, and meshB is very fine,
    then setting **extrap_from_A = False** is usually the best. 