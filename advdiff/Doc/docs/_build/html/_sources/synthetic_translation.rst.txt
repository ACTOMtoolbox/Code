Synthetic translation
===========================================


Synthetic translation
----------------------

The **allow_synthetic_translation** variable in :doc:`setup.ini <./setup>` allows us to run simulations on almost any kind of velocity file, 
regardless of the coordinates provided.
For example, let us say that our riskmap/sources are taken from the Gulf of Mexico. 
Ideally, we should provide the AdvDiff module velocity data from the exact same location. 
But sometimes we may not have access to this data. If we provided data from another location, 
then the UTM coordinates of the riskmap and the velocity files will likely not overlap at all.
And as such, we would not be able to extract any local velocity grids around each of the sources from the global velocity data. 


Synthetic translation allows us to utilize velocity data at a different location regardless of where.
If **allow_synthetic_translation = True** in :doc:`setup.ini <./setup>`, AND if the
UTM coordinates of the global velocity data does not overlap with the source locations, 
then the UTM coordinates of the velocity data will be translated and fit onto the sources. As such,
a simulation can be run without any issue. 
A warning will be printed to the console if this happens, so that the user is aware of this synthetic change.


If **allow_synthetic_translation = False** in :doc:`setup.ini <./setup>` however, then no translation of the velocity data will take place.
Any source which does not overlap with the velocity data, will be ignored and not run.

The user also has the ability to set if synthetic scaling should take place in addition to translation. This can allow velocities defined in a smaller regions
to fit over a large area.

Ultimately, if you want to run realistic scenatio, make sure that the UTM coordinates of the sources and the velocity data do overlap.

.. note::

    If the velocities actually do overlap with the sources just fine, then no translation will take place as it is not needed.
    This is regardless of whether or not synthetic translation is enabled or not. So synthetic translation should not affect 
    your results if velocities and sources do align properly.
