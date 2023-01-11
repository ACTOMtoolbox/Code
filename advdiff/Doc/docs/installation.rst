.. _installation:

Installation
=============

To install AdvDiff, simply clone the repository from the `AdvDiff repository <https://github.com/KetilFIversen/AdvDiff>`_.

To clone the repository, open up Git Bash and change the current working directory in Git Bash to the location where you want the cloned directory.
Then type the following:

.. code-block:: console

   (your-dir) $ git clone https://github.com/KetilFIversen/AdvDiff.git

Alternatively, you can use `GitHub Desktop <https://desktop.github.com/>`_ and clone using the GUI.

.. figure:: .//images//Clone_repo.png
   :width: 50 %
   :align: center

   How to clone a repository using GitHub Desktop.

The directory should now contain the necessary files to run the AdvDiff module.



Pre-requisite packages
------------------------

The AdvDiff module requires several python packages to be installed for it to run correctly. 
It is recommended to create a virtual environment for which the packages can be installed
in isolation without affecting other modules (see `Anaconda <https://www.anaconda.com/>`_).

.. note::

   The AdvDiff module will attempt to use the `Cartopy module <https://scitools.org.uk/cartopy/docs/latest/>`_ for plotting maps, if it is installed. 
   However Cartopy is not listed in either the :download:`requirements.txt file <../../requirements.txt>`or the :download:`environment.yml file <../../environment.yml>`,
   as it is not essential.
   This is because I could not get an install of Cartopy to work in docker. Cartopy is not essential for the module to run, only for making a map.
   If Cartopy is not installed, it will simply give a warning that Cartopy cannot be imported
   and move on...





Manual installation of packages with Anaconda:
---------------------------------------------------

The required packages
are listed in the :download:`requirements.txt file <../../requirements.txt>`.

.. include:: ..//..//requirements.txt
   :literal:


To manually install the packages, open up an instance of the Anaconda Prompt, and create a conda environment with python version 3.8 (IMPORTANT).

.. code-block:: console

   (conda) $ conda create -n myenvname python=3.8.12

Then activate the conda environment

.. code-block:: console

   (conda) $ conda activate myenvname

Then install the packages as listed in requirements.txt with the specific version numbers given.

.. code-block:: console

   (myenvname) $ conda -c conda-forge install package=x.xx.x


Once all packages are installed, you should now have a virtual environment ready to run AdvDiff. 




Installation of packages using environment.yml with Anaconda:
------------------------------------------------------------------

Alternatively, you can install the entire required environment from a .yml file. The :download:`environment.yml file <../../environment.yml>` file contains the following.

.. include:: ..//..//environment.yml
   :literal:

To install the packages, open up an instance of the Anaconda Prompt, naviage to the relevant directory where environment.yml is located, then type the following:

.. code-block:: console
   
   (conda) $ conda env create -f environment.yml

This should create a Conda environment with the name ACTOM-Advection which contain all the necessary pre-requisite packages.
If this does not work, try the manual installation instead.


.. note::

   The python version number is important. I have tried other versions but FiPY and/or gmsh have given me a lot of trouble.

   However, not all packages need to be installed with the exact version number given above. The version numbers are there to show a configuration that is known to work.
   The only packages that are important to get the correct version number for are FiPY and gmsh. This is because FiPy seems to be very picky with what version
   of gmsh it likes to cooperate with. Without the proper version of gmsh, you will not be able to generate unstructured triangle meshes, and there will be an 
   error if you try to do so.

   Packages such as tqdm, numpy, scipy, sklearn and so on may be more flexible but this has not been tested.
