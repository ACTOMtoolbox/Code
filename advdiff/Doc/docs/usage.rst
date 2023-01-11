Usage
=====


.. _quickstart:

Quick-Start
----------------


Setup.ini
~~~~~~~~~~

The setup.ini config file determines the parameteres to be used when running the AdvDiff module. Here the user can for example configure how the module should be run, 
what outputs should be computed, what solvers should be used, grid specifications and how outputs should be handled.
For now, the config file can be left unchanged and left at default, however when running your own scenarios, you will have to configure 
what files to read etc in setup.ini. Read more :doc:`here <./setup>`.



Running a single scenario
~~~~~~~~~~~~~~~~~~~~~~~~~~

The cloned repository should come prepackaged with the necessary inputs and a preconfigured setup.ini, ready to run a test scenario. 
It may be helpful to run a test to get accustomed to how the module is run, and what types of inputs and outputs are expected. 

To run the module, leave everything in the setup.ini config file at default values.
Then open a python shell with the virtual environment containing all necessary dependencies.
Then navigate to the directory where the repository was cloned and type the following:

.. code-block:: console

   (venv) $ python main.py

This will run the AdvDiff module using the configurations as defined in setup.ini. The outputs will be stored in the Outdata/ folder.


To visualize the results generated in the Outdata/ folder, open a python shell with the virtual environment containing all necessary dependencies, and navigate to the directory where the repository was cloned. 
Type the following:

.. code-block:: console

   (venv) $ python visualize.py

This will run the visualizer module using the condigurations as defined in setup.ini. The outputs will be stored in the Figures/ folder.


.. note::

   It is not necessary to run the visualizer if you have set **plot = True** in the config file. 
   If True, then the visualizer will be automatically be run after a simulation terminates.



Running multiple scenarios in an ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also run multiple different scenarios by using the ensemble module. It allows you to run several simulations with different;

* Velocity files

* Starting dates

* Time-scales

* Diffusion rates

* Boundary condtition types

Without having to rewrite the config file and manually calling for main.py over and over. The output data will then be stored in Outdata/Scenario_number/.
Similarly, any figures will be stored in Figures/Scenario_number/.

To do this, you need to change some of the variables in the config file then run the enseble module. 
For example, if you want to run two simulations with the same velocity data,
but with different starting dates, you can;

1. Set **velocity_file** to the desired velocity data to read.

2. Set **time_start = Random** so the module picks a random starting date for each simulation.

3. Set **time_seed = 111, 222**. This means that first simulation's starting date will be picked with the random seed 111, 
   while the second simulation will be set with random seed 222. The two different scenarios **MUST BE COMMA SEPARATED**.

Then run the ensemble module by typing the following:

.. code-block:: console

   (venv) $ python ensemble.py

The ensemble module will now call on main.py first with **time_seed = 111** then **time_seed = 222**, giving us two different scenarios.

.. note::
   
   You can also combine multiple comma separated variables. Say, if you also wanted to vary the diffusion rate between the different scenarios,
   simply set **D = 0.0, 1.0**. Now the module will run 4 scenarios: (**time_seed = 111**, **D = 0.0**), (**time_seed = 111**, **D = 1.0**),
   (**time_seed = 222**, **D = 0.0**) and (**time_seed = 222**, **D = 1.0**).
   When constructing your own custom ensemble sets, it is important to note that the different variables are comma separated.

.. note::

   Attempting to run main.py with comma separated variables will lead to an error! If you want to run single scenarios, ensure that
   the variables do not contain multiple different options! The exception is **stat_variables** and **step_N_stats**.

.. _important_commands:

Important commands
~~~~~~~~~~~~~~~~~~~

Have a look at the :download:`Important_Commands.txt file <../../Important_Commands.txt>` to see other important commands used to run and maintain the code.

.. include:: ..//..//Important_Commands.txt
   :literal:



.. _docker:

Build a docker image
----------------------

To build an AdvDiff docker container, simply run the Dockerfile in the main directory. This can be done in VScode by right-clicking and selecting build image as shown in the figure below.

.. figure:: .//images//Docker_Build.png
   :width: 50 %
   :align: center

   How to build a docker image in VScode.

Everything in the dockerfile is preconfigured to install an instance of the correct python version and the correct dependencies into the image.
The python dependencies are installed via use of the requirements.txt file. 

However, it is important to note that gmsh requires additional packages on linux distros. That is why libgl1-mesa-glx, ffmpeg, libsm6, libxext6 and libglu1 are set to be installed.
It is critically important that these packages are installed in order. Without them, the AdvDiff will crash whenever triangle meshes are used.

The docker container is set up such that it will run independently without any inputs hooked to it. Therefore you should be able to run a test example from within the docker container without 
having to do anything. However, if you want the docker image to hook to content outside the container, you will have to do some manual work.
The most important folders are the Indata, Outdata and Figures folders. Indata contain all the necessary user input that the module needs, while the Outdata and Figures folder contain
any and all potential outputs of the module. So some work may have to be done to hook the three folders into a directory contained outside the docker container. 
Moreover, the setup.ini file located within the Indata folder contains all the userdefined settings. If you want to change what files to load when running the module, you will have
to modify the filenames within.