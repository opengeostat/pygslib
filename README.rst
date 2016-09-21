PyGSLIB
=======

This is two things: 

- a GSLIB FORTRAN code wrapped into python
- a set of Python/Cython Function and Classes for drillhole processing,
  block model creation and manipulation, search neighborhood, VTK
  data visualization and exporting and non-linear geostatistical 
  applications 


Current version
----------
version = '0.0.0.3.8'


This means that we are in a very early developing stage and the package 
is experimental!


Need some help? 
------ 
User manual at https://opengeostat.github.io/pygslib/
Youtube chanel https://www.youtube.com/c/opengeostat


Ipython notebook templates and examples.
--------------------------
The easiest way to use PyGSLIB is to modify the Ipython notebook 
provided as template and examples. Just change some input and enjoy 
the results. 

Keep in mind that new functionality is not implemented in the examples yet. 

You may work with a vtk viewer (e.j. Paraview) for 3D visualization. 

Notes
-----
If you are planning to use or modify this library you may understand 
the code organization. 

The code is organized in two separated folders 

- cython_code
- for_code

The fortran code is in the folder ``for_code`` and has 2 levels of 
implementation: 


1. A "**low level**" python module(s) generated automatically from 
   FORTRAN 90 code (modified from the original FORTRAN 77 code and 
   some non standard GSLIB programs). To compile the fortran code 
   into python module we use f2py. This modules are 
   "hided" (named __f*module name*.so) 
2. The python module(s) interfacing the module auto-generated with f2py. 
   These are for high end users only and to develop algorithms. 
   The examples (Ipython notebooks) that use this code are named with 
   prefix *_raw*.


Installation in Anaconda distribution (Linux/Window/ {OS not implemented yet})
------------
The easiest way to install and work with PyGSLIB is using an Anaconda 
(conda) distribution. To install PyGSLIB in the root environment of 
your anaconda distribution simply type in a terminal:  


``conda install -c opengeostat pygslib``



Installation from source (from github.com)
--------------------
This is the most update but unstable development version. You may manually 
install all the dependencies and make sure you have gfortran available:: 


    git clone https://github.com/opengeostat/pygslib.git
    cd pygslib
    python setup.py install 
    
    or 

    python setup.py develop
    

To update this module as contributor, make changes and the update git (requesting a pull).

To update the pypi repository::

    python setup.py sdist upload -r pypi

To update conda repository(Linux)::

    conda skeleton pypi pygslib
    conda build pygslib
    anaconda upload /home/adrian/anaconda/conda-bld/linux-64/pygslib-0.0.0.3.#-nppy27_0.tar.bz2


Usage
-----
See the Ipython noteebooks provided in the folder ``pygslib/Ipython_templates``. 


License 
-------
Copyright 2016, Adrian Martinez Vargas

Supported by Opengeostat Consulting @ http://opengeostat.com/

                                                                 
This software may be modified and distributed under the terms of the 
MIT license.  See the LICENSE.txt file for details.

Monday 12 August 2016


