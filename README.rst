PyGSLIB
=======

This is an open source python module for mineral resource estimation and geostatistics. It consists of four main sub modules:  
- ``gslib``. This is for geostatistics and interpolation. It was built with [GSLIB Fortran 77 code] ( http://www.statios.com/Quick/gslib.html) enhanced and linked to Python. 
- ``drillhole``. This is for basic drillhole operations, such as compositing and desurveying.
- ``blockmodel``. This is for block modelling, it has functions to fill wireframes with blocks, reblocking, among others.
- ``vtktools``. This is for 3D computational geometry based on VTK, for example, to select samples within wireframes. It also handles VTK files.
- ``nonlinear``. This module is under construction! It is an experimental module for nonlinear geostatistics based on the Discrete Gaussian Model.

Current version
----------
version = '0.0.0.3.8.4'


This means that we are in a very early developing stage and the package is experimental!


Need some help? 
------ 
Visit [the PyGSLIB online help](https://opengeostat.github.io/pygslib/)
There is also a [Youtube channel](https://www.youtube.com/c/opengeostat) with some demonstrations

You can also find some examples at the [Ipython notebook templates and examples]( https://github.com/opengeostat/pygslib/tree/master/pygslib/Ipython_templates), but keep in mind that this may not include some new functionality.
Notes
-----
If you are planning to use or modify this library you may understand the code organization. 
The code is organized in two separated folders
- one with Cython code named [cython_code]( https://github.com/opengeostat/pygslib/tree/master/cython_code)
- one with Fortran 90 code named [for_code](https://github.com/opengeostat/pygslib/tree/master/for_code)


Installation in Anaconda distribution (Linux, Window and OS)
------------
The easiest way to install and work with PyGSLIB is using an Anaconda 
(conda) distribution. To install PyGSLIB in the root environment of 
your anaconda distribution simply type in a terminal:  

``conda install -c opengeostat pygslib``


Installation from source (from github.com)
--------------------
This is the most update but unstable development version. You may manually 
install all the dependencies and make sure you have gfortran available.  For development use the following shell script/commands. 
```
$ git clone https://github.com/opengeostat/pygslib.git
$ cd pygslib
$ python setup.py develop
```
  
    
To update this module as contributor, make changes and the update git (requesting a pull).


Usage
-----
See  [tutorial] (https://opengeostat.github.io/pygslib/Tutorial.html),  [video demonstrations]( https://youtu.be/SEwKy6wJbLE) and the [Ipython notebook examples] ( https://github.com/opengeostat/pygslib/tree/master/pygslib/Ipython_templates). 


License 
-------
Copyright 2016, Adrian Martinez Vargas

Supported by Opengeostat Consulting @ http://opengeostat.com/
                                                                 
This software may be modified and distributed under the terms of the 
MIT and GPL licenses.  See the LICENSE.txt file for details.

Monday 12 August 2016


