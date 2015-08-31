PyGSLIB
=======

This is GSLIB fortran code wrapped into python

What is implemented? 

* function to import GSLIB/Geoeas files into pandas DataFrames
* function to calculate directional varigrams (using gslib gamv)
* modified fortran version of gslib gamv function to implement  downhole variograms and variograms with lithocodes 
* modified fortran version of gslib gamv function to implement variogram cloud 

The implement the rest of the GSLIB programs is in process


Algorithms
----------
PyGSLIB implements algoritms, those are pure python functions than calls GSLIB fortran code. The algorithms functions are easy to use and include some graphical output implemented. 


Ipython notebook templates 
--------------------------
The easies way to use PyGSLIB is to modify the Ipython notebook  provided as template. Just change some input and enjoy the results. 

Notes
-----
If you are planning to use or modify this library you may understand the code organization. This python package has 3 levels of implementation: 

1. A "**low level**" python module(s) generated automatically from fortran 90 code (modified from the original fortran 77 code and some non standard GSLIB programs). To compile the fortran code into python module we use f2py. This modules are "hided" (named __f*module name*.so) 
2. The python module(s) interfacing the module auto-generated with f2py and fortran code. These are for high end users only and to develop algorithms. The examples (Ipython notebooks) that use this code are named with prefix *_raw*.
3. The algorithms modules, which are intended to simplify the use of pygslib.  

Installation
------------
The easiest way to install most Python packages is via ``easy_install`` or ``pip``:

    $ pip install opengeostat

You may need access to gfortran compiler to compile the fortran code. 

TODO: generate a conda binary distribution


Usage
-----
See the Ipython noteebooks provided in the folder `pygslib/Ipython_templates`. 


Copyright 2015, Adrian Martinez Vargas
                                                                        
This software may be modified and distributed under the terms  of the MIT license.  See the LICENSE.txt file for details.  

Friday, 28. August 2015 03:32PM 
