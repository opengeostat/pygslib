PyGSLIB
=======

This is GSLIB fortran code wrapped into python

What is implemented? 

* function to import GSLIB/Geoeas files into pandas DataFrames
* function to calculate directional varigrams (using gslib gamv)
* modified fortran version of gslib gamv funtion to implement  downhole variograms and variograms with lithocodes 
* modified fortran version of gslib gamv funtion to implement variogram cloud 

The implement the rest of the GSLIB programs is in process


Algorithms
----------
PyGSLIB implements algoritms, those are pure python functions thant calls GSLIB fortran code. The algorithms fuctions are easy to use and include some graphical output implemented  with matplotlib. 


Ipython notebook templates 
--------------------------
The easies way to use PyGSLIB is to modify the Ipython notebook  provided as template. Just change some input and enjoy the results. 

Notes
-----
If you are planning to use or modify this library you may understand  the code organization.  This python package has 3 levels of implementation: 

1. A "**low level**" python module(s) generated automatically from   fortran 90 code (modified from the original fortran 77 code and  some non standard GSLIB programs). For this we use f2py. 
2. The python module(s) interfacing the module auto-generated   with f2py and fortran code.
3. The algorithms module, which are intended to simplify more   the use of pygslib with simple an general applications. 

Installation
------------
The easiest way to install most Python packages is via ``easy_install`` or ``pip``:

    $ pip install opengeostat

You may need access to gfortran compiler to compile the fortran code. 

TODO: generate conda binary distribution


Usage
-----
See the Ipython noteebooks provided in the folder `pygslib/Ipython_templates`. 


Copyright 2015, Adrian Martinez Vargas
                                                                        
This software may be modified and distributed under the terms  of the MIT license.  See the LICENSE.txt file for details.  

Friday, 28. August 2015 03:32PM 
