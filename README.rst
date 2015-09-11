PyGSLIB
=======

This is GSLIB FORTRAN code wrapped into python

What is implemented? 

* function to import GSLIB/Geoeas files into pandas DataFrames
* function to calculate directional varigrams (using gslib gamv)
* modified FORTRAN version of GSLIB gamv function to implement  downhole variograms and variograms with lithocodes 
* modified FORTRAN version of GSLIB gamv function to implement variogram cloud 

The implement the rest of the GSLIB programs is in process


Algorithms
----------
PyGSLIB implements algorithms, those are pure python functions than calls GSLIB FORTRAN code. The algorithms functions are easy to use and include some graphical output implemented. 


Ipython notebook templates 
--------------------------
The easies way to use PyGSLIB is to modify the Ipython notebook  provided as template. Just change some input and enjoy the results. 

Notes
-----
If you are planning to use or modify this library you may understand the code organization. This python package has 3 levels of implementation: 

1. A "**low level**" python module(s) generated automatically from FORTRAN 90 code (modified from the original FORTRAN 77 code and some non standard GSLIB programs). To compile the fortran code into python module we use f2py. This modules are "hided" (named __f*module name*.so) 
2. The python module(s) interfacing the module auto-generated with f2py and FORTRAN code. These are for high end users only and to develop algorithms. The examples (Ipython notebooks) that use this code are named with prefix *_raw*.
3. The algorithms modules, which are intended to simplify the use of pygslib.

Installation in Anaconda distribution (Linux/Window/OS)
------------
The easiest way to install and work with PyGSLIB is to use Anaconda (conda) distributions. The binary packages of PyGSLIB are not in the Anaconda server yet (https://binstar.org/). This work is in process but for now you will have to compile the code.

TODO: generate a conda binary distribution

To install PyGSLIB in the root environment of your anaconda distribution follow the instructions below. 


Binary Installation in Anaconda 64 bits distribution   (Linux)
------------
There is a binary distribution at binstar.org. First install conda or anaconda distribution 
if you don't have one and run the command. 

``conda install -c https://conda.binstar.org/opengeostat pygslib``


Installation from sources (pypi.python.org) in Anaconda 32/64 bits distribution   (Linux)
------------
Install dependencies: 

 
``$ conda install numpy pandas matplotlib``



Install PyGSLIB with  ``easy_install`` or ``pip``:



``$ pip install pygslib``



You may need access to gfortran compiler to compile the FORTRAN code. This is usually available in Linux most linux distributions. 


Installation from sources (pypi.python.org) in Anaconda 32 bits distribution (Windows)
______________________________
Install dependencies, including mingw which comes with gfortran: 


``C:\>conda install mingw numpy pandas matplotlib``


Install PyGSLIB with  ``easy_install`` or ``pip`` using gfortran 32 bits compiler


``C:\>pip install --global-option build_ext --global-option --compiler=mingw32 pygslib``



Installation from sources (pypi.python.org) in Anaconda 64 bits distribution  (Windows)
______________________________
Install dependencies: 

 

``C:\>conda install numpy pandas matplotlib`` 



Install mingw with 64 bit compiler



``C:\>conda install -c https://conda.binstar.org/omnia mingwpy ``



Install PyGSLIB with  `easy_install` or `pip` using gfortran 64 bits compiler:


``C:\>pip install --global-option build_ext --global-option --compiler=mingw32 pygslib``

If you get an error like this::

    File "C:\Users\Your_Path_Here\Anaconda\envs\test3\lib\site-packages\numpy\distutils\fcompiler\gnu.py", 
    line 337, in get_libraries raise NotImplementedError("Only MS compiler supported with gfortran on win64")
    NotImplementedError: Only MS compiler supported with gfortran on win64



Don't worry, this is a known issue in numpys distutils. Go to the file 

``C:\Users\YYOUR_USER_NAME\Anaconda\lib\site-packages\numpy\distutils\fcompiler\gnu.py``

or this file, if you are installing PyGSLIB in an environment

``C:\Users\YYOUR_USER_NAME\Anaconda\envs\YOUR_ENVIRONMENT\lib\site-packages\numpy\distutils\fcompiler\gnu.py``

around the line 337 you will see::

    if is_win64():
        c_compiler = self.c_compiler
        if c_compiler and c_compiler.compiler_type == "msvc":
            return []
        else:
            raise NotImplementedError("Only MS compiler supported with gfortran on win64")



rewrite the code like this::

	if is_win64():
		c_compiler = self.c_compiler
		if c_compiler and c_compiler.compiler_type == "msvc":
		    return []
		else:
		    return [] #raise NotImplementedError("Only MS compiler supported with gfortran on win64")



and rerun


``C:\>pip install --global-option build_ext --global-option --compiler=mingw32 pygslib``


This may fix the problem


Installation from source (from github.com)
--------------------
This is the most update but unstable development version. You may install all the dependencies 
manually and make sure you have a gfortran available. 


	git clone https://github.com/opengeostat/pygslib.git
	cd pygslib
	python setup.py install 



Usage
-----
See the Ipython noteebooks provided in the folder ``pygslib/Ipython_templates``. 



License 
-------
Copyright 2015, Adrian Martinez Vargas

Supported by Opengeostat Consulting @ http://opengeostat.com/

                                                                 
This software may be modified and distributed under the terms  of the MIT license.  See the LICENSE.txt file for details.

Wed 02 Sep 2015 

