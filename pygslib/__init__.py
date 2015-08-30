# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""PyGSLIB A module that wraps GSLIB fortran code into python.

This module provides an interface to FGSLIB (a python module 
generated automatically with f2p and modified gslib fortran code)


Notes
-----
    The module is organized in raw functions, linked directly directly 
    to `FGSLIB`, and algorithm functions. For example, the function 
    `gamv` call directly the gslib gamv function but the output is 
    a complicated set of variables. To simplify the work we implemented 
    `gamv_alg_directional_1V` for directional variogram of a single 
    variable.

    Algorithm functions are named with the following convention. 
    `{parentfunction}_alg_{algoritm explanation}`. 


Copyright 2015, Adrian Martinez Vargas
                                                                        
This software may be modified and distributed under the terms        
of the MIT license.  See the LICENSE.txt file for details.  
         
"""                                                                      

import os.path
import pandas as pd
import __fgslib
import platform
import warnings
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------------------------------------------
#
#    General functions
#
#-----------------------------------------------------------------------------------------------------------------

def version():
    """
    Version of the gslib fortran code and the python wrapper 

    Returns
    -------
    dict 
        {'python':python wrap version, 'fortran': gslib fortran version}, [operating system] 

        each version have 

       {  'major':major , 
          'minor':minor , 
          'maintenance': maintenance, 
          'build': build, 
          'month':month, 
          'year': year}
    """

    major, minor , maintenance, build, month, year= __fgslib.version()
    
    foversion={'major':major, 
               'minor':minor , 
               'maintenance': maintenance, 
               'build': build, 
               'month':month, 
               'year': year}
    pyversion={'major':0, 
               'minor':0 , 
               'maintenance':0, 
               'build':2, 
               'month':7, 
               'year':2016}

    osplatform=platform.platform()

    return {'fortran version': foversion, 
            'python version': pyversion,
            'platform': osplatform }


#read GSLIB file
def read_gslib_file(fname, maxnvar=500):
    """ Read a gslib file
    
    The native fortran code is used to read geoeas/gslib files

    Parameters
    ----------    
    fname : str 
        file name path (absolute or relative)
    
    maxnvar : Optional[int]
        maximum number of variables in file
    
    Returns
    -------
    data : DataFrame 
        The file is returned in memory as a pandas dataframe


    Examples
    --------
    >>> import pygslib as gslib     
    >>> mydata= gslib.read_gslib_file('dataset/cluster.dat') 
    >>> print mydata.head(n=5)
              Xlocation  Ylocation  Primary  Secondary  Declustering Weight
        0       39.5       18.5     0.06       0.22                1.619
        1        5.5        1.5     0.06       0.27                1.619
        2       38.5        5.5     0.08       0.40                1.416
        3       20.5        1.5     0.09       0.39                1.821
        4       27.5       14.5     0.09       0.24                1.349
    """

    # make sure the file exists 
    assert os.path.isfile(fname), "invalid file name"
        
    comment_line,varnames,nvar,error = __fgslib.read_header(fname ,maxnvar)
    #the var names requires this work around in f2py
    vnames=[]
    for i in varnames[0:nvar]:
        vnames.append(str.strip(i.tostring()))

    #now read data 
    maxdat,error=__fgslib.read_ndata(fname, nvar)
    table,error=__fgslib.read_data(fname, nvar, maxdat)
    
    assert error==0, "there was an error= %r importing the data" % error
    
    return pd.DataFrame(table, columns=vnames)


#generate coordinates for gridded files
def addcoord(nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, grid):
    """ Insert X, Y, Z coordinates to a grid file 
    
    The code use modified function addcoord from gslib

    Parameters
    ----------    
    nx,ny,nz : int 
        number of rows

    xmn,ymn,zmn : double 
        origin of coordinates (block centroid)

    xsiz,ysiz,zsiz : double 
        block sizes

    grid : DataFrame 
        Pandas dataframe containing variables
    
    Returns
    -------
    data : DataFrame 
        Pandas dataframe (drid) with tree extra columns x, y, z

    Notes 
    -------
    Grid is expected to have xmn*ymn*zmn rows with geoeas grid format order

    """

    assert  type(grid)==pd.DataFrame , 'the parameter grid is not a Pandas Dataframe, incorrect data type (it is %s)' % type(grid)
    assert  len(grid) ==  nx*ny*nz , 'len(grid) !=  nx*ny*nz, unexpected grid number of raws'
    #make sure you don't overwrite x,y,z
    assert  'x' not in  grid, 'x already exist in grid, use grid.drop("x", axis=1, inplace=True) to remove x' 
    assert  'y' not in  grid, 'y already exist in grid, use grid.drop("y", axis=1, inplace=True) to remove y'
    assert  'z' not in  grid, 'z already exist in grid, use grid.drop("z", axis=1, inplace=True) to remove z'


    #adding columns (not inplace)
    x,y,z = __fgslib.addcoord(nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz)   
    tmpgrid= grid.copy(deep=True)
    tmpgrid.insert(loc=0, column='z', value=z, allow_duplicates=False)
    tmpgrid.insert(loc=0, column='y', value=y, allow_duplicates=False)
    tmpgrid.insert(loc=0, column='x', value=x, allow_duplicates=False)
    
    return pd.DataFrame(tmpgrid)


#-----------------------------------------------------------------------------------------------------------------
#
#    Variograms GAMV low level functions
#
#-----------------------------------------------------------------------------------------------------------------
    
def gamv(parameters):
    """Calculate variograms with the modified GSLIB gamv function
    
    This function is similar to the GSLIB gamv function but with some minor differences:
    
        a) the indicator automatic transformation is not applied, you may define 
           indicator variables externally
        b) a bhid variable is always required, to avoid downhole variogram use 
           `bhid=None` in the parameter file or `bhid = 0` in the file 
            **Warning:** bhid may be an integer
           
        c) Z variable is required, you can use constant coordinate to reduce 
           dimensions, for example use Z=None in the parameter file or 
           Z = 0 in the file 
        d) The variogram cloud was implemented but it only works for the first 
           variogram and first direction. Some variogram types will ignore this
    
    the parameter file here is a dictionary with the following keys.
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, for example:



            parameters = { 
                    'x'      :  X,                      # X coordinates, array('f') with bounds (nd), nd is number of data points
                    'y'      :  Y,                      # Y coordinates, array('f') with bounds (nd)
                    'z'      :  Z,                      # Z coordinates, array('f') with bounds (nd)
                    'bhid'   :  bhid,                   # bhid for downhole variogram, array('i') with bounds (nd)    
                    'vr'     :  VR,                     # Variables, array('f') with bounds (nd,nv), nv is number of variables
                    'tmin'   : -1.0e21,                 # trimming limits, float
                    'tmax'   :  1.0e21,                 # trimming limits, float
                    'nlag'   :  10,                     # number of lags, int
                    'xlag'   :  4,                      # lag separation distance, float                
                    'xltol'  :  2,                      # lag tolerance, float
                    'azm'    : [0,0,90],                # azimut, array('f') with bounds (ndir)
                    'atol'   : [90,22.5,22.5],          # azimut tolerance, array('f') with bounds (ndir)
                    'bandwh' : [50,10,10],              # bandwith h, array('f') with bounds (ndir)
                    'dip'    : [0,0,0],                 # dip, array('f') with bounds (ndir)
                    'dtol'   : [10,10,10],              # dip tolerance, array('f') with bounds (ndir)
                    'bandwd' : [10,10,10],              # bandwith dit, array('f') with bounds (ndir)
                    'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100],                   # variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1,1,1,1,1,1,1],         # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1,1,1,1,1,1,1],         # head var., array('i') with bounds (nvarg)
                    'ivtype' : [1,3,4,5,6,7,8],         # variogram type, array('i') with bounds (nvarg)
                    'maxclp' : 50000}                   # maximum numver of variogram point cloud to use, input int
                    
                   

    the equivalent GSLIB parameter file may be something like this:

                              Parameters for GAMV
                          *******************
        
        START OF PARAMETERS:
        cluster.dat                       -file with data
        1   2   0   4                     -   columns for X, Y, Z coordinates, bhid
        1   3                             -   number of variables,col numbers
        -1.0e21     1.0e21                -   trimming limits
        gamv.out                          -file for variogram output
        10                                -number of lags
        4.0                               -lag separation distance
        2.0                               -lag tolerance
        3                                 -number of directions
        0.0  90.0 50.0   0.0  10.0  10.0  -azm,atol,bandh,dip,dtol,bandv
        0.0  22.5 10.0   0.0  10.0  10.0  -azm,atol,bandh,dip,dtol,bandv
        90.  22.5 10.0   0.0  10.0  10.0  -azm,atol,bandh,dip,dtol,bandv
        0                                 -standardize sills? (0=no, 1=yes)
        7                                 -number of variograms
        1   1   1                         -tail var., head var., variogram type
        1   1   3                         -tail var., head var., variogram type
        1   1   4                         -tail var., head var., variogram type
        1   1   5                         -tail var., head var., variogram type
        1   1   6                         -tail var., head var., variogram type
        1   1   7                         -tail var., head var., variogram type
        1   1   8                         -tail var., head var., variogram type
        50000                             -maxclp, maximum number of point cloud
    
    Returns
    -------   
    out : (ndarray,ndarray, ndarray,ndarray,ndarray,ndarray,ndarray, ndarray, ndarray, ndarray, ndarray)

      for directional variograms we have:
       pdis   Distance of pairs falling into this lag 
       pgam   Semivariogram, covariance, correlogram,... value
       phm    Mean of the tail data
       ptm    Mean of the head data
       phv    Variance of the tail data
       ptv    Variance of the head data
       pnump  Number of pairs

      for variogram cloud points:
       cldi   data index of the head 
       cldj   data index of the tail
       cldg   Semivariogram, covariance, ... value
       cldh   Distance of each pair
      
    Notes
    -----
    The output is a tuple of numpy 3D ndarrays (pdis,pgam, phm,ptm,phv,ptv,pnump) with dimensions 
    (nvarg, ndir, nlag+2), representing the  experimental variograms output, and 1D array 
    (cldi, cldj, cldg, cldh) representing variogram cloud for 


    The variables with prefix `p` are for directional variogram and are generated 
    by GSLIB fortran function  `__fgslib.writeout`. This is similar to the output of 
    the GSLIB standalone program `gamv`.

    The variables with prefix `cld` are for variogram cloud. The variogram cloud 
    was implemented by modifying the original gamv fortran funtion. This only 
    works for the first variogram/direction and and only if the variogram type is 1, 2, 6, 7 and 8. 
    That is: 
        traditional semivariogram (1), traditional cross semivariogram (2), 
        pairwise relative semivariogram (6), semivariogram of logarithms (7), 
        semimadogram(8)
    
    """

    # ERROR (multivariate vr need to be extended in 1D array)
    # TODO: fix this 


    np,dis, gam, hm, tm, hv, tv, cldi, cldj, cldg, cldh, l = __fgslib.gamv(**parameters)
    
    if l==parameters['maxclp']:
        warnings.warn( 'Warning: l == maxclp; maximum number ' + \
                       'of point clouds reached, increase maxclp' + \
                       'or modify variogram parameters')
    
    #remove crap data from variogram cloud
    cldi=cldi[:l]
    cldj=cldj[:l]
    cldg=cldg[:l]
    cldh=cldh[:l]

    # get output like in gslib
    
    ndir = len(parameters['azm']) 
    nvarg = len(parameters['ivtype']) 
    nlag = parameters['nlag']

    
    pdis,pgam, phm,ptm,phv,ptv,pnump = __fgslib.writeout(nvarg,ndir,nlag,np,dis,gam,hm,tm,hv,tv)
      
    
    
    return pdis,pgam, phm,ptm,phv,ptv,pnump, cldi, cldj, cldg, cldh    

#define a dummy  parameter for variogram 
gamv_parameter_template = { 
                    'x'      :  [1,2,3,4,5,6,7,8,9],    # X coordinates, array('f') with bounds (nd), nd is number of data points
                    'y'      :  [0,0,0,0,0,0,0,0,0],    # Y coordinates, array('f') with bounds (nd)
                    'z'      :  [0,0,0,0,0,0,0,0,0],    # Z coordinates, array('f') with bounds (nd)
                    'bhid'   :  [0,0,0,0,0,0,0,0,0],    # bhid for downhole variogram, array('i') with bounds (nd)    
                    'vr'     :  [0,1,1,0,0,1,1,0,0],    # Variables, array('f') with bounds (nd,nv), nv is number of variables
                    'tmin'   : -1.0e21,                 # trimming limits, float
                    'tmax'   :  1.0e21,                 # trimming limits, float
                    'nlag'   :  3,                      # number of lags, int
                    'xlag'   :  1,                      # lag separation distance, float                
                    'xltol'  :  -5,                     # lag tolerance, float
                    'azm'    : [0],                     # azimut, array('f') with bounds (ndir)
                    'atol'   : [90],                    # azimut tolerance, array('f') with bounds (ndir)
                    'bandwh' : [50],                    # bandwith h, array('f') with bounds (ndir)
                    'dip'    : [0],                     # dip, array('f') with bounds (ndir)
                    'dtol'   : [22.5],                  # dip tolerance, array('f') with bounds (ndir)
                    'bandwd' : [50],                    # bandwith dit, array('f') with bounds (ndir)
                    'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100],                   # variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1],                     # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1],                     # head var., array('i') with bounds (nvarg)
                    'ivtype' : [1],                     # variogram type, array('i') with bounds (nvarg)
                    'maxclp' : 50000}                   # maximum numver of variogram point cloud to use, input int
"""
This is a dummy gamv parameter dict. 
{ 'x'      :  [1,2,3,4,5,6,7,8,9],    # X coordinates, array('f') with bounds (nd), nd is number of data points
  'y'      :  [0,0,0,0,0,0,0,0,0],    # Y coordinates, array('f') with bounds (nd)
  'z'      :  [0,0,0,0,0,0,0,0,0],    # Z coordinates, array('f') with bounds (nd)
  'bhid'   :  [0,0,0,0,0,0,0,0,0],    # bhid for downhole variogram, array('i') with bounds (nd)    
  'vr'     :  [0,1,1,0,0,1,1,0,0],    # Variables, array('f') with bounds (nd,nv), nv is number of variables
  'tmin'   : -1.0e21,                 # trimming limits, float
  'tmax'   :  1.0e21,                 # trimming limits, float
  'nlag'   :  3,                      # number of lags, int
  'xlag'   :  1,                      # lag separation distance, float                
  'xltol'  :  -5,                     # lag tolerance, float
  'azm'    : [0],                     # azimut, array('f') with bounds (ndir)
  'atol'   : [90],                    # azimut tolerance, array('f') with bounds (ndir)
  'bandwh' : [50],                    # bandwith h, array('f') with bounds (ndir)
  'dip'    : [0],                     # dip, array('f') with bounds (ndir)
  'dtol'   : [22.5],                  # dip tolerance, array('f') with bounds (ndir)
  'bandwd' : [50],                    # bandwith dit, array('f') with bounds (ndir)
  'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
  'sills'  : [100],                   # variance used to std the sills, array('f') with bounds (nv)
  'ivtail' : [1],                     # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
  'ivhead' : [1],                     # head var., array('i') with bounds (nvarg)
  'ivtype' : [1],                     # variogram type, array('i') with bounds (nvarg)
  'maxclp' : 50000}                   # maximum numver of variogram point cloud to use, input int
}

Notes
-----
Note that data are passed as array like (and not as file path, like in the standalone gslib function gamv).
Any numpy or list may work as input data but make sure that `x,y,z,bhid,vr` have the same dimentions. 
For multivariate data use ndarray/lists to define `vr`. 

"""

    

#review and fix this
def __check_vario__(ivtail,ivhead,ivtype,vrmin,vrmax):
    """ This is the python version of the f77 function check in gamv.f
    This function is for internal use only. 

    

    """
    
    nvarg=len(ivtype)
   


    for iv in range(nvarg):
                
        it = int(abs(ivtype[iv]))

        assert it>=1 and it<=8 , "ivtype is not valid: check parameters... it>=1 and it<=8"


        
        if(it == 1): title='Semivariogram          :'
        if(it == 2): title='Cross Semivariogram    :'
        if(it == 3): title='Covariance             :'
        if(it == 4): title='Correlogram            :'
        if(it == 5): title='General Relative       :'
        if(it == 6): title='Pairwise Relative      :'
        if(it == 7): title='Variogram of Logarithms:'
        if(it == 8): title='Semimadogram           :'
        if(it == 9): title='Indicator 1/2 Variogram:'
        if(it == 10):title='Indicator 1/2 Variogram:'

        # Check for possible errors or inconsistencies:
            
        if it == 2:
            if ivtail[iv] == ivhead[iv]:
                warnings.warn('  WARNING: cross variogram with the same variable!')
        elif it == 5:
            if ivtail[iv] != ivhead[iv]:
                warnings.warn('  WARNING: cross variogram with the same variable!')
            if vrmin[ivtail[iv]-1] < 0.0 and vrmax[ivtail[iv]-1] > 0.0:   
                warnings.warn('  WARNING: cross general relative variogram are')
            if vrmin[ivhead[iv]-1] < 0.0 and vrmax[ivhead[iv]-1] > 0.0: 
                warnings.warn('  WARNING: there are both positive and negative values - lag mean could be zero!')
        elif it == 6:
            if ivtail[iv] != ivhead[iv]:
                warnings.warn('  WARNING: cross pairwise relative variogram are difficult to interpret!')
            if vrmin[ivtail[iv]-1] < 0.0 and vrmax[ivtail[iv]-1] > 0.0:
                warnings.warn('  WARNING: there are both positive and negative values - pair means could be zero!')
            if vrmin[ivhead[iv]-1] < 0.0 and vrmax[ivhead[iv]-1] > 0.0:
                warnings.warn('  WARNING: there are both positive and negative values - pair means could be zero!')
        elif it == 7:
            if ivtail[iv] != ivhead[iv]: 
                warnings.warn('  WARNING: cross logarithmic variograms may be difficult to interpret!')
            if vrmin[ivtail[iv]-1] < 0.0 or vrmin[ivhead[iv]-1] < 0.0:  
                warnings.warn('  WARNING: there are zero or negative values - logarithm undefined!')
        
        return 1




#check parameters
def check_gamv_par(parameters):
    """Check the variogram parameter dictionary
    This function verify the validity of the parameter file and check for possible errors or inconsistencies

    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, for example:



            parameters = { 
                    'x'      :  X,                      # X coordinates, array('f') with bounds (nd), nd is number of data points
                    'y'      :  Y,                      # Y coordinates, array('f') with bounds (nd)
                    'z'      :  Z,                      # Z coordinates, array('f') with bounds (nd)
                    'bhid'   :  bhid,                   # bhid for downhole variogram, array('i') with bounds (nd)    
                    'vr'     :  VR,                     # Variables, array('f') with bounds (nd,nv), nv is number of variables
                    'tmin'   : -1.0e21,                 # trimming limits, float
                    'tmax'   :  1.0e21,                 # trimming limits, float
                    'nlag'   :  10,                     # number of lags, int
                    'xlag'   :  4,                      # lag separation distance, float                
                    'xltol'  :  2,                      # lag tolerance, float
                    'azm'    : [0,0,90],                # azimut, array('f') with bounds (ndir)
                    'atol'   : [90,22.5,22.5],          # azimut tolerance, array('f') with bounds (ndir)
                    'bandwh' : [50,10,10],              # bandwith h, array('f') with bounds (ndir)
                    'dip'    : [0,0,0],                 # dip, array('f') with bounds (ndir)
                    'dtol'   : [10,10,10],              # dip tolerance, array('f') with bounds (ndir)
                    'bandwd' : [10,10,10],              # bandwith dit, array('f') with bounds (ndir)
                    'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100],                   # variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1,1,1,1,1,1,1],         # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1,1,1,1,1,1,1],         # head var., array('i') with bounds (nvarg)
                    'ivtype' : [1,3,4,5,6,7,8],         # variogram type, array('i') with bounds (nvarg)
                    'maxclp' : 50000}                   # maximum numver of variogram point cloud to use, input int


    Returns
    -------   
    out : bool 
      return True (or 1) if the parameter dictionary is valid

    """

    #check the parameters are correct 
    assert set(gamv_parameter_template)==set(parameters)

    nvar= len(parameters['vr'])/len(parameters['x'])
    ndir= len(parameters['azm'])
    nvarg=len(parameters['bandwd'])
    assert nvar >=1 , "nvar is too small: check parameters"
    assert len(parameters['sills'])==nvar , "incorrect sill parameters; parameters['sills'].shape[1]!=nvar"
    assert parameters['nlag']>=1, "nlag is too small: check parameters"
    assert  parameters['xlag'] > 0, "xlag is too small: check parameter file"
    if  parameters['xltol'] <= 0:
        warnings.warn('xltol is too small: resetting to xlag/2')
        parameters['xltol']=parameters['xlag']*0.5

    assert  ndir > 0, "ndir is too small: check parameter file"

    for i in range(ndir):
        assert parameters['bandwh']>=0, "Horizontal bandwidth is too small!"
        assert parameters['bandwd']>=0, "Vertical bandwidth is too small!"
    
    assert nvarg>0, "nvarg is too small: check parameters"

    vrmin=np.array(np.min(parameters['vr'], axis=0), ndmin=1)   #  this is to force array if np.min is scalar
    vrmax=np.array(np.max(parameters['vr'], axis=0), ndmin=1)
    #now check if the variogram is ok
    assert __check_vario__(parameters['ivtail'],
                       parameters['ivhead'],
                       parameters['ivtype'],
                       vrmin,
                       vrmax)==1, "error in the experimental variogram parameters"
    

    return 1




#-----------------------------------------------------------------------------------------------------------------
#
#    Variograms GAM, for gridded data
#
#-----------------------------------------------------------------------------------------------------------------
def gam(parameters):
    """Calculate variograms with the modified GSLIB gam function (for gridded data)
       
    The parameter file here is a dictionary with the following keys.
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, for example:


            parameters = { 
                    'nx'     :  50,                     # number of rows in the gridded data
                    'ny'     :  50,                     # number of columns in the gridded data
                    'nz'     :  1,                      # number of levels in the gridded data
                    'xsiz'   :  1,                      # size of the cell in x direction 
                    'ysiz'   :  1,                      # size of the cell in y direction
                    'zsiz'   :  1,                      # size of the cell in z direction
                    'bhid'   :  bhid,                   # bhid for downhole variogram, array('i') with bounds (nd)    
                    'vr'     :  VR,                     # Variables, array('f') with bounds (nd,nv), nv is number of variables
                    'tmin'   : -1.0e21,                 # trimming limits, float
                    'tmax'   :  1.0e21,                 # trimming limits, float
                    'nlag'   :  10,                     # number of lags, int
                    'ixd'    : [1,0],                   # directiom x 
                    'iyd'    : [1,0],                   # directiom y 
                    'izd'    : [1,0],                   # directiom z 
                    'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100, 200],              # variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1,1,2,2],               # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1,1,2,2],               # head var., array('i') with bounds (nvarg)
                    'ivtype' : [1,3,1,3]}               # variogram type, array('i') with bounds (nvarg)

                    
                  
    the equivalent GSLIB parameter file may be something like this:

                                      Parameters for GAM
                                      ******************

                    START OF PARAMETERS:
                    ../data/true.dat      -file with data
                    2   1   2             -   number of variables, column numbers
                    -1.0e21     1.0e21    -   trimming limits
                    gam.out               -file for variogram output
                    1                     -grid or realization number
                    50   0.5   1.0        -nx, xmn, xsiz
                    50   0.5   1.0        -ny, ymn, ysiz
                     1   0.5   1.0        -nz, zmn, zsiz
                    2  10                 -number of directions, number of lags
                     1  0  0              -ixd(1),iyd(1),izd(1)
                     0  1  0              -ixd(2),iyd(2),izd(2)
                    1                     -standardize sill? (0=no, 1=yes)
                    4                     -number of variograms
                    1   1   1             -tail variable, head variable, variogram type
                    1   1   3             -tail variable, head variable, variogram type
                    2   2   1             -tail variable, head variable, variogram type
                    2   2   3             -tail variable, head variable, variogram type

    
    Returns
    -------   
    out : (ndarray,ndarray, ndarray,ndarray,ndarray,ndarray,ndarray, ndarray, ndarray, ndarray, ndarray)

      for directional variograms we have:
       pdis   Distance of pairs falling into this lag 
       pgam   Semivariogram, covariance, correlogram,... value
       phm    Mean of the tail data
       ptm    Mean of the head data
       phv    Variance of the tail data
       ptv    Variance of the head data
       pnump  Number of pairs

      
    Notes
    -----
    The output is a tuple of numpy 3D ndarrays (pdis,pgam, phm,ptm,phv,ptv,pnump) with dimensions 
    (nvarg, ndir, nlag+2), representing the  experimental variograms output

    The variables with prefix `p` are for directional variogram and are generated 
    by GSLIB fortran function  `__fgslib.writeout`. This is similar to the output of 
    the GSLIB standalone program `gam`.
    """

    # ERROR (multivariate vr need to be extended in 1D array)
    # TODO: fix this 
    np, gam, hm, tm, hv, tv = __fgslib.gamma(**parameters)
    

    # get output like in gslib 
    nvarg = len(parameters['ivtype']) 
    nlag = parameters['nlag']
    ixd = parameters['ixd']
    xsiz = parameters['xsiz']
    iyd = parameters['iyd']
    ysiz = parameters['ysiz']
    izd = parameters['izd']
    zsiz = parameters['zsiz']

    
    pdis,pgam,phm,ptm,phv,ptv,pnump = __fgslib.writeout_gam(nvarg,nlag,ixd,xsiz,iyd,ysiz,izd,zsiz,np,gam,hm,tm,hv,tv)
      
    
    
    return pdis,pgam, phm,ptm,phv,ptv,pnump 



