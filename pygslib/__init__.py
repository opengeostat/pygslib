# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""PyGSLIB A module that wraps GSLIB fortran code into python.

This module provides an interface to FGSLIB (a python module 
generated automatically with f2p and modified gslib fortran code)

Copyright 2015, Adrian Martinez Vargas
                                                                        
This software may be modified and distributed under the terms        
of the MIT license.  See the LICENSE.txt file for details.  
         
"""                                                                      

import os.path
import pandas as pd
import __rotscale
import __block_covariance
import __read_gslib
import __addcoord
import __kt3d
import __plot
import __declus
import __variograms
import __dist_transf
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
    Return version information (only python code)

    Returns
    -------
    dict 
        {'python':python wrap version}, [operating system] 

        each version have 

       {  'major':major , 
          'minor':minor , 
          'maintenance': maintenance, 
          'build': build, 
          'month':month, 
          'year': year}
    """

    pyversion={'major':0, 
               'minor':0 , 
               'maintenance':0, 
               'build':3, 
               'month':9, 
               'year':2015}

    osplatform=platform.platform()

    return {'python version': pyversion,
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
        
    comment_line,varnames,nvar,error = __read_gslib.read_header(fname ,maxnvar)
    #the var names requires this work around in f2py
    vnames=[]
    for i in varnames[0:nvar]:
        vnames.append(str.strip(i.tostring()))

    #now read data 
    maxdat,error=__read_gslib.read_ndata(fname, nvar)
    table,error=__read_gslib.read_data(fname, nvar, maxdat)
    
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

    assert  type(grid)==pd.DataFrame , 'the parameter grid is not a pandas DataFrame, incorrect data type (it is %s)' % type(grid)
    assert  len(grid) ==  nx*ny*nz , 'len(grid) !=  nx*ny*nz, unexpected grid number of rows'
    #make sure you don't overwrite x,y,z
    assert  'x' not in  grid, 'x already exist in grid, use grid.drop("x", axis=1, inplace=True) to remove x' 
    assert  'y' not in  grid, 'y already exist in grid, use grid.drop("y", axis=1, inplace=True) to remove y'
    assert  'z' not in  grid, 'z already exist in grid, use grid.drop("z", axis=1, inplace=True) to remove z'


    #adding columns (not inplace)
    x,y,z = __addcoord.addcoord(nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz)   
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
                    'bandwh' : [50,10,10],              # bandwith 'horizontal', array('f') with bounds (ndir)
                    'dip'    : [0,0,0],                 # dip, array('f') with bounds (ndir)
                    'dtol'   : [10,10,10],              # dip tolerance, array('f') with bounds (ndir)
                    'bandwd' : [10,10,10],              # bandwith 'vertical', array('f') with bounds (ndir)
                    'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100],                   # variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1,1,1,1,1,1,1],         # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1,1,1,1,1,1,1],         # head var., array('i') with bounds (nvarg)
                    'ivtype' : [1,3,4,5,6,7,8],         # variogram type, array('i') with bounds (nvarg)
                    'maxclp' : 50000}                   # maximum number of variogram point cloud to use, input int
                    
                   

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
    by GSLIB fortran function  `writeout`. This is similar to the output of 
    the GSLIB standalone program `gamv`.

    The variables with prefix `cld` are for variogram cloud. The variogram cloud 
    was implemented by modifying the original gamv fortran funtion. This only 
    works for the first variogram/direction and and only if the variogram type is 1, 2, 6, 7 and 8. 
    That is: 
        traditional semivariogram (1), traditional cross semivariogram (2), 
        pairwise relative semivariogram (6), semivariogram of logarithms (7), 
        semimadogram(8)
    
    """

    np,dis, gam, hm, tm, hv, tv, cldi, cldj, cldg, cldh, l = __variograms.gamv(**parameters)
    
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

    
    pdis,pgam, phm,ptm,phv,ptv,pnump = __variograms.writeout(nvarg,ndir,nlag,np,dis,gam,hm,tm,hv,tv)
      
    
    
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
                    'bandwd' : [50],                    # bandwith d, array('f') with bounds (ndir)
                    'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100],                   # variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1],                     # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1],                     # head var., array('i') with bounds (nvarg)
                    'ivtype' : [1],                     # variogram type, array('i') with bounds (nvarg)
                    'maxclp' : 50000}                   # maximum number of variogram point cloud to use, input int
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
  'bandwd' : [50],                    # bandwith d, array('f') with bounds (ndir)
  'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
  'sills'  : [100],                   # variance used to std the sills, array('f') with bounds (nv)
  'ivtail' : [1],                     # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
  'ivhead' : [1],                     # head var., array('i') with bounds (nvarg)
  'ivtype' : [1],                     # variogram type, array('i') with bounds (nvarg)
  'maxclp' : 50000}                   # maximum number of variogram point cloud to use, input int
}

Notes
-----
Note that data are passed as array like (and not as file path, like in the standalone gslib function gamv).
Any numpy or list may work as input data but make sure that `x,y,z,bhid,vr` have the same dimensions. 
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
                    'bandwd' : [10,10,10],              # bandwith d, array('f') with bounds (ndir)
                    'isill'  : 0,                       # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100],                   # variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1,1,1,1,1,1,1],         # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1,1,1,1,1,1,1],         # head var., array('i') with bounds (nvarg)
                    'ivtype' : [1,3,4,5,6,7,8],         # variogram type, array('i') with bounds (nvarg)
                    'maxclp' : 50000}                   # maximum number of variogram point cloud to use, input int


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
                    'vr'     :  VR,                     # Variables, array('f') with bounds (nd*nv), nv is number of variables
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

    For two variables use flatten array with order C, for example, if mydata is a Pandas DataFrame use 
    
	>> vr=mydata[['U', 'V']].values.flatten(order='FORTRAN')
    >> parameters['vr']=vr


    The output is a tuple of numpy 3D ndarrays (pdis,pgam, phm,ptm,phv,ptv,pnump) with dimensions 
    (nvarg, ndir, nlag+2), representing the  experimental variograms output

    The variables with prefix `p` are for directional variogram and are generated 
    by GSLIB fortran function  `writeout`. This is similar to the output of 
    the GSLIB standalone program `gam`.
    """

    np, gam, hm, tm, hv, tv = __variograms.gamma(**parameters)
    

    # get output like in gslib 
    nvarg = len(parameters['ivtype']) 
    nlag = parameters['nlag']
    ixd = parameters['ixd']
    xsiz = parameters['xsiz']
    iyd = parameters['iyd']
    ysiz = parameters['ysiz']
    izd = parameters['izd']
    zsiz = parameters['zsiz']

    
    pdis,pgam,phm,ptm,phv,ptv,pnump = __variograms.writeout_gam(nvarg,nlag,ixd,xsiz,iyd,ysiz,izd,zsiz,np,gam,hm,tm,hv,tv)
      
    
    
    return pdis,pgam, phm,ptm,phv,ptv,pnump 



# ----------------------------------------------------------------------------------------------------------------
#
#    Rotations
#
#-----------------------------------------------------------------------------------------------------------------
def setrot(ang1=0,ang2=0,ang3=0,anis1=1,anis2=1,ind=1,maxrot=1):
    """Update a rotation and scaling matrix given a set of angles
    
    This is necessary for some interval GSLIB function (ex. cova3)
    
    The rotations are odd, dont use this to rotate data...
    
    TODO: implement code from http://www.ccgalberta.com/ccgresources/report06/2004-403-angle_rotations.pdf
    
     
    Parameters
    ----------
        ang1 : input float, Azimuth angle  (rotation along Z)
        ang2 : input float, Dip angle (downdip positive) (rotation along X)
        ang3 : input float, Plunge angle   (rotation along Z)
        anis1 : input float, First anisotropy ratio (semimayor direction/mayor direction)
        anis2 : input float, Second anisotropy ratio (minor direction/mayor direction)
        ind : input int, index 
        maxrot : input int, number of matrices (for example use 3 for a ZXZ, rotation)
                               
    
    Returns
    -------   
        rotmat : rank-3 array('d') with bounds (maxrot,3,3)
        
    
    """
    return __kt3d.setrot(ang1,ang2,ang3,anis1,anis2,ind,maxrot)


def rotcoord(X,Y,Z,ang1=0,ang2=0,ang3=0,anis1=1,anis2=1,ind=1,invert=0):
    """Rotate and revert rotation of data with 3d Coordinates
    
    This is implemented as in http://www.ccgalberta.com/ccgresources/report06/2004-403-angle_rotations.pdf
    
    Note
    ---------
        This is a pure python function
    
     
    Parameters
    ----------
        ang1 : input float, Azimuth angle  (rotation along Z)
        ang2 : input float, Dip angle (downdip positive) (rotation along X)
        ang3 : input float, Plunge angle   (rotation along Z)
        anis1 : input float, First anisotropy ratio (semimayor direction/mayor direction)
        anis2 : input float, Second anisotropy ratio (minor direction/mayor direction)
        ind : input int, index 
        maxrot : input int, number of matrices (for example use 3 for a ZXZ, rotation)
                               
    
    Returns
    -------   
        rotmat : rank-3 array('d') with bounds (maxrot,3,3)
        
    
    """
    return __kt3d.setrot(ang1,ang2,ang3,anis1,anis2,ind,maxrot)




# ---------------------------------------------------------------------------------------------------------------
#
#    Variograms cova3
#
#-----------------------------------------------------------------------------------------------------------------
def cova3(parameters):
    """Calculate variogram values given a variogram function and direction
       
    The parameters here is a dictionary with the following keys.
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, for example:


             parameters = { 
                'x1'     :  0,            # X coordinates, point 1
                'y1'     :  0,            # Y coordinates, point 1
                'z1'     :  0,            # Z coordinates, point 1
                'x2'     :  0,            # X coordinates, point 2
                'y2'     :  0,            # Y coordinates, point 2
                'z2'     :  0,            # Z coordinates, point 2
                'nst'    :  2,            # number of nested structures, array('i') with bounds (ivarg), 
                                          # ivarg is variogram number
                'it'     :  [3, 3],       # structure type,  array('i') with bounds (ivarg)        
                'c0'     :  [0.1],        # nugget,  array('f') with bounds (ivarg)        
                'cc'     :  [0.4, 0.5],   # variance, array('f') with bounds (nvarg*nst[0])
                'aa'     :  [8, 16],       # parameter a (or range), array('f') with bounds (nvarg*nst[0])
                'irot'   :  0,            # index of the rotation matrix for the first nested structure
                                          # the second nested structure will use irot+1, the third irot+2, and so on
                'rotmat' :  rotmat}       # rotation matrices (output of the funtion setrot)
                               
    
    Returns
    -------   
    out : (cmax,cova)
       
       cmax is the maximum covariance value (required to invert the variogram into covariance function)
       cova is the actual covariance value. 
       
       Note: The variogram can be calculated as v= cmax - cova
       
      
    Notes
    -----
       rotmat may be optained with setrot

    """
      
    # calculate this:  rotmat = gslib.__kt3d.setrot(ang1=0,ang2=0,ang3=0,anis1=1,anis2=1,ind=1,maxrot=1)
      
    cmax,cova = __kt3d.cova3(**parameters)
    
    return cmax,cova



# ----------------------------------------------------------------------------------------------------------------
#
#    Variograms block_covariance
#
#-----------------------------------------------------------------------------------------------------------------
def block_covariance(parameters): 
    """Calculate the block covariance
       
    The parameters here is a dictionary with the following keys.
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, for example:


             parameters = { 
                'xdb'  :  [0 , 0, 1, 1],  # X coordinates of discretazation points
                'ydb'  :  [1 , 0, 1, 0],  # Y coordinates of discretazation points
                'zdb'  :  [0 , 0, 0, 0],  # Z coordinates of discretazation points
                'it'     :  [3, 3],        # structure type,  array('i') with bounds (ivarg)        
                'c0'     :  [0.1],        # nugget,  array('f') with bounds (ivarg)        
                'cc'     :  [0.4, 0.5],   # variance, array('f') with bounds (nvarg*nst[0])
                'aa'     :  [8, 16],       # parameter a (or range), array('f') with bounds (nst)
                'aa1'    :  [8, 16],       # parameter a (or range), array('f') with bounds (nst)
                'aa2'    :  [8, 16],       # parameter a (or range), array('f') with bounds (nst)
                'irot'   :  0,            # index of the rotation matrix for the first nested structure
                                          # the second nested structure will use irot+1, the third irot+2, and so on
                'rotmat' :  rotmat,       # rotation matrices (output of the funtion setrot)
                
                'ang1'   : [0, 0],              # input rank-1 array('d') with bounds (nst)
                'ang2'   : [0, 0],              # input rank-1 array('d') with bounds (nst)
                'ang3'   : [0, 0]}                # input rank-1 array('d') with bounds (nst)
    
    Returns
    -------   
    out : cbb
       
       cbb: float,  is the actual block covariance. 

    """

    unbias,cbb = __block_covariance.block_covariance(**parameters)

    return cbb

# ----------------------------------------------------------------------------------------------------------------
#
#    Kriging kt3d get kriging matrix size
#
# ----------------------------------------------------------------------------------------------------------------
def kt3d_getmatrix_size(ktype,idrif,na):
    """Calculate the number of kriging equations (univariate)
       
    The number of kriging equations for univariate kriging is calculated 
    as number of data (na) number of unbias conditions (mdt)
    
    mdt is calculated depending on the kriging type (ktype) and the 
    number of drift terms (defined in idrif)
    
    Parameters
    ----------
    ktype : input, int. Kriging type, 0=SK,1=OK,2=non-st SK,3=exdrift
    idrif : in/output array('i/boolean') with bounds (9). 
                Indicator of drift term (used if idrif(i)==1)
                The following drift terms are defined: 
                    x,y,z,xx,yy,zz,xy,xz,zy 
                    0,1,2, 3, 4, 5, 6, 7, 8
     
    na : input int. Number of data points used to estimate
    
    Returns
    -------   
    mdt : int.  Number of unbias elements in the kriging matrix
    kneq : int. Number of kriging equations
    error : int. If > 1 there was an error 
                error = 1, (idrif(i) < 0 or  idrif(i) > 1)
                error = 100, if(na >= 1 and na <= mdt)
                     no enough samples to estimate all drift terms

    """
    
    mdt,kneq,error = __kt3d.kt3d_getmatrix_size(ktype,idrif,na)
    
    if (error>0):
        warnings.warn('Error > 0, check your parameters')
    
    return mdt,kneq,error
    
    
# ----------------------------------------------------------------------------------------------------------------
#
#    Kriging kt3d
#
# ----------------------------------------------------------------------------------------------------------------
def kt3d(parameters):
    """Estimate with univariate kriging in a single block/polygon/point
    
    This function is similar to the GSLIB Kt3D function but with some differences:
    
        a) The estimate is in a single point/block or polygon
        b) The the target is difined as a set of discretization points.
           If the number of discretization points is greater than one
           point kriging will be performed, otherwise block kriging 
           will be performed. The size of the block is defined by the 
           location and number of discretization points. The location 
           can be arbitrary, this allows to do polygon kriging or block 
           kriging with random discretization points            
        c) The block covariance is an input parameter
        d) The kriging equations an other debug results are
           returned as arrays
    
    the parameter file here is a dictionary with the following keys.
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, for example:



            parameters = { 
                    'xa'      :  [ 0.,1.],    # data x coordinates, array('f') with bounds (na), na is number of data points
                    'ya'      :  [ 1.,0.],    # data y coordinates, array('f') with bounds (na)
                    'za'      :  [ 0.,0.],    # data z coordinates, array('f') with bounds (na)
                    'vra'     :  [ 1.,9.],    # variable, array('f') with bounds (na)
                    'vea'     :  [ 2.,7.],    # external drift, array('f') with bounds (na)
                    'xdb'     :  [0.],        # discretization points x coordinates, array('f') with bounds (ndb), ndb is number of discretization points
                    'ydb'     :  [0.],        # discretization points y coordinates, array('f') with bounds (ndb)
                    'zdb'     :  [0.],        # discretization points z coordinates, array('f') with bounds (ndb)
                    'extest'  :  5.,          # external drift at target block, 'f' 
                    'cbb'     :  .5,          # block covariance , 'f'. Use block_covariance(parameters) to calculate its value
                    'radius'  :  1.,          # search radius, 'f'. Only used to rescale some values
                    'c0'      :  [0.25],      # nugget,  array('f') with bounds (1)  
                    'it'      :  [2, 3],      # structure type,  array('i') with bounds (nst)          
                    'cc'      :  [0.25, 0.5], # variance, array('f') with bounds (nst)
                    'aa'      :  [5, 30],     # parameter a (or range mayor), array('f') with bounds (nst)
                    'aa1'     :  [5, 30],     # parameter a1 (or range semimayor), array('f') with bounds (nst)
                    'aa2'     :  [5, 30],     # parameter a2 (or range minor), array('f') with bounds (nst)
                    'ang1'    :  [0.,0.],     # rotation angle 1, array('f') with bounds (nst)
                    'ang2'    :  [0.,0.],     # rotation angle 2, array('f') with bounds (nst)
                    'ang3'    :  [0.,0.],     # rotation angle 3, array('f') with bounds (nst)
                    'ktype'   :  1,           # kriging type, 'i' (-0=SK,1=OK,2=non-st SK,3=exdrift)
                    'skmean'  :  0.,          # mean for simple kriging, 'f'
                    'unest'   :  numpy.nan    # value for unestimated, 'f'
                    'idrift'   :  [0,0,0,0,0,0,0,0,0]  # drift terms,  array('i') with bounds (9)     
                                                      # the following drift terms are used
                                                      # x,y,z,xx,yy,zz,xy,xz,zy  
                    'kneq'    :  3}                   # number of kriging equations, 'f'
  
    
    Returns
    -------   
        est : float, estimated value
        estv : float, kriging variance
        estt : float, estimated trend 
        estvt : float, kriging variance of the trend estimate
        w : rank-1 array('d') with bounds (na), weight for the estimate
        wt : rank-1 array('d') with bounds (na), weight for the trend estimate
        error : int, error = 0 if all ok. 
                error = 700   ! rescale factor
                error = 100   ! no enough samples to estimate the drift terms
                error = 900000  ! warning, estimate with one sample
                error = 1     ! allocation error
                error = 10    ! wrong size for the kriging matrix, use kt3d_getmatrix_size to calculate right size
                error = 2     ! deallocation error
                error = 20    ! singular matrix
        kmatrix : rank-2 array('d') with bounds (kneq,kneq)
        kvector : rank-2 array('d') with bounds (1,kneq)
        ksolution : rank-2 array('d') with bounds (1,kneq)
        
    
    """

    est,estv,estt,estvt,w,wt,error,kmatrix,kvector,ksolution = __kt3d.kt3d(**parameters)
    
    if (error>0):
        warnings.warn('Error > 0, check your parameters')

    return est,estv,estt,estvt,w,wt,error,kmatrix,kvector,ksolution


# ----------------------------------------------------------------------------------------------------------------
#
#    Declustering
#
# ----------------------------------------------------------------------------------------------------------------
def declus(parameters):
    """Decluster data and run test with different declustering sizes
    
    The parameter file here is a dictionary with the following keys
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, for example:

            parameters = { 
                    'x'      :  mydata.x,     # data x coordinates, array('f') with bounds (na), na is number of data points
                    'y'      :  mydata.y,     # data y coordinates, array('f') with bounds (na)
                    'z'      :  mydata.z,     # data z coordinates, array('f') with bounds (na)
                    'vr'     :  mydata.vr,    # variable, array('f') with bounds (na)
                    'anisy'  :  5.,           # Y cell anisotropy (Ysize=size*Yanis), 'f' 
                    'anisz'  :  3.,           # Z cell anisotropy (Zsize=size*Zanis), 'f' 
                    'minmax' :  1,            # 0=look for minimum declustered mean (1=max), 'i' 
                    'ncell'  :  10,           # number of cell sizes, 'i' 
                    'cmin'   :  5,            # minimum cell sizes, 'i' 
                    'cmax'   :  50,           # maximum cell sizes, 'i'. Will be update to cmin if ncell == 1
                    'noff'   :  8,            # number of origin offsets, 'i'. This is to avoid local minima/maxima
                    'maxcel' :  100000}       # maximum number of cells, 'i'. This is to avoid large calculations, if MAXCEL<1 this check will be ignored
     
    Returns
    -------     
        wtopt : rank-1 array('d') with bounds (nd), weight value
        vrop  : float, declustered mean
        wtmin : float, weight minimum
        wtmax : float, weight maximum
        error : int,   runtime error: 
                error = 10   ! ncellt > MAXCEL ' check for outliers - increase cmin and/or MAXCEL'
                error = 1    ! allocation error
                error = 2    ! deallocation error
        xinc  : float, cell size increment
        yinc  : float, cell size increment
        zinc  : float, cell size increment
        rxcs  : rank-1 array('d') with bounds (ncell + 1), xcell size 
        rycs  : rank-1 array('d') with bounds (ncell + 1), ycell size 
        rzcs  : rank-1 array('d') with bounds (ncell + 1), zcell size 
        rvrcr : rank-1 array('d') with bounds (ncell + 1), declustering weight  
    
    Note
    -------
    
    Minimun and maximum valid data values are not tested. Filter out vr 
    tmin, tmax values on x,y,z,vr before using this function.
    
    """

    wtopt,vrop,wtmin,wtmax,error,xinc,yinc,zinc,rxcs,rycs,rzcs,rvrcr = __declus.declus(**parameters)
    
    if (error>0):
        warnings.warn('Error > 0, check your parameters')

    return wtopt,vrop,wtmin,wtmax,error,xinc,yinc,zinc,rxcs,rycs,rzcs,rvrcr


#***********************************************************************
# 
# 
#     New functions (non standard GSLIB)
# 
# 
#***********************************************************************
def rotscale(parameters):
    """ Shift, Rotate and rescale a set of 3D coordinates
    
    The parameter file here is a dictionary with the following keys
    
    Parameters
    ----------
        parameters  :  dict
            This is a dictionary with key parameter (case sensitive) and values, for example:

            parameters = { 
                    'x'      :  mydata.x,     # data x coordinates, array('f') with bounds (na), na is number of data points
                    'y'      :  mydata.y,     # data y coordinates, array('f') with bounds (na)
                    'z'      :  mydata.z,     # data z coordinates, array('f') with bounds (na)
                    'x0'     :  0             # new X origin of coordinate , 'f' 
                    'y0'     :  0             # new Y origin of coordinate , 'f'
                    'z0'     :  0             # new Z origin of coordinate , 'f'
                    'ang1'   :  45.,          # Z  Rotation angle, 'f' 
                    'ang2'   :  0.,           # X  Rotation angle, 'f' 
                    'ang2'   :  0.,           # Y  Rotation angle, 'f' 
                    'anis1'  :  1.,           # Y cell anisotropy, 'f' 
                    'anis2'  :  1.,           # Z cell anisotropy, 'f' 
                    'invert' :  0}            # 0 do rotation, <> 0 invert rotation, 'i' 

     
    Returns
    -------     
        xr : rank-1 array('d') with bounds (nd), new X coordinate
        yr : rank-1 array('d') with bounds (nd), new X coordinate
        zr : rank-1 array('d') with bounds (nd), new X coordinate
 
    
    Note
    -------
    This is non standard gslib function and is based on the paper:
        http://www.ccgalberta.com/ccgresources/report06/2004-403-angle_rotations.pdf
        
    The rotation is  {Z counter clockwise ; X clockwise; Y counter clockwise} [-ZX-Y]
    
    """
    xr,yr,zr = __rotscale.rotscale(**parameters)

    return xr,yr,zr


