# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
PyGSLIB.gslib, Module to interface GSLIB programs. This includes the
standard programs available at
http://www.statios.com/software/gslib77_ls.tar.gz and non-standard GSLIB
programs provided by contributors.


Note
----
The original Fortran 77 code was modified. That was required to convert
the existing standalone executables into shared libraries. The original
code was also modified to introduce extra functionality, for example,
maximum number of samples per drillhole in the program ``kt3D`` and
variogram cloud in the program ``gamv``. The main changes consisted of:

   - the old Fortran 77 as converted to Fortran 90.
   - data and parameters are directly transferred between
     Fortran and Python functions and not through data files
     and parameter files.

Copyright 2016, Adrian Martinez Vargas

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE.txt file for details.


"""

# general python modules
import os.path
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import vtk
import sys

#gslib fortran code
from pygslib.gslib import __gslib__kt3d
from pygslib.gslib import __gslib__postik
from pygslib.gslib import __read_gslib
from pygslib.gslib import __addcoord
from pygslib.gslib import __rotscale
from pygslib.gslib import __block_covariance
from pygslib.gslib import __general
from pygslib.gslib import __plot
from pygslib.gslib import __declus
from pygslib.gslib import __variograms
from pygslib.gslib import __dist_transf
from pygslib.gslib import __bigaus
from pygslib.gslib import __bicalib
from pygslib.gslib import __trans
from pygslib.gslib import __draw
from pygslib.gslib import __dm2csv

#pygslib modules
import pygslib.vtktools as vtktools
import pygslib




#-----------------------------------------------------------------------------------------------------------------
#
#    Initialization
#
#-----------------------------------------------------------------------------------------------------------------

# Set default parameters in the Fortran common module
__gslib__kt3d.set_unest (np.nan)
__gslib__postik.set_unest (np.nan)




#-----------------------------------------------------------------------------------------------------------------
#
#    Set nan values
#
#-----------------------------------------------------------------------------------------------------------------
def set_nan_value(value=np.nan):
    """Set non-estimated value in the modules KT3D and POSTIK

    This will change the value assigned to non-estimated values

    Parameters
    ----------
        value  :  numeric (optional, default np.nan)

    Note
    -----
    This will change the behavior of KT3D in the actual python section.


    to see the actual non-estimated value you may call ``__gslib__kt3d.UNEST``

    """

    __gslib__kt3d.set_unest (value)
    __gslib__postik.set_unest (value)




#-----------------------------------------------------------------------------------------------------------------
#
#    General functions
#
#-----------------------------------------------------------------------------------------------------------------

# save drillhole data as gslib file
def df2gslib(data, fname, comment = 'gslib file', decimals = 4):
    '''df2gslib(data, filename, comment = 'gslib file', decimals = 4)

    Save pandas DataFrame to gslib format.

    Parameters
    ----------
    fname : pandas DataFrame
        Dataframe to output
    fname : str
        output file name/path (absolute or relative)

    comment : str dafault ('gslib file')
        commnet line

    decimals : int dafault (4)
        number of decimals

    Warning
    -------
    Non numerical data may produce unexpected results
    The dataframe index is not saved

    '''

    nvar = len(data.columns)
    varnames = ''

    for c in data.columns[0:-1]:
        varnames = varnames + c + '\n'
    varnames = varnames + data.columns[-1]
    tex_data= data.to_string(float_format='%.{}f'.format(decimals), header=False, index = False)

    gslib_data = '{}\n{}\n{}\n{}\n'.format(comment,nvar,varnames,tex_data)

    with open(fname, 'w') as the_file:
        the_file.write(gslib_data)


#read GSLIB file
def read_gslib_file(fname, maxnvar=500):
    """ Reads a gslib file

    The native Fortran code is used to read GEOEAS (gslib) files

    Parameters
    ----------
    fname : str
        file name path (absolute or relative)

    maxnvar : Optional[int]
        maximum number of variables in file

    Returns
    -------
    data : Pandas DataFrame
        The file is returned in memory as a Pandas DataFrame


    Examples
    --------
    >>> import pygslib as gslib
    >>> mydata= gslib.gslib.read_gslib_file('dataset/cluster.dat')
    >>> print (mydata.head(n=5))
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

    assert error == 0, "gslib error = {} at __read_gslib.read_header".format(error)

    #the var names requires this work around in f2py
    vnames=[]
    if sys.version_info[0]==2:                     # python 2.X
        for i in varnames[0:nvar]:
            vnames.append(str.strip(i.tostring()))
    else:                                          # python 3.X
        for i in varnames[0:nvar]:
            a=''
            for j in i:
                a= a+j.decode("utf-8")
            vnames.append(a.strip())


    #now read data
    maxdat,error=__read_gslib.read_ndata(fname, nvar)
    table,error=__read_gslib.read_data(fname, nvar, maxdat)

    assert error==0, "there was an error= %r importing the data" % error

    return pd.DataFrame(table, columns=vnames)


#read datamine table to csv
def dm2csv_ep(fname, outfile, format='F15.4'):
    """ Convert extended precision datamine file to .CSV

    This is an experimental function an may work only if
    the datamine file is in extended precision format.


    Parameters
    ----------
    fname : str
        datamine file name/path (absolute or relative)

    outfile : str
        csv file name/path (absolute or relative)

    format : str dafault ('F15.4')
        string format for numeric values


    Notes
    -----
    All numeric outputs are formatted with the string format.
    This is a limitation and may create non-readable ascii
    outputs. Non-readable outputs are automatically replaced
    with the string '********'.

    Examples
    --------
    >>> import pygslib as gslib
    >>> mydata= gslib.gslib.dm2csv_ep('points.dm', 'points.csv', format='F10.2')
    >>>
    """

    # make sure the file exists
    assert os.path.isfile(fname), "invalid file name"

    __dm2csv.dm2csv_ep(fname, outfile, format)

def dm2csv_sp(fname, outfile, format='F15.4'):
    """ Convert a single precision datamine file to .CSV

    This is an experimental function an may work only if
    the datamine file is in single precision format.


    Parameters
    ----------
    fname : str
        datamine file name/path (absolute or relative)

    outfile : str
        csv file name/path (absolute or relative)

    format : str dafault ('F15.4')
        string format for numeric values


    Notes
    -----
    All numeric outputs are formatted with the string format.
    This is a limitation and may create non-readable ascii
    outputs. Non-readable outputs are automatically replaced
    with the string '********'.

    Examples
    --------
    >>> import pygslib as gslib
    >>> mydata= gslib.gslib.dm2csv_sp('points.dm', 'points.csv', format='F10.2')
    >>>
    """


    # make sure the file exists
    assert os.path.isfile(fname), "invalid file name"

    __dm2csv.dm2csv_sp(fname, outfile, format)

#generate coordinates for gridded files
def addcoord(nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, grid):
    """ Insert X, Y, Z coordinates into a grid file


    Parameters
    ----------
    nx,ny,nz : int
        number of rows.

    xmn,ymn,zmn : double
        origin of coordinates (block centroid).

    xsiz,ysiz,zsiz : double
        block sizes.

    grid : DataFrame
        Pandas dataframe with the input grid variables.

    Returns
    -------
    data : Pandas DataFrame
        Pandas dataframe (the grid) with extra columns x, y, z.

    Note
    ----
    The input grid is expected to have xmn*ymn*zmn rows
    with GEOEAS grid order

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
    """Calculate experimental variograms on sparse data

    This function is similar to the GSLIB gamv function but with some minor differences:

     - the indicator automatic transformation is not applied; you may define
       indicator variables externally.
     - a bhid variable is always required, to avoid downhole variogram use
       ``bhid=None`` or ``bhid = 0``.
     - Z variable is required, you can use constant coordinate to reduce
       dimensions, for example use Z=None in the parameter file or
       Z = 0 in the file.
     - The variogram cloud was implemented but it only works for the first
       variogram and first direction. It only works for some variogram types.


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::


            parameters = {
                    'x'      :  X,             # X coordinates, array('f') with bounds (nd), nd is number of data points
                    'y'      :  Y,             # Y coordinates, array('f') with bounds (nd)
                    'z'      :  Z,             # Z coordinates, array('f') with bounds (nd)
                    'bhid'   :  bhid,          # bhid for downhole variogram, array('i') with bounds (nd)
                    'vr'     :  VR,            # Variables, array('f') with bounds (nd,nv), nv is number of variables
                    'tmin'   : -1.0e21,        # trimming limits, float
                    'tmax'   :  1.0e21,        # trimming limits, float
                    'nlag'   :  10,            # number of lags, int
                    'xlag'   :  4,             # lag separation distance, float
                    'xltol'  :  2,             # lag tolerance, float
                    'azm'    : [0,0,90],       # azimuth, array('f') with bounds (ndir)
                    'atol'   : [90,22.5,22.5], # azimuth tolerance, array('f') with bounds (ndir)
                    'bandwh' : [50,10,10],     # bandwidth 'horizontal', array('f') with bounds (ndir)
                    'dip'    : [0,0,0],        # dip, array('f') with bounds (ndir)
                    'dtol'   : [10,10,10],     # dip tolerance, array('f') with bounds (ndir)
                    'bandwd' : [10,10,10],     # bandwidth 'vertical', array('f') with bounds (ndir)
                    'isill'  : 0,              # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100],          # variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1,1,1,1,1,1,1],# tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1,1,1,1,1,1,1],# head var., array('i') with bounds (nvarg)
                    'ivtype' : [1,3,4,5,6,7,8],# variogram type, array('i') with bounds (nvarg)
                    'maxclp' : 50000}          # maximum number of variogram point cloud to use, input int

    Warnings
    --------
    bhid must be an array of integers

    Returns
    -------
       pdis :  Distance of pairs falling into this lag
       pgam :  Semivariogram, covariance, correlogram,... value
       phm  :  Mean of the tail data
       ptm  :  Mean of the head data
       phv  :  Variance of the tail data
       ptv  :  Variance of the head data
       pnump:  Number of pairs
       cldi :  data index of the head
       cldj :  data index of the tail
       cldg :  Semivariogram, covariance, ... value
       cldh :  Distance of each pair


    Note
    -----
    The output variables with prefix *cld* are for variogram cloud and
    with prefix *p* are for directional variograms

    The variables with prefix *p* are similar to the output generated
    by the GSLIB standalone program **gamv**.

    The variogram cloud only works for the first variogram/direction
    and only if the variogram of type 1, 2, 6, 7 and 8.

    The variogram types are:

     - traditional semivariogram (1)
     - traditional cross semivariogram (2)
     - pairwise relative semivariogram (6)
     - semivariogram of logarithms (7)
     - semimadogram(8)

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
    """Check the variogram calculation parameters

    Verifies the parameters of the ``gamv`` function
    and checks for possible errors or inconsistencies

    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters


    The dictionary with parameters may be as follows::


             parameters = {
               # Coordinates of the data points
                'x'      :  X,              # array('f') with bounds (nd), nd is number of data points
                'y'      :  Y,              # array('f') with bounds (nd)
                'z'      :  Z,              # array('f') with bounds (nd)
                'bhid'   :  bhid,           # drillhole id for downhole variogram, array('i') with bounds (nd)
                'vr'     :  VR,             # Variables, array('f') with bounds (nd,nv), nv is number of variables
                'tmin'   : -1.0e21,         # trimming limits, float
                'tmax'   :  1.0e21,         # trimming limits, float
                'nlag'   :  10,             # number of lags, int
                'xlag'   :  4,              # lag separation distance, float
                'xltol'  :  2,              # lag tolerance, float
                'azm'    : [0,0,90],        # azimuth, array('f') with bounds (ndir)
                'atol'   : [90,22.5,22.5],  # azimuth tolerance, array('f') with bounds (ndir)
                'bandwh' : [50,10,10],      # bandwidth h, array('f') with bounds (ndir)
                'dip'    : [0,0,0],         # dip, array('f') with bounds (ndir)
                'dtol'   : [10,10,10],      # dip tolerance, array('f') with bounds (ndir)
                'bandwd' : [10,10,10],      # bandwidth d, array('f') with bounds (ndir)
                'isill'  : 0,               # standardize sills? (0=no, 1=yes), int
                'sills'  : [100],           # variance used to std the sills, array('f') with bounds (nv)
                'ivtail' : [1,1,1,1,1,1,1], # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                'ivhead' : [1,1,1,1,1,1,1], # head var., array('i') with bounds (nvarg)
                'ivtype' : [1,3,4,5,6,7,8], # variogram type, array('i') with bounds (nvarg)
                'maxclp' : 50000}           # maximum number of variogram point cloud to use, input int


    Returns
    -------
    out : bool
      return True (or 1) if the parameter dictionary is valid

    Example
    -------
    >>> assert pygslib.gslib.check_gamv_par(parameters)==1 , 'sorry this parameter file is wrong'
    >>>

    """

    #check the parameters are correct
    #assert set(gamv_parameter_template)==set(parameters)

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
        assert parameters['bandwh'][i]>=0, "Horizontal bandwidth is too small!"
        assert parameters['bandwd'][i]>=0, "Vertical bandwidth is too small!"

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
    """Calculates experimental variograms on gridded data


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::


            parameters = {
                    'nx'     :  50,       # number of rows in the gridded data
                    'ny'     :  50,       # number of columns in the gridded data
                    'nz'     :  1,        # number of levels in the gridded data
                    'xsiz'   :  1,        # size of the cell in x direction
                    'ysiz'   :  1,        # size of the cell in y direction
                    'zsiz'   :  1,        # size of the cell in z direction
                    'bhid'   :  bhid,     # bhid for downhole variogram, array('i') with bounds (nd)
                    'vr'     :  VR,       # Variables, array('f') with bounds (nd*nv), nv is number of variables
                    'tmin'   : -1.0e21,   # trimming limits, float
                    'tmax'   :  1.0e21,   # trimming limits, float
                    'nlag'   :  10,       # number of lags, int
                    'ixd'    : [1,0],     # direction x
                    'iyd'    : [1,0],     # direction y
                    'izd'    : [1,0],     # direction z
                    'isill'  : 0,         # standardize sills? (0=no, 1=yes), int
                    'sills'  : [100, 200],# variance used to std the sills, array('f') with bounds (nv)
                    'ivtail' : [1,1,2,2], # tail var., array('i') with bounds (nvarg), nvarg is number of variograms
                    'ivhead' : [1,1,2,2], # head var., array('i') with bounds (nvarg)
                    'ivtype' : [1,3,1,3]} # variogram type, array('i') with bounds (nvarg)


    Returns
    -------
       pdis :  Distance of pairs falling into this lag
       pgam :  Semivariogram, covariance, correlogram,... value
       phm  :  Mean of the tail data
       ptm  :  Mean of the head data
       phv  :  Variance of the tail data
       ptv  :  Variance of the head data
       pnump:  Number of pairs


    Note
    -----

    For two variables use flatten array with order F, for example::

        vr=mydata[['U', 'V']].values.flatten(order='FORTRAN')
        parameters['vr']=vr


    The output is a tuple of numpy 3D ndarrays (pdis,pgam, phm,ptm,phv,
    ptv,pnump) with dimensions (nvarg, ndir, nlag+2), representing the
    experimental variograms output.

    The output is similar to the one generated by the
    GSLIB standalone program `gam`.
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

#-----------------------------------------------------------------------------------------------------------------
#
#    Variograms in all 3D directions
#
#-----------------------------------------------------------------------------------------------------------------
def gamv3D(parameters):
    """Calculates experimental variograms in all 3D directions


    Parameters
        ----------
        parameters = {
            'x' : x,                   # X coordinates, array('f')
            'y' : y,                   # Y coordinates, array('f')
            'z' : z,                   # Z coordinates, array('f')
            'bhid' : zone,             # bhid or zone for downhole variogram or fault areas, array('i')
            'vr' :   variables,        # variables, array('f','f') with bounds (x.shape[0], variables.shape[1])
            'tmin  : -1.0e21,          # trimming limits, float  (optional, default -1.0e21)
            'tmax' : 1.0e21,           # trimming limits, float  (optional, default  1.0e21)
            'nlag' : 10,               # number of lags, int
            'xlag' : 4,                # lag separation distance, float
            'ndir' : 18,               # number of directions, int
            'ndip' : 18,               # number of dip directions, int
            'orgdir' :                 # clockwise rotation added to dir
            'orgdip' :                 # clockwise rotation added down dip
            'isill' : 1,               # standardize sills? (0=no, 1=yes), int 
            'sills' : [100],           # variance used to std the sills, array('f') with bounds (variables.shape[1])
            'ivtail' : [1,1,1,1,1,1,1],# tail varariable for each variogram, array('i') 
            'ivhead' : [1,1,1,1,1,1,1],# head varariable for each variogram, array('i') 
            'ivtype' : [1,3,4,5,6,7,8],# variogram type for each variogram, array('i')
        }

        Returns
        -------
        npair : number of pair, rank-4 array('d') with bounds (nlag,ndir,ndip,nvarg)
        dis   : average distance, rank-4 array('d') with bounds (nlag,ndir,ndip,nvarg)
        gam   : variogram value, rank-4 array('d') with bounds (nlag,ndir,ndip,nvarg)
        hm    : head mean, rank-4 array('d') with bounds (nlag,ndir,ndip,nvarg)
        tm    : tail mean, rank-4 array('d') with bounds (nlag,ndir,ndip,nvarg)
        hv    : head variance, rank-4 array('d') with bounds (nlag,ndir,ndip,nvarg)
        tv    : tail variance, rank-4 array('d') with bounds (nlag,ndir,ndip,nvarg)

        these are 4D arrays with index (lag,dir,dip,variogram). It may fit in a i,j,k structured vtk grid with 7*nvariograms data fields

    Note
    -----


    """

    if 'tmin' not in parameters.keys():
        parameters['tmin'] = -1.0e21

    if 'tmax' not in parameters.keys():
        parameters['tmin'] =  1.0e21

    #get 4D arrays with dimesions (nlag,ndir,ndip,nvarg)
    npair,dis,gam,hm,tm,hv,tv = __variograms.gamv3d(**parameters)



    return np,dis,gam,hm,tm,hv,tv

# ----------------------------------------------------------------------------------------------------------------
#
#    Rotations
#
#-----------------------------------------------------------------------------------------------------------------
def setrot(ang1=0,ang2=0,ang3=0,anis1=1,anis2=1,ind=1,maxrot=1):
    """Updates a rotation and scaling matrix given a set of angles

    This is necessary for some interval GSLIB function (ex. cova3)

    Warning
    -------
    The rotations are peculiar, do not use this to rotate data.


    Parameters
    ----------
        ang1 : input float, Azimuth angle  (rotation along Z)
        ang2 : input float, Dip angle (downdip positive) (rotation along X)
        ang3 : input float, Plunge angle   (rotation along Z)
        anis1 : input float, First anisotropy ratio (semi mayor direction/mayor direction)
        anis2 : input float, Second anisotropy ratio (minor direction/mayor direction)
        ind : input int, index
        maxrot : input int, number of matrices (for example use 3 for a ZXZ, rotation)


    Returns
    -------
        rotmat : rank-3 array('d') with bounds (maxrot,3,3)


    """
    return __general.setrot(ang1,ang2,ang3,anis1,anis2,ind,maxrot)


def rotcoord(ang1=0,ang2=0,ang3=0,anis1=1,anis2=1,ind=1,invert=0):
    """Generates rotation matrix

    This is the rotation matrix used internally in standard
    GSLIB programs


    Parameters
    ----------
        ang1 : input float, Azimuth angle  (rotation along Z)
        ang2 : input float, Dip angle (downdip positive) (rotation along X)
        ang3 : input float, Plunge angle   (rotation along Z)
        anis1 : input float, First anisotropy ratio (semi mayor direction/mayor direction)
        anis2 : input float, Second anisotropy ratio (minor direction/mayor direction)
        ind : input int, index
        maxrot : input int, number of matrices (for example use 3 for a ZXZ, rotation)


    Returns
    -------
        rotmat : rank-3 array('d') with bounds (maxrot,3,3)


    """
    return __general.setrot(ang1,ang2,ang3,anis1,anis2,ind,maxrot)




# ---------------------------------------------------------------------------------------------------------------
#
#    Variograms cova3
#
#-----------------------------------------------------------------------------------------------------------------
def cova3(parameters):
    """Calculates variogram values given a variogram function and direction


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::

             parameters = {
               # Coordinates of point 1 and 2
                'x1'    :  0,         # X coordinates, point 1, float
                'y1'    :  0,         # Y coordinates, point 1, float
                'z1'    :  0,         # Z coordinates, point 1, float
                'x2'    :  1,         # X coordinates, point 2, float
                'y2'    :  0,         # Y coordinates, point 2, float
                'z2'    :  0,         # Z coordinates, point 2, float
               # Variogram model
                'nst'   :  2,         # number of nested structures
                'it'    :  [3, 3],    # structure type,  array('i') with bounds (ivarg)
                'c0'    :  [0.1],     # nugget,  array('f') with bounds (ivarg)
                'cc'    :  [1, 1.4],  # variance, array('f') with bounds (nvarg*nst[0])
                'aa'    :  [8, 22],   # parameter a (or range), array('f') with bounds (nvarg*nst[0])
                'irot'  :  1,         # index of the rotation matrix for the first nested structure
                                      # the second nested structure will use irot+1, the third irot+2, and so on
                'rotmat':  rotmat}    # rotation matrices (output of the function setrot)


    Returns
    -------
    cmax,cova : float64, float64

       cmax is the maximum covariance value (required to invert the variogram into covariance function)
       cova is the actual covariance value.

       Note: The variogram can be calculated as v= cmax - cova


    Note
    -----
       rotmat may be obtained with setrot

    Todo:
        calculate irot and rotmat internally and remove it from parameter dict.

    Example
    --------
    >>> # this is the covariance between the points x1, x2
    >>> cmax,cova=pygslib.gslib.cova3(parameters_mod)
    >>> print (cmax, cova)
    >>> 2.40999984741 2.24302244186

    """

    # calculate this:  rotmat = gslib.__kt3d.setrot(ang1=0,ang2=0,ang3=0,anis1=1,anis2=1,ind=1,maxrot=1)

    cmax,cova = __general.cova3(**parameters)

    return cmax,cova



# ----------------------------------------------------------------------------------------------------------------
#
#    Variograms block_covariance
#
#-----------------------------------------------------------------------------------------------------------------
def block_covariance(parameters):
    """Calculates block covariance

    This function calculates the block covariance given a variogram
    model.
    The anisotropy can be different in each nested structure.
    The block size is defined by input discretization points. Note
    that discretization points can have any arbitrary locations,
    for example

     - regular discretization
     - random  discretization
     - or discretization in an irregular polygon


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters


    Parameters are E-Type in a dictionary as follows::

             parameters = {
               # Coordinates of the discretization points
                'xdb'  :  [0, 0, 1, 1], # array('f')
                'ydb'  :  [1, 0, 1, 0], # array('f')
                'zdb'  :  [0, 0, 0, 0], # array('f')
               # Variogram model
                'it'   :  [3, 2],       # structure type, array('i')  with bounds (nst)
                'c0'   :  [0.1],        # nugget, array('f')
                'cc'   :  [0.4, 0.5],   # variance, array('f')  with bounds (nst)
                'aa'   :  [8, 16],      # ranges in main direction, array('f')   with bounds (nst)
                'aa1'  :  [8, 16],      # ranges in second direction, array('f') with bounds (nst)
                'aa2'  :  [8, 16],      # ranges in third direction, array('f')  with bounds (nst)
                'ang1' : [0, 0],        # input rank-1 array('d') with bounds (nst)
                'ang2' : [0, 0],        # input rank-1 array('d') with bounds (nst)
                'ang3' : [0, 0]}        # input rank-1 array('d') with bounds (nst)

    Returns
    -------
    cbb : float
       block covariance.

    Example
    -------
    >>> print (pygslib.gslib.block_covariance(parameters_blk))
    0.803182760643

    """

    unbias,cbb = __block_covariance.block_covariance(**parameters)

    return cbb


# ----------------------------------------------------------------------------------------------------------------
#
#    Declustering
#
# ----------------------------------------------------------------------------------------------------------------
def declus(parameters):
    """Decluster data and run declustering test with different
    declustering cell sizes


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::

            parameters = {
                   # database
                    'x'      :  mydata.x,     # data x coordinates, array('f') with bounds (na), na is number of data points
                    'y'      :  mydata.y,     # data y coordinates, array('f') with bounds (na)
                    'z'      :  mydata.z,     # data z coordinates, array('f') with bounds (na)
                    'vr'     :  mydata.vr,    # variable, array('f') with bounds (na)
                   # cell size
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
        rvrcr : rank-1 array('d') with bounds (ncell + 1), declustering mean

    Note
    -------

    Minimun and maximum valid data values are not tested. Filter out ``vr``
    values before using this function.

    TODO
    ----
    Create nice output for errors


    """

    assert parameters['ncell']>0, 'Error, parameter ncell<=0'


    wtopt,vrop,wtmin,wtmax,error,xinc,yinc,zinc,rxcs,rycs,rzcs,rvrcr = __declus.declus(**parameters)

    if (error>0):
        print ('There was an Error = {}, check your parameters and see __doc__'.format(error))
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
    """ Rotates and rescales a set of 3D coordinates

    The new rotated and rescaled system of coordinates will have
    origin at [x = 0, y = 0, z = 0]. This point corresponds to
    [x0,y0,z0] (the pivot point) in the original system of coordinates.

    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters


    Parameters are pased in a dictionary as follows::

            parameters = {
                    'x'      :  mydata.x,# data x coordinates, array('f') with bounds (na), na is number of data points
                    'y'      :  mydata.y,# data y coordinates, array('f') with bounds (na)
                    'z'      :  mydata.z,# data z coordinates, array('f') with bounds (na)
                    'x0'     :  0,       # pivot point coordinate X, 'f'
                    'y0'     :  0,       # pivot point coordinate Y, 'f'
                    'z0'     :  0,       # pivot point coordinate Z, 'f'
                    'ang1'   :  45.,     # Z  Rotation angle, 'f'
                    'ang2'   :  0.,      # X  Rotation angle, 'f'
                    'ang3'   :  0.,      # Y  Rotation angle, 'f'
                    'anis1'  :  1.,      # Y cell anisotropy, 'f'
                    'anis2'  :  1.,      # Z cell anisotropy, 'f'
                    'invert' :  0}       # 0 do rotation, <> 0 invert rotation, 'i'


    Returns
    -------
        xr : rank-1 array('d') with bounds (nd), new X coordinate
        yr : rank-1 array('d') with bounds (nd), new X coordinate
        zr : rank-1 array('d') with bounds (nd), new X coordinate


    Note
    -------
    This is nonstandard gslib function and is based on the paper:

        http://www.ccgalberta.com/ccgresources/report06/2004-403-angle_rotations.pdf

    The rotation is  {Z counter clockwise ; X clockwise; Y counter clockwise} [-ZX-Y]

    """
    xr,yr,zr = __rotscale.rotscale(**parameters)

    return xr,yr,zr


#-----------------------------------------------------------------------------------------------------------------
#
#    Kriging with Kt3D
#
#-----------------------------------------------------------------------------------------------------------------

def kt3d(parameters):
    """Estimates with the GSLIB program KT3D

    This is a wrap for the GSLIB Fortran code of the program KT3D Version 2.0,
    originally in Fortran 77. Only minor changes were included, the most
    relevant are:

     - support for maximum number of samples per drillhole was
       implemented
     - the direct file output was redirected to numpy arrays
     - the input (for grid estimate) is now as numpy arrays.
       The grid definition was retained because is an input
       of the GSLIB search super-block algorithm.
     - the trimming limits were removed; you may filter out
       undesired values before estimating

    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters


    The dictionary with parameters may be as follows::


        kt3d_parameters = {
            # Input Data
            # ----------
            'x' : ,   # 1D array('f'), x coordinates of input data
            'y' : ,   # 1D array('f'), y coordinates of input data
            'z' : ,   # 1D array('f'), z coordinates of input data
            'vr' : ,   # 1D array('f'), primary variable
            've' : ,   # (optional) 1D array('f'), extra variable, external drift
            'bhidint': , #  (optional) 1D array('i'), drillhole ID, or fault zone ID
            # Output grid definition (only used for superblock search, can use dummy values)
            # ----------
            'nx' : ,   # int, size of the grid
            'ny' : ,   # int,
            'nz' : ,   # int,
            'xmn' : ,   # float, origing of coordinate
            'ymn' : ,   # float,
            'zmn' : ,   # float,
            # Output data
            # ----------
            'xsiz' : ,   # float, block size in x
            'ysiz' : ,   # float, block size in x
            'zsiz' : ,   # float, block size in x
            'nxdis' : ,   # int, number of discretization points in x
            'nydis' : ,   # int, number of discretization points in y
            'nzdis' : ,   # int, number of discretization points in z
            'outx' : ,   # 1D array('f'), x coordinates of output data
            'outy' : ,   # 1D array('f'), y coordinates of output data
            'outz' : ,   # 1D array('f'), z coordinates of output data
            'outextve' : ,   # 1D array('f') external drift in output data
            # Search parameters
            # ----------
            'radius'     : ,   # float, search distance in direction 1
            'radius1'    : ,   # float, search distance in direction 2
            'radius2'    : ,   # float, search distance in direction 3
            'ndmax'      : ,   # int, maximum number of points
            'ndmin'      : ,   # int, minimum number of points
            'noct'       : ,   # (optional) int, maximum number of samples per octant, if zero octans not used, default == 0
            'nbhid'      : ,   # (optional) int, maximum number of samples from each drillhole, default not used == 0
            'sang1'      : ,   # (optional) float, rotation angle 1
            'sang2'      : ,   # (optional) float, rotation angle 1
            'sang3'      : ,   # (optional) float, rotation angle 1
            # Kriging parameters and options
            # ----------
            'idrif'      : ,   # (optional) array('i'), array of 9 indicators to use drift models
            'itrend'     : ,   # (optional) int, if == 1 the trend will be estimated
            'ktype'      : ,   # (optional) int, kriging type: 0 skmean, 1 ordinary kriging, 2 simple kriging with local means, 3 kriging with an external drift
            'skmean'     : ,   # (optional) float, simple kriging mean for ktype = 0
            'koption'    : ,   # (optional) int, set to 0 for kriging a grid of points or blocks, to 1 for cross validation with the data in datafl and to 2 for jackknifing
            'iktype'     : ,   # (optional) int, median indicator kriging, default False = 0
            'cut'        : ,   # (optional) 1D array('f'), thresholds for median indicator kriging, default = []
            'idbg'       : ,   # (optional) int, debug ?, defaul False == 0
            # Inverse of the power of the distance parameter
            'id2power'   : ,   # (optional) float, inverse of the distance powe, defaul 2
            # Variogram parameters
            # ----------
            'c0'         : ,   # float, nugget value
            'it'         : ,   # array('i'), structures type, on for each structure: 1 Spherical, 2 Exponential, 3 Gaussian, 4 Power, 5 Cosine hole effect
            'cc'         : ,   # array('f'), structures variance, one for each structure
            'aa'         : ,   # array('f'), structures range/practical range in direction 1, one for each structure
            'aa1'        : ,   # array('f'), structures range/practical range in direction 2, one for each structure
            'aa2'        : ,   # array('f'), structures range/practical range in direction 3, one for each structure
            'ang1'       : ,   # (optional) array('f'), rotation angle 1, one for each structure, defaul array of zeros
            'ang2'       : ,   # (optional) array('f'), rotation angle 2, one for each structure, defaul array of zeros
            'ang3'       : }   # (optional) array('f'), rotation angle 3, one for each structure, defaul array of zeros


    the nine idrif[1,...,8] indicators terms are
        idrif[0] linear drift in x
        idrif[1] linear drift in y
        idrif[2] linear drift in z
        idrif[3] quadratic drift in x
        idrif[4] quadratic drift in y
        idrif[5] quadratic drift in z
        idrif[6] cross quadratic drift in xy
        idrif[7] cross quadratic drift in xz
        idrif[8] cross quadratic drift in yz


    Returns
    -------
    output :  dict, estimation results at target points/blocks
    debug  :  debug output for the last block estimated
    estimate : estimation summary


    The estimation results ``output`` will be as follows:

    If iktype=0 returns variable estimate as::

        {'outest'    : ,  # kriging estimate
        'outkvar'    : ,  # kriging variance
        'outidpower' : ,  # inverse of the power of the distance estimate
        'outnn'      : }  # nearest neightbor estimate

        all these are rank-1 array('f') with bounds (nout)

    If iktype=1 returns median indicator kriging estimate as::

       {'outcdf'     : }  #rank-2 array('f') with bounds (nout,ncut)


    The debug output ``debug`` will be only for the last block estimated
    but calculated for all the blocks. Make sure to set ``idbg``
    to zero for estimation in large models.

    If idbg = 0 ``debug`` will be an empty dictionary.
    If debug > 0 ``debug`` will be as follows::

        {'cbb'       : ,  #float
        'neq'        : ,  #int
        'na'         : ,  #int
        'dbgxdat'    : ,  #rank-1 array('f') with bounds (ndmax)
        'dbgydat'    : ,  #rank-1 array('f') with bounds (ndmax)
        'dbgzdat'    : ,  #rank-1 array('f') with bounds (ndmax)
        'dbgvrdat'   : ,  #rank-1 array('f') with bounds (ndmax)
        'dbgwt'      : ,  #rank-1 array('f') with bounds (ndmax + 11)
        'dbgxtg'     : ,  #float
        'dbgytg'     : ,  #float
        'dbgztg'     : ,  #float
        'dbgkvector' : ,  #rank-1 array('f') with bounds (ndmax + 11)
        'dbgkmatrix' : ,
        'ellipsoid'  :  } #vtkpolydata object with ellipsoid

    The output ``estimate`` is to provide an estimation summary
    but this is not implemented yet.


    Note
    -----
    Not all the parameters are used for calculation, it depends on the
    kriging type and options

    Optional parameters can be removed, knowing that kt3D will create
    internal arrays/variable initialized to zero value

    If using nbhid > 0 the hole id number (bhid) is required. Hole
    IDs may be integers from one to total number of drillholes.
    Use function pygslib.Drillhole.txt2intID(table_name) to get
    a correct bhid number.

    Some important stats and run progress are only available in the
    stdout (the terminal) and will not be available in Ipython Notebooks.

    """

    # add dummy cut if not in parameters

    if 'cut' not in parameters:
        parameters['cut'] =[0]

    if 'cut' in parameters and parameters['cut']==None:
        parameters['cut'] =[0]

    if 'id2power' not in parameters:
        parameters['id2power']=2.0

    assert parameters['id2power'] >= 0, 'Error: parameter id2power <0 '

    assert parameters['x'].shape[0]>0, 'Error: input array x with size zero'

    assert parameters['outx'].shape[0]>0, 'Error: output array x with zero size'


    # check that bhid is provided if nbhid > 0
    if 'nbhid' in parameters:
        if parameters['nbhid']> 0 :
            # make sure exists
            assert parameters['bhidint'] is not None, 'Error: BHID required if nbhid > 0'

            # make sure there is no drillhole number equal to zero
            assert 0 not in parameters['bhidint'], 'Error: BHID == 0 detected, BHIDs may be > 0 '

            # make sure the is bhid are correct
            assert max(parameters['bhidint'])<=len(parameters['bhidint']), 'Error: bhidint > ndata detected, bhidint must be in interval [1,ndata] '


        # check not using octants and min
        if (parameters['nbhid']> 0 and parameters['noct']> 0):
            #pass
            warnings.warn('\nWarning: !!!!! Using octants and maximum ' + \
                          ' number of samples per drillholes at the same '  + \
                          ' time may produce unexpected results !!!!!!')


    # prepare the output
    output = {}
    debug  = {}
    estimate={}

    try:

        (output['outest'],
        output['outkvar'],
        output['outcdf'],
        output['cbb'],
        output['neq'],
        output['na'],
        output['dbgxdat'],
        output['dbgydat'],
        output['dbgzdat'],
        output['dbgvrdat'],
        output['dbgwt'],
        output['dbgxtg'],
        output['dbgytg'],
        output['dbgztg'],
        output['dbgkvector'],
        output['dbgkmatrix'],
        error,
        fwarnings,
        output['outidpower'],
        output['outnn'],
        output['outlagr'],
        output['outwmean']) =__gslib__kt3d.pykt3d(**parameters)

    except:

        raise NameError('unexpected error in FORTRAN kt3d code')


    if len(error.strip())>0:
        raise NameError(str(error).strip())

    if len(fwarnings.strip())>0:
        warnings.warn(str(fwarnings).strip())

    estimate['x']=parameters['outx']
    if 'outy' in parameters:
        estimate['y']=parameters['outy']
    if 'outz' in parameters:
        estimate['z']=parameters['outz']

    estimate['outidpower']=output['outidpower']
    estimate['outnn']=output['outnn']

    # clean a bit the output
    if 'iktype' in parameters:
        if parameters['iktype']==1:
            estimate['outcdf'] = output['outcdf']
        else:
            estimate['outest']=output['outest']
            estimate['outkvar']=output['outkvar']
    else:
            estimate['outest']=output['outest']
            estimate['outkvar']=output['outkvar']


    if parameters['ktype']==1:
        estimate['outlagrange'] = output['outlagr']

    if parameters['ktype']==0 or parameters['ktype']==2:
        estimate['outweightofmean'] = output['outwmean']

    debug ={}
    if 'idbg' in parameters:
        if parameters['idbg']>0:
            for key in ['cbb','neq', 'na','dbgxtg', 'dbgytg', 'dbgztg']:
                debug[key]=output[key]
            debug['dbgxdat'] = output['dbgxdat'][:output['na']]
            debug['dbgydat'] = output['dbgydat'][:output['na']]
            debug['dbgzdat'] = output['dbgzdat'][:output['na']]
            debug['dbgvrdat'] = output['dbgvrdat'][:output['na']]
            debug['dbgwt'] = output['dbgwt'][:output['na']]
            debug['dbgkvector'] = output['dbgkvector'][:output['neq']]
            debug['dbgkmatrix'] = output['dbgkmatrix'][:output['neq'],:output['neq']]
            debug['lagrange'] = output['dbgwt'][output['na']:output['neq']]



            #plot data and weight

            fig, ax= plt.subplots()

            plt.axis((debug['dbgxtg'] - parameters['radius'],
                          debug['dbgxtg'] + parameters['radius'],
                          debug['dbgytg'] - parameters['radius'],
                          debug['dbgytg'] + parameters['radius']))

            plt.plot(parameters['x'],parameters['y'], '.', color = '0.25')
            plt.scatter(debug['dbgxdat'], debug['dbgydat'], s=100+debug['dbgwt'], c=debug['dbgwt'], alpha=0.6)
            plt.plot (debug['dbgxtg'],debug['dbgytg'], 'r*')

            plt.colorbar()

            debug['plotxy'] = ax.get_figure()


            # get Ellipsoid in vtk format
            # a) we generate an sphere with center at (0,0,0) and radius a
            source = vtk.vtkSphereSource()
            source.SetThetaResolution (32)
            source.SetPhiResolution (32)
            source.SetRadius(parameters['radius'])
            source.Update()
            ellipsoid = source.GetOutput()

            # now we rotate the points... gslib rotations
            # a) get points and rotate + rescale

            p = vtktools.GetPointsInPolydata(ellipsoid)

            # rescall and rotate

            anis1 = float(parameters['radius1'])/parameters['radius']
            anis2 = float(parameters['radius2'])/parameters['radius']

            p[:,1] = p[:,1]*anis1
            p[:,2] = p[:,2]*anis2

            ang1 = parameters['sang1']
            ang2 = parameters['sang2']
            ang3 = parameters['sang3']

            # this emulates the GSLIB set rottion matrix function
            rotmat = __setrotmatrix__(ang1,ang2,ang3, 1., 1.)


            # Transposing the rotation matrix produces back transformation
            pp = np.dot(rotmat.T, p.T).T


            # b) translate
            xr = pp[:,0] + debug['dbgxtg']
            yr = pp[:,1] + debug['dbgytg']
            zr = pp[:,2] + debug['dbgztg']

            p[:,0] = xr
            p[:,1] = yr
            p[:,2] = zr

            # update ellipsoid
            vtktools.SetPointsInPolydata(ellipsoid, p)

            debug['ellipsoid'] = ellipsoid

    #TODO: Implement the estimation summary

    summary = {}

    return estimate, debug, summary


def __setrotmatrix__(ang1,ang2,ang3, anis1, anis2):
    """
    Internal function

    This emulates the setrot.f90 function

    """
    rotmat = np.zeros([3,3])

    # constants
    DEG2RAD=3.141592654/180.0
    EPSLON=1.e-20

    # set angles in radians
    if ang1 >= 0.0 and ang1 < 270.0:
        alpha = (90.0   - ang1) * DEG2RAD
    else:
        alpha = (450.0  - ang1) * DEG2RAD

    beta  = -1.0 * ang2 * DEG2RAD
    theta =        ang3 * DEG2RAD


    # Get the required sines and cosines:

    sina  = np.sin(alpha)
    sinb  = np.sin(beta)
    sint  = np.sin(theta)
    cosa  = np.cos(alpha)
    cosb  = np.cos(beta)
    cost  = np.cos(theta)

    # Construct the rotation matrix in the required memory:

    afac1 = 1.0 / max(anis1,EPSLON)
    afac2 = 1.0 / max(anis2,EPSLON)
    rotmat[0,0] =       (cosb * cosa)
    rotmat[0,1] =       (cosb * sina)
    rotmat[0,2] =       (-sinb)
    rotmat[1,0] = afac1*(-cost*sina + sint*sinb*cosa)
    rotmat[1,1] = afac1*(cost*cosa + sint*sinb*sina)
    rotmat[1,2] = afac1*( sint * cosb)
    rotmat[2,0] = afac2*(sint*sina + cost*sinb*cosa)
    rotmat[2,1] = afac2*(-sint*cosa + cost*sinb*sina)
    rotmat[2,2] = afac2*(cost * cosb)

    return rotmat

#-----------------------------------------------------------------------------------------------------------------
#
#    PostIK
#
#-----------------------------------------------------------------------------------------------------------------

def postik(parameters):
    """Postprocess indicator kriging output (cdf)


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::

        postik_parameters = {
            # output option, output parameter
            'iout'   : ,   # int. 1 E-type,2 P and means outpar,3 p-quantile for outpar=p, and 4 conditional variance
            'outpar' : ,   # float. Parameter for iout (ie cutoff to calculate P)
            # the thresholds
            'ccut1'  : ,   # 1D array of floats. Cutoff used in MIK
            # volume support?, type, varred
            'ivol'   : ,   # int. If 1 the support correction is applied
            'ivtyp'  : ,   # int. 1 for affine correction and indirect lognormal correction
            'varred' : ,   # float. Volumen correction, ussually r~ Block variance/Point variance
            # minimum and maximum Z value
            'zmin'   : ,   # float. Minimum value in local CDF
            'zmax'   : ,   # float. Maximum value in local CDF
            # lower,middle and upper tail: option, parameter
            'ltail'  : ,   # int. Lower tail interpolation function, 1 linear, 2  power model, 3 tabulated quantiles
            'ltpar'  : ,   # float. Lower tail function parameter
            'middle'  : ,  # int. Middle CDF interpolation function 1 linear, 2  power model, 3 tabulated quantiles
            'mpar'  : ,    # float.  Middle CDF segment function parameter
            'utail'  : ,   # int. Upper tail interpolation function, 1 linear, 2  power model, 3 tabulated quantiles, 4 hyperbolic
            'utpar'  : ,   # float. Uper tail function parameter
            # maximum discretization
            'maxdis' : ,   # int. Discretization of the local CDF.
            # 1D arrays with global distribution
            'vr'     : ,   # 1D array of floats for table look-up if tabulated quantiles are used as interpolation function
            'wt'     : ,   # 1D array of floats with wights on table look-up
            # 2D array with IK3D output (continuous)
            'p'      : }   # 2D array of floats. input data

    Returns
    -------
    out1 : rank-1 array('f') with bounds (na)
    out2 : rank-1 array('f') with bounds (na)
    out3 : rank-1 array('f') with bounds (na)


    If ``iout == 1`` (E-Type estimate) ``out1`` will be the E-Type
    estimate. ``out2`` and ``out3`` will be NaN

    If ``iout == 2`` (prob > cutoff) ``out1`` will be the probability
    above the cutoff, ``out2`` the mean above the cutoff and ``out3``
    the mean below the cutoff

    If ``iout == 3`` (p-quantile estimate) ``out1`` will be the
    Z value corresponding to the CDF == p. ``out2`` and ``out3`` will be
    NaN

    If ``iout == 4`` (conditional variance) ``out1`` will be the
    conditional variance. ``out2`` and ``out3`` will be NaN



    Example
    -------

    >>>     out1,out2,out3,error = postik(parameters)

    TODO
    ----
    Fix bug. The algorithm fails if there are nans in probabilities.

    """

    if parameters['vr'] is None:
        assert parameters['ltail']!=3 and parameters['middle']!=3 and parameters['utail']!=3
        parameters['vr'] = [0.,0.,0.]
        parameters['wt'] = [1.,1.,1.]

    # prepare the output

    out1,out2,out3,error = __gslib__postik.postik(**parameters)


    if error>0:
        raise NameError('Error = {}'.format(error))

    return out1,out2,out3


#-----------------------------------------------------------------------------------------------------------------
#
#    Plot fuctions (declustering supported)
#
#-----------------------------------------------------------------------------------------------------------------
def histgplt(parameters, fig = None, ax = None, title = 'Bin probability', label = 'Data', alpha= 0.7, grid=True, loc = 2, color='b', c = 0.001, edgecolor= 'black'):
    """Plot histogram. It uses declustering weight.


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::

        histgplt_parameters = {
            # histogram limits and definition.
            'hmin' : , # input float (Optional: set as minimum value in dataset). Minimun value in the histogram.
            'hmax' : , # input float (Optional: set as maximum value in dataset). Maximum value in the histogram.
            'ncl'  : , # input int (Optional: set as 10). Number of bins/clases.
            'iwt'  : , # input boolean (Optional: set True). Use weight variable?
            'ilog' : , # input boolean (Optional: set False). If true uses log scale, otherwise uses arithmetic
            'icum' : , # input boolean (Optional: set False). If true uses cumulative histogram, otherwise plots frequency histograms
            'va'   : , # input rank-1 array('d') with bounds (nd). Variable
            'wt'   : } # input rank-1 array('d') with bounds (nd) (Optional, set to array of ones). Declustering weight.



    Returns
    -------
    dict: a dictionary with some data parameters

    {'binval' : , # rank-1 array('d') with bounds (ncl). Class frequency.
     'nincls' : , # rank-1 array('d') with bounds (ncl). Count of values in the class (not using declustering weight).
     'cl'     : , # rank-1 array('d') with bounds (ncl). Class.
     'clwidth': , # rank-1 array('d') with bounds (ncl). Class width
     'xpt025' : , # float. Percentile 0.025
     'xlqt'   : , # float. Percentile 0.25
     'xmed'   : , # float. Mediam (Percentile 0.50)
     'xuqt'   : , # float. Percentile 0.75
     'xpt975' : , # float. Percentile 0.975
     'xmin'   : , # float. Minimum value.
     'xmax'   : , # float. Maximum value.
     'xcvr'   : , # float. Coeficient of variation.
     'xmen'   : , # float. Mean
     'xvar'   : , # float. Variance
     'xfrmx'  : , # float. Maximum Class Frequency
     'dcl'    : } # float. Class width (thmax-thmin)/real(ncl), only usefull if using arithmetic scale.

    fig: a matplotlib figure

    Example
    -------

    >>>

    """

    # set weight if not included in the input
    if 'wt' not in parameters:
        parameters['wt']= np.ones(parameters['va'].shape[0], dtype = parameters['va'].dtype)

    # set undefine parameters

    if 'ilog' not in parameters:
        parameters['ilog']= 1
    if parameters['ilog'] is None:
        parameters['ilog']= 1
    if parameters['ilog'] < 0:
        parameters['ilog']= 1

    if 'hmin' not in parameters:
        parameters['hmin']= parameters['va'].min()
    if parameters['hmin'] is None:
        parameters['hmin']= parameters['va'].min()
    if parameters['hmin']<=0 and parameters['ilog']==1:
        parameters['hmin']= c

    if 'hmax' not in parameters:
        parameters['hmax']= parameters['va'].max()
    if parameters['hmax'] is None:
        parameters['hmax']= parameters['va'].max()
    if parameters['hmax']<=0 and parameters['ilog']==1:
        raise NameError('All data in parameters["va"]<=0 but you are using logarithmic sale. Make your data positive and try again')

    if 'ncl' not in parameters:
        parameters['ncl']= 10
    if parameters['ncl'] is None:
        parameters['ncl']= 10
    if parameters['ncl'] <= 0:
        parameters['ncl']= 10

    if 'iwt' not in parameters:
        parameters['iwt']= 0
    if parameters['iwt'] is None:
        parameters['iwt']= 0

    if 'icum' not in parameters:
        parameters['icum']= 0
    if parameters['icum'] is None:
        parameters['icum']= 0



    # prepare the output

    binval,nincls,cl, clwidth,xpt025,xlqt,xmed,xuqt,xpt975, \
    xmin,xmax,xcvr,xmen,xvar,xfrmx,dcl,error = __plot.histplt(**parameters)

    out1 = {'binval' : binval, # rank-1 array('d') with bounds (ncl). Class frequency.
     'nincls' : nincls , # rank-1 array('d') with bounds (ncl). Count of values in the class (not using declustering weight).
     'cl'     : cl , # rank-1 array('d') with bounds (ncl). Class.
     'clwidth': clwidth, # rank-1 array('d') with bounds (ncl). Class width
     'xpt025' : xpt025, # float. Percentile 0.025
     'xlqt'   : xlqt, # float. Percentile 0.25
     'xmed'   : xmed, # float. Mediam (Percentile 0.50)
     'xuqt'   : xuqt, # float. Percentile 0.75
     'xpt975' : xpt975, # float. Percentile 0.975
     'xmin'   : xmin, # float. Minimum value.
     'xmax'   : xmax, # float. Maximum value.
     'xcvr'   : xcvr, # float. Coeficient of variation.
     'xmen'   : xmen, # float. Mean
     'xvar'   : xvar, # float. Variance
     'xfrmx'  : xfrmx, # float. Maximum Class Frequency
     'dcl'    : dcl} # float. Class width (thmax-thmin)/real(ncl), only usefull if using arithmetic scale.

    # create a figure or update one

    if ax is None and fig is None:
        fig, ax = plt.subplots(1,1)
        ax.set_title(title)
        if parameters['ilog']>0:
            ax.set_xscale('log')
    else:
        if ax is not None and fig is not None:
            pass
        else:
            raise NameError('Parameter error on ax or fig, both may be None or not None')


    ax.bar (cl-clwidth/2., binval, width=-clwidth, align='center', alpha=alpha, color=color, label = label, edgecolor= edgecolor)
    ax.grid(grid)
    ax.legend(loc=loc)

    return out1, ax, fig


def cdfplt(parameters, fig = None, ax = None, title = 'Bin probability', label = 'Data',  grid=True, loc = 2, color='b', xlog = True, ylog=True, style = '+'):
    """Plot cdf. It uses declustering weight.


    Parameters
    ----------
        parameters  :  dict
            dictionary with calculation parameters

    The dictionary with parameters may be as follows::

        cdfplt_parameters = {
            # CDF limits and definition.
            'iwt'  : , # input boolean (Optional: set True). Use weight variable?
            'va'   : , # input rank-1 array('d') with bounds (nd). Variable
            'wt'   : } # input rank-1 array('d') with bounds (nd) (Optional, set to array of ones). Declustering weight.



    Returns
    -------
    dict: a dictionary with some data parameters

    {'binval' : , # rank-1 array('d') with bounds (ncl). Class frequency.
     'cl'     : , # rank-1 array('d') with bounds (ncl). Class.
     'xpt025' : , # float. Percentile 0.025
     'xlqt'   : , # float. Percentile 0.25
     'xmed'   : , # float. Mediam (Percentile 0.50)
     'xuqt'   : , # float. Percentile 0.75
     'xpt975' : , # float. Percentile 0.975
     'xmin'   : , # float. Minimum value.
     'xmax'   : , # float. Maximum value.
     'xcvr'   : , # float. Coeficient of variation.
     'xmen'   : , # float. Mean
     'xvar'   : } # float. Class width (thmax-thmin)/real(ncl), only usefull if using arithmetic scale.

    fig: a matplotlib figure

    Example
    -------

    >>>

    """

    # set weight if not included in the input
    if 'wt' not in parameters:
        parameters['wt']= np.ones(parameters['va'].shape[0], dtype = parameters['va'].dtype)

    if 'iwt' not in parameters:
        parameters['iwt']= 0
    if parameters['iwt'] is None:
        parameters['iwt']= 0


    # prepare the output

    binval,cl,xpt025,xlqt,xmed,xuqt,xpt975,xmin,xmax, \
    xcvr,xmen,xvar,error = __plot.probplt(**parameters)

    out1 = {'binval' : binval, # rank-1 array('d') with bounds (ncl). Class frequency.
     'cl'     : cl , # rank-1 array('d') with bounds (ncl). Class.
     'xpt025' : xpt025, # float. Percentile 0.025
     'xlqt'   : xlqt, # float. Percentile 0.25
     'xmed'   : xmed, # float. Mediam (Percentile 0.50)
     'xuqt'   : xuqt, # float. Percentile 0.75
     'xpt975' : xpt975, # float. Percentile 0.975
     'xmin'   : xmin, # float. Minimum value.
     'xmax'   : xmax, # float. Maximum value.
     'xcvr'   : xcvr, # float. Coeficient of variation.
     'xmen'   : xmen, # float. Mean
     'xvar'   : xvar} # float. Variance

    # create a figure or update one

    if ax is None and fig is None:
        fig, ax = plt.subplots(1,1)
        ax.set_title(title)
        if xlog:
            ax.set_xscale('log')
        if ylog:
            ax.set_yscale('log')
    else:
        if ax is not None and fig is not None:
            pass
        else:
            raise NameError('Parameter error on ax or fig, both may be None or not None')


    ax.plot (cl, binval, color=color, label = label)
    ax.grid(grid)
    ax.legend(loc=loc)

    return out1, ax, fig
