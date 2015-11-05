'''
PyGSLIB Blockmodel, Module to handle blockmodel data.  

Copyright (C) 2015 Adrian Martinez Vargas 

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
any later version.
   
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
'''


cimport numpy as np
import numpy as np
import pandas as pd
import warnings

#-------------------------------------------------------------------
#  General functions
#-------------------------------------------------------------------
cpdef x2ix(np.ndarray [double, ndim=1] x,
           double xorg,
           double dx):
    """
    x2ix(x, xorg,dx)
    
    Calculates the block index (ix, iy or iz) 
    
    Return array of ix indices representing the raw, column or level 
    in which a point with coordinate x is located. To get the three
    indexes ix, iy, iz you may run this functions three times, changing
    coordinates to x, y and z, respectively.    
    
    You can use this function to find the block which contain a points 
    with arbitrary coordinates.  
    
    Parameters
    ----------
    x    : 1D array of floats 
           arbitrary coordinates x (or y, or z)
    xorg : float
           origin of coordinate (lower left corner of the block) 
    dx   : float 
           size of the block
    
    Returns
    -------
    ix : 1D array of integers
        Index of the blocks where points with coordinets x are located
    
    See Also
    --------
    ind2ijk, ijk2ind
    
    Notes
    -----
    The index start at zero (first block)
    If there is a point i with coordinate x[i] < xorg  and error will be 
    raised
    
    
    Examples
    --------
    >>> a=np.array([2.1,7.8,15.9,20,30,31])
    >>> x2ix(a,xorg=2., dx=1)
    array([ 0,  5, 13, 18, 28, 29])

    """
    
    assert  xorg < x.min(), '\nError:\n x.min < xorg, values out of grid. Redefine xorg<= %f ' %  x.min()
    
    # output 
    cdef np.ndarray [long, ndim=1, cast=True] ix= np.empty([x.shape[0]], dtype= int) 
    
    ix[:]= (x-xorg)/dx # numpy do casting authomatically
                       # we start from 0 
                       # 

    return ix


cpdef ind2ijk(np.ndarray [long, ndim=1] ix,
              np.ndarray [long, ndim=1] iy,
              np.ndarray [long, ndim=1] iz,
              unsigned int nx,
              unsigned int ny, 
              unsigned int nz):
    """
    ind2ijk(ix,iy,iz, nx,ny,nz)
    
    Calculates the IJK block index
    
    The IJK block index is an unique identifier of each block position.
    This is equivalent to the position of a block in a gslib grid file. 
    All the points within a block will have the same IJK number. 
    
    
    Parameters
    ----------
    ix, iy, iz : 1D array of integers 
           arbitrary raw, level and column indices
    nx,ny,nz   : integers
           number of blocks per row, column, level
    
    Returns
    -------
    ijk : 1D array of integers
        Unique identifier of block location
    
    See Also
    --------
    x2ix, ijk2ind
    
    Notes
    -----
    The index ijk start at zero (first block) and ends at nx*ny*nz
    
    
    Examples
    --------
    >>> # a 2D grid with 2x2x1 cells 
    >>> ix=np.array([0,1,0,1])
    >>> iy=np.array([1,1,0,0])
    >>> iz=np.array([0,0,0,0])
    >>> ind2ijk(ix,iy,iz,2,2,1)
    array([2, 3, 0, 1])

    """
    
    assert ix.shape[0]==iy.shape[0]==iz.shape[0], 'Error: wrong shape ix, iy and iz may have same shape'
    assert nx>0 and ny>0 and nz>0, 'Error: nx, ny and nz may be >=1 '
    assert ix.min()>=0, 'Error: Negative index in ix'
    assert iy.min()>=0, 'Error: Negative index in iy'
    assert iz.min()>=0, 'Error: Negative index in iz'
    assert ix.max()<nx or iy.max()<ny or iz.max()<nz , 'Error: Index out or range, ex. ix.max>=nx. Review nx,ny,nz!'
    
    
    cdef np.ndarray [long, ndim=1, cast=True] ijk= np.empty([ix.shape[0]], dtype= int) 
    
    
    # output based on gslib page 21, loc formula but index starting at zero 
    ijk[:]= iz*nx*ny+ iy*nx + ix
                           
    return ijk

cpdef ijk2ind(np.ndarray [long, ndim=1] ijk,
              unsigned int nx,
              unsigned int ny, 
              unsigned int nz):
    """
    ijk2ind(ijk,nx,ny,nz)
    
    Calculates the raw, column, level indices ``ix, iy, iz`` from IJK.
    
    The IJK block index is an unique identifier of each block position.
    This is equivalent to the position of a block in a gslib grid file. 
    All the points within a block will have the same IJK number. 
    
    From IJK you can calculate ix,iy and iz as:
    
    iz =  ijk / nx*ny
    iy = (ijk-iz*nx*ny)/nx
    ix =  ijk-iz*nx*ny - iy*nx
    
    
    Parameters
    ----------
    ijk        : 1D array of integers 
           arbitrary raw, level and column indices
    nx,ny,nz   : integers
           number of blocks per row, column, level
    
    Returns
    -------
    ix,iy,iz : 1D array of integers
        The raw, column and level indices
    
    See Also
    --------
    x2ix, ijk2ind
    
    Notes
    -----
    The indices ix, iy and iz start at zero
    
    
    Examples
    --------
    >>> ijk= np.array([0, 1, 2, 3, 4, 5, 6, 7])
    >>> ix,iy,iz = ijk2ind(ijk,2,2,2)
    >>> print ix
    >>> print iy
    >>> print iz
    [0 1 0 1 0 1 0 1]
    [0 0 1 1 0 0 1 1]
    [0 0 0 0 1 1 1 1]

    """ 
    cdef float fnx, fny,fnz
    
    assert nx>0 and ny>0 and nz>0, 'Error: nx, ny and nz may be >=1 '
    assert nx*ny*nz>=ijk.shape[0], 'Error: out of range, nx*ny*nz<ijk.shape[0] '
    
    if  ijk.min()<0:
        warnings.warn('\nWarning:\n Negative indices. Review ijk!')
    if  ijk.max()>nx*ny*nz:
        warnings.warn('\nWarning:\n Index out or range, ex. ijk.max()>nx*ny*nz. Review ijk!')
    
    cdef np.ndarray [long, ndim=1, cast=True] ix= np.empty([ijk.shape[0]], dtype= int) 
    cdef np.ndarray [long, ndim=1, cast=True] iy= np.empty([ijk.shape[0]], dtype= int) 
    cdef np.ndarray [long, ndim=1, cast=True] iz= np.empty([ijk.shape[0]], dtype= int) 
    
    # to ensure float division
    fnx=nx
    fny=ny
    fnz=nz
    
    # output based on gslib page 21, loc formula but index starting at zero 
    # make sure to put float in parentesis ex. (fnx*fny)
    # <int>(ijk[i]/(fnx*fny)) != <int>(ijk[i]/fnx*fny)
    for i in range(ijk.shape[0]):
        iz[i] = <int>(ijk[i]/(fnx*fny))
        iy[i] = <int>((ijk[i]-iz[i]*(fnx*fny))/fnx)
        ix[i] = <int>(ijk[i]-iz[i]*(fnx*fny) - iy[i]*fnx)
        
    return ix,iy,iz

#-------------------------------------------------------------------
#  Drillhole class
#-------------------------------------------------------------------
cdef class Blockmodel:
    """
    Blockmodel(nx,ny,nz,xorg,yorg,zorg,dx,dy,dz)
    
    Blockmodel working database object with functions to handle 
    block models.
    
    Parameters
    ----------
    nx,ny,nz : integer
        Number of rows columns and levels  
    xorg,yorg,zorg : float
        Coordinates of the lower left corner (not centroid) 
    dx,dy,dz : float
        Size of the parent block
    
    Attributes
    ----------

    
    Notes
    -----
 
    
    """ 
    cdef readonly int nx,ny,nz
    cdef readonly double xorg,yorg,zorg,dx,dy,dz
        
    
    def __cinit__(self,nx,ny,nz,xorg,yorg,zorg,dx,dy,dz):
        assert nx>0 and ny>0 and nz>0 
        assert xorg>0 and yorg>0 and zorg
        assert dx>0 and dy>0 and dz>0
        
        self.nx=nx
        self.ny=ny
        self.nz=nz
        
        self.xorg=xorg
        self.yorg=yorg
        self.zorg=zorg

        self.dx=dx
        self.dy=dy
        self.dz=dz
            
