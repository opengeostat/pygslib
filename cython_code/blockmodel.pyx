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

import vtk
import vtk.util.numpy_support as vtknumpy
cimport numpy as np
import numpy as np
import pandas as pd
import warnings
import pygslib

#-------------------------------------------------------------------
#  General functions
#-------------------------------------------------------------------
cpdef x2ix(np.ndarray [double, ndim=1] x,
           double xorg,
           double dx):
    """x2ix(np.ndarray [double, ndim=1] x, double xorg, double dx)

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
        Index of the blocks where points with coordinates x are located

    Example
    -------
    >>>
    >>> a=np.array([2.1,7.8,15.9,20,30,31])
    >>> x2ix(a,xorg=2., dx=1)
    array([ 0,  5, 13, 18, 28, 29])
    >>>

    See Also
    --------
    ind2ijk, ijk2ind

    Note
    ----
    The index starts at zero (first block)
    If there is a point i with coordinate x[i] < xorg and error will be
    raised

    """

    assert  xorg < x.min(), '\nError:\n x.min ={0} < xorg = {1}, values out of grid. Redefine xorg<= {0}'.format(x.min(),xorg)


    # output
    cdef np.ndarray [long, ndim=1, cast=True] ix= np.empty([x.shape[0]], dtype= int)

    ix[:]= (x-xorg)/dx # numpy do casting authomatically
                       # we start from 0
                       #

    return ix

cpdef ix2x(np.ndarray [long, ndim=1] ix,
           double xorg,
           double dx):
    """ix2x(np.ndarray [long, ndim=1] ix, double xorg, double dx)

    Calculates the block coordinate (x, y or z)

    Return array of coordinates x calculated from indices ix.


    Parameters
    ----------
    ix   : 1D array of positive integers
           arbitrary index ix (or iy, or iz)
    xorg : float
           origin of coordinate (lower left corner of the block)
    dx   : float
           size of the block

    Returns
    -------
    x  : 1D array of floats
        Coordinates x corresponding to index ix

    See Also
    --------
    x2ix

    Notes
    -----
    The index starts at zero (first block)
    If there is a point i with coordinate x[i] < xorg and error will be
    raised

    """

    assert  ix.min()>=0 , '\nError: ix<0'

    # output
    cdef np.ndarray [double, ndim=1, cast=True] x= np.empty([ix.shape[0]], dtype= float)

    x = ix[:]*dx+xorg+dx/2


    return x


cpdef ind2ijk(np.ndarray [long, ndim=1] ix,
              np.ndarray [long, ndim=1] iy,
              np.ndarray [long, ndim=1] iz,
              unsigned int nx,
              unsigned int ny,
              unsigned int nz):
    """ind2ijk(np.ndarray [long, ndim=1] ix, np.ndarray [long, ndim=1] iy, np.ndarray [long, ndim=1] iz, unsigned int nx, unsigned int ny, unsigned int nz)

    Calculates the IJK block index

    The IJK block index is a unique identifier of each block position.
    This is equivalent to the position of a block in a `gslib` grid file.
    All the points within a block will have the same IJK number.

    IJK is calculated as follow:

    ``ijk = iz*nx*ny+ iy*nx + ix ``


    Parameters
    ----------
    ix,iy,iz : 1D array of integers
           arbitrary raw, level and column indices
    nx,ny,nz   : integers
           number of blocks per row, column, level

    Returns
    -------
    ijk : 1D array of integers
        Unique identifier of block location

    Example
    --------
    >>>
    >>> # a 2D grid with 2x2x1 cells
    >>> ix=np.array([0,1,0,1])
    >>> iy=np.array([1,1,0,0])
    >>> iz=np.array([0,0,0,0])
    >>> ind2ijk(ix,iy,iz,2,2,1)
    array([2, 3, 0, 1])
    >>>

    See Also
    --------
    x2ix, ijk2ind

    Note
    -----
    The index IJK start at zero (first block) and ends at nx*ny*nz

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
    """ijk2ind(np.ndarray [long, ndim=1] ijk, unsigned int nx, unsigned int ny, unsigned int nz)

    Calculates the raw, column, level indices ``ix, iy, iz`` from IJK.

    The IJK block index is a unique identifier of each block position.
    This is equivalent to the position of a block in a `gslib` grid file.
    All the points within a block will have the same IJK number.

    From IJK you can calculate ix, iy and iz as::

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

    Example
    -------
    >>>
    >>> ijk= np.array([0, 1, 2, 3, 4, 5, 6, 7])
    >>> ix,iy,iz = ijk2ind(ijk,2,2,2)
    >>> print ix
    >>> print iy
    >>> print iz
    [0 1 0 1 0 1 0 1]
    [0 0 1 1 0 0 1 1]
    [0 0 0 0 1 1 1 1]
    >>>

    See Also
    --------
    x2ix, ijk2ind

    Note
    ----
    The indices ix, iy and iz start at zero

    """
    cdef double fnx, fny,fnz

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
    """ Blockmodel(nx,ny,nz,xorg,yorg,zorg,dx,dy,dz)
    Blockmodel object.



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
    nx, ny, nz : int
        number of rows, columns, levels

    xorg, yorg, zorg : float
        coordinates of the left lower corner of the first block

    dx, dy, dz : float
        sizes of the blocks

    bmtable : Pandas DataFrame
        Table with blocks

    gridtable : Pandas Dataframe


    """
    cdef readonly int nx,ny,nz  # nrows,cols,levels
    cdef readonly int ndx,ndy,ndz  # ndiscretization
    #cdef readonly bint randisc     # is the discretization random?
    #cdef readonly bint subcell     # has subcells?
    #cdef readonly bint percent     # has percent?  Can be both percent+subcell
    cdef readonly double xorg,yorg,zorg,dx,dy,dz # origin of coordinate + block size
    cdef readonly object bmtable    # the actual block model
    cdef readonly object grid2DXY  # the grid table
    cdef readonly object grid2DXZ  # the grid table
    cdef readonly object grid2DYZ  # the grid table

    def __cinit__(self,nx,ny,nz,xorg,yorg,zorg,dx,dy,dz):
        assert nx>0 and ny>0 and nz>0
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

        # self.bmtable=None

    cpdef set_block_size(self, float dx,float dy, float dz):
        """set_block_size(float dx,float dy, float dz)

        Set block sizes

        Parameters
        ----------
        dx,dy,dz : float
            sizes of the blocks

        Examples
        --------
        >>>
        >>> myblockmodel.set_block_size(dx=10,dy=10,dz=5)
        >>>

        """



        assert dx>0 and dy>0 and dz>0
        self.dx=dx
        self.dy=dy
        self.dz=dz


    cpdef set_origin(self, float xorg,float yorg, float zorg):
        """set_origin(float xorg,float yorg, float zorg)

        Set the block model origin of coordinate. This is the lower
        left corner of the lower left block (not-the centroid).

        Parameters
        ----------
        xorg,yorg,zorg : float
            coordinates of the left lower corner of the first block

        Examples
        --------
        >>>
        >>> myblockmodel.set_origin(xorg=-5,yorg=-5, zorg=0)
        >>>

        """
        assert xorg>0 and yorg>0 and zorg>0
        self.xorg=xorg
        self.yorg=yorg
        self.zorg=zorg

    cpdef set_rcl(self, int nx,int ny, int nz):
        """set_rcl(int nx,int ny, int nz)

        Set number of blocks at row, column, level

        Two conditions apply

        - nx,ny and nz may be greater than zero
        - nx*ny*nz may be equal to the actual number of blocks

        Examples
        --------
        >>>
        >>> myblockmodel.set_rcl(nx=10,ny=10,nz=10)
        >>> myblockmodel.set_rcl(nx=10,ny=100,nz=1)
        >>>

        Note
        ----
        This function basically reorders the blocks along rows, columns
        and levels without changing the total number of blocks.


        """
        assert nx>0 and ny>0 and nz>0, 'nx<=0 or ny<=0 or nz<=0'
        assert nx*ny*nz==self.nx*self.ny*self.nz, 'nx*ny*nz may be equal to {}'.format(self.nx*self.ny*self.nz)

        self.nx=nx
        self.ny=ny
        self.nz=nz


    cpdef set_blocks(self, object bmtable):
        """set_blocks(object bmtable)

        Assigns an external block model stored in a Pandas DataFrame
        table.

        One or many conditions apply

         - the block model has IJK field
         - or/and has IX, IY and IZ fields
         - or/and has XC, YC and ZC fields


        Parameters
        ----------
        bmtable : Pandas DataFrame object
            block model with special coordinate or index fields


        Examples
        --------
        >>>
        >>> bmtable=pd.DataFrame({'IJK':np.array([1,2,3,4], dtype=np.int32)})
        >>> myblockmodel.set_blocks(bmtable)
        >>>

        Note
        ----
        Make sure IJK, IX,IY and IZ have dtype int32. XC,YC and ZC may
        have dtype float32

        One of the above conditions is required. In case that two or
        more of the special fields mentioned above are defined the
        validity of the fields is not verified, for example the
        valid IJK representation of fields IX, IY and IZ is not tested.

        """
        cdef bint has_ijk, has_ixyz, has_cxyz

        assert isinstance(bmtable, pd.DataFrame), "bmtable is not a pandas DataFrame"

        # test the type of block model
        has_ijk = 'IJK' in bmtable.columns
        has_ixyz = set(('IX', 'IY', 'IZ')).issubset(bmtable.columns)
        has_cxyz = set(('XC', 'YC', 'ZC')).issubset(bmtable.columns)

        assert has_ijk or has_ixyz or has_cxyz , "Fields IJK or IX,IY,IZ or XC,YC,ZC defined in the input bmtable"

        self.bmtable=bmtable

    cpdef delete_blocks(self):
        """delete_blocks()

        Deletes the block model table

        Examples
        --------
        >>>
        >>> myblockmodel.delete_blocks()
        >>> print myblockmodel.bmtable
        None
        >>>

        Note
        ----
        This functions makes Blockmodel.bmtable=None. The data will
        be preserved in any external instance of bmtable.

        """

        self.bmtable=None

    cpdef calc_ixyz_fromxyz(self, bint overwrite=False):
        """calc_ixyz_fromxyz(bint overwrite=False)

        Calculates the IX, IY, IZ fields from XC, YC, ZC coordinates
        If IX, IY, IZ exist and overwrite=True the existing values
        will be overwritten.


        Parameters
        ----------
        overwrite : Boolean, default False

        Example
        -------
        >>>
        >>> myblockmodel.calc_ixyz_fromxyz()
        >>> myblockmodel.calc_ixyz_fromxyz(overwrite=True)
        >>>

        """

        cdef bint has_ixyz

        assert isinstance(self.bmtable, pd.DataFrame), 'Error: No bmtable loaded or created yet'
        assert set(('XC', 'YC', 'ZC')).issubset(self.bmtable.columns), 'Error: No XC,YC,ZC coordinates in bmtable'
        if overwrite==False:
            has_ixyz = set(('IX', 'IY', 'IZ')).issubset(self.bmtable.columns)
            assert has_ixyz==False, 'IX,IY,IZ already exist in bmtable, set overwrite=True to overwrite'

        self.bmtable['IX']= x2ix (self.bmtable['XC'].values.astype('float'), self.xorg, self.dx)
        self.bmtable['IY']= x2ix (self.bmtable['YC'].values.astype('float'), self.yorg, self.dy)
        self.bmtable['IZ']= x2ix (self.bmtable['ZC'].values.astype('float'), self.zorg, self.dz)

    cpdef calc_xyz_fromixyz(self, bint overwrite=False):
        """calc_xyz_fromixyz(bint overwrite=False)

        Calculates the XC, YC, ZC coordinates from IX, IY, IZ indices
        If XC, YC, ZC exist and overwrite=True the existing values
        will be overwritten.


        Parameters
        ----------
        overwrite : Boolean, default False

        Examples
        --------
        >>>
        >>> myblockmodel.calc_xyz_fromixyz()
        >>> myblockmodel.calc_xyz_fromixyz(overwrite=True)
        >>>

        """
        cdef bint has_xyz

        assert isinstance(self.bmtable, pd.DataFrame), 'Error: No bmtable loaded or created yet'
        assert set(('IX', 'IY', 'IZ')).issubset(self.bmtable.columns), 'Error: No IX,IY,IZ coordinates in bmtable'
        if overwrite==False:
            has_ixyz = set(('XC', 'YC', 'ZC')).issubset(self.bmtable.columns)
            assert has_ixyz==False, 'XC,YC,ZC already exist in bmtable, set overwrite=True to overwrite'

        self.bmtable['XC']= ix2x (self.bmtable['IX'].values.astype('int'), self.xorg, self.dx)
        self.bmtable['YC']= ix2x (self.bmtable['IY'].values.astype('int'), self.yorg, self.dy)
        self.bmtable['ZC']= ix2x (self.bmtable['IZ'].values.astype('int'), self.zorg, self.dz)

    cpdef calc_ijk(self, bint overwrite=False):
        """calc_ijk(bint overwrite=False)

        Calculates the IJK field from IX, IY, IZ indices
        If IJK exist and overwrite=True the existing values
        will be overwritten.


        Parameters
        ----------
        overwrite : Boolean, default False
                If True overwrites any existing IJK field

        Example
        -------
        >>>
        >>> myblockmodel.calc_ijk()
        >>> myblockmodel.calc_ijk(overwrite=True)
        >>>
        """

        assert isinstance(self.bmtable, pd.DataFrame), 'Error: No bmtable loaded or created yet'
        assert set(('IX', 'IY', 'IZ')).issubset(self.bmtable.columns),  'Error: No IX,IY,IZ indices in bmtable'
        if overwrite==False:
            assert 'IJK' not in self.bmtable.columns, 'IJK already exist in bmtable, set overwrite=True to overwrite'

        self.bmtable['IJK'] =  ind2ijk(self.bmtable['IX'].values,
                                       self.bmtable['IY'].values,
                                       self.bmtable['IZ'].values,
                                       self.nx,
                                       self.ny,
                                       self.nz)

    cpdef calc_ixyz_fromijk(self, bint overwrite=False):
        """calc_ixyz_fromijk(bint overwrite=False)

        Calculates the IX, IY, IZ fields from IJK index
        If IX, IY, IZ exist and overwrite=True the existing values
        will be overwritten.


        Parameters
        ----------
        overwrite : Boolean, default False

        Example
        -------
        >>>
        >>> myblockmodel.calc_ixyz_fromijk()
        >>> myblockmodel.calc_ixyz_fromijk(overwrite=True)
        >>>

        """

        cdef bint has_ixyz

        assert isinstance(self.bmtable, pd.DataFrame), 'Error: No bmtable loaded or created yet'
        assert 'IJK' in self.bmtable.columns, 'Error: No IJK index in bmtable'
        if overwrite==False:
            has_ixyz = set(('IX', 'IY', 'IZ')).issubset(self.bmtable.columns)
            assert has_ixyz==False, 'IX,IY,IZ already exist in bmtable, set overwrite=True to overwrite'

        self.bmtable['IX'],self.bmtable['IY'],self.bmtable['IZ'] =  ijk2ind(self.bmtable['IJK'].values,self.nx,self.ny,self.nz)


    cpdef reblock(self, int nx, int ny, int nz, double xorg, double yorg, double zorg, double dx, double dy, double dz):
        """reblock(int nx, int ny, int nz, double xorg, double yorg, double zorg, double dx, double dy, double dz)

        Reblocks a block model to new block model with arbitrary definition

        The model is reblocked using as index the IJK corresponding
        to the input model definition, which is calculated from current
        XC, YC, ZC coordinates. The funtion computes mean and counts of
        all points with same IJK.

        Note that this is like a point to block conversion and subblock
        overlapping are not taken into account.

        Parameters
        ----------
        nx,ny,nz : integer
            Number of rows columns and levels in the new model
        xorg,yorg,zorg : float
            Coordinates of the lower left corner of the new model
        dx,dy,dz: floats
            New block dimensions. Note that new dimentions may be
            greater or equal than old dimensions in all directions.


        Returns
        -------
        block: pygslib.blockmodel.Blockmodel instance
            model reblocked

        err: list
            list of variables excluded (ej no numeric)

        Warning
        -------
        The existing model will be overwritten

        """

        assert dx>=self.dx  and dy>=self.dy and dz>=self.dz, 'Error: dx<self.dx, only reblocking to larger blocks supported'

        assert isinstance(self.bmtable, pd.DataFrame), 'Error: No bmtable loaded or created yet'
        assert set(('XC', 'YC', 'ZC')).issubset(self.bmtable.columns), 'Error: No XC,YC,ZC in bmtable'

        tmpmod= pygslib.blockmodel.Blockmodel(nx, ny, nz,
                                              xorg,yorg,zorg,
                                              dx,dy,dz)


        # add tables as a deep copy
        tmpmod.set_blocks(self.bmtable.copy(deep =True))

        # recalculate ix iy and iz
        tmpmod.calc_ixyz_fromxyz(overwrite=True)

        #recalculate IJK
        tmpmod.calc_ijk(overwrite=True)

        # this is the count of blocks
        tmpmod.bmtable['NBLK'] = 1


        npoints = int(dx/self.dx*dy/self.dy*dz/self.dz)

        # now create a group by ijk
        np = tmpmod.bmtable[['NBLK','IJK']].groupby('IJK').sum()
        np['RD'] = np['NBLK']/float(npoints)

        err = []
        for i in tmpmod.bmtable.columns:
            if i in ['IJK','RD','NBLK','IX','IY','IZ','XC','YC','ZC']:
                continue
            try:
                np[i] = tmpmod.bmtable[[i,'IJK']].groupby('IJK').mean()
                np[i+'_sum'] = tmpmod.bmtable[[i,'IJK']].groupby('IJK').sum()
                np[i+'_count'] = tmpmod.bmtable[[i,'IJK']].groupby('IJK').count()
            except:
                err.append(i)


        if len(err)>0:
            warnings.warn('\nWarning:\n Non numeric columns were ignored!')


        tmpmod.set_blocks(np.reset_index())

        # set index as integre 32 bits
        tmpmod.bmtable['IJK'] = tmpmod.bmtable['IJK'].values.astype(int)

        # recalculate ix iy and iz
        tmpmod.calc_ixyz_fromijk(overwrite=True)

        #recalculate XC,YC,ZC
        tmpmod.calc_xyz_fromixyz(overwrite=True)



        return tmpmod, err


    cpdef create_IJK(self, bint overwrite=False):

        """create_IJK(bint overwrite=False)

        Creates a new block model table with IJK indices.

        Examples
        --------
        >>>
        >>> create_IJK(overwrite=True)
        >>>

        Note
        ----
        A new set of blocks will be created if there is not block
        defined. If there are blocks in the model and overwrite==True
        the blocks will be removed first.

        """

        if overwrite==False:
             assert self.bmtable is None, 'Error: bmtable already exist, set overwrite=True to overwrite'

        self.bmtable=pd.DataFrame({'IJK': np.arange(self.nx*self.ny*self.nz, dtype=int)})


    cpdef create_3Dgrid(self, bint overwrite=False):
        """create_3Dgrid(bint overwrite=False)

        Creates a full block 3D grid/block model


        Parameters
        ----------
        overwrite : boolean
               overwrite flag, if true the entire grid will be
               overwritten

        Example
        -------
        >>>
        >>> mymodel.create_3Dgrid(overwrite=False)
        >>>


        """

        #check the model do not exist
        if overwrite==False:
             assert self.bmtable is None, 'Error: bmtable already exist, set overwrite=True to overwrite'


        # reusing existing functions
        self.create_IJK(overwrite)
        self.calc_ixyz_fromijk()
        self.calc_xyz_fromixyz()



    cpdef blockinsurface(self,
                    object surface,
                    str field,
                    double azm=0,
                    double dip =90,
                    int test=1,
                    bint overwrite=False):

        """blockinsurface(object surface, str field, double azm=0, double dip =90, int test=1, bint overwrite=False)

        Creates blocks given a VTK surface (polydata) depending on
        test criteria:

        Parameters
        ----------
        surface : VTK polydata
               this may work with any 3D object..., no necessarily triangles
        field : str
               Name of the new field with selection results
        azm, dip: float, default 0, 90
               rotation defining the direction we will use to test the points
               azm 0 will point north and dip positive means downward direction
               (like surface drillholes)
        test    : integer, default 1
               ``test == 1`` test inside closed surface. Here we use
                 vtkOBBTree::InsideOrOutside. Closed surface are required
               ``test == 2`` test 'above' surface
               ``test == 3`` test 'below' surface
               ``test == 4`` test 'inside' surface (the surface can be open)
         overwrite : boolean
               overwrite flag, if true and field exist in bmtable the
               values will be overwritten

        Example
        -------
        >>>
        >>> fillblocks(surface, azm, dip, test, overwrite=False)
        >>>

        Note
        ----
        This function calls vtktools.pointquering for all the points
        existing in the block model. The function is not optimized
        for large model.

        A new set of blocks will be created if there is not block
        defined. If there are blocks in the model and overwrite==True
        the blocks will be removed first.

        """

        assert isinstance(surface, vtk.vtkPolyData), 'Error: No bmtable loaded or created yet'
        assert isinstance(self.bmtable, pd.DataFrame), 'Error: No bmtable loaded or created yet'
        assert set(('XC', 'YC', 'ZC')).issubset(self.bmtable.columns), 'Error: No XC,YC,ZC coordinates in bmtable'

        if overwrite==False:
             assert field not in self.bmtable.columns, 'Error: field {} already exist in bmtable, set overwrite=True to overwrite'.format(field)

        self.bmtable[field], p1=pygslib.vtktools.pointquering(surface, azm, dip, self.bmtable['XC'].values, self.bmtable['YC'].values, self.bmtable['ZC'].values, test)

    cpdef blocks2vtkImageData(self, str path=''):
        """blocks2vtkImageData(str path)

        Returns blocks of a full grid as a vtkImageData and optionaly save
        it to a *.vti file.

        Parameters
        ----------
        path : string (default '')
               file name and path. The file extension (*.vtr) will be
               added automatically. The file will not be saved if path is ''

        Returns
        --------
        vtkImageData with full grid model

        Examples
        --------
        >>>
        >>> blocks2vtkImageData('myfile')
        >>>

        Note
        ----
        This will only work for full grid, in other words, if all the
        nx*ny*nz are defined.

        All the fields defined in the block model will be exported

        """

        cdef np.ndarray [double, ndim=1] x,y,z

        assert  self.bmtable.shape[0]==self.nx*self.ny*self.nz, 'Error: this work only with full grid, ex. bmtable.shape[0]==nx*ny*nz'

        data = {}

        for i in self.bmtable.columns:
            data[i]=self.bmtable[i].values


        vtkgrid = pygslib.vtktools.grid2vtkImageData(nx=self.nx,
                                                     ny=self.ny,
                                                     nz=self.nz,
                                                     xorg=self.xorg,
                                                     yorg=self.yorg,
                                                     zorg=self.zorg,
                                                     dx=self.dx,
                                                     dy=self.dy,
                                                     dz=self.dz,
                                                     cell_data=data)

        if path!='':
          pygslib.vtktools.SaveImageData(vtkgrid, path)

        return vtkgrid


    cpdef blocks2vtkUnstructuredGrid(self, str path, str varname= None):
        """blocks2vtkUnstructuredGrid(str path, str varname = None)

        Exports blocks of a partial grid to a vtkUnestructuredGrid cells.

        Parameters
        ----------
        path : string
               file name and path, without extension. The file extension
               (*.vtr) will be added automatically.

        varname : string
                If None export all columns, otherwise export only
                the column varname


        Examples
        --------
        >>>
        >>> blocks2vtkUnstructuredGrid('myfile', 'AU')
        >>>

        Note
        ----
        Require XC, YC and ZC.

        Only numeric variables can be exported
        TODO: export all

        This call pygslib.vtktools.partialgrid2vtkfile

        """


        cdef np.ndarray [double, ndim=1] x,y,z

        x= self.bmtable['XC'].values
        y= self.bmtable['YC'].values
        z= self.bmtable['ZC'].values
        var = []
        prop_name = []

        if varname is None:
            # prepare colums as array of numpy data
            for i in self.bmtable.columns:
                var.append(self.bmtable[i].values)
                prop_name.append(i)
        else:
            var = [self.bmtable[varname]]
            prop_name = [varname]

        if ('DX' in self.bmtable.columns):
          pygslib.vtktools.partialgrid2vtkfile(path=path, x=x, y=y, z=z,
                               DX=self.dx, DY=self.dy, DZ=self.dz,
                               dx=self.bmtable['DX'].values,
                               dy=self.bmtable['DY'].values,
                               dz=self.bmtable['DZ'].values,
                               var=var, varname=prop_name)
        else:
            pygslib.vtktools.partialgrid2vtkfile(path=path, x=x, y=y, z=z,
                                 DX=self.dx, DY=self.dy, DZ=self.dz,
                                 var=var, varname=prop_name)

    cpdef blocks2vtkUnstructuredGrid_p(self, str path, str varname= None):
        """blocks2vtkUnstructuredGrid(str path, str varname = None)

        Exports block centroids of a partial grid to a vtkUnistructuredGrid points.

        Parameters
        ----------
        path : string
               file name and path, without extension. The file extension
               (*.vtr) will be added automatically.

        varname : string
                If None export all columns, otherwise export only
                the column varname


        Examples
        --------
        >>>
        >>> blocks2vtkUnstructuredGrid('myfile', 'AU')
        >>>

        Note
        ----
        Require XC, YC and ZC.

        Only numeric variables can be exported
        TODO: export all

        This call pygslib.vtktools.partialgrid2vtkfile

        """


        cdef np.ndarray [double, ndim=1] x,y,z

        x= self.bmtable['XC'].values
        y= self.bmtable['YC'].values
        z= self.bmtable['ZC'].values
        var = []
        prop_name = []

        if varname is None:
            # prepare colums as array of numpy data
            for i in self.bmtable.columns:
                var.append(self.bmtable[i].values)
                prop_name.append(i)
        else:
            var = [self.bmtable[varname]]
            prop_name = [varname]


        pygslib.vtktools.partialgrid2vtkfile_p(path, x, y, z,
                             self.dx, self.dy, self.dz, var, prop_name)

    cpdef block2point(self, np.ndarray [double, ndim=1] x,
                            np.ndarray [double, ndim=1] y,
                            np.ndarray [double, ndim=1] z,
                            str prop_name):
        """block2point(np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, str prop_name)

        Assigns a block property to the points inside a block


        Parameters
        ----------
        x,y,z : float 1D numpy arrays
            coordinates of the points. Array shapes may be equal
        prop_name : str
            name of the property at block model


        Example
        -------
        >>>
        >>> Au_at_point= point2block( x, y, z, prop_name= 'Au')
        >>>


        Note
        ----
        Points not laying in an existing block will be returned with
        value numpy.nan
        """

        cdef long i, n
        cdef np.ndarray [double, ndim=1] var

        # note that prop can be an array of array with any type
        assert x.shape[0] == y.shape[0] == z.shape[0]


        var = np.empty(x.shape[0], dtype = float)
        var[:]= np.nan


        # make sure there is ijk field in model

        assert 'IJK' in self.bmtable.columns, 'Error: no IJK field in block model'


        # get index ix, iy, iz
        ix = x2ix (x, self.xorg, self.dx)
        iy = x2ix (y, self.yorg, self.dy)
        iz = x2ix (z, self.zorg, self.dz)

        # calculate ijk
        ijk = ind2ijk(ix,iy,iz, self.nx,self.ny,self.nz)

        # set ijk as key for indexing
        n = x.shape[0]
        self.bmtable.set_index('IJK', inplace= True)


        for i in range(n):
            if ijk[i] in self.bmtable.index:
                var[i]=self.bmtable.at[ijk[i],prop_name]


        self.bmtable.reset_index('IJK', inplace= True)

        return var


    cpdef point2block(self, np.ndarray [double, ndim=1] x,
                            np.ndarray [double, ndim=1] y,
                            np.ndarray [double, ndim=1] z,
                            object prop,
                            str prop_name,
                            bint average = False,
                            bint overwrite = False):
        """point2block(np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z,object prop, str prop_name, bint average = False, bint overwrite = False)

        Assigns data in a point to a block

        Parameters
        ----------
        x,y,z : Float 1D numpy arrays
            coordinates of the points. Array shapes may be equal
        prop : array like
            array with property shape[0] == x.shape[0]
            note that prop can be an array of array, for example, an
            array of indicators arrays
        prop_name: str
            name of the property at block model
        average : bool
            If there are more than one point at block
            average the points. If true the average will be created,
            this may fail if prop is not numeric. If False the last
            point in the block will be used to assign the value.
        overwrite : bool
            If False an existing property will not be overwritten

        Example
        -------
        >>>
        >>> point2block( x, y, z, prop, prop_name= 'Au',  average = False, overwrite = False)
        >>>

        TODO: pass array of properties np.ndarray [double, ndim=2] prop and [prop_names]

        Note
        ----
        Points not laying in an existing block will be ignored.

        """

        cdef long i, n
        cdef np.ndarray [long, ndim=1] bad

        # note that prop can be an array of array with any type
        assert x.shape[0] == y.shape[0] == z.shape[0] == prop.shape[0]

        # check that the column exist and is not overwrite
        if overwrite==False:
            assert prop_name not in self.bmtable.columns, 'The field {} already exist, use overwrite = True to rewrite'.format(prop_name)

        bad = np.zeros(x.shape[0], dtype = int)


        # make sure there is ijk field in model

        assert 'IJK' in self.bmtable.columns, 'Error: no IJK field in block model'


        # get index ix, iy, iz
        ix = x2ix (x, self.xorg, self.dx)
        iy = x2ix (y, self.yorg, self.dy)
        iz = x2ix (z, self.zorg, self.dz)

        # calculate ijk
        ijk = ind2ijk(ix,iy,iz, self.nx,self.ny,self.nz)

        # set ijk as key for indexing
        n = x.shape[0]
        self.bmtable.set_index('IJK', inplace= True)

        if average:
            self.bmtable[prop_name] = 0.
        else:
            self.bmtable[prop_name] = np.nan

        for i in range(n):
            if (ijk[i] in self.bmtable.index) and (np.isfinite(prop[i])):
                #TODO: loop each property
                if average:
                    self.bmtable.at[ijk[i],prop_name] = (self.bmtable.at[ijk[i],prop_name] + prop[i])/2.
                else:
                    self.bmtable.at[ijk[i],prop_name] =  prop[i]
            else:
                bad[i]=1


        self.bmtable.reset_index('IJK', inplace= True)

        return bad


    cpdef fillwireframe(self,
                    object surface,
                    float toll=0,
                    bint overwrite=False):

        """fillwireframe(object surface, float toll=0, bint overwrite=False)

        Creates a full block model given a VTK surface (polydata) using
        vtkPolyDataToImageStencil. The results consist of blocks with
        parameter ``in``  with value between 0 and 1.

        ``in = 0`` represent blocks completely outside the wireframe
        ``in = 1`` represent blocks completely inside the wireframe
        ``1 >in > 0`` represent blocks cut by the wireframe

        The parameter in is calculated as the sum of corners of the
        block, each corner is equal to 1 if is inside the block or equal
        to 0 if is outside.

        The parameter ``in`` can be used as percentage of inclusion in
        the block with an error proportional to the size of the block
        and the resolution/complexity of the wireframe.

        The method may fail if the wireframe is narrower than a block,
        e.j. larg blocks in narrow veins. In this case you can set
        a tolerance value 1>=toll>=0. The output may not be used as
        volume, only as selection of blocks within a wireframe
        (with ``in>0``).


        Parameters
        ----------
        surface : VTK polydata
               this may work with any 3D object..., no necessarily triangles
        tol    : float, default 0.0
               number in interval [0. , 1.]
               you may use tol>0 if the blocks are large and the
               wireframe is narrower than the block size.
        overwrite : boolean
               overwrite flag, if true the entire block model will be
               overwritten

        Example
        -------
        >>>
        >>> mymodel.fillwireframe(surface, toll=0, overwrite=False)
        >>>

        Note
        ----
        The tolerance only works on x, y plane, not in z.

        This function works well with closed surfaces

        This function works well with large and complicated wireframes,
        for example, Leapfrog grade shells.

        """

        #check the model do not exist
        if overwrite==False:
            assert self.bmtable is None, 'Error: bmtable already exist, set overwrite=True to overwrite'


        #create and empty grid (the background)
        grid= vtk.vtkImageData()
        grid.SetSpacing([self.dx,self.dy,self.dz])
        grid.SetOrigin([self.xorg,self.yorg,self.zorg])
        grid.SetDimensions([self.nx+1,self.ny+1,self.nz+1])

        pscalars = vtk.vtkFloatArray()
        pscalars.SetName('__in')

        for i in range ((self.nx+1)*(self.ny+1)*(self.nz+1)):
            pscalars.InsertTuple1(i, 1.0) # use here 1.0 or 100.0

        grid.GetPointData().SetScalars(pscalars)

        # create the vtkfilter to select point corners of image voxels
        # in solid
        pol2stenc = vtk.vtkPolyDataToImageStencil()
        pol2stenc.SetInputData(surface)
        pol2stenc.SetOutputSpacing([self.dx,self.dy,self.dz])
        pol2stenc.SetOutputOrigin([self.xorg,self.yorg,self.zorg])
        pol2stenc.SetOutputWholeExtent(grid.GetExtent())

        #SetTolerance in order to include more blocks 50/50...
        # not work at elevation z
        pol2stenc.SetTolerance(toll)
        pol2stenc.Update()

        #create the image (regular grid)
        imgstenc = vtk.vtkImageStencil()
        imgstenc.SetInputData(grid)
        imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
        imgstenc.ReverseStencilOff();
        imgstenc.SetBackgroundValue(0.0);
        imgstenc.Update();

        #create the volume percentage (average of 8 '__in' point data)
        p2c = vtk.vtkPointDataToCellData()
        p2c.SetInputConnection(0, imgstenc.GetOutputPort(0))
        p2c.PassPointDataOn()
        p2c.Update()

        # populate the block model

        #first we create an empty full model
        self.create_IJK(overwrite=True)
        self.calc_ixyz_fromijk()
        self.calc_xyz_fromixyz()

        #then we get the data from vtk using vtk.util._numpy_support
        self.bmtable['__in']= vtknumpy.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray(0))

        # this return the actual model in VTK format, in case you what to use it
        return p2c.GetOutput()



    # -----------------------------------------
    #      GRID functions
    # -----------------------------------------

    cpdef create_gridxy(self, bint overwrite=False):
        """create_gridxy(bint overwrite=False)

        Creates a full block 2D grid in defined in the XY direction.

        The grid definition is taken from the block model definition in
        the X, and Y directions.


        Parameters
        ----------
        overwrite : boolean
               overwrite flag, if true the entire grid will be
               overwritten

        Example
        -------
        >>>
        >>> mymodel.create_gridxy(overwrite=False)
        >>>


        """

        #check the model do not exist
        if overwrite==False:
            assert self.grid2DXY is None, 'Error: grid2DXY already exist, set overwrite=True to overwrite'


		# create one layer model
        self.grid2DXY=pd.DataFrame({'IJK': np.arange(self.nx*self.ny, dtype=int)})
        self.grid2DXY['IX'],self.grid2DXY['IY'],self.grid2DXY['IZ'] =  ijk2ind(self.grid2DXY['IJK'].values,self.nx,self.ny,1)

        self.grid2DXY['XC']= ix2x (self.grid2DXY['IX'].values.astype('int'), self.xorg, self.dx)
        self.grid2DXY['YC']= ix2x (self.grid2DXY['IY'].values.astype('int'), self.yorg, self.dy)
        self.grid2DXY['ZC']= self.zorg





    cpdef create_gridxz(self, bint overwrite=False):
        """create_gridxz(bint overwrite=False)

        Creates a full block 2D grid in defined in the XZ direction.

        The grid definition is taken from the block model definition in
        the X, and Z directions.


        Parameters
        ----------
        overwrite : boolean
               overwrite flag, if true the entire grid will be
               overwritten

        Example
        -------
        >>>
        >>> mymodel.create_gridxz(overwrite=False)
        >>>


        """


        #check the model do not exist
        if overwrite==False:
            assert self.grid2DXZ is None, 'Error: grid2DXZ already exist, set overwrite=True to overwrite'

		# create one layer model
        self.grid2DXZ=pd.DataFrame({'IJK': np.arange(self.nx*self.nz, dtype=int)})
        self.grid2DXZ['IX'],self.grid2DXY['IY'],self.grid2DXY['IZ'] =  ijk2ind(self.grid2DXY['IJK'].values,self.nx,1,self.nz)

        self.grid2DXZ['XC']= ix2x (self.grid2DXY['IX'].values.astype('int'), self.xorg, self.dx)
        self.grid2DXZ['YC']= self.yorg
        self.grid2DXZ['ZC']= ix2x (self.grid2DXY['IZ'].values.astype('int'), self.zorg, self.dz)


    cpdef create_gridyz(self, bint overwrite=False):
        """create_gridyz(bint overwrite=False)

        Creates a full block 2D grid in defined in the YZ direction.

        The grid definition is taken from the block model definition in
        the Y, and Z directions.


        Parameters
        ----------
        overwrite : boolean
               overwrite flag, if true the entire grid will be
               overwritten

        Example
        -------
        >>>
        >>> mymodel.create_gridyz(overwrite=False)
        >>>


        """

        #check the model do not exist
        if overwrite==False:
            assert self.grid2DYZ is None, 'Error: grid2DYZ already exist, set overwrite=True to overwrite'

		# create one layer model
        self.grid2DYZ=pd.DataFrame({'IJK': np.arange(self.ny*self.nz, dtype=int)})
        self.grid2DYZ['IX'],self.grid2DXY['IY'],self.grid2DXY['IZ'] =  ijk2ind(self.grid2DXY['IJK'].values,1, self.ny,self.nz)

        self.grid2DYZ['XC']= self.xorg
        self.grid2DYZ['YC']= ix2x (self.grid2DXY['IY'].values.astype('int'), self.yorg, self.dy)
        self.grid2DYZ['ZC']= ix2x (self.grid2DXY['IZ'].values.astype('int'), self.zorg, self.dz)
