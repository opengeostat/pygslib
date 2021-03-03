'''
PyGSLIB vtktools, Module with tools to work with VTK.

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

import os
import sys
import vtk
import warnings
from IPython.display import Image
from libc.math cimport sin
from libc.math cimport cos
cimport numpy as np
import numpy as np
import vtk.util.numpy_support as vtknumpy
from vtk.numpy_interface import dataset_adapter as vtkdsa
from scipy.interpolate import Rbf
import scipy.spatial.distance as scipy_dist
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import ezdxf

# ----------------------------------------------------------------------
#   Functions for volume calculation
# ----------------------------------------------------------------------
cpdef calculate_volume(object mesh):
    """calculate_volume(object mesh)

    Takes a single closed vtkPolydata Mesh, composed by triangles,
    and calculate its volume.

    This function may produce wrong results if surfaces are not closed or
    poorly defined. Use results with caution.

    Parameters
    ----------
    mesh : vtkPolyData
        mesh surface (not tested with polylines...)


    Returns
    -------
    float with valume

    Note
    ----
    The volume is computed three times, with input as it is and with Normals
    recalculated in two different directions (flipped and no fliepped). Results
    are unreliable if the volumes are different than, let say 1%.

    """

    # normal flip
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(mesh)
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.ComputePointNormalsOn()
    normals.FlipNormalsOn()
    normals.Update()
    new_mesh= normals.GetOutput()

    # normal non flip
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(mesh)
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.ComputePointNormalsOn()
    normals.FlipNormalsOff()
    normals.Update()
    new_mesh2= normals.GetOutput()

    vol = vtk.vtkMassProperties()
    vol.SetInputData(mesh)
    vol0 = vol.GetVolume()

    vol = vtk.vtkMassProperties()
    vol.SetInputData(new_mesh)
    vol1 = vol.GetVolume()

    vol = vtk.vtkMassProperties()
    vol.SetInputData(new_mesh2)
    vol2 = vol.GetVolume()

    return vol0,vol1,vol2

# ----------------------------------------------------------------------
#   Functions for point querying
# ----------------------------------------------------------------------
cpdef vtk_raycasting(object surface, object pSource, object pTarget, object obbTree = None):
    """vtk_raycasting(object surface, object pSource, object pTarget, object obbTree)

    Intersects a line defined by two points with a polydata vtk object,
    for example a closed surface.

    This function is to find intersection and inclusion of points
    within, below or over a surface. This is handy to select points
    or define blocks inside solids.

    Parameters
    ----------
    surface : VTK polydata
           this may work with any 3D object..., no necessarily triangles
    pSource, pTarget: tuples with 3 float
           point location defining a ray
    obbTree: vtk.vtkOBBTree()
           tree generated in previous test


    Returns
    -------
    intersect : integer, with values -1, 0, 1
        0 means that no intersection points were found
        1 means that some intersection points were found
       -1 means the ray source lies within the surface
    pointsIntersection : tuple of tuples with 3 float
        this is a tuple with coordinates tuples defining
        intersection points

    pointsVTKIntersectionData : an Vtk polydata object
        similar to pointsIntersection but in Vtk polydata format

    obbTree : a vtk.vtkOBBTree object


    See Also
    --------
    vtk_show, loadSTL

    Note
    ----
    This Code was modified from the original code by Adamos Kyriakou
    published in https://pyscience.wordpress.com/

    """
    cdef int intersect, idx

    if obbTree is None:
      obbTree = vtk.vtkOBBTree()
      obbTree.SetDataSet(surface)
      obbTree.BuildLocator()

    pointsVTKintersection = vtk.vtkPoints()
    pointsVTKintersection.SetDataTypeToDouble()


    #Perform the ray-casting. The IntersectWithLine method returns an int
    #code which if 0 means that no intersection points were found.
    #If equal to 1 then some points were found.
    #If equal to -1 then the ray source lies within the surface
    intersect=obbTree.IntersectWithLine(pSource, pTarget, pointsVTKintersection, None)

    #'Convert' the intersection points into coordinate tuples and store them
    #under the list named pointsIntersection
    pointsVTKIntersectionData = pointsVTKintersection.GetData()
    noPointsVTKIntersection = pointsVTKIntersectionData.GetNumberOfTuples()
    pointsIntersection = []

    # >>>>>>>> optimize this
    for idx in range(noPointsVTKIntersection):
        _tup = pointsVTKIntersectionData.GetTuple3(idx)
        pointsIntersection.append(_tup)

    return intersect, pointsIntersection, pointsVTKIntersectionData, obbTree


cpdef pointquering(object surface,
                   double azm,
                   double dip,
                   np.ndarray [double, ndim=1] x,
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   int test,
                   object obbTree = None):
    """pointquering(object surface, double azm, double dip, np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, int test, object obbTree = None)

    Find points inside, over and below surface/solid


    Parameters
    ----------
    surface : VTK polydata
           this may work with any 3D object..., no necessarily triangles
    azm, dip: float
           rotation defining the direction we will use to test the points
           azm 0 will point north and dip positive meas downward direction
           (like surface drillholes)
    x,y,z   : 1D array of floats
           coordinates of the points to be tested
    test    : integer
           If ``test==1`` it queries points inside closed surface. It
           only works with closed solids.
           If ``test==2`` it queries point above surface.
           If ``test==3`` it queries point below surface.
           If ``test==4`` it queries point above and below the surface.
           ``test==4`` is similar to ``test==1`` but it works with
           open surfaces.
    obbTree: vtk.vtkOBBTree()
           obb tree. If not avaliable wil be generated.
           Building the tree is time consuming. You can generated a obbtree 
           with `pygslib.vtktools.vtk_raycasting`,  and then use it as 
           input if you plan to run multiple iterations of the pointquering function. 

    Returns
    -------
    inside : 1D array of integers
        Indicator of point inclusion with values [0,1]
        0 means that the point is not inside, above or below surface
        1 means that the point is inside, above or below surface
    p1 : 1D array size 3
        rays to show where are pointing in the output


    See Also
    --------
    vtk_raycasting

    Note
    ----
    The test 1 requires the surface to be close, to find points between
    two surfaces you can use test=4
    The test, 2,3 and 4 use raycasting and rays orientation are defined
    by azm and dip values. The results will depend on this direction,
    for example, to find the points between two open surfaces you may
    use ray direction perpendicular to surfaces
    If a surface is closed the points over and below a surface will
    include points inside surface

    """


    assert x.shape[0]==y.shape[0]==z.shape[0]

    cdef int intersect, intersect2
    cdef float razm
    cdef float rdip
    cdef float DEG2RAD
    cdef np.ndarray [long, ndim=1] inside = np.zeros([x.shape[0]], dtype= int)
    cdef np.ndarray [double, ndim=1] p0 = np.empty([3], dtype= float)
    cdef np.ndarray [double, ndim=1] p1 = np.empty([3], dtype= float)
    cdef np.ndarray [double, ndim=1] p2 = np.empty([3], dtype= float)


    DEG2RAD=3.141592654/180.0
    razm = azm * DEG2RAD
    rdip = -dip * DEG2RAD


    # construct obbTree: see oriented bounding box (OBB) trees
    # see http://gamma.cs.unc.edu/SSV/obb.pdf

    if obbTree is None:
      obbTree = vtk.vtkOBBTree()
      obbTree.SetDataSet(surface)
      obbTree.BuildLocator()

    pointsVTKintersection = vtk.vtkPoints()
    pointsVTKintersection.SetDataTypeToDouble()

    # test each point
    for i in range(x.shape[0]):

        #create vector of the point tested (target)
        p0[0]=x[i]
        p0[1]=y[i]
        p0[2]=z[i]

        if test == 1:
            # The return value is +1 if outside, -1 if inside, and 0 if undecided.
            intersect=obbTree.InsideOrOutside(p0)   # see also vtkSelectEnclosedPoints
            if intersect==-1:
                inside[i]=1

            continue

        if test > 1:
            # in this case we use the ray intersection from point tested
            # and a source point located at infinit in the azm, dip dir

            # convert degree to rad and correct sign of dip


            # calculate the coordinate of an unit vector
            p1[0] = sin(razm) * cos(rdip)  # x
            p1[1] = cos(razm) * cos(rdip)  # y
            p1[2] = sin(rdip)              # z

            if test == 2:
                # rescale to infinite
                p1[:]=p1[:]*1e+15
                # and translate
                p1[0]=p1[0]+x[i]
                p1[1]=p1[1]+y[i]
                p1[2]=p1[2]+z[i]

                # 0 if no intersections, -1 if point 'p0' lies inside
                # closed surface, +1 if point 'a0' lies outside the closed surface
                intersect=obbTree.IntersectWithLine(p0, p1, None, None)
                inside[i]=abs(intersect)

                continue

            if test == 3:
                # invert and rescale to infinite
                p1[:]=-p1[:]*1e+15
                # and translate
                p1[0]=p1[0]+x[i]
                p1[1]=p1[1]+y[i]
                p1[2]=p1[2]+z[i]

                # 0 if no intersections, -1 if point 'p0' lies inside
                # closed surface, +1 if point 'a0' lies outside the closed surface
                intersect=obbTree.IntersectWithLine(p0, p1, None, None)
                inside[i]=abs(intersect)

                continue

            if test == 4:
                # invert and rescale to infinite
                p1[:]=p1[:]*1e+15
                p2[:]=-p1[:]
                # and translate
                p1[0]=p1[0]+x[i]
                p1[1]=p1[1]+y[i]
                p1[2]=p1[2]+z[i]

                p2[0]=p2[0]+x[i]
                p2[1]=p2[1]+y[i]
                p2[2]=p2[2]+z[i]

                # 0 if no intersections, -1 if point 'p0' lies inside
                # closed surface, +1 if point 'a0' lies outside the closed surface
                intersect = obbTree.IntersectWithLine(p0, p1, None, None)
                intersect2 = obbTree.IntersectWithLine(p0, p2, None, None)
                if intersect!=0 and intersect2!=0:
                    inside[i]=1
                else:
                    inside[i]=0

                continue


    # Prepare a ray to see where our rays are pointing in the output
    p1[0] = sin(razm) * cos(rdip)  # x
    p1[1] = cos(razm) * cos(rdip)  # y
    p1[2] = sin(rdip)              # z
    p1[:]=p1[:]*1e+2


    return inside, p1



cpdef pointinsolid(object surface,
                   np.ndarray [double, ndim=1] x,
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   double tolerance = .000001):
    """pointinsolid(object surface, np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, double tolerance = .000001)

    Finds points inside a closed surface using vtkSelectEnclosedPoints


    Parameters
    ----------
    surface : VTK polydata
           this may work with any 3D object..., no necessarily triangles
    x,y,z   : 1D array of floats
           coordinates of the points to be tested


    Returns
    -------
    inside : 1D array of integers
        Indicator of point inclusion with values [0,1]
        0 means that the point is not inside
        1 means that the point is inside

    See Also
    --------
    vtk_raycasting, pointquering

    Note
    ----
    It is assumed that the surface is closed and manifold.
    Note that if this check is not performed and the surface is not
    closed, the results are undefined.

    TODO: add check to verify that a surface is closed

    """


    cdef int i

    assert x.shape[0]==y.shape[0]==z.shape[0]


    points = vtk.vtkPoints();
    points.SetDataTypeToDouble()

    for i in range(x.shape[0]):
        points.InsertNextPoint(x[i],y[i],z[i]);

    pointsPolydata = vtk.vtkPolyData();
    pointsPolydata.SetPoints(points);

    #Points inside test
    selectEnclosedPoints = vtk.vtkSelectEnclosedPoints();
    selectEnclosedPoints.SetTolerance(tolerance)
    selectEnclosedPoints.SetInputData(pointsPolydata);
    selectEnclosedPoints.SetSurfaceData(surface);
    selectEnclosedPoints.Update();

    insideArray = vtk.vtkDataArray.SafeDownCast(
        selectEnclosedPoints.GetOutput().GetPointData().GetArray("SelectedPoints"));

    inside = vtknumpy.vtk_to_numpy(insideArray)


    return inside


# ----------------------------------------------------------------------
#   Functions to modify polydata
# ----------------------------------------------------------------------
cpdef GetPointsInPolydata(object polydata):
    """GetPointsInPolydata(object polydata)

    Returns the point coordinates in a polydata (e.j. wireframe)
    as a numpy array.


    Parameters
    ----------
    polydata : VTK polydata
        vtk object with data, ej. wireframe


    Returns
    -------
    p : 2D array of float
        Points coordinates in polydata object

    See also
    --------
    SetPointsInPolydata, SavePolydata

    """
    npoints = polydata.GetNumberOfPoints()

    p = np.empty([npoints,3])

    for i in range(npoints):
        p[i,0],p[i,1],p[i,2] = polydata.GetPoint(i)

    return p


cpdef SetPointsInPolydata(object polydata, object points):
    """SetPointsInPolydata(object polydata, object points)

    Set the points coordinates in a VTK polydata object (e.j. wireframe).


    Parameters
    ----------
    polydata : VTK polydata
        vtk object with data, ej. wireframe
    points   : 2D array
        New points coordinates to be set in the polydata
        points shape may be equal to [polydata.GetNumberOfPoints(),3]


    See also
    --------
    GetPointsInPolydata, SavePolydata

    """
    npoints = polydata.GetNumberOfPoints()

    assert npoints == points.shape[0], 'Error: Wrong shape[0] in points array'
    assert points.shape[1]==3, 'Error: Wrong shape[1] in points array'

    p = vtk.vtkPoints()
    p.SetDataTypeToDouble()
    p.SetNumberOfPoints(npoints)


    for i in range(npoints):
        p.InsertPoint(i,(points[i,0],points[i,1],points[i,2]))

    polydata.SetPoints(p)


cpdef getbounds(object polydata):
    """getbounds(object polydata)

    Returns the bounding box limits x,y,z minimum and maximum


    Parameters
    ----------
    polydata : VTK polydata
        vtk object with data, ej. wireframe


    Returns
    -------
    xmin,xmax,ymin,ymax,zmin,zmax : float
        The geometry bounding box

    """

    return polydata.GetBounds()

# ----------------------------------------------------------------------
#   Functions to generate vtk 3D block
# ----------------------------------------------------------------------
cpdef grid2vtkImageData(
               int nx, int ny, int nz,
               double xorg, double yorg, double zorg,
               double dx, double dy, double dz,
               object cell_data={},
               object point_data={}):
    """grid2vtkImageData(nx, ny, nz, xorg, yorg, zorg, dx, dy, dz, cell_data, point_data)

    Utility function to creates a vtkImageData object from full grid dataself.
    It is designed to be used from `Blockmodel.blocks2vtkImageData`

    Data must be ordered in a ix, iy, iz sequence.


    Parameters
    ----------
    nx,ny,nz : int
        number of rows, columns and leves of the grid
    xorg,yorg,zorg : double
        origing of coordinates of the grid
    dx,dy,dz : double
        block dimensions
    cell_data : dictionary
        dictionary with {'dataname':datavalues} with righ size for cell
    point_data: dictionary
        dictionary with {'dataname':datavalues} with righ size for points

    See also
    --------
    blockmodel.Blockmodel.blocks2vtkImageData
    """

    ufgrid = vtk.vtkImageData()
    ufgrid.SetOrigin(xorg,yorg,zorg)
    ufgrid.SetSpacing(dx,dy,dz)
    ufgrid.SetDimensions(nx+1,ny+1,nz+1)

    for i in cell_data:
      ufgrid.GetCellData().AddArray(vtkdsa.numpyTovtkDataArray(cell_data[i], name=i, array_type=None))

    for i in point_data:
      ufgrid.GetPointData().AddArray(vtkdsa.numpyTovtkDataArray(point_data[i], name=i, array_type=None))

    return ufgrid


cpdef partialgrid2vtkfile(
                   np.ndarray [double, ndim=1] x,
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   double DX,
                   double DY,
                   double DZ,
                   object var,
                   object varname,
                   object dx = None,
                   object dy = None,
                   object dz = None,
                   str path = None):
    """partialgrid2vtkfile(np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, double DX, double DY, double DZ, object var, object varname, np.ndarray [double, ndim=1] dx, np.ndarray [double, ndim=1] dy, np.ndarray [double, ndim=1] dz, str path)

    Saves data in the cells of a VTK Unstructured Grid file


    Parameters
    ----------
    x,y,z : np.ndarray
        coordinates of the points
    DX,DY,DZ : float
        block sizes
    var : list of 1D arrays
        array with variable values.
    varname : list of strings
        variables names.
    dx, dy, dz: None or array, default None
        variable block size in each direction. If None DX,DY,DZ are used
    path : str (optional)
        file path (relative or absolute) and name. If None file will not be saved
    Returns
    ------
    vtk unstructured grid
    """

    # number of cells/blocks
    nc = x.shape[0]

    if dx is None:
      dx = np.ones(nc)*DX

    if dy is None:
      dy = np.ones(nc)*DY

    if dz is None:
      dz = np.ones(nc)*DZ

    #number of variables
    nvar = len(var)

    # number of points (block vertex)
    npt = nc*8

    # Create array of the points and ID
    pcoords = vtk.vtkFloatArray()
    pcoords.SetNumberOfComponents(3)
    pcoords.SetNumberOfTuples(npt)

    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    voxelArray = vtk.vtkCellArray()

    # create vertex (points)
    id=0
    #for each block
    for i in range (nc):
        # for each vertex
        pcoords.SetTuple3(id, x[i]+dx[i]/2., y[i]-dy[i]/2., z[i]-dz[i]/2.)
        pcoords.SetTuple3(id+1, x[i]-dx[i]/2., y[i]-dy[i]/2., z[i]-dz[i]/2.)
        pcoords.SetTuple3(id+2, x[i]+dx[i]/2., y[i]+dy[i]/2., z[i]-dz[i]/2.)
        pcoords.SetTuple3(id+3, x[i]-dx[i]/2., y[i]+dy[i]/2., z[i]-dz[i]/2.)
        pcoords.SetTuple3(id+4, x[i]+dx[i]/2., y[i]-dy[i]/2., z[i]+dz[i]/2.)
        pcoords.SetTuple3(id+5, x[i]-dx[i]/2., y[i]-dy[i]/2., z[i]+dz[i]/2.)
        pcoords.SetTuple3(id+6, x[i]+dx[i]/2., y[i]+dy[i]/2., z[i]+dz[i]/2.)
        pcoords.SetTuple3(id+7, x[i]-dx[i]/2., y[i]+dy[i]/2., z[i]+dz[i]/2.)

        id+=8

    # add points to the cell
    points.SetData(pcoords)

    # Create the cells.
    #for each block
    id=-1
    for i in range (nc):
        # add next cell
        voxelArray.InsertNextCell(8)


        # for each vertex
        for j in range (8):
            id+=1
            voxelArray.InsertCellPoint(id)

    # create the unestructured grid
    ug = vtk.vtkUnstructuredGrid()
    # Assign points and cells
    ug.SetPoints(points)
    ug.SetCells(vtk.VTK_VOXEL, voxelArray)

    # asign scalar
    for i in range(nvar):
        # implement string variables here
        dtype = var[i].dtype

        if  dtype==np.int8 or dtype==np.int16 or dtype==np.int32 or dtype==np.int64 or dtype==np.float16 or dtype==np.float32 or dtype==np.float64:
            cscalars = vtknumpy.numpy_to_vtk(var[i])
            cscalars.SetName(varname[i])
            ug.GetCellData().AddArray(cscalars)
        else:
            # this is fos string array. Not optimized...
            cscalars = vtk.vtkStringArray()
            cscalars.SetName(varname[i])
            cscalars.SetNumberOfComponents(1)
            cscalars.SetNumberOfTuples(nc)
            for l in range(nc):
                cscalars.SetValue(l,str(var[i][l]))
            ug.GetCellData().AddArray(cscalars)

    # Clean before saving...
    # this will remove duplicated points
    extractGrid = vtk.vtkExtractUnstructuredGrid()
    extractGrid.SetInputData(ug)
    extractGrid.PointClippingOff()
    extractGrid.ExtentClippingOff()
    extractGrid.CellClippingOn()
    extractGrid.MergingOn()
    extractGrid.SetCellMinimum(0)
    extractGrid.SetCellMaximum(nc)
    extractGrid.Update()

    if path is not None:
      # add extension to path
      if not path.lower().endswith('.vtu'):
          path = path + '.vtu'

      writer = vtk.vtkXMLUnstructuredGridWriter();
      writer.SetFileName(path);
      writer.SetInputData(extractGrid.GetOutput())
      writer.Write()

    return extractGrid.GetOutput()


cpdef rotate_vtk(float xorg, float yorg, float zorg, float ang1, float ang2, float ang3, object obj):
    """rotate_vtk(float xorg, float yorg, float zorg, float ang1, float ang2, float ang3, object obj)

    rotate a vtk object


    Parameters
    ----------
    xorg , yorg, zorg : float
        pivot point coordinates used for rotation
    ang1, ang2, ang3 : float
        angles of ZXZ rotation
    polyData : object
        vtk polyData object to be rotated

    Returns
    ---
    VTK object with rotations applied
    """
    # translate
    transL1 = vtk.vtkTransform()
    transL1.Translate( -xorg, -yorg, -zorg)

    tf = vtk.vtkTransformFilter()
    tf.SetInputData(obj)
    tf.SetTransform(transL1)
    tf.Update()

    # rotate 1
    tmp = tf.GetOutput()
    transL1 = vtk.vtkTransform()
    transL1.RotateZ(ang1)

    tf = vtk.vtkTransformFilter()
    tf.SetInputData(tmp)
    tf.SetTransform(transL1)
    tf.Update()

    # rotate 2
    tmp = tf.GetOutput()
    transL1 = vtk.vtkTransform()
    transL1.RotateX(ang2)

    tf = vtk.vtkTransformFilter()
    tf.SetInputData(tmp)
    tf.SetTransform(transL1)
    tf.Update()

    # rotate 3
    tmp = tf.GetOutput()
    transL1 = vtk.vtkTransform()
    transL1.RotateX(ang3)

    tf = vtk.vtkTransformFilter()
    tf.SetInputData(tmp)
    tf.SetTransform(transL1)
    tf.Update()

    return tf.GetOutput()


cpdef back_rotate_vtk(float xorg, float yorg, float zorg, float ang1, float ang2, float ang3, object obj):
    """rotate_vtk(float xorg, float yorg, float zorg, float ang1, float ang2, float ang3, object obj)

    back rotate/transform a vtk object


    Parameters
    ----------
    xorg , yorg, zorg : float
        pivot point coordinates used for rotation
    ang1, ang2, ang3 : float
        angles of ZXZ rotation
    polyData : object
        vtk polyData object to be rotated

    Returns
    ---
    VTK object with rotations applied
    """

    # rotate 3
    transL1 = vtk.vtkTransform()
    transL1.RotateX(-ang3)

    tf = vtk.vtkTransformFilter()
    tf.SetInputData(obj)
    tf.SetTransform(transL1)
    tf.Update()

    # rotate 2
    tmp = tf.GetOutput()
    transL1 = vtk.vtkTransform()
    transL1.RotateX(-ang2)

    tf = vtk.vtkTransformFilter()
    tf.SetInputData(tmp)
    tf.SetTransform(transL1)
    tf.Update()

    # rotate 1
    tmp = tf.GetOutput()
    transL1 = vtk.vtkTransform()
    transL1.RotateZ(-ang1)

    tf = vtk.vtkTransformFilter()
    tf.SetInputData(tmp)
    tf.SetTransform(transL1)
    tf.Update()

    # translate
    tmp = tf.GetOutput()
    transL1 = vtk.vtkTransform()
    transL1.Translate( xorg, yorg, zorg)

    tf = vtk.vtkTransformFilter()
    tf.SetInputData(tmp)
    tf.SetTransform(transL1)
    tf.Update()

    return tf.GetOutput()


cpdef partialgrid2vtkfile_p(str path,
                   np.ndarray [double, ndim=1] x,
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   double DX,
                   double DY,
                   double DZ,
                   object var,
                   object varname):
    """partialgrid2vtkfile_p(str path, np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, double DX, double DY, double DZ, object var, object varname)

    Saves data in the points of a VTK Unstructured Grid file


    Parameters
    ----------
    path : str
        file path (relative or absolute) and name
    x,y,z : np.ndarray
        coordinates of the points
    DX,DY,DZ : float
        block sizes
    data : array like
        array with variable values.
    """

    # add extension to path
    if not path.lower().endswith('.vtu'):
        path = path + '.vtu'


    # number of points/blocks cetroids
    np = x.shape[0]

    #number of variables
    nvar = len(var)

    # Create array of the points and ID
    pcoords = vtk.vtkFloatArray()
    pcoords.SetNumberOfComponents(3)
    pcoords.SetNumberOfTuples(np)

    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    pointArray = vtk.vtkCellArray()


    #for each block
    for i in range (np):
        # for each vertex
        pcoords.SetTuple3(i, x[i], y[i], z[i])

    # add points
    points.SetData(pcoords)

    # Create the cells.
    #for each block
    for i in range (np):
        # add next cell
        pointArray.InsertNextCell(1)
        pointArray.InsertCellPoint(i)

    # create the unestructured grid
    ug = vtk.vtkUnstructuredGrid()
    # Assign points and cells
    ug.SetPoints(points)
    ug.SetCells(vtk.VTK_VERTEX, pointArray)

    # asign scalar
    for i in range(nvar):
        cscalars = vtknumpy.numpy_to_vtk(var[i])
        cscalars.SetName(varname[i])
        ug.GetPointData().AddArray(cscalars)

    # save results
    writer = vtk.vtkXMLUnstructuredGridWriter();
    writer.SetFileName(path);
    writer.SetInputData(ug)
    writer.Write()


# ----------------------------------------------------------------------
#   triangilations
# ----------------------------------------------------------------------
def delaunay2D (   np.ndarray [double, ndim=1] x,
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   constraints = None):
    """delaunay2D (np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, constraints = None):
    Creates a triangulated Surface


    Parameters
    ----------
    x,y,z : np.ndarray [double, ndim=1]
        Coordinates of the input points
    constraints: vtkPolydata or None
        constraint polygons, lines or polylines

    Returns
    -------
    vtkPolyData with wireframe

    """

    # number of cells/blocks
    np = x.shape[0]


    # Create array of points
    pcoords = vtk.vtkFloatArray()
    pcoords.SetNumberOfComponents(3)
    pcoords.SetNumberOfTuples(np)

    for i in range (np):
        pcoords.SetTuple3(i, x[i], y[i], z[i])

    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    points.SetData(pcoords)


    # create input polydata
    myPolyData = vtk.vtkPolyData()
    myPolyData.SetPoints(points)


    # Triangulate
    delaunay = vtk.vtkDelaunay2D()
    delaunay.SetInputData(myPolyData)

    if constraints != None:
        delaunay.SetSourceData(constraints)

    delaunay.Update()

    return delaunay.GetOutput()


def delaunay2D_polylines (object polydata, constraints = None):
    """delaunay2D (polydata, constraints = None):

    Creates a triangulated Surface from vtkPolyData object


    Parameters
    ----------
    polydata : vtkPolyData
        polylines or any other vtkPolyData with points
    constraints: vtkPolydata or None
        constraint polygons, lines or polylines
    Returns
    -------
    vtkPolyData with wireframe

    """

    # Triangulate
    delaunay = vtk.vtkDelaunay2D()
    delaunay.SetInputData(polydata)

    if constraints != None:
        delaunay.SetSourceData(constraints)

    delaunay.Update()

    return delaunay.GetOutput()

# ----------------------------------------------------------------------
#   Functions to interpolate surfaces
# ----------------------------------------------------------------------
cpdef rbfinterpolate(np.ndarray [double, ndim=1] x,
                 np.ndarray [double, ndim=1] y,
                 np.ndarray [double, ndim=1] z,
                 object xg = None,
                 object yg = None,
                 object xorg = None,
                 object yorg = None,
                 object dx = None,
                 object dy = None,
                 object nx = None,
                 object ny = None,
                 double tol=0.01,
                 double tolg = 0.1,
                 str method = 'linear',
                 double epsilon=100,
                 constraints = None,
                 bint snap = True,
                 bint remove_duplicates = True,
                 bint triangulate = True):
    """
    rbfinterpolate(  np.ndarray [double, ndim=1] x,
                     np.ndarray [double, ndim=1] y,
                     np.ndarray [double, ndim=1] z,
                     object xg = None,
                     object yg = None,
                     object xorg = None,
                     object yorg = None,
                     object dx = None,
                     object dy = None,
                     object nx = None,
                     object ny = None,
                     double tol=0.01,
                     double tolg = 0.1,
                     str method = 'linear',
                     double epsilon=100,
                     constraints = None,
                     bint snap = True,
                     bint remove_duplicates = True,
                     bint triangulate = True)
    Creates a grid surface interpolated with rbf. Optionaly will include (snap)
    input points.


    Parameters
    ----------
    x,y,z : np.ndarray [double, ndim=1]
        Coordinates of the input points
    xg,yg : np.ndarray [double, ndim=1]
        Coordinates of the grid (or irregular) target points
    xorg, yorg, dx, dy: floats
        coordinates of the lower left point in the 2D grid and grid spacing
    nx, ny: int
        number of rows and cols of the 2D grid
    tol: double (default 0.01)
        an error/warn will be raised if any pair of input points is within `tol` distance
    tolg: double (default 0.1)
        points xg,yg will be removed if distance from closes x,y point is within `tolg` distance
    method: str (default 'linear')
        any of ['multiquadric', 'inverse', 'gaussian', 'linear', 'cubic', 'quintic', 'thin_plate']
    epsilon: double (default 100)
        Adjustable constant for gaussian or multiquadrics functions - defaults
        to approximate average distance between nodes (which is a good start).
    constraints: vtkPolydata or None
        constraint polygons, lines or polylines
    snap: boolean (default True)
        if True input points [x,y,z] will be append to target points [xg,yg,zg]
    remove_duplicates: boolean (default True)
        remove the second duplicated points located within distance = tol

    Returns
    -------
    vtkPolyData with wireframe

    """

    assert method in ['multiquadric', 'inverse', 'gaussian', 'linear', 'cubic', 'quintic', 'thin_plate']
    assert (xg is None)==(yg is None)

    if xg is None:
      xg,yg = np.meshgrid(np.arange(xorg,xorg+dx*(nx+1),dx),np.arange(yorg,yorg+dy*(ny+1),dy))
      xg =xg.ravel()
      yg =yg.ravel()

    #Check for duplicates around m metres for computational stability
    f = np.ones(shape=(len(x)), dtype=bool)

    if remove_duplicates:
      ktree = cKDTree(np.stack((x, y), axis=-1))
      duplicates = ktree.query_pairs(r=tol)

      if len(duplicates)>0:
        warnings.warn("There are {} duplicated points and second set of points was removed:\n {}".format(len(duplicates)>0,duplicates))
        f[np.array(tuple(duplicates))[:,1]]=False


    # generate rbf and interpolate in grid
    rbfi = Rbf(x[f], y[f], z[f], epsilon=epsilon, function = method)
    zg = rbfi(xg, yg)

    # Remove points from grid close to data
    if remove_duplicates and snap:
      mask = np.ones(xg.shape[0], dtype = bool)
      if tolg is None:
        tolg = dx/2.2
      for i in range(xg.shape[0]):
        if ktree.query([xg[i],yg[i]])[0]<tolg:
          mask[i] = False

      xg = xg[mask]
      yg = yg[mask]
      zg = zg[mask]

    # merge data
    if snap:
      # get max to ignore points out of the grid limits
      xmax = xorg+dx*nx
      ymax = yorg+dy*ny
      ff = (x>=xorg) & (x<=xmax) & (y>=yorg ) & (y<=ymax) & (f)
      xa = np.concatenate((x[ff],xg))
      ya = np.concatenate((y[ff],yg))
      za = np.concatenate((z[ff],zg))
      # triangulate and calculate normals
      if triangulate:
        mesh = delaunay2D (xa, ya, za, constraints = constraints)
        return calculate_normals(mesh),xa, ya, za
      else:
        return None,xa, ya, za
    else:
      if triangulate:
        mesh = delaunay2D (xg, yg, zg, constraints = constraints)
        return calculate_normals(mesh),xg, yg, zg
      else:
        return None,xg, yg, zg


# ----------------------------------------------------------------------
#   Functions to read write vtk and other formats
# ----------------------------------------------------------------------
cpdef loadVTP(str filenameVTP):
    """loadVTP(str filenameVTP)

    Load an XML VTK Polydata file

    Parameters
    ----------
    filenameVTP : file path

    Returns
    -------
    polydata : VTK Polydata object



    """

    # check file exists
    assert os.path.isfile(filenameVTP), 'Error: file path'

    readerVTP = vtk.vtkXMLPolyDataReader()
    readerVTP.SetFileName(filenameVTP)
    # 'update' the reader i.e. read the .VTP file
    readerVTP.Update()

    polydata = readerVTP.GetOutput()

    # If there are no points in 'vtkPolyData' something went wrong
    if polydata.GetNumberOfPoints() == 0:
        raise ValueError(
            "No point data could be loaded from '" + filenameVTP)
        return None

    return polydata


cpdef loadSTL(str filenameSTL):
    """loadSTL(str filenameSTL)

    Load a STL wireframe file

    Parameters
    ----------
    filenameSTL : file path

    Returns
    -------
    polydata : VTK Polydata object

    """

    # check file exists
    assert os.path.isfile(filenameSTL), 'Error: file path'

    readerSTL = vtk.vtkSTLReader()
    readerSTL.SetFileName(filenameSTL)
    # 'update' the reader i.e. read the .stl file
    readerSTL.Update()

    polydata = readerSTL.GetOutput()

    # If there are no points in 'vtkPolyData' something went wrong
    if polydata.GetNumberOfPoints() == 0:
        raise ValueError(
            "No point data could be loaded from '" + filenameSTL)
        return None

    return polydata

cpdef dmtable2wireframe(
                   x,
                   y,
                   z,
                   pid1,
                   pid2,
                   pid3,
                   bint indexone = False,
                   str filename = None):
    """dmtable2wireframe(np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, np.ndarray [long, ndim=1] pid1, np.ndarray [long, ndim=1] pid2, np.ndarray [long, ndim=1] pid3, bint indexone = False, str filename = None)

    Takes a wireframe defined by two tables and creates a VTK polydata wireframe object.

    The input tables are as follow:

    - x, y, z : This is the points table
    - pid1,pid2,pid3: This is another table defining triangles with
      three point IDs.

    Parameters
    ----------
    x,y,z          : 1D array of floats
        coordinates of the points to be tested
    pid1,pid2,pid3 : 1D array of integers
        triangles defined by three existing point IDs
    indexone:  boolean (Default False)
        If false pid index start at zero, otherwise pid index
        start at one
    filename   : Str (Default None)
        file name and path. If provided the file wireframe is saved to
        this file.

    Returns
    -------
    surface : VTK polydata
           wireframe imported

    TODO
    ----
    Add variables to triangles (vtk cell) and to points (vtk points)


    Note
    ----
    ``pid1,pid2,pid3`` are the row number in the point table of the
    three points that form a triangle.

    ``pid`` indices start at zero of indexone == False,
    otherwise ``pid`` indices start at one

    For example:

    a point table

    +----+----+----+
    | x  | y  | z  |
    +====+====+====+
    | .0 | .0 | .0 |
    +----+----+----+
    | 1. | .0 | .0 |
    +----+----+----+
    | 0. | 1. | .0 |
    +----+----+----+

    a triangle table

    +----+----+----+
    |pid1|pid2|pid3|
    +====+====+====+
    | 0  | 1  | 2  |
    +----+----+----+

    """


    cdef int i

    assert x.shape[0]==y.shape[0]==z.shape[0],          'Error: x, y, z or pid with different dimension'
    assert pid1.shape[0]==pid1.shape[0]==pid1.shape[0], 'Error: pid1,pid2 or pid3 with different dimension'

    # reset index to start from zero if they start at one
    if indexone == True:
        pid1 = pid1 - 1
        pid2 = pid2 - 1
        pid3 = pid3 - 1

    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()

    for i in range(x.shape[0]):
        points.InsertNextPoint(x[i],y[i],z[i]);


    #create triangles
    triangle = vtk.vtkTriangle()
    triangles = vtk.vtkCellArray()

    for i in range(pid1.shape[0]):
        triangle.GetPointIds().InsertId ( 0, int(pid1[i]));
        triangle.GetPointIds().InsertId ( 1, int(pid2[i]));
        triangle.GetPointIds().InsertId ( 2, int(pid3[i]));

        triangles.InsertNextCell ( triangle )

    # Create a polydata object
    trianglePolyData = vtk.vtkPolyData()
    # Add the geometry and topology to the polydata
    trianglePolyData.SetPoints(points)
    trianglePolyData.SetPolys(triangles)

    # if required save the file
    if filename!=None:
        # check extension
        if not filename.endswith('.mp3'):
            filename = filename + '.vtp'

        writer =  vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(trianglePolyData)
        writer.SetDataModeToBinary()
        writer.Write();


    return trianglePolyData


cpdef SavePolydata(object polydata, str path):
    """SavePolydata(object polydata, str path)

    Saves polydata into a VTK XML polydata file ('*.vtp')


    Parameters
    ----------
    polydata : VTK polydata
        vtk object with data, ej. wireframe
    path : str
        Extension (*.vtp) will be added if not provided


    """

    assert polydata.GetClassName()=='vtkPolyData', 'error input vtk object is of type {}, a vtkPolyData was expected'.format(polydata.GetClassName())


    # add extension to path
    if not path.lower().endswith('.vtp'):
        path = path + '.vtp'

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(polydata)
    writer.Write()

cpdef SaveSTL(object polydata, str path):
    """SaveSTL(object polydata, str path)

    Saves polydata into a STL file ('*.stl')


    Parameters
    ----------
    polydata : VTK polydata
        vtk object with data, ej. wireframe
    path : str
        Extension (*.stl) will be added if not provided


    """

    assert polydata.GetClassName()=='vtkPolyData', 'error input vtk object is of type {}, a vtkPolyData was expected'.format(polydata.GetClassName())


    # add extension to path
    if not path.lower().endswith('.stl'):
        path = path + '.stl'

    writer = vtk.vtkSTLWriter()
    writer.SetFileName(path)
    writer.SetInputData(polydata)
    writer.Write()

cpdef LoadImageData(str path):
    """ReadImageData(object grid, str path)

    Read a vtkImageData ('*.vti')

    Parameters
    ----------
    path : str
        Extension (*.vti) will be added if not provided

    Returns
    -------
    object : vtkImageData
    """

    # add extension to path
    if not path.lower().endswith('.vti'):
        path = path + '.vti'

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(path)
    reader.Update()

    return reader.GetOutput()


cpdef SaveImageData(object grid, str path):
    """SaveImageData(object grid, str path)

    Saves a vtkImageData into a VTK XML grid file ('*.vti')

    Parameters
    ----------
    grid : vtkImageData
        vtk object with image (regular) grid
    path : str
        Extension (*.vti) will be added if not provided
    """

    assert grid.GetClassName()=='vtkImageData', 'error input vtk object is of type {}, a vtkImageData was expected'.format(grid.GetClassName())

    # add extension to path
    if not path.lower().endswith('.vti'):
        path = path + '.vti'

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(path)
    writer.SetDataModeToBinary()
    writer.SetInputData(grid)
    writer.Write()

cpdef SaveRectilinearGrid(object grid, str path):
    """SaveRectilinearGrid(object grid, str path)

    Saves a vtkRectilinearGrid into a VTK XML file ('*.vtr')

    Parameters
    ----------
    grid : vtkRectilinearGrid
        vtk object with grid (rectilinear) grid
    path : str
        Extension (*.vtr) will be added if not provided

    """

    assert grid.GetClassName()=='vtkRectilinearGrid', 'error input vtk object is of type {}, a vtkRectilinearGrid was expected'.format(grid.GetClassName())

    # add extension to path
    if not path.lower().endswith('.vtr'):
        path = path + '.vtr'

    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetFileName(path)
    writer.SetDataModeToBinary()
    writer.SetInputData(grid)
    writer.Write()

cpdef SaveUnstructuredGrid(object grid, str path):
    """SaveUnstructuredGrid(object grid, str path)

    Saves a vtkUnstructuredGrid into a VTK XML file ('*.vtu')

    Parameters
    ----------
    grid : vtkUnstructuredGrid
        vtk object with grid (rectilinear) grid
    path : str
        Extension (*.vtu) will be added if not provided

    """
    assert grid.GetClassName()=='vtkUnstructuredGrid', 'error input vtk object is of type {}, a vtkUnstructuredGrid was expected'.format(grid.GetClassName())

    # add extension to path
    if not path.lower().endswith('.vtu'):
        path = path + '.vtu'

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(path)
    writer.SetDataModeToBinary()
    writer.SetInputData(grid)
    writer.Write()


cpdef PolyData2dxf(object mesh, str path):
    """PolyData2dxf(object mesh, str path)

    Saves a mesh vtkPolydataData into a dxf file ('*.dxf')
    Only vtkTriangle cells will be exported.

    Parameters
    ----------
    mesh : vtkPolyData
        vtk object with vtkTriangle cells
    path : str
        File name and path. Extension (*.dxf) will be added if not provided

    Note
    -----
    Requires the python modele ezdxf

    """
    try:
      import ezdxf
    except:
      print("This function requires the python module ezdxf")
      raise

    # add extension to path
    if not path.lower().endswith('.dxf'):
        path = path + '.dxf'

    dxfdw = ezdxf.new()
    dxfmd = dxfdw.modelspace()

    for i in range(mesh.GetNumberOfCells()):
        #check the type is 5 (a tringle)
        if mesh.GetCellType(i)==5:
            dxfmd.add_3dface((mesh.GetCell(i).GetPoints().GetPoint(0),
                              mesh.GetCell(i).GetPoints().GetPoint(1),
                              mesh.GetCell(i).GetPoints().GetPoint(2)))

    dxfdw.saveas(path)

cpdef PolyData2dxf_Quad(object mesh, str path):
    """PolyData2dxf_Quad(object mesh, str path)

    Saves a mesh vtkPolydataData into a dxf file ('*.dxf')
    Only vtkQuad cells will be exported.

    Parameters
    ----------
    mesh : vtkPolyData
        vtk object with vtkQuad cells
    path : str
        File name and path. Extension (*.dxf) will be added if not provided

    Note
    -----
    Requires the python modele ezdxf

    """
    try:
      import ezdxf
    except:
      print("This function requires the python module ezdxf")
      raise

    # add extension to path
    if not path.lower().endswith('.dxf'):
        path = path + '.dxf'

    dxfdw = ezdxf.new()
    dxfmd = dxfdw.modelspace()

    for i in range(mesh.GetNumberOfCells()):
        #check the type is 5 (a tringle)
        if mesh.GetCellType(i)==9:
            dxfmd.add_3dface((mesh.GetCell(i).GetPoints().GetPoint(0),
                              mesh.GetCell(i).GetPoints().GetPoint(1),
                              mesh.GetCell(i).GetPoints().GetPoint(2),
                              mesh.GetCell(i).GetPoints().GetPoint(3)))

    dxfdw.saveas(path)

cpdef dxf2PolyData(str path):
    """dxf2PolyData(path)

    Reads a mesh in dxf format ('*.dxf') and returns a Polydata Mesh
    Only 3DFACE from dxf will be imported.

    Parameters
    ----------
    path : str
        File name and path of dxf file.

    Returns
    -------
    vtkPolyData with vtkTriangle cells

    Note
    -----
    Requires the python modele ezdxf

    """

    # read the dxf
    dwg = ezdxf.readfile(path)

    # get the dxf model space
    msp = dwg.modelspace()

    # extract all the triangles
    triangles = []
    for e in msp.query('3DFACE'):
        # get point coordinates of the triangle (one per row)
        triangles.append([e.dxf.vtx0, e.dxf.vtx1, e.dxf.vtx2])

    # Put trisngles in vtkPolydata
    meshpoints = vtk.vtkPoints()
    meshpoints.SetDataTypeToDouble()
    meshtriangle = vtk.vtkTriangle()
    meshtriangles = vtk.vtkCellArray()
    id = 0
    for i in triangles:
        meshpoints.InsertNextPoint(i[0])
        meshpoints.InsertNextPoint(i[1])
        meshpoints.InsertNextPoint(i[2])
        meshtriangle = vtk.vtkTriangle()
        meshtriangle.GetPointIds().InsertId(0,id)
        meshtriangle.GetPointIds().InsertId(1,id+1)
        meshtriangle.GetPointIds().InsertId(2,id+2)
        meshtriangles.InsertNextCell(meshtriangle)
        id = id+3

    mesh = vtk.vtkPolyData()
    mesh.SetPoints(meshpoints)
    mesh.SetPolys(meshtriangles)

    # doing a bit of cleanup (remove duplicated points)
    cleanPolyData = vtk.vtkCleanPolyData()
    cleanPolyData.SetInputData(mesh)
    cleanPolyData.Update()

    return cleanPolyData.GetOutput()

cpdef dxf2PolyData_line3d(str path, bint clean = False):
    """dxf2PolyData_line3d(path, clean = False)

    Reads 3D lines (POLYLINE) in dxf format ('*.dxf') and returns a Polydata.

    Parameters
    ----------
    path : str
        File name and path of dxf file.
    clean : boolean, default False
        Clean duplicted points

    Returns
    -------
    vtkPolyData, points : a vtkPolyData, and list of point coordinates tuples

    Note
    -----
    Requires the python modele ezdxf

    """

    # read the dxf
    dwg = ezdxf.readfile(path)

    # get the dxf model space
    msp = dwg.modelspace()

    # extract lines and points and save it in a dict
    lines = msp.query('POLYLINE')  # implement also LWPOLYLINE
    lins = {}
    lid = 0
    pid = 0
    pnts = []
    for l in range(len(lines)):
        lins[l] = []
        for v in lines[l]:
            pnts.append(v.dxf.location)
            lins[l].append(pid)
            pid = pid + 1

    # create poinst
    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    for p in pnts:
        points.InsertNextPoint(p)

    # create lines and cell array
    cells = vtk.vtkCellArray()
    for l in range(len(lines)):
        polyLine = vtk.vtkPolyLine()
        polyLine.GetPointIds().SetNumberOfIds(len(lins[l]))
        for i in range(len(lins[l])):
            polyLine.GetPointIds().SetId(i, lins[l][i])
        cells.InsertNextCell(polyLine)

    # create polydata
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(cells)

    if clean:
        # doing a bit of cleanup (remove duplicated points)
        cleanPolyData = vtk.vtkCleanPolyData()
        cleanPolyData.SetInputData(polyData)
        cleanPolyData.Update()
        return cleanPolyData.GetOutput(), pnts
    else:
        return polyData, pnts

cpdef dxf2PolyData_lwpoly(str path, bint clean = False):
    """dxf2PolyData_lwpoly(path, clean = False)

    Reads LWPOLYLINE in dxf format ('*.dxf') and returns a Polydata.

    Parameters
    ----------
    path : str
        File name and path of dxf file.
    clean : boolean
        Clean duplicted points

    Returns
    -------
    vtkPolyData, points : a vtkPolyData, and list of point coordinates tuples

    Note
    -----
    Requires the python modele ezdxf

    """

    # read the dxf
    dwg = ezdxf.readfile(path)

    # get the dxf model space
    msp = dwg.modelspace()

    # extract lines and points and save it in a dict
    lines = msp.query('LWPOLYLINE')  # implement also LWPOLYLINE
    lins = {}
    lid = 0
    pid = 0
    pnts = []
    for l in range(len(lines)):
        zi = lines[l].dxf.elevation
        lins[l] = []
        for li in lines[l]:
            pnts.append([float(li[0]),float(li[1]),float(zi)])
            lins[l].append(pid)
            pid = pid + 1

    # create poinst
    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    for p in pnts:
        points.InsertNextPoint(p)

    # create lines and cell array
    cells = vtk.vtkCellArray()
    for l in range(len(lines)):
        polyLine = vtk.vtkPolyLine()
        polyLine.GetPointIds().SetNumberOfIds(len(lins[l]))
        for i in range(len(lins[l])):
            polyLine.GetPointIds().SetId(i, lins[l][i])
        cells.InsertNextCell(polyLine)

    # create polydata
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(cells)

    if clean:
        # doing a bit of cleanup (remove duplicated points)
        cleanPolyData = vtk.vtkCleanPolyData()
        cleanPolyData.SetInputData(polyData)
        cleanPolyData.Update()
        return cleanPolyData.GetOutput(), pnts
    else:
        return polyData, pnts

cpdef dxf2PolyData_Quad(str path):
    """dxf2PolyData_Quad(path)

    Reads a mesh in dxf format ('*.dxf') and returns a Polydata Mesh
    Only 3DFACE of 4 side quadrilaters will be imported.

    Parameters
    ----------
    path : str
        File name and path of dxf file.

    Returns
    -------
    vtkPolyData with vtkQuad cells

    Note
    -----
    Requires the python modele ezdxf

    """

    # read the dxf
    dwg = ezdxf.readfile(path)

    # get the dxf model space
    msp = dwg.modelspace()

    # extract all the triangles
    quads = []
    for e in msp.query('3DFACE'):
        # get point coordinates of the triangle (one per row)
        quads.append([e.dxf.vtx0, e.dxf.vtx1, e.dxf.vtx2,e.dxf.vtx3])

    # Put quads in vtkPolydata
    meshpoints = vtk.vtkPoints()
    meshpoints.SetDataTypeToDouble()
    meshquad = vtk.vtkQuad()
    meshquads = vtk.vtkCellArray()
    id = 0
    for i in quads:
        meshpoints.InsertNextPoint(i[0])
        meshpoints.InsertNextPoint(i[1])
        meshpoints.InsertNextPoint(i[2])
        meshpoints.InsertNextPoint(i[3])
        meshquad= vtk.vtkQuad()
        meshquad.GetPointIds().InsertId(0,id)
        meshquad.GetPointIds().InsertId(1,id+1)
        meshquad.GetPointIds().InsertId(2,id+2)
        meshquad.GetPointIds().InsertId(3,id+3)
        meshquads.InsertNextCell(meshquad)
        id = id+4

    mesh = vtk.vtkPolyData()
    mesh.SetPoints(meshpoints)
    mesh.SetPolys(meshquads)

    # doing a bit of cleanup (remove duplicated points)
    cleanPolyData = vtk.vtkCleanPolyData()
    cleanPolyData.SetInputData(mesh)
    cleanPolyData.Update()

    return cleanPolyData.GetOutput()


def surpac_TRI2polydata(triangles):
    """surpac_TRI2polydata(triangles)

    Creates vtk polydata from triangled with coordinates:

        ['x1','y1','z1','x2','y2','z2','x3','y3','z3']

    Parameters
    ----------
    triangles : numpy 2d array
        array of triangles with shape [n,9]

    Returns
    -------
    vtkPolyData with vtkTriangle cells

    """

    assert triangles.shape[1]==9

    n = triangles.shape[0]

    # Put trisngles in vtkPolydata
    meshpoints = vtk.vtkPoints()
    meshpoints.SetDataTypeToDouble()
    meshtriangle = vtk.vtkTriangle()
    meshtriangles = vtk.vtkCellArray()

    id = 0
    for i in range(n):
        meshpoints.InsertNextPoint(triangles[i,0:3])
        meshpoints.InsertNextPoint(triangles[i,3:6])
        meshpoints.InsertNextPoint(triangles[i,6:9])
        meshtriangle = vtk.vtkTriangle()
        meshtriangle.GetPointIds().InsertId(0,id)
        meshtriangle.GetPointIds().InsertId(1,id+1)
        meshtriangle.GetPointIds().InsertId(2,id+2)
        meshtriangles.InsertNextCell(meshtriangle)
        id = id+3

    mesh = vtk.vtkPolyData()
    mesh.SetPoints(meshpoints)
    mesh.SetPolys(meshtriangles)

    # doing a bit of cleanup (remove duplicated points)
    cleanPolyData = vtk.vtkCleanPolyData()
    cleanPolyData.SetInputData(mesh)
    cleanPolyData.Update()

    return cleanPolyData.GetOutput()

def line2polydata(line):
    """line2polydata(line)

    Creates vtk polydata from lines with coordinates:

        ['x1','x2','...','xn']

    Parameters
    ----------
    line : numpy 2d array
        array of points with shape [n,3]

    Returns
    -------
    vtkPolyData with vtkLine cells

    """

    points= vtk.vtkPoints()
    points.SetDataTypeToDouble()
    vline = vtk.vtkLine()
    lines = vtk.vtkCellArray()

    for l in range(line.shape[0]):
        points.InsertNextPoint(line[l][0], line[l][1], line[l][2])
        vline.GetPointIds().InsertId(l,l)

    lines.InsertNextCell(vline)

    linesPolyData = vtk.vtkPolyData()
    linesPolyData.SetPoints(points)
    linesPolyData.SetLines(lines)

    return linesPolyData

cpdef decimate_polyline(object polyline, float target_reduction = 0.9, float maximum_error=0.01):
    """decimate_polyline(polyline, target_reduction = 0.9, maximum_error=0.01)

    Reduce the number o pints in polylines

    Parameters
    ----------
    polyline : vtkPolyData
        polylines
    
    target_reduction : float, default  0.9
        Desired reduction in the total number of polygons (e.g., if TargetReduction is set 
        to 0.9, this filter will try to reduce the data set to 10% of its original size).
    maximum_error: float, default 0.01
        largest decimation error that is allowed during the decimation process.
        This may limit the maximum reduction that may be achieved. The maximum 
        error is specified as a fraction of the maximum length of the input data bounding box.
    

    Returns
    -------
    vvtkPolyData with decimated polylines

    """

    # decimate 
    decimate = vtk.vtkDecimatePolylineFilter()
    decimate.SetInputData(polyline)
    decimate.SetTargetReduction(target_reduction)
    decimate.SetMaximumError(maximum_error)
    decimate.Update()
    
    return decimate.GetOutput()



# ----------------------------------------------------------------------
#   Functions for imlicit modeling
# ----------------------------------------------------------------------
cpdef calculate_normals(object mesh, bint flip_normals=False):
    """calculate_normals(mesh, flip_normals=False)

    Takes a vtkPolydata Mesh and calculate its normal vectors.

    Implicit functions use vtkPolydata normals to define the sign of
    the implicit distances. The results may be spurios if normals are not
    calculated or inconsistent (for example pointing downward in a topografic
    surface). In this case you may use this fuction to calculate normals.

    Parameters
    ----------
    mesh : vtkPolyData
        mesh surface (not tested with polylines...)

    Returns
    -------
    vvtkPolyData with normals

    """
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(mesh)
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.ComputePointNormalsOn()
    if flip_normals:
      normals.FlipNormalsOn()
    normals.Update()
    return normals.GetOutput()


cpdef implicit_surface(object mesh, bint update_normals=True):
    """implicit_surface(mesh)

    Takes a vtkPolydata Mesh (a surface or a solid) and generates an implicit
    vtk function that can be used to evaluate the signed distance d to the
    nearest point in the mesh.
    The sign of the distance d can be used to evaluate if a point x is
    inside - outside a solid, or above - below a surface, or in the surface.
    The distance d will be zero for points in a surface, negative for interior
    points and positive for exterior points.

    The vtkImplicitPolyDataDistance objects generated by this functions can
    be used as imput of vtkImplicitBoolean.

    Parameters
    ----------
    mesh : vtkPolyData
        mesh surface (not tested with polylines...)
    update_normals: boolean (defaul True)
        if true normal vectors will be recalculated

    Returns
    -------
    vtkImplicitPolyDataDistance (this is an implicit function that can be evaluated)

    See Also
    ------
    evaluate_implicit_function, boolean_implicit

    """
    # generate instance of tkImplicitPolyDataDistance and set input
    if update_normals:
      mesh = calculate_normals(mesh)
    implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
    implicitPolyDataDistance.SetInput(mesh)

    return implicitPolyDataDistance


cpdef evaluate_implicit_points( object implicit_mesh,
                                np.ndarray [double, ndim=1] x,
                                np.ndarray [double, ndim=1] y,
                                np.ndarray [double, ndim=1] z,
                                cap_dist=None,
                                bint normalize=False):
    """evaluate_implicit_points(implicit_mesh, x, y, z, cap_dist=None, normalize=False)

    Calculates the signed distances d between input points,
    with corrdinates x,y,z, and the nearest point in a vtkImplicitPolyDataDistance
    object (usually generated with the function implicit_surface)
    Optionally, the distance d will be truncated to a treshold `cap_dist` and
    normalized to a value between [-1 and 1] if `normalize==True`

    The sign of the distance d can be used to evaluate if a point x is
    inside - outside a solid, or above - below a surface, or in the surface.
    The distance d will be zero for points in a surface, negative for interior
    points and positive for exterior points.

    Parameters
    ----------
    implicit_mesh : vtkImplicitPolyDataDistance
        implicit mesh surface (not tested with polylines...)
    x, y, z: array of floats
        coordinates of the points
    cap_dist : unsigned float expected (default None)
        thresuld to truncate istance d to cap_dist if d is positive or to
        -cap_dist if d is negative. No truncation will be applied if cap_dist is
        None
    normalize: boolean (default False)
        If normalize is True and cap_dist is not None the distance will be
        normalize with a vaue in the interval [-1, 1] by applying the
        fuction d = d/cap_dist

    Returns
    -------
    float with distances d

    See Also
    ------
    evaluate_implicit_grid, implicit_surface

    """
    assert x.shape[0]==y.shape[0]==z.shape[0]
    d = np.empty([x.shape[0]])
    # loop on each point # TODO: optimize sending vtk array with coordinates: see virtual void 	EvaluateFunction (vtkDataArray *input, vtkDataArray *output)
    for i in range(x.shape[0]):
      d[i] = implicit_mesh.FunctionValue(x[i],y[i],z[i])

    if cap_dist is not None:
      cap_dist = abs(cap_dist)
      d[d>cap_dist] = cap_dist
      d[d<-cap_dist] = -cap_dist
    if cap_dist is not None and normalize:
      d = d/cap_dist
    if cap_dist is None and normalize:
      warnings.warn("Normalization was not applied because cap_dist is None")

    return d


cpdef define_region_grid( float xorg, float yorg, float zorg,
                          float dx, float dy, float dz,
                          int nx, int ny, int nz,
                          snapping_points = None,
                          float tol= 0.01):
    """define_region_grid(xorg, yorg, zorg, dx, dy,  dz, nx, ny, nz, snapping_points = None, tol= 0.01)

    The region grid is an irregular grid of tetras (a vtk.vtkUnstructuredGrid)
    within the limits of the working area or project, generated with vtkDelaunay3D
    of regularly spaced points. Alternativelly you can append points, defined as
    vtk geometries, to ensure snapping of implicit functions.

    The region grid is required as input to generate closed wireframes from open
    and closed surfaces. Note that a grid define corner points of a standar
    block model, if your block model has 2 block along direction x the grid has
    nx+1 == 3 points. This fuction internally redefine nx=nx+1 to match block
    model definitions. However, you may define region grid with at least one
    block wider because boundary points are set to an arbitrary negative float
    number to ensure that surfaces are closed.


    Parameters
    ----------
    xorg, yorg, zorg, dx, dy,  dz: floats
        coordinates of the lower left point in the grid and grid spacing
    nx, ny, nz: int
        number of rows, cols and levels in the grid
    snapping_points: array like of vtk geometries
        each element may contain snapping points, for example unit surfaces
        or contact points

    Returns
    -------
    (array, vtkUnstructuredGrid) array of floats and vtk grid with distances `d`

    See Also
    ------
    evaluate_implicit_points, implicit_surface, blockmodel.Blockmodel.fillwireframe

    """

    nx=  nx + 1
    ny = ny + 1
    nz = nz + 1

    ufgrid = vtk.vtkImageData()
    ufgrid.SetOrigin(xorg,yorg,zorg)
    ufgrid.SetSpacing(dx,dy,dz)
    ufgrid.SetDimensions(nx,ny,nz)

    if snapping_points is not None:
        appfil = vtk.vtkAppendFilter()
        appfil.AddInputData(ufgrid)
        for i in snapping_points:
          appfil.AddInputData(i)
        appfil.Update()

        ufgrid2=appfil.GetOutput()

        delny = vtk.vtkDelaunay3D()
        delny.SetInputData(ufgrid2)
        delny.SetTolerance(tol)
        delny.SetAlpha(max(dx*2,dy*2,dz*2))
        delny.BoundingTriangulationOff()
        delny.Update()

        result = vtk.vtkUnstructuredGrid()
        result.DeepCopy(delny.GetOutput())

        return result

    else:
        ufgrid.AllocateScalars(vtk.VTK_INT, 1)
        ufgrid.GetPointData().GetAttribute(0).Fill(0)
        thresh = vtk.vtkThreshold()
        thresh.SetInputData(ufgrid)
        thresh.ThresholdByUpper(-1)
        thresh.Update()
        return thresh.GetOutput()



cpdef evaluate_region(object region,
                      object implicit_func,
                      str func_name,
                      bint invert=False,
                      float capt = -1000000,
                      bint closed=True):
    """evaluate_region(region, implicit_func, func_name, invert=False, capt = -1000000, closed = True)

    Evaluate an implicit function within a region.

    The results is a distance scalar with name func_name. The sign of distances
    are defined by geometry normals. Positive distances are within the object,
    the sign can be inverted by setting invert=True. Cells in the outer surface
    are set to capt to ensure closed surfaces are retrieved from contour
    functions.


    Parameters
    ----------
    region : vtkUnstructuredGrid
        region to evaluate
    implicit_func : vtkImplicitFuntion or vtkImplicitPolyDataDistance
        implicit function to evaluate
    func_name: str
        name of the scalar asigned to points in the region
    invert: boolean (default False)
        if True the sign of the distance to the implicit function will be inverted
    cap: float (default -1000000)
        any distance over this value will be set equal to this value.
        points in the outer surface of the region will be set to this value

    Returns
    -------
    (vtkUnstructuredGrid, ndarray)
        vtk object with the region and numpy array with distances

    See Also
    ------
    evaluate_implicit_points, implicit_surface, define_region_grid

    """

    xmin,xmax, ymin,ymax, zmin,zmax = region.GetBounds()
    n = region.GetNumberOfPoints()

    # make sure we get the sign right
    capt = -np.abs(capt)

    # generate points and scalar
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    d = np.zeros(n)

    for i in range(n):
      x[i] ,y[i] , z[i] = region.GetPoint(i)
      d[i] = implicit_func.FunctionValue(x[i],y[i],z[i])
      if invert:
        d[i] = -d[i]

    # make outer surface equal to capt
    if capt is not None:
      d[d<capt]= capt
      d[d>-capt]= -capt
    if closed:
      d[x<=xmin] = capt
      d[y<=ymin] = capt
      d[z<=zmin] = capt
      d[x>=xmax] = capt
      d[y>=ymax] = capt
      d[z>=zmax] = capt

    # pudate implisit function value
    tmp = vtk.vtkUnstructuredGrid()
    tmp.DeepCopy(region)
    tmp.GetPointData().AddArray(vtkdsa.numpyTovtkDataArray(d, name=func_name, array_type=None))

    return tmp, d

cpdef set_region_field(object region, d, str func_name):
        """set_region_field(region, d, func_name)

        Set region value.

        This function is usefull to asign values to a region (a vtkUnstructuredGrid).

        You can uses it to define externally defined implisit functions. For example,
        Given n implicit distances d1, d2, ..., dn the intersection is
        d_intersect = min (d1, d2, ..., dn)

        Parameters
        ----------
        region : vtkUnstructuredGrid
            region to evaluate
        d : expects a numpy array
            field value
        func_name: str
            name asigned to d in the vtk object

        Returns
        -------
        vtk Object

        """

        assert (region.GetNumberOfPoints()==d.shape[0])

        tmp = vtk.vtkUnstructuredGrid()
        tmp.DeepCopy(region)
        tmp.GetPointData().AddArray(vtkdsa.numpyTovtkDataArray(d, name=func_name, array_type=None))

        return tmp

cpdef extract_surface(object region, str func_name, float threshold = 0.0):
        """extract_surface(object region, str func_name, float threshold = 0.0)

        Extract surface defined by an implisit function evaluated in the region

        Parameters
        ----------
        region : vtkUnstructuredGrid
            region with implicit function evaluated
        func_name: str
            name of the implicit function
        threshold: float (default 0.0)
            value to contour the implisit surface

        Returns
        -------
        vtkPolydata with surface

        """

        surface = vtk.vtkContourGrid()
        surface.SetInputData(region)
        surface.SetInputArrayToProcess(0, 0, 0, vtk.vtkAssignAttribute.POINT_DATA, func_name)
        surface.SetValue(0,threshold)
        surface.Update()

        return surface.GetOutput()

cpdef evaluate_implicit_grid( object implicit_mesh,
                              float xorg, float yorg, float zorg,
                              float dx, float dy, float dz,
                              int nx, int ny, int nz,
                              x=None,y=None,z=None,
                              cap_dist=None,
                              bint normalize=False,
                              bint closed = False,
                              str name = ''):
    """evaluate_implicit_grid(mplicit_mesh, xorg, yorg, zorg, dx, dy,  dz, nx, ny, nz,x = None,y = None,z = None, cap_dist=None, normalize=False, closed = False, name = '')

    Calculates the signed distances `d` between points in a regular 3Dgrid and
    the nearest point in a vtkImplicitPolyDataDistance
    object (usually generated with the function implicit_surface)
    Optionally, the distance d will be truncated to a treshold `cap_dist` and
    normalized to a value between [-1 and 1] if `normalize==True`. The
    outer points of the grid will be set to cap_dist if `closed` is True,
    that ensure that any implicit surface generated with the grid is closed.
    The points in the grid are calculated internally if x,y, or z are None.
    Note that a grid define corner points of a standar block model, if your
    block mode has 2 block along direction x the grid has nx+1 == 3 points.

    The sign of the distance d can be used to evaluate if a point x is
    inside - outside a solid, or above - below a surface, or in the surface.
    Evaluating the four points defining a block may allow identify blocks
    cut by a surface. This may fail for closed surfaces thinner than the
    block sizes, such as narrow gold veins, this issue can be resolved using
    two separated sufaces to define a vein.

    Parameters
    ----------
    implicit_mesh : vtkImplicitPolyDataDistance
        implicit mesh surface (not tested with polylines...)
    xorg, yorg, zorg, dx, dy,  dz: floats
        coordinates of the lower left point in the grid and grid spacing
    nx, ny, nz: int
        number of rows, cols and levels in the grid
    x, y, z: numpy arrays or None
        (optional) coordinates of grid poins, if None coordinates are calculated
    cap_dist : unsigned float expected (default None)
        thresuld to truncate istance d to cap_dist if d is positive or to
        -cap_dist if d is negative. No truncation will be applied if cap_dist is
        None
    normalize: boolean (default False)
        If normalize is True and cap_dist is not None the distance will be
        normalize with a vaue in the interval [-1, 1] by applying the
        fuction d = d/cap_dist
    closed: boolean (defaul False)
        If closed is True and cap_dist is not None then points in the outer
        wall will be set to cap_dist. The sign of cap_dist will depermine
        if the surface closes above or below the surface, or inside vs outside
        the solid.
    name: string (default '')
        name subfix to be applied to the distance d array in vtk file,
        for example, if name is 'topo', the output array will be 'dist_topo'

    Returns
    -------
    (array, vtkImageData) array of floats and vtk grid with distances `d`

    See Also
    ------
    evaluate_implicit_points, implicit_surface, blockmodel.Blockmodel.fillwireframe

    """
    cdef int ijk
    assert (x is None)==(y is None)==(z is None)

    nx=  nx + 1
    ny = ny + 1
    nz = nz + 1

    mask = np.empty(nx*ny*nz, dtype=bool)
    mask[:] = False

    if x is not None:
      assert x.shape[0]==y.shape[0]==z.shape[0]
      assert x.shape[0] == nx*ny*nz
      ijk = -1
      for i in range(nz):
        for j in range(ny):
          for k in range(nx):
            ijk = ijk+1
            if i==0 or j==0 or k == 0 or i==nz-1 or j==ny-1 or k == nx-1:
              mask[ijk] = True

    else:
      # calculate x,y,z coordinates
      x = np.empty(nx*ny*nz, dtype=float)
      y = np.empty(nx*ny*nz, dtype=float)
      z = np.empty(nx*ny*nz, dtype=float)

      ijk = -1
      for i in range(nz):
        for j in range(ny):
          for k in range(nx):
            ijk = ijk+1
            if i==0 or j==0 or k == 0 or i==nz-1 or j==ny-1 or k == nx-1:
              mask[ijk] = True

            x[ijk]= k*dx+xorg
            y[ijk]= j*dy+yorg
            z[ijk]= i*dz+zorg

    # apply evaluation
    d=evaluate_implicit_points(implicit_mesh, x, y, z, cap_dist, normalize)

    if closed and cap_dist is not None:
      d[mask] = -cap_dist

    ufgrid = vtk.vtkImageData()
    ufgrid.SetOrigin(xorg,yorg,zorg)
    ufgrid.SetSpacing(dx,dy,dz)
    ufgrid.SetDimensions(nx,ny,nz)
    ufgrid.GetPointData().AddArray(vtkdsa.numpyTovtkDataArray(d, name='dist_{}'.format(name), array_type=None))


    return d, ufgrid


cpdef clip_with_surface( object region, object implicit_surface, str how= 'inside'):
    """clip_with_surface(region, surface, how= 'inside')

    Clip a region defined as a regular grid (a vtkImageData) or an irregular
    grid (a vtkUnstructuredGrid) with a surface. The clipped mode
    `inside` or `outside` the implicit_surface is defined by surface normals.
    For example, clipping a region with a topografic surface with normals
    pointing upward, using `how == inside` will remove the air (above topo)
    portion of the region.

    You can model individual lithology units defined by a set of open surfaces
    by using as input the output vtkUnstructuredGrid resulting from previous
    clipping operations.


    Parameters
    ----------
    region: a vtkImageData or vtkUnstructuredGrid
        3D vtk object defining a region to be clipped
    implicit_surface: a vtkImplicitPolyDataDistance
        implicit surface with consitent normals. Can be obtained with the
        function implicit_surface()

    Returns
    -------
    (vtkUnstructuredGrid, vtkPolydata) 3D object and outer surface of the region clipped

    See Also
    ------
    evaluate_implicit_points, implicit_surface, blockmodel.Blockmodel.fillwireframe

    """
    # clip region
    clip = vtk.vtkClipDataSet()
    clip.SetInputData(region)
    clip.SetClipFunction(implicit_surface)
    clip.InsideOutOn()
    if how=='outside':
      clip.InsideOutOff()
    clip.Update()
    clip_region = clip.GetOutput()

    # get outer surface
    # extract surface
    gfilter=vtk.vtkGeometryFilter()
    gfilter.SetInputData(clip_region)
    gfilter.Update()

    # Triangle filter to remove quads and ensure only triangles are retrived
    cutTriangles = vtk.vtkTriangleFilter()
    cutTriangles.SetInputData(gfilter.GetOutput())
    cutTriangles.Update()

    # Calc Normals
    clip_surf=calculate_normals(cutTriangles.GetOutput())

    return clip_region,clip_surf

# ----------------------------------------------------------------------
#   Functions for polygon operations (inside/outside test, area, etc.)
# ----------------------------------------------------------------------
cpdef inside_polygon(object points, object polygon_pts = None, object polygon = None, object normal = None, object bounds = None):
    """inside_polygon(object points, object polygon_pts = None, object polygon = None, object normal = None, object bounds = None)

    Test is foints are inside a polygon.

    The polygon can be passed as array of x,y,z coordinates or as vtkpolygon
    object. The polygons cannot have any internal holes, and cannot
    self-intersect. Define the polygon with n-points ordered in the
    counter-clockwise direction; do not repeat the last point.

    The function will calculate the normal and bounding box of the polygon.
    However, you can pass it by reference to speeds repated calculations with
    the same polygon. The normal is used to define the direction of the
    computation in 3D.


    Parameters
    ----------
    points: array of floats[3]
        3D points to be tested for inclusion and exclusion inside the polygon
    polygon_pts: array of floats[3] defining the polygon (optional)
        if polygon_pts is None, then you may provide polygon
    polygon: vtkPolygon
        if polygon is None, then you may provide polygon_pts
    normal: floats[3]
        normal of the polygon. If none it will be calculated internally
    bounds: floats[6]
        bounding coordinates of the polygon. If none it will be calculated internally

    Returns
    -------
    (array of int, vtkPolygon, floats[3], floats[6]) these are the:
      results with 1 if the point is inside or 0 if outside,
      polygon as a vtkPolygon object,
      normal floats[3] with polygon normal,
      and bounds floats[6]

    you can use polygon, normal, bounds as input parameters in subsequent call
    to this function.

    See Also
    ------
    polygon_area

    """
    result = np.zeros(len(points), dtype = int)

    if polygon is None:
        assert polygon_pts is not None

        polygon = vtk.vtkPolygon()

        for i in polygon_pts:
            polygon.GetPoints().InsertNextPoint(i)

    # compute normal
    if normal is None:
        normal=[0.0,0.0,0.0]
        polygon.ComputeNormal(polygon.GetPoints(), normal)
    # compute bounds
    if bounds is None:
        bounds = [0.,0.,0.,0.,0.,0.];
        polygon.GetPoints().GetBounds(bounds)

    for i in range(len(points)):
        result[i] = polygon.PointInPolygon(points[i],
                       polygon.GetPoints().GetNumberOfPoints(),
                       vtknumpy.vtk_to_numpy(polygon.GetPoints().GetData()).ravel(),
                       bounds,
                       normal)

    return result, polygon, normal, bounds


def polygon_area(object polygon_pts = None, object polygon = None):
    """polygon_area(object polygon_pts = None, object polygon = None)

    Calculates the area of a polygon

    The polygon can be passed as array of x,y,z coordinates or as vtkpolygon
    object. The polygons cannot have any internal holes, and cannot
    self-intersect. Define the polygon with n-points ordered in the
    counter-clockwise direction; do not repeat the last point.

    Parameters
    ----------
    polygon_pts: array of floats[3] defining the polygon (optional)
        if polygon_pts is None, then you may provide polygon
    polygon: vtkPolygon
        if polygon is None, then you may provide polygon_pts

    Returns
    -------
    floats
      this is the planar area of the polygon (perpendicular to its normal (?))

    Note
    ----
    the polygon normal is calculated internally

    See Also
    ------
    inside_polygon

    """

    if polygon is None:
        assert polygon_pts is not None

        polygon = vtk.vtkPolygon()

        for i in polygon_pts:
            polygon.GetPoints().InsertNextPoint(i)

    # compute normal
    normal=[0.0,0.0,0.0]
    polygon.ComputeNormal(polygon.GetPoints(), normal)

    # compute area
    return polygon.ComputeArea(polygon.GetPoints(),
                        polygon.GetPoints().GetNumberOfPoints(),
                        list(range(polygon.GetPoints().GetNumberOfPoints())),
                        normal)
