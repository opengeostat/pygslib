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



cpdef vtk_raycasting(object surface, object pSource, object pTarget):
    """vtk_raycasting(object surface, object pSource, object pTarget)

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


    See Also
    --------
    vtk_show, loadSTL

    Note
    ----
    This Code was modified from the original code by Adamos Kyriakou
    published in https://pyscience.wordpress.com/

    """
    cdef int intersect, idx

    obbTree = vtk.vtkOBBTree()
    obbTree.SetDataSet(surface)
    obbTree.BuildLocator()
    pointsVTKintersection = vtk.vtkPoints()

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

    return intersect, pointsIntersection, pointsVTKIntersectionData



# ----------------------------------------------------------------------
#   Functions for point querying
# ----------------------------------------------------------------------
cpdef pointquering(object surface,
                   double azm,
                   double dip,
                   np.ndarray [double, ndim=1] x,
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   int test):
    """pointquering(object surface, double azm, double dip, np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, int test)

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
    obbTree = vtk.vtkOBBTree()
    obbTree.SetDataSet(surface)
    obbTree.BuildLocator()
    pointsVTKintersection = vtk.vtkPoints()

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



# ----------------------------------------------------------------------
#   Functions for point querying
# ----------------------------------------------------------------------
cpdef pointinsolid(object surface,
                   np.ndarray [double, ndim=1] x,
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   double tolerance = .000001 ):
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
#   Inport triangles into wireframe
# ----------------------------------------------------------------------
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

    points = vtk.vtkPoints();


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
    p.SetNumberOfPoints(npoints)


    for i in range(npoints):
        p.InsertPoint(i,(points[i,0],points[i,1],points[i,2]))

    polydata.SetPoints(p)


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

    # add extension to path
    if not path.lower().endswith('.vtp'):
        path = path + '.vtp'

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(polydata)
    writer.Write()

cpdef SaveImageData(object image, str path):
    """SaveImageData(object image, str path)

    Saves a vtkImageData into a VTK XML image file ('*.vti')

    Parameters
    ----------
    image : vtkImageData
        vtk object with image (regular) grid
    path : str
        Extension (*.vti) will be added if not provided


    """

    # add extension to path
    if not path.lower().endswith('.vti'):
        path = path + '.vti'

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(path)
    writer.SetDataModeToBinary()
    writer.SetInputData(image)
    writer.Write()


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


cpdef grid2vtkImageData(
               int nx, int ny, int nz,
               double xorg, double yorg, double zorg,
               double dx, double dy, double dz,
               object data):
    """grid2vtkImageData(int nx, int ny, int nz, double xorg, double yorg, double zorg, double dx, double dy, double dz, np.ndarray [double, ndim=1] data)

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
    data : dictionary
        dictionary with {'dataname':datavalues}

    See also
    --------
    blockmodel.Blockmodel.blocks2vtkImageData
    """

    ufgrid = vtk.vtkImageData()
    ufgrid.SetOrigin(xorg+dx/2,yorg+dy/2,zorg+dz/2)
    ufgrid.SetSpacing(dx,dy,dz)
    ufgrid.SetDimensions(nx+1,ny+1,nz+1)

    for i in data:
      ufgrid.GetCellData().AddArray(vtkdsa.numpyTovtkDataArray(data[i], name=i, array_type=None))

    return ufgrid


cpdef partialgrid2vtkfile(str path,
                   np.ndarray [double, ndim=1] x,
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   double DX,
                   double DY,
                   double DZ,
                   object var,
                   object varname):
    """partialgrid2vtkfile(str path, np.ndarray [double, ndim=1] x, np.ndarray [double, ndim=1] y, np.ndarray [double, ndim=1] z, double DX, double DY, double DZ, object var, object varname)

    Saves data in the cells of a VTK Unstructured Grid file


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


    # number of cells/blocks
    nc = x.shape[0]

    #number of variables
    nvar = len(var)

    # number of points (block vertex)
    np = nc*8

    # Create array of the points and ID
    pcoords = vtk.vtkFloatArray()
    pcoords.SetNumberOfComponents(3)
    pcoords.SetNumberOfTuples(np)

    points = vtk.vtkPoints()
    voxelArray = vtk.vtkCellArray()

    # create vertex (points)
    id=0
    #for each block
    for i in range (nc):
        # for each vertex
        pcoords.SetTuple3(id, x[i]+DX/2., y[i]-DY/2., z[i]-DZ/2.)
        pcoords.SetTuple3(id+1, x[i]-DX/2., y[i]-DY/2., z[i]-DZ/2.)
        pcoords.SetTuple3(id+2, x[i]+DX/2., y[i]+DY/2., z[i]-DZ/2.)
        pcoords.SetTuple3(id+3, x[i]-DX/2., y[i]+DY/2., z[i]-DZ/2.)
        pcoords.SetTuple3(id+4, x[i]+DX/2., y[i]-DY/2., z[i]+DZ/2.)
        pcoords.SetTuple3(id+5, x[i]-DX/2., y[i]-DY/2., z[i]+DZ/2.)
        pcoords.SetTuple3(id+6, x[i]+DX/2., y[i]+DY/2., z[i]+DZ/2.)
        pcoords.SetTuple3(id+7, x[i]-DX/2., y[i]+DY/2., z[i]+DZ/2.)

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
        cscalars = vtknumpy.numpy_to_vtk(var[i])
        cscalars.SetName(varname[i])
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


    # save results
    writer = vtk.vtkXMLUnstructuredGridWriter();
    writer.SetFileName(path);
    writer.SetInputData(extractGrid.GetOutput())
    writer.Write()


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
                   Alpha = 0.0,
                   Tolerance = 0.001,
                   constraints = None):
    """
    Creates a triangulated Surface


    Parameters
    ----------
    x,y,z : np.ndarray [double, ndim=1]
        Coordinates of the input points
    Alpha: double, default(0.0)
        For a non-zero alpha value, only edges or triangles contained
        within a sphere centered at mesh vertices will be output.
        Otherwise, only triangles will be output.
    Tolerance : double, default(0.001)
        Specify a tolerance to control discarding of closely spaced points.
        This tolerance is specified as a fraction of the diagonal length of
        the bounding box of the points.
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
