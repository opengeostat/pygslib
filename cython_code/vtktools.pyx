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

import vtk
from IPython.display import Image
from libc.math cimport sin
from libc.math cimport cos
cimport numpy as np
import numpy as np
import pyevtk.hl
import vtk.util.numpy_support as vtknumpy


cpdef vtk_show(renderer, width=400, height=300, 
             camera_position=None, 
             camera_focalpoint=None):
    """
    
    vtk_show(renderer, 
             width=400, 
             height=300,
             camera_position=None,
             camera_focalpoint=None)
    
    Display a vtk renderer in Ipython Image
    
    Parameters
    ----------
    renderer : VTK renderer
           renderer with vtk objects and properties 
    width, height: float
           Size of the image, default (400,300)
    camera_position, camera_focalpoint: tuples with 3 float 
           camera position, camera focal point (ex. center of mass)
           default None. If not None the the camera will be overwrite
    Returns
    -------
    Image : an IPython display Image object
      
    See Also
    --------
    polydata2renderer, loadSTL
       
    Notes
    -----
    This Code was modified fron the original code by Adamos Kyriakou
    published in https://pyscience.wordpress.com/
    
    """
    
    #update the camera if provided
    if camera_position!=None and camera_focalpoint!=None:
        camera =vtk.vtkCamera()
        camera.SetPosition(*camera_position)
        camera.SetFocalPoint(*camera_focalpoint)
        renderer.SetActiveCamera(camera)
    
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetOffScreenRendering(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(width, height)
    renderWindow.Render()
     
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.Update()
     
    writer = vtk.vtkPNGWriter()
    writer.SetWriteToMemory(1)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()
    data = str(buffer(writer.GetResult()))
    
    return Image(data)


cpdef loadSTL(filenameSTL):
    """
    
    loadSTL(filenameSTL)
    
    Load a STL wireframe file
    
    Parameters
    ----------
    filenameSTL : file path
    
    Returns
    -------
    polydata : VTK Polydata object 
           
       
    Notes
    -----
    Code by Adamos Kyriakou
    published in https://pyscience.wordpress.com/
    
    """
    
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
    
cpdef polydata2renderer(polydata, color=(1.,0.,0.), 
                      opacity=1, background=(1.,1.,1.)):    
    """
    
    polydata2renderer(polydata, 
             color=None, 
             opacity=None, 
             background=None)
    
    Creates vtk renderer from vtk polydata
    
    Parameters
    ----------
    polydata : VTK polydata
    
    color: tuple with 3 floats (RGB)
        default color=(1.,0.,0.)
    opacity: float 
        default 1 (opaque)
    background: tuple with 3 floats (RGB)
        default background=(1.,1.,1.)
        
        
    Returns
    -------
    VtkRenderer : an Vtk renderer containing VTK polydata 
       
    See Also
    --------
    vtk_show, loadSTL
    
    Notes
    -----
    This Code was modified fron the original code by Adamos Kyriakou
    published in https://pyscience.wordpress.com/
    
    """
    
    VtkMapper = vtk.vtkPolyDataMapper()
    VtkMapper.SetInputData(polydata)

    VtkActor = vtk.vtkActor()
    VtkActor.SetMapper(VtkMapper)
    VtkActor.GetProperty().SetColor(*color)
    VtkActor.GetProperty().SetOpacity(opacity)

    VtkRenderer = vtk.vtkRenderer()
    VtkRenderer.SetBackground(*background)
    VtkRenderer.AddActor(VtkActor)
    
    return VtkRenderer


cpdef addPoint(renderer, p, radius=1.0, color=(0.0, 0.0, 0.0)):
    """
    addPoint(renderer, 
             p, 
             radius=1.0, 
             color=(0.0, 0.0, 0.0))
    
    Adds a point into an existing VTK renderer 
    
    Parameters
    ----------
    renderer : VTK renderer
           renderer with vtk objects and properties 
    p: tuple with 3 float
           point location
    radius: radius of the point  
           Default 1.
    color: tuple with 3 float
            Default (0.0, 0.0, 0.0)
    
    Returns
    -------
    VtkRenderer : an Vtk renderer containing VTK polydata 
       
    See Also
    --------
    vtk_show, loadSTL
    
    Notes
    -----
    This Code was modified fron the original code by Adamos Kyriakou
    published in https://pyscience.wordpress.com/
    
    """
    point = vtk.vtkSphereSource()
    point.SetCenter(*p)
    point.SetRadius(radius)
    point.SetPhiResolution(100)
    point.SetThetaResolution(100)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(point.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)

    renderer.AddActor(actor)
    
    return renderer



    
    
cpdef addLine(renderer, p1, p2, color=(0.0, 0.0, 1.0)):
    """
    addLine(renderer, 
             p1,
             p2,  
             color=(0.0, 0.0, 1.0))
    
    Adds a line into an existing VTK renderer 
    
    Parameters
    ----------
    renderer : VTK renderer
           renderer with vtk objects and properties 
    p1,p2: tuple with 3 float
           point location
    radius: radius of the point  
           Default 1.
    color: tuple with 3 float
            Default (0.0, 0.0, 0.0)
    
    Returns
    -------
    VtkRenderer : an Vtk renderer containing VTK polydata 
       
    See Also
    --------
    vtk_show, loadSTL
    
    Notes
    -----
    This Code was modified fron the original code by Adamos Kyriakou
    published in https://pyscience.wordpress.com/
    
    """
    line = vtk.vtkLineSource()
    line.SetPoint1(*p1)
    line.SetPoint2(*p2)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(line.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)

    renderer.AddActor(actor)
    
    return renderer 

cpdef vtk_raycasting(surface, pSource, pTarget):
    """
    vtk_raycasting(surface, pSource, pTarget)
    
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
    
    Notes
    -----
    This Code was modified fron the original code by Adamos Kyriakou
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
cpdef pointquering(surface, 
                   double azm,
                   double dip,
                   np.ndarray [double, ndim=1] x, 
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   int test):
    """
    pointquering(surface, azm, dip, x, y, z, test)
    
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
           1 test inside closed surface. Here we use 
             vtkOBBTree::InsideOrOutside. Closed surface are required
           2 test 'above' surface 
           3 test 'below' surface 
           4 test 'inside' surface (the surface can be open)
    
    Returns
    -------
    inside : 1D array of integers 
        Indicator of point inclusion with values [0,1]  
        0 means that the point is not inside, above or below surface
        1 means that the point is inside, above or below surface
    
    See Also
    --------
    vtk_raycasting
    
    Notes
    -----
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
cpdef pointinsolid(surface, 
                   np.ndarray [double, ndim=1] x, 
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   double tolerance = .000001 ):
    """
    inside = pointquering(surface, x, y, z)
    
    Find points inside a closed surface using vtkSelectEnclosedPoints
        
    
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
    
    Notes
    -----
    It is assumed that the surface is closed and manifold.  
    Note that if this check is not performed and the surface is not 
    closed, the results are undefined.
    
    """
    
    #TODO: add check to verify that the surface is closed 
    
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
cpdef tables2wireframe( 
                   np.ndarray [double, ndim=1] x, 
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   np.ndarray [long, ndim=1] pid,
                   np.ndarray [long, ndim=1] pid1,
                   np.ndarray [long, ndim=1] pid2,
                   np.ndarray [long, ndim=1] pid3,
                   str filename = None):
    """
    mywireframe = tables2wireframe(x, y, z, pid1,pid2,pid3) 
    
    Takes a wireframe defined by two table: 
    
    - x, y, z : This is the points table
    - pid1,pid2,pid3: This is another table defining triangles with 
      three point IDs. 
        
    Parameters
    ---------- 
    x,y,z          : 1D array of floats
        coordinates of the points to be tested        
    pid1,pid2,pid3 : 1D array of integers
        triangles defined by three existing point IDs
    str filename   : Str (Default None)
        file name and path. If provided the file wireframe is saved to 
        this file.   
    
    Returns
    -------
    surface : VTK polydata
           wireframe imported
    
    
    Notes
    -----
    ``pid`` is the row number in the point table.
    
    ``pid`` indices start at zero
    
    For example: 
    
    a point table
    
    | x | y | z |
    -------------
    | 0 | 0 | 0 |
    | 1 | 0 | 0 | 
    | 0 | 1 | 0 |
    
    a triangle table
    
    |pid1|pid2|pid3|
    -------------
    | 0  | 1  | 2  |
    
    """
    
    # TODO: test this 
    
    cdef int i
    
    assert x.shape[0]==y.shape[0]==z.shape[0],          'Error: x, y, z or pid with different dimension'
    assert pid1.shape[0]==pid1.shape[0]==pid1.shape[0], 'Error: pid1,pid2 or pid3 with different dimension'



    points = vtk.vtkPoints();
    
    
    for i in range(x.shape[0]):
        points.InsertNextPoint(x[i],y[i],z[i]);


    #create triangles
    triangle = vtk.vtkTriangle()
    triangles = vtk.vtkCellArray()

    for i in range(pid1.shape[0]):
        triangle.GetPointIds().InsertId ( 0, pid1[i]);
        triangle.GetPointIds().InsertId ( 1, pid2[i] );
        triangle.GetPointIds().InsertId ( 2, pid3[i] );
        
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



















cpdef getbounds(polydata):    
    """
    
    boundingbox(polydata)
    
    Returns the bounding box limits x,y,z minimum and maximum
    
    
    Parameters
    ----------
    polydata : VTK polydata
    
        
    Returns
    -------
    xmin,xmax, ymin,ymax, zmin,zmax : float
        The geometry bounding box 
    
    """
    
    return polydata.GetBounds()


cpdef points2vtkfile(str path,                    
                   np.ndarray [double, ndim=1] x, 
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   data):
    """
    
    bpoints2vtkfile(path, x,y,z, data)
    
    Save points in vtk file
    
    
    Parameters
    ----------
    
    
        
    Returns
    -------
    
    
    """
    
    pyevtk.hl.pointsToVTK(path, x, y, z, data = data) 
    
    
cpdef grid2vtkfile(str path,                    
                   np.ndarray [double, ndim=1] x, 
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   data):  
    """
    
    grid2vtkfile(path, x,y,z, data)
    
    save grid in vtk file if x,y,z are 1D it saves VtkRectilinearGrid,
    if are 3D it saves VtkStructuredGrid
    
    
    Parameters
    ----------
    
    
        
    Returns
    -------
    
    
    """

    pyevtk.hl.gridToVTK(path, x, y, z, cellData=data)


cpdef partialgrid2vtkfile(str path,                    
                   np.ndarray [double, ndim=1] x, 
                   np.ndarray [double, ndim=1] y,
                   np.ndarray [double, ndim=1] z,
                   double DX, 
                   double DY,
                   double DZ,
                   object var,
                   object varname):
    """
    
    partialgrid2vtkfile(path, x,y,z, dx,dy,dz,var, varname)
    
    save unstructurw vtk grid of parent blocks into a file
    
    
    Parameters
    ----------
    
    
        
    Returns
    -------
    
    
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
