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


def vtk_show(renderer, width=400, height=300, 
             camera_position=None, 
             camera_focalpoint=None):
    """
    
    vtk_show(renderer, 
             width=400, from libc.math cimport sin
from libc.math cimport cos
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


def loadSTL(filenameSTL):
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
    
def polydata2renderer(polydata, color=(1.,0.,0.), 
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


def addPoint(renderer, p, radius=1.0, color=(0.0, 0.0, 0.0)):
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
    
    
def addLine(renderer, p1, p2, color=(0.0, 0.0, 1.0)):
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
            intersect=obbTree.InsideOrOutside(p0)
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
