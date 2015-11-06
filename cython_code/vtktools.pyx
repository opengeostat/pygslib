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


def vtk_show(renderer, width=400, height=300, 
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
