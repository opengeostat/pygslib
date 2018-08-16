TODO list:  
=====

Priority
----

- [ ] (conda package) Update conda module
- [ ] (nonlinear) Implement Uniform conditioning
- [ ] (nonlinear) Implement localization
- [ ] (nonlinear) Implement calc cdf function using pygslib.gslib.__plot.probplt
- [ ] (nonlinear) Add upper and lower tail CDF model fitting tool
- [ ] (blockmodel) implement grade tonnage report and block cdf
- [ ] (gslib) complete variogram3d and python interface
- [ ] (gslib) create an Ipython interface with widgets for friendly variogram modeling 
- [ ] (blockmodel) implement grade tonnage report and block cdf
- [ ] Add bokeh to dependencies or remove fuctions using this package  


Do this later...
---
- [x] (user manual) Update user manual and find a way to update faster the user manual (probably a separate git repo for doc sphinx?). Using shinx page on git\docs
- [ ] (user manual) Simplify examples
- [ ] (general) Review interface and improve user experience
- [ ] (general) Add GUI for ipython
- [ ] (general) Add default plot parameters accessible via **kwargs arguments. 
- [ ] (drillhole) Implement the prototype *Drillhole generator _v3.ipynb*. This is to drill drilholes in simulated block models. 
- [ ] (drillhole) Implement bench composite.
- [ ] (drillhole) Add a dictionary with drillholes names, traces and properties for faster programming in class members and to create user's reports. Add summary of drillholes in Drillhole class.
- [ ] (drillhole) Create export in other formats, example inventor (see vtkIVWriter).
- [ ] (drillhole) Create a dictionary with warning and errors in validations and do only one warning. 
- [ ] (drillhole) Remove one interval survey error, fix it automatically and generate warning.
- [ ] (blockmodel) Reimplement the block fill with vtkUniformGrid, this is similar to vtkImage but is the base for AMR objects 
- [ ] (blockmodel) Implement block subcells (AMR or Datamine style?)
- [ ] (blockmodel) Implement implicit modeler functions using vtkImplicitModeller
- [ ] (blockmodel) Implement block split   
- [ ] (gslib,blockmodel) Add VTK export of debug data in KT3D function, including search ellipse and target block/point
- [ ] (vtk, general) Send data directly to Paraview using Collaboration
- [ ] (gslib) implement kriging with proportions. This is to support paper
- [ ] (gslib) implement a simplified gslib kriging function for better maintenance and implementation. Separate from search and other functions. 
- [ ] (gslib) consider rewriting new non gslib kriging module in cython or fortran
- [ ] (gslib) implement simulation programs
- [ ] (gslib) implement cokriging 
- [ ] (general) Create script to install Paraview python module. Add file "paraview.pth" with text "C:\{ParaView path}\bin\Lib\site-packages"
- [ ] (vtk) Replace save polydata (for points) in vtk legacy and add save/read polydata (for surfaces) in xml
- [ ] (vtk) Add boolean operations (find alternative to vtkbooleanoperationpolydatafilter, e.j. vtkImplicitBoolean)
- [ ] (vtk) Implement Unfolding with vtkwarpvector (hint, vector interpolation (angle + scallar) required)
- [ ] (vtk) Add slice operation
- [ ] (vtk) Replace pyevtk with vtk.util.numpy_support
- [ ] (vtk) Add export 3D grid (centroids points / no blocks)
- [ ] (drillhole) Generate reposts with reportlab??? for example for log files. Or use pygslib.charttable to creat nice plots in Ipython
 
Avoid bugs and improve compilation
----
- [ ] Add explicit interface (blocks) for external or dummy procedures in FORTRAN code.  This is creating this warning ``Could not found the body of interfaced routine``
- [ ] Check compilation warning ``getarrdims:warning: assumed shape array, using 0 instead of '*'``
- [ ] Check compilation warnings ``Warning: Possible change of value in conversion`` and ``Type mismatch in argument``. This already created some issues in some of the FORTRAN programs by transferring Real*8 to Real arrays.
- [ ] Clean a bit the code ``Warning: Label 4 at (1) defined but not used`` and  ``Warning: Unused dummy argument``.
- [ ] Initialize arrays ``Warning: '####' may be used uninitialized in this function``
- [ ] test dtype of dict parameters at python level to avoid dtype complains and error from fortran at gslib.kt3d function
- [ ] review cython c declarations to avoid rounding error. Use ``cdef double`` instead of ``cdef float``.
- [ ] Implement constants at python level to define integer inputs like variogram types, type of desurvey, kriging options, etc. 