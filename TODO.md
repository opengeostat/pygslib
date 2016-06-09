TODO list:  
=====

Avoid bugs and improve compilation
----
- [ ] Add explicit interface (blocks) for external or dummy procedures in FORTRAN code.  This is creating this warning ``Could not found the body of interfaced routine``
- [ ] Check compilation warning ``getarrdims:warning: assumed shape array, using 0 instead of '*'``
- [ ] Check compilation warnings ``Warning: Possible change of value in conversion`` and ``Type mismatch in argument``. This already created some issues in some of the FORTRAN programs by transferring Real*8 to Real arrays.
- [ ] Clean a bit the code ``Warning: Label 4 at (1) defined but not used`` and  ``Warning: Unused dummy argument``.
- [ ] Initialize arrays ``Warning: ‘####’ may be used uninitialized in this function``
- [ ] Initialize arrays

Improve algorithms and testing
----
- [ ] Complete the estimation part for n node, blocks or polygons (kt3D program)
- [ ] Optimize speed/Identify slow code and bottleneck 


User manual and help
----
- [ ] Simplify examples
- [ ] Do all in one example
- [ ] Review doc string
- [ ] Write user manual

TODO at drillhole module 
----
- [ ] Implement the prototype *Drillhole generator _v3.ipynb*
- [ ] Validate merge tables.
- [ ] Implement bench composite.
- [ ] Add HTMLbar tables.
- [ ] Update Ipython templates.
- [ ] Optimize code for desurvey and composite
- [x] Do more testing on desurvey drillhole
- [x] Make numeric but non float float before importing
- [x] Validate and test composite.
- [x] In VTK export remove NaNs in coordinates internally and produce warning if there are NaNs.
- [x] Create warning in desurvey if there are intervals undefined (not created)
- [ ] Add a dictionary with drillholes names, traces and properties for faster programming in class members and to create user's reports. Add summary of drillholes in Drillhole class.
- [ ] Generate reposts with reportlab??? for example for log files. 
- [ ] VTK export is too slow, optimize. Consider to use directly VTK code in C+ or Python. VTK write XML polydata format, in other words: PolyData (.vtp) — Serial vtkPolyData (unstructured).
- [ ] Create export in other formats, example inventor (see vtkIVWriter).
- [ ] Create a dictionary with warning and errors in validations and do only one warning. 
- [ ] Remove the need for one interval error ?

TODO at block model module 
----
- [x] The block fill is too slow with large/complicated wireframes, optimize (testing vtkPolyDataToImageStencil)
- [x] Implement block percentage in solid. **Note the actual solution is an approximation** 
- [ ] Reimplement the block fill with vtkUniformGrid, this is similar to vtkImage but is the base for AMR objects 
- [ ] Implement block subcells (AMR or Datamine style?)
- [ ] Implement implicit modeler functions using vtkImplicitModeller
- [ ] Implement block split
- [ ] Implement reblock
- [x] Export partial block models as unstructured VTK grid
- [x] Modify blocks2vtkUnstructuredGrid to export all variables in the model  
- [ ] Update examples. Add new functionality to examples 
- [ ] Add function to report summary with minimum maximum coordinates and other properties. Include limits at parent centroids and corner
- [ ] Add a simplified function to create a full parent model (create_IJK, then calc_ixyz_fromijk, then calc_xyz_fromixyz) 
- [ ] Add validation of IJK vs IX,IY,IZ, vs XC,YC and ZC. Hint: Use calc_ixyz_fromijk, calc_xyz_fromixyz and calc_ijk to create temporary variables and compare with numpy.isclose.
- [ ] Add VTK export of debug data in KT3D function, including search ellipse and target block/point


TODO at interpolators module 
----
- [ ] Review or Reimplement 

TODO at neighborhood module 
----
- [ ] Review or Reimplement 


TODO at nonlinear module 
----
- [ ] Fix the issues with gaussian anamorphosis and do more testing. 

TODO at vtk tools module 
----
- [x] pointquering with vtkOBBTree is too slow with some wireframes, try vtkCellLocator or optimize with grid/image locator or vtkSelectEnclosedPoints
- [ ] Replace save polydata (for points) in vtk legacy and add save/read polydata (for surfaces) in xml
- [ ] Add boolean operations (find alternative to vtkbooleanoperationpolydatafilter, e.j. vtkImplicitBoolean)
- [ ] Implement Unfolding with vtkwarpvector (hint, vector interpolation (angle + scallar) required)
- [ ] Add slice operation
- [x] Do more testing on VTK selection
- [ ] Replace pyevtk with vtk.util.numpy_support
- [ ] Update the code to properly export grids, use this code:

    ``` python

    writer = vtk.vtkStructuredPointsWriter();
    writer.SetFileName("myfilename.vtk");
    writer.SetInputData(grid)
    writer.Write()

    ```
