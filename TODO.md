TODO list:  
=====

Avoid bugs and improve compilation
----
- [] Add explicit interface (blocks) for external or dummy procedures in FORTRAN code.  This is creating this warning ``Could not found the body of interfaced routine``
- [] Check compilation warning ``getarrdims:warning: assumed shape array, using 0 instead of '*'``
- [] Check compilation warnings ``Warning: Possible change of value in conversion`` and ``Type mismatch in argument``. This already created some issues in some of the FORTRAN programs by transferring Real*8 to Real arrays.
- [] Clean a bit the code ``Warning: Label 4 at (1) defined but not used`` and  ``Warning: Unused dummy argument``.
- [] Initialize arrays ``Warning: ‘####’ may be used uninitialized in this function``
- [] Initialize arrays

Improve algorithms and testing
----
- [] Complete the estimation part for n node, blocks or polygons (kt3D program)
- [] Optimize speed/Identify slow code and bottleneck 

User manual and help
----
- [] Simplify examples
- [] Do all in one example
- [] Review doc string
- [] Write user manual

TODO at drillhole module 
----
- [] Validate merge tables.
- [] Implement bench composite.
- [] Add HTMLbar tables.
- [] Update Ipython templates.
- [x] Do more testing on desurvey drillhole
- [x] Make numeric but non float float before importing
- [x] Validate and test composite.
- [x] In VTK export remove NaNs in coordinates internally and produce warning if there are NaNs.
- [x] Create warning in desurvey if there are intervals undefined (not created)
- [] Add a dictionary with drillholes names, traces and properties for faster programming in class members and to create user's reports. Add summary of drillholes in Drillhole class.
- [] Generate reposts with reportlab??? for example for log files. 
- [] VTK export is too slow, optimize. Consider to use directly VTK code in C+ or Python. VTK write XML polydata format, in other words: PolyData (.vtp) — Serial vtkPolyData (unstructured).
- [] Create export in other formats, example inventor.
- [] Create a dictionary with warning and errors in validations and do only one warning. 
- [] Remove the need for one interval error ?

TODO at block model module 
----
- [] The block fill is too slow with large/complicated wireframes, optimize
- [] Implement block percentage in solid
- [] Implement block split
- [] Implement reblock
- [] Implement block subcells (AMR or Datamine style?)
- [] Export partial block models


TODO at interpolators module 
----

TODO at neighborhood module 
----


TODO at nonlinear module 
----
- [] Fix the issues with gaussian anamorphosis and do more testing. 

TODO at vtk tools module 
----
- [x] Do more testing on VTK selection
