TODO list:  
=====

Avoid bugs and improve compilation
----
- [] Add explicit interface (blocks) for external or dummy procedures in FORTRAN code.  This is creating this warning ``Could not found the body of interfaced routine``
- [] Check compilation warning ``getarrdims:warning: assumed shape array, using 0 instead of '*'``
- [] Check compilation warnings ``Warning: Possible change of value in conversion`` and ``Type mismatch in argument``. This already created some issues in some of the FORTRAN programs by transferring Real*8 to Real arrays.
- [] Clean a bit the code ``Warning: Label 4 at (1) defined but not used`` and  ``Warning: Unused dummy argument``.
- [] Initialize arrays ``Warning: ‘####’ may be used uninitialized in this function``

Improve algorithms and testing
----
- [] Complete the estimation part for n node, blocks or polygons (kt3D program)
- [] Fix the issuaes with gaussian anamorphosis and do more testing. 
- [x] Do more testing on VTK selection
- [x] Do more testing on desurvey drillhole


User interface and help
----
- [] Simplify examples
- [] Do all in one example
- [] Review doc string
- [] Write user manual

TODO at drillhole module 
----
- [] Add summary of drillholes in Drillhole class. 
- [] Add a dictionary with drillholes names, traces and properties for faster programming in class members and to create user's reports.
- [] Create a dictionary with warning and errors in validations and do only one warning. 
- [] Generate reposts with reportlab??? for example for log files. 
- [] Remove the need for one interval error (holes with survey at collar only).  This, is really slow and not user friendly.
- [] Create warning in desurvey if there are intervals undefined (not created)
- [x] In VTK export remove NaNs in coordinates internally and produce warning if there are NaNs.
- [] VTK export is too slow, optimize. Consider to use directly VTK code in C+ or Python
- [] Create export in other formats, example inventor.
- [] VTK write XML polydata format, in other words: PolyData (.vtp) — Serial vtkPolyData (unstructured).
- [] Add lithocode and more options to Downhole composite.
- [] Validate and test composite.
- [] Implement bench composite.
- [] Add HTMLbar tables.

TODO at block model module 
----

TODO at interpolators module 
----

TODO at neighborhood module 
----


TODO at nonlinear module 
----


TODO at vtk tools module 
----

