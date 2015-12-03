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
- [] Do more testing on GTK selection
- [] Do more testing on desurvey drillhole


User interface and help
----
- [] Simplify examples
- [] Do all in one example
- [] Review doc string
- [] Write user manual
