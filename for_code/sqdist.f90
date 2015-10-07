!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 2015 Adrian Martinez Vargas                            %
!                                                                      %
! This software may be modified and distributed under the terms        %
! of the MIT license.  See the LICENSE.txt file for details.           %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
! The functions and subroutines below were modified from the
! version 2.0 of the gslib code written in fortran 77
! 
! The objective is to add functionality to GSLIB and to link 
! this code with python using f2py. It uses f90 code convention. 
! The arrays are dynamic and externally declared (from python)
!
! the code was converted from Fortran 77 to Fortran 90 using F2F.pl
! 
! for more information please refer to:
! - gslib77 source code: http://www.statios.com/software/gslib77_ls.tar.gz
! - GSLIB: Geostatistical Software Library and User's Guide. Second edition 
!   by Clayton V. Deutsch, Andre G. Journel, 1997.
! - F2PY Users Guide: http://docs.scipy.org/doc/numpy-dev/f2py/
! - F2F https://bitbucket.org/lemonlab/f2f
!-----------------------------------------------------------------------

real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,maxrot,rotmat)
    !-----------------------------------------------------------------------
    !
    !    Squared Anisotropic Distance Calculation Given Matrix Indicator
    !    ***************************************************************
    !
    ! This routine calculates the anisotropic distance between two points
    !  given the coordinates of each point and a definition of the
    !  anisotropy.
    !
    !
    ! INPUT VARIABLES:
    !
    !   x1,y1,z1         Coordinates of first point
    !   x2,y2,z2         Coordinates of second point
    !   ind              The rotation matrix to use
    !   maxrot           The maximum number of rotation matrices dimensioned
    !   rotmat           The rotation matrices
    !
    !
    !
    ! OUTPUT VARIABLES:
    !
    !   sqdis           The squared distance accounting for the anisotropy
    !                      and the rotation of coordinates (if any).
    !
    !
    ! NO EXTERNAL REFERENCES
    !
    !
    !-----------------------------------------------------------------------
     
    implicit none
    
    ! input
    integer, intent(in) ::  maxrot, ind
    real*8, intent(in):: x1,y1,z1,x2,y2,z2
    real*8, intent(in), dimension(maxrot,3,3) ::  rotmat

    
    ! output
    ! real*8 :: sqdis

    ! Internal 
    real*8 :: cont,dx,dy,dz
    integer  :: i

    !
    ! Compute component distance vectors and the squared distance:
    !
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx &
                   + rotmat(ind,i,2) * dy &
                   + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
      end do
    return
end function sqdist
