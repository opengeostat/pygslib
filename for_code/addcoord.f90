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

subroutine addcoord(nx, ny, nz, xmn, ymn, zmn, xsiz, ysiz, zsiz, xx, yy, zz) 


    !
    !              Calculate coordinates of an ordered grid 
    !              ****************************************
    ! INPUT VARIABLES:
    !
    !   nx, ny, nz           Number nodes in the grid
    !   xmn, ymn, zmn        Coordinate of the centroid (lower left corner)
    !   xsiz, ysiz, zsiz     Sizes of the cell/block 
    !
    ! OUTPUT VARIABLES:
    !
    !   xx(n), yy(n), zz(n)  coordinate arrays with size n= nx*ny*nz
    !-----------------------------------------------------------------------
    
    real*8,  intent(in) :: xmn, ymn, zmn, xsiz, ysiz, zsiz
    integer, intent(in) :: nx, ny, nz
    real*8,    intent(out), dimension(nx*ny*nz) :: xx, yy, zz
    
    real*8 :: x, y, z
    integer :: i

    !
    ! Loop over grid
    !
    i=0
    do iz=1,nz
          z = zmn + dble(iz-1)*zsiz
          do iy=1,ny
                y = ymn + dble(iy-1)*ysiz
                do ix=1,nx
                      x = xmn + dble(ix-1)*xsiz
                      i=i+1
                      xx(i)=x
                      yy(i)=y
                      zz(i)=z
                end do
          end do
    end do

    return

end subroutine addcoord
