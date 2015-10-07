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

subroutine getindx(n,min,siz,loc,index,inflag)
    !-----------------------------------------------------------------------

    !     Gets the coordinate index location of a point within a grid
    !     ***********************************************************


    ! n       number of "nodes" or "cells" in this coordinate direction
    ! min     origin at the center of the first cell
    ! siz     size of the cells
    ! loc     location of the point being considered
    ! index   output index within [1,n]
    ! inflag  true if the location is actually in the grid (false otherwise
    !         e.g., if the location is outside then index will be set to
    !         nearest boundary



    !-----------------------------------------------------------------------
    implicit none
    

    !input
    integer, intent(in) :: n
    real*8, intent(in) ::  min, siz, loc

    !output
    integer, intent(out) :: index 
    logical, intent(out) :: inflag



    ! Compute the index of "loc":

    index = int( (loc-min)/siz + 1.5 )

    ! Check to see if in or out:

    if(index < 1) then
        index  = 1
        inflag = .FALSE. 
    else if(index > n) then
        index  = n
        inflag = .FALSE. 
    else
        inflag = .TRUE. 
    end if

    ! Return to calling program:

    return

end subroutine getindx
