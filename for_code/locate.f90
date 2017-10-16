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
subroutine locate(xx,n,is,ie,x,j)
    ! -----------------------------------------------------------------------
    ! 
    ! Given an array "xx" of length "n", and given a value "x", this routine
    ! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
    ! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
    ! returned to indicate that x is out of range.
    ! 
    ! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
    ! -----------------------------------------------------------------------

    real, intent(in), dimension(n) :: xx
    real :: x
    integer, intent (out) :: j
    !f2py intent(out) j  
    
    ! Initialize lower and upper methods:
 
    if(is.le.0) is = 1
    jl = is-1
    ju = ie
    if(xx(n).le.x) then
        j = ie
        return
    end if

    ! If we are not done then compute a midpoint: 
    10 if(ju-jl.gt.1) then
        jm = (ju+jl)/2

        ! Replace the lower or upper limit with the midpoint:
        if((xx(ie).gt.xx(is)).eqv.(x.gt.xx(jm))) then
              jl = jm
        else
              ju = jm
        endif
        go to 10
    endif

    ! Return with the array index:

    j = jl

    return

end subroutine locate


subroutine dlocate(xx,n,is,ie,x,j)
    !-----------------------------------------------------------------------

    ! Given an array "xx" of length "n", and given a value "x", this routine
    ! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
    ! must be monotonic, either increasing or decreasing.  j=0 or j=n is
    ! returned to indicate that x is out of range.

    ! Modified to set the start and end points by "is" and "ie"

    ! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
    !-----------------------------------------------------------------------
    
    real*8, intent (in) :: x
    real*8, intent (in), dimension(n) :: xx
    integer, intent (in) :: n,is,ie
    
    integer, intent (out) :: j
    !f2py intent(out) j
    
    integer :: jl,ju,jm, iss
    
    iss = is
    
    ! Initialize lower and upper methods:
    if(iss.le.0) iss = 1
    jl = iss-1
    ju = ie
    if(xx(n).le.x) then
        j = ie
        return
    end if
    
    ! If we are not done then compute a midpoint:

    10 if(ju-jl > 1) then
    
        jm = (ju+jl)/2
    
        ! Replace the lower or upper limit with the midpoint:
    
        if((xx(ie) > xx(iss)).eqv.(x > xx(jm))) then
            jl = jm
        else
            ju = jm
        endif
        go to 10
    endif

    ! Return with the array index:

    j = jl
    
    return
end subroutine dlocate
