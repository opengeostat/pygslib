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

subroutine gauinv(p,xp,ierr)
    !-----------------------------------------------------------------------

    ! Computes the inverse of the standard normal cumulative distribution
    ! function with a numerical approximation from : Statistical Computing,
    ! by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95.



    ! INPUT/OUTPUT:

    !   p    = double precision cumulative probability value: dble(psingle)
    !   xp   = G^-1 (p) in single precision
    !   ierr = 1 - then error situation (p out of range), 0 - OK


    !-----------------------------------------------------------------------
    real :: xp
    real*8 :: p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim,p
    save   p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim

    ! Coefficients of approximation:

    data lim/1.0e-10/
    data p0/-0.322232431088/,p1/-1.0/,p2/-0.342242088547/, &
    p3/-0.0204231210245/,p4/-0.0000453642210148/
    data q0/0.0993484626060/,q1/0.588581570495/,q2/0.531103462366/, &
    q3/0.103537752850/,q4/0.0038560700634/

    ! Check for an error situation:

    ierr = 1
    if(p < lim) then
        xp = -1.0e10
        return
    end if
    if(p > (1.0-lim)) then
        xp =  1.0e10
        return
    end if
    ierr = 0

    ! Get k for an error situation:

    pp   = p
    if(p > 0.5) pp = 1 - pp
    xp   = 0.0
    if(p == 0.5) return

    ! Approximate the function:

    y  = dsqrt(dlog(1.0/(pp*pp)))
    xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) / &
    ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
    if(real(p) == real(pp)) xp = -xp

    ! Return with G^-1(p):

    return
end subroutine gauinv
