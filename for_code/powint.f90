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

real function powint(xlow,xhigh,ylow,yhigh,xval,pow)
    !-----------------------------------------------------------------------

    ! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
    !                 for a value of x and a power pow.

    !-----------------------------------------------------------------------
    real :: EPSLON, xlow,xhigh,ylow,yhigh,xval, pow

    EPSLON=1.0e-20

    if((xhigh-xlow) < EPSLON) then
        powint = (yhigh+ylow)/2.0
    else
        powint = ylow + (yhigh-ylow)* &
        (((xval-xlow)/(xhigh-xlow))**pow)
    end if

    return

end function powint

real*8 function dpowint(xlow,xhigh,ylow,yhigh,xval,pow)
    !-----------------------------------------------------------------------

    ! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
    !                 for a value of x and a power pow.

    !-----------------------------------------------------------------------
    implicit real*8 (a-h,o-z)

    parameter(EPSLON=1.0e-20)

    if((xhigh-xlow) < EPSLON) then
        dpowint = (yhigh+ylow)/2.0
    else
        dpowint = ylow + (yhigh-ylow)* &
        (((xval-xlow)/(xhigh-xlow))**pow)
    end if

    return
end function
