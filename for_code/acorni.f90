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

real*8 function acorni(ixv, KORDEI)
    !-----------------------------------------------------------------------

    ! Fortran implementation of ACORN random number generator of order less
    ! than or equal to 12 (higher orders can be obtained by increasing the
    ! parameter value MAXORD).


    ! NOTES: 1. The common block
    !           IACO is used to transfer data into the function.

    !        2. Before the first call to ACORN the common block IACO must
    !           be initialised by the user, as follows. The values of
    !           variables in the common block must not subsequently be
    !           changed by the user.

    !             KORDEI - order of generator required ( must be =< MAXORD)

    !             ixv(1) - seed for random number generator
    !                      require 0 < ixv(1) < MAXINT

    !             (ixv(I+1),I=1,KORDEI)
    !                    - KORDEI initial values for generator
    !                      require 0 =< ixv(I+1) < MAXINT

    !        3. After initialisation, each call to ACORN generates a single
    !           random number between 0 and 1.

    !        4. An example of suitable values for parameters is

    !             KORDEI   = 10
    !             MAXINT   = 2**30
    !             ixv(1)   = an odd integer in the (approximate) range
    !                        (0.001 * MAXINT) to (0.999 * MAXINT)
    !             ixv(I+1) = 0, I=1,KORDEI


    ! Author: R.S.Wikramaratna,                           Date: October 1990
    
    ! Note: 
    
    ! This function was modified in 2015 by Adrian Martinez (opengeostat)
    ! * the common block data transference was replaced by parameter transference
    ! * the function definition is different, before it was 
    !   double precision function acorni(idum), with idum as dummy variable. 
    ! 
    ! Warning: you may update the way you call this function in preexisting gslib code
    !-----------------------------------------------------------------------

    implicit none


    ! input 
    integer, intent(in) :: KORDEI

    !inout
    real*8, intent(inout), dimension (KORDEI+1) :: ixv

    ! internal 
    integer i, MAXINT

    MAXINT = 2**30

    do i=1,KORDEI
        ixv(i+1)=(ixv(i+1)+ixv(i))
        if(ixv(i+1) >= MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
    end do

    acorni=dble(ixv(KORDEI+1))/MAXINT

    return


end function acorni
