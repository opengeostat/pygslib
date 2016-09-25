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

subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,maxrot,rotmat)
    !-----------------------------------------------------------------------

    !              Sets up an Anisotropic Rotation Matrix
    !              **************************************

    ! Sets up the matrix to transform Cartesian coordinates to coordinates
    ! accounting for angles and anisotropy (see manual for a detailed
    ! definition):


    ! INPUT PARAMETERS:

    !   ang1             Azimuth angle for principal direction
    !   ang2             Dip angle for principal direction
    !   ang3             Third rotation angle
    !   anis1            First anisotropy ratio
    !   anis2            Second anisotropy ratio
    !   ind              matrix indicator to initialize
    !   maxrot           maximum number of rotation matrices dimension 
    !                    for example maxrot = number of structures in  variogram + 1
    !   rotmat           rotation matrices


    ! This code was modified from original f77 GSLIB code (v.2)
    ! Mayor changes
    ! rotmat is dynamically defined
    ! DEG2RAD,EPSLON redefined as variables
    ! maxrot is now variable defined externally... 


    !-----------------------------------------------------------------------

    implicit none

    
    
    ! input
    real, intent(in) :: ang1,ang2,ang3,anis1,anis2
    integer, intent(in) :: ind, maxrot

    ! output
    real*8, intent(out), dimension(maxrot,3,3) :: rotmat

    ! internal variables
    real::  afac1,afac2,sina,sinb,sint, cosa,cosb,cost, alpha, beta, theta
    
    !parameters
    real :: DEG2RAD,EPSLON

    DEG2RAD=3.141592654/180.0
    EPSLON=1.e-20

    ! Converts the input angles to three angles which make more
    !  mathematical sense:

    !         alpha   angle between the major axis of anisotropy and the
    !                 E-W axis. Note: Counter clockwise is positive.
    !         beta    angle between major axis and the horizontal plane.
    !                 (The dip of the ellipsoid measured positive down)
    !         theta   Angle of rotation of minor axis about the major axis
    !                 of the ellipsoid.

    if(ang1 >= 0.0 .AND. ang1 < 270.0) then
        alpha = (90.0   - ang1) * DEG2RAD
    else
        alpha = (450.0  - ang1) * DEG2RAD
    endif
    beta  = -1.0 * ang2 * DEG2RAD
    theta =        ang3 * DEG2RAD

    ! Get the required sines and cosines:

    sina  = sin(alpha)
    sinb  = sin(beta)
    sint  = sin(theta)
    cosa  = cos(alpha)
    cosb  = cos(beta)
    cost  = cos(theta)

    ! Construct the rotation matrix in the required memory:

    afac1 = 1.0 / max(anis1,EPSLON)
    afac2 = 1.0 / max(anis2,EPSLON)
    rotmat(ind,1,1) =       (cosb * cosa)
    rotmat(ind,1,2) =       (cosb * sina)
    rotmat(ind,1,3) =       (-sinb)
    rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
    rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
    rotmat(ind,2,3) = afac1*( sint * cosb)
    rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
    rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
    rotmat(ind,3,3) = afac2*(cost * cosb)

    ! Return to calling program:

    return
end subroutine setrot

subroutine dsetrot(ang1,ang2,ang3,anis1,anis2,ind,maxrot,rotmat)
    !-----------------------------------------------------------------------

    !              Sets up an Anisotropic Rotation Matrix
    !              **************************************

    ! Sets up the matrix to transform Cartesian coordinates to coordinates
    ! accounting for angles and anisotropy (see manual for a detailed
    ! definition):


    ! INPUT PARAMETERS:

    !   ang1             Azimuth angle for principal direction
    !   ang2             Dip angle for principal direction
    !   ang3             Third rotation angle
    !   anis1            First anisotropy ratio
    !   anis2            Second anisotropy ratio
    !   ind              matrix indicator to initialize
    !   maxrot           maximum number of rotation matrices dimension 
    !                    for example maxrot = number of structures in  variogram + 1
    !   rotmat           rotation matrices


    ! This code was modified from original f77 GSLIB code (v.2)
    ! Mayor changes
    ! rotmat is dynamically defined
    ! DEG2RAD,EPSLON redefined as variables
    ! maxrot is now variable defined externally... 


    !-----------------------------------------------------------------------

    implicit none

    
    
    ! input
    real*8, intent(in) :: ang1,ang2,ang3,anis1,anis2
    integer, intent(in) :: ind, maxrot

    ! output
    real*8, intent(out), dimension(maxrot,3,3) :: rotmat

    ! internal variables
    real*8 ::  afac1,afac2,sina,sinb,sint, cosa,cosb,cost, alpha, beta, theta
    
    !parameters
    real*8 :: DEG2RAD,EPSLON

    DEG2RAD=3.141592654/180.0
    EPSLON=1.e-20

    ! Converts the input angles to three angles which make more
    !  mathematical sense:

    !         alpha   angle between the major axis of anisotropy and the
    !                 E-W axis. Note: Counter clockwise is positive.
    !         beta    angle between major axis and the horizontal plane.
    !                 (The dip of the ellipsoid measured positive down)
    !         theta   Angle of rotation of minor axis about the major axis
    !                 of the ellipsoid.

    if(ang1 >= 0.0 .AND. ang1 < 270.0) then
        alpha = (90.0   - ang1) * DEG2RAD
    else
        alpha = (450.0  - ang1) * DEG2RAD
    endif
    beta  = -1.0 * ang2 * DEG2RAD
    theta =        ang3 * DEG2RAD

    ! Get the required sines and cosines:

    sina  = dble(sin(alpha))
    sinb  = dble(sin(beta))
    sint  = dble(sin(theta))
    cosa  = dble(cos(alpha))
    cosb  = dble(cos(beta))
    cost  = dble(cos(theta))

    ! Construct the rotation matrix in the required memory:

    afac1 = 1.0 / dble(max(anis1,EPSLON))
    afac2 = 1.0 / dble(max(anis2,EPSLON))
    rotmat(ind,1,1) =       (cosb * cosa)
    rotmat(ind,1,2) =       (cosb * sina)
    rotmat(ind,1,3) =       (-sinb)
    rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
    rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
    rotmat(ind,2,3) = afac1*( sint * cosb)
    rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
    rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
    rotmat(ind,3,3) = afac2*(cost * cosb)

    ! Return to calling program:

    return
end subroutine dsetrot
