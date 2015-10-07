!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 2015 Adrian Martinez Vargas                            %
!                                                                      %
! This software may be modified and distributed under the terms        %
! of the MIT license.  See the LICENSE.txt file for details.           %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine rotscale(X,Y,Z,nd,X0,Y0,Z0,ang1,ang2,ang3,anis1,anis2,invert,Xr,Yr,Zr)
    !-----------------------------------------------------------------------

    !              Rotate coordinates
    !              ******************

    ! This is implemented as in http://www.ccgalberta.com/ccgresources/report06/2004-403-angle_rotations.pdf



    ! INPUT PARAMETERS:
    
    !   X,Y, Z           X,Y,Z arrays with coordinates
    !   nd               number of data points
    !   X0,Y0,Z0         New origin of coordinate
    !   ang1             Azimuth angle for principal direction
    !   ang2             Dip angle for principal direction
    !   ang3             Third rotation angle
    !   anis1            First anisotropy ratio
    !   anis2            Second anisotropy ratio
    !   invert           If 0 do rotation, if 1 invert rotation


    ! Rerurn
    !   Xr,Yr,Zr         New rotated and scaled X,Y,Z arrays with coordinates


    !-----------------------------------------------------------------------

    implicit none

    ! input
    real*8, intent(in), dimension(nd) ::X,Y,Z
    real*8, intent(in) :: ang1,ang2,ang3,anis1,anis2, X0,Y0,Z0
    integer, intent(in) :: invert, nd

    ! output
    real*8, intent(out), dimension(nd) ::Xr,Yr,Zr

    ! internal variables
    real*8 ::   alpha, beta, theta, sina,sinb,sint, &
                cosa,cosb,cost,afac1,afac2, &
                X1,Y1,Z1, X2,Y2,Z2
    integer :: i
         
    
    !parameters
    real*8 :: DEG2RAD,EPSLON

    DEG2RAD=3.141592654/180.0
    EPSLON=1.e-20

    alpha = ang1 * DEG2RAD
    beta  = ang2 * DEG2RAD
    theta = ang3 * DEG2RAD

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
    
    if (invert == 0) then
        do i=1,nd
            !shift 
            X1 = X(i)-X0
            Y1 = Y(i)-Y0
            Z1 = Z(i)-Z0
            ! rotate using equation 4 on http://www.ccgalberta.com/ccgresources/report06/2004-403-angle_rotations.pdf
            X2= (cosa*cost+sina*sinb*sint)*X1 + (-sina*cost+cosa*sinb*sint)*Y1 +(-cosb*sint)*Z1
            Y2=                (sina*cosb)*X1 +                 (cosa*cosb)*Y1 +      (sinb)*Z1
            Z2= (cosa*sint-sina*sinb*cost)*X1 + (-sina*sint-cosa*sinb*cost)*Y1 + (cosb*cost)*Z1
            !rescale
            Xr(i)= X2
            Yr(i)= Y2*afac1
            Zr(i)= Z2*afac2
        end do
    else
        do i=1,nd
            !shift 
            X1 = X(i)+X0
            Y1 = Y(i)+Y0
            Z1 = Z(i)+Z0
            ! rotate using equation 5 on http://www.ccgalberta.com/ccgresources/report06/2004-403-angle_rotations.pdf
            X2= (cosa*cost+sina*sinb*sint)*X1 + (sina*cosb)*Y1 + (cosa*sint-sina*sinb*cost)*Z1
            Y2=(-sina*cost+cosa*sinb*sint)*X1 + (cosa*cosb)*Y1 +(-sina*sint-cosa*sinb*cost)*Z1
            Z2=               (-cosb*sint)*X1 +      (sinb)*Y1 +                (cosb*cost)*Z1
            !rescale
            Xr(i)= X2
            Yr(i)= Y2/afac1
            Zr(i)= Z2/afac2
        end do
    end if 
    
    ! Return to calling program:
    return
    
end subroutine rotscale
