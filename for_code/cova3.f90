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

subroutine cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,c0,it,cc,aa, &
    irot,maxrot,rotmat,cmax,cova)
    !-----------------------------------------------------------------------

    !                    Covariance Between Two Points
    !                    *****************************

    ! This subroutine calculated the covariance associated with a variogram
    ! model specified by a nugget effect and nested varigoram structures.
    ! The anisotropy definition can be different for each nested structure.



    ! INPUT VARIABLES:

    !   x1,y1,z1         coordinates of first point
    !   x2,y2,z2         coordinates of second point
    !   nst              number of nested structures (same number for all variograms)
    !   ivarg            variogram number (set to 1 unless doing cokriging
    !                       or indicator kriging and in this case use same number of structures in all variograms)
    !   c0(ivarg)        isotropic nugget constant
    !   it(i)            type of each nested structure:
    !                      1. spherical model of range a;
    !                      2. exponential model of parameter a;
    !                           i.e. practical range is 3a
    !                      3. gaussian model of parameter a;
    !                           i.e. practical range is a*sqrt(3)
    !                      4. power model of power a (a must be gt. 0  and
    !                           lt. 2).  if linear model, a=1,c=slope.
    !                      5. hole effect model
    !   cc(i)            multiplicative factor of each nested structure.
    !                      (sill-c0) for spherical, exponential,and gaussian
    !                      slope for linear model.
    !   aa(i)            parameter "a" of each nested structure.
    !   irot             index of the rotation matrix for the first nested
    !                    structure (the second nested structure will use
    !                    irot+1, the third irot+2, and so on)
    !   maxrot           size of rotation matrix arrays
    !   rotmat           rotation matrices
    ! 
    !  Note that aa, it and cc are 1D arrays with size (MXVARG*MAXNST).     
    !  MAXNST was removed from this code and recalculated as MAXROT=MAXNST+1
    !  MXVARG is equal to ivarg
    ! 
    ! OUTPUT VARIABLES:

    !   cmax             maximum covariance
    !   cova             covariance between (x1,y1,z1) and (x2,y2,z2)



    ! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
    !                      rotmat    computes rotation matrix for distance
    !-----------------------------------------------------------------------

    implicit none
    !external references
    real*8 :: sqdist


    ! input
    real, intent(in) :: x1,y1,z1,x2,y2,z2
    integer, intent(in) :: ivarg, irot,maxrot
    integer, intent(in) :: nst
    real, intent(in), dimension(ivarg) :: c0
    integer, intent(in),dimension(nst*ivarg) :: it
    real, intent(in), dimension(nst*ivarg) :: cc, aa
    real*8, intent(in), dimension(maxrot,3,3) :: rotmat

    ! output
    real, intent(out) :: cmax, cova

    ! internal variables
    real ::    hsqd, h, hr
    integer :: ir, is, ist, istart
    
    !parameters
    real :: DEG2RAD,EPSLON,PI,PMX
 
    DEG2RAD=3.141592654/180.0
    EPSLON=1.e-20
    PI=3.14159265
    PMX=999.
    EPSLON=1.e-10


 
    ! Calculate the maximum covariance value (used for zero distances and
    ! for power model covariance):


    istart = 1 + (ivarg-1)*nst
    cmax   = c0(ivarg)
    do is=1,nst
        ist = istart + is - 1
        if(it(ist) == 4) then
            cmax = cmax + PMX
        else
            cmax = cmax + cc(ist)
        endif
    end do

    ! Check for "zero" distance, return with cmax if so:

    hsqd = real(sqdist(x1,y1,z1,x2,y2,z2,irot,maxrot,rotmat))
    if(dble(hsqd) < EPSLON) then
        cova = cmax
        return
    endif

    ! Loop over all the structures:

    cova = 0.0
    do is=1,nst
        ist = istart + is - 1
    
    ! Compute the appropriate distance:
    
        if(ist /= 1) then
            ir = min((irot+is-1),maxrot)
            hsqd=real(sqdist(x1,y1,z1,x2,y2,z2,ir,maxrot,rotmat))
        end if
        h = real(sqrt(hsqd))
    
    ! Spherical Variogram Model?
    
        if(it(ist) == 1) then
            hr = h/aa(ist)
            if(hr < 1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
        
        ! Exponential Variogram Model?
        
        else if(it(ist) == 2) then
            cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
        
        ! Gaussian Variogram Model?
        
        else if(it(ist) == 3) then
            cova = cova + cc(ist)*exp(-(3.0*h/aa(ist)) &
            *(3.0*h/aa(ist)))
        
        ! Power Variogram Model?
        
        else if(it(ist) == 4) then
            cova = cova + cmax - cc(ist)*(h**aa(ist))
        
        ! Hole Effect Model?
        
        else if(it(ist) == 5) then
        !                 d = 10.0 * aa(ist)
        !                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
            cova = cova + cc(ist)*cos(h/aa(ist)*PI)
        endif
    end do

    ! Finished:

    return

end subroutine cova3


subroutine dcova3(x1,y1,z1,x2,y2,z2,ivarg,nst,c0,it,cc,aa, &
    irot,maxrot,rotmat,cmax,cova)
    !-----------------------------------------------------------------------

    !                    Covariance Between Two Points
    !                    *****************************

    ! This subroutine calculated the covariance associated with a variogram
    ! model specified by a nugget effect and nested varigoram structures.
    ! The anisotropy definition can be different for each nested structure.



    ! INPUT VARIABLES:

    !   x1,y1,z1         coordinates of first point
    !   x2,y2,z2         coordinates of second point
    !   nst              number of nested structures (same number for all variograms)
    !   ivarg            variogram number (set to 1 unless doing cokriging
    !                       or indicator kriging and in this case use same number of structures in all variograms)
    !   c0(ivarg)        isotropic nugget constant
    !   it(i)            type of each nested structure:
    !                      1. spherical model of range a;
    !                      2. exponential model of parameter a;
    !                           i.e. practical range is 3a
    !                      3. gaussian model of parameter a;
    !                           i.e. practical range is a*sqrt(3)
    !                      4. power model of power a (a must be gt. 0  and
    !                           lt. 2).  if linear model, a=1,c=slope.
    !                      5. hole effect model
    !   cc(i)            multiplicative factor of each nested structure.
    !                      (sill-c0) for spherical, exponential,and gaussian
    !                      slope for linear model.
    !   aa(i)            parameter "a" of each nested structure.
    !   irot             index of the rotation matrix for the first nested
    !                    structure (the second nested structure will use
    !                    irot+1, the third irot+2, and so on)
    !   maxrot           size of rotation matrix arrays
    !   rotmat           rotation matrices
    ! 
    !  Note that aa, it and cc are 1D arrays with size (MXVARG*MAXNST).     
    !  MAXNST was removed from this code and recalculated as MAXROT=MAXNST+1
    !  MXVARG is equal to ivarg
    ! 
    ! OUTPUT VARIABLES:

    !   cmax             maximum covariance
    !   cova             covariance between (x1,y1,z1) and (x2,y2,z2)



    ! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
    !                      rotmat    computes rotation matrix for distance
    !-----------------------------------------------------------------------

    implicit none
    !external references
    real*8:: dsqdist


    ! input
    real*8, intent(in) :: x1,y1,z1,x2,y2,z2
    integer, intent(in) :: ivarg, irot,maxrot
    integer, intent(in) :: nst
    real*8, intent(in), dimension(ivarg) :: c0
    integer, intent(in),dimension(nst*ivarg) :: it
    real*8, intent(in), dimension(nst*ivarg) :: cc, aa
    real*8, intent(in), dimension(maxrot,3,3) :: rotmat

    ! output
    real*8, intent(out) :: cmax, cova

    ! internal variables
    real*8 ::    hsqd, h, hr
    integer :: ir, is, ist, istart
    
    !parameters
    real*8 :: DEG2RAD,EPSLON,PI,PMX
 
    DEG2RAD=3.141592654/180.0
    EPSLON=1.e-20
    PI=3.14159265
    PMX=999.
    EPSLON=1.e-10


 
    ! Calculate the maximum covariance value (used for zero distances and
    ! for power model covariance):


    istart = 1 + (ivarg-1)*nst
    cmax   = c0(ivarg)
    do is=1,nst
        ist = istart + is - 1
        if(it(ist) == 4) then
            cmax = cmax + PMX
        else
            cmax = cmax + cc(ist)
        endif
    end do

    ! Check for "zero" distance, return with cmax if so:

    hsqd = dsqdist(x1,y1,z1,x2,y2,z2,irot,maxrot,rotmat)
    if(dble(hsqd) < EPSLON) then
        cova = cmax
        return
    endif

    ! Loop over all the structures:

    cova = 0.0
    do is=1,nst
        ist = istart + is - 1
    
    ! Compute the appropriate distance:
    
        if(ist /= 1) then
            ir = min((irot+is-1),maxrot)
            hsqd=dsqdist(x1,y1,z1,x2,y2,z2,ir,maxrot,rotmat)
        end if
        h = dble(dsqrt(hsqd))
    
    ! Spherical Variogram Model?
    
        if(it(ist) == 1) then
            hr = h/aa(ist)
            if(hr < 1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
        
        ! Exponential Variogram Model?
        
        else if(it(ist) == 2) then
            cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
        
        ! Gaussian Variogram Model?
        
        else if(it(ist) == 3) then
            cova = cova + cc(ist)*exp(-(3.0*h/aa(ist)) &
            *(3.0*h/aa(ist)))
        
        ! Power Variogram Model?
        
        else if(it(ist) == 4) then
            cova = cova + cmax - cc(ist)*(h**aa(ist))
        
        ! Hole Effect Model?
        
        else if(it(ist) == 5) then
        !                 d = 10.0 * aa(ist)
        !                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
            cova = cova + cc(ist)*cos(h/aa(ist)*PI)
        endif
    end do

    ! Finished:

    return

end subroutine dcova3
