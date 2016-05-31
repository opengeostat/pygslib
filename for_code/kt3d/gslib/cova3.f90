    subroutine cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa, &
    irot,MAXROT,rotmat,cmax,cova)
!-----------------------------------------------------------------------

!                    Covariance Between Two Points
!                    *****************************

! This subroutine calculated the covariance associated with a variogram
! model specified by a nugget effect and nested varigoram structures.
! The anisotropy definition can be different for each nested structure.



! INPUT VARIABLES:

!   x1,y1,z1         coordinates of first point
!   x2,y2,z2         coordinates of second point
!   nst(ivarg)       number of nested structures (maximum of 4)
!   ivarg            variogram number (set to 1 unless doing cokriging
!                       or indicator kriging)
!   MAXNST           size of variogram parameter arrays
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
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices


! OUTPUT VARIABLES:

!   cmax             maximum covariance
!   cova             covariance between (x1,y1,z1) and (x2,y2,z2)



! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
!                      rotmat    computes rotation matrix for distance
!-----------------------------------------------------------------------
    parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-10)
    integer ::   nst(*),it(*)
    real ::      c0(*),cc(*),aa(*)
    real*8 ::    rotmat(MAXROT,3,3),hsqd,sqdist

! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):

    istart = 1 + (ivarg-1)*MAXNST
    cmax   = c0(ivarg)
    do is=1,nst(ivarg)
        ist = istart + is - 1
        if(it(ist) == 4) then
            cmax = cmax + PMX
        else
            cmax = cmax + cc(ist)
        endif
    end do

! Check for "zero" distance, return with cmax if so:

    hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
    if(real(hsqd) < EPSLON) then
        cova = cmax
        return
    endif

! Loop over all the structures:

    cova = 0.0
    do is=1,nst(ivarg)
        ist = istart + is - 1
    
    ! Compute the appropriate distance:
    
        if(ist /= 1) then
            ir = min((irot+is-1),MAXROT)
            hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
        end if
        h = real(dsqrt(hsqd))
    
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
