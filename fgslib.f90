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
! version 2.0 of the gslib code written in fortran 77, unless explained
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


! compilation instructions
! requires f2py
!  use f2py -c -m fgslib fgslib.f90

! usage example in python
! >> import pygslib as gs
! >> [2]: dir(gs)
! >> [3]: print gs.read_data.__doc__


!*********************************************************************************
!     Subroutine for version control
!*********************************************************************************
subroutine version(major, minor , maintenance, build, month, year)
    !output the actual library version
    ! we use here the version convention major.minor[.maintenance[.build]]
    ! see http://en.wikipedia.org/wiki/Software_versioning

    integer, intent(out) ::    major, minor , maintenance, build, month, year
    
    major=0
    minor=0
    maintenance=0
    build=3
    month=9
    year=2015

    return

end subroutine version


!*********************************************************************************
!     Subroutines in GSLIB (auxiliary functions)
!*********************************************************************************
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
end subroutine setrot



real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,maxrot,rotmat)
    !-----------------------------------------------------------------------
    !
    !    Squared Anisotropic Distance Calculation Given Matrix Indicator
    !    ***************************************************************
    !
    ! This routine calculates the anisotropic distance between two points
    !  given the coordinates of each point and a definition of the
    !  anisotropy.
    !
    !
    ! INPUT VARIABLES:
    !
    !   x1,y1,z1         Coordinates of first point
    !   x2,y2,z2         Coordinates of second point
    !   ind              The rotation matrix to use
    !   maxrot           The maximum number of rotation matrices dimensioned
    !   rotmat           The rotation matrices
    !
    !
    !
    ! OUTPUT VARIABLES:
    !
    !   sqdis           The squared distance accounting for the anisotropy
    !                      and the rotation of coordinates (if any).
    !
    !
    ! NO EXTERNAL REFERENCES
    !
    !
    !-----------------------------------------------------------------------
     
    implicit none
    
    ! input
    integer, intent(in) ::  maxrot, ind
    real*8, intent(in):: x1,y1,z1,x2,y2,z2
    real*8, intent(in), dimension(maxrot,3,3) ::  rotmat

    
    ! output
    ! real*8 :: sqdis

    ! Internal 
    real*8 :: cont,dx,dy,dz
    integer  :: i

    !
    ! Compute component distance vectors and the squared distance:
    !
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx &
                   + rotmat(ind,i,2) * dy &
                   + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
      end do
    return
end function sqdist


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
    real*8,external :: sqdist


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

    hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,maxrot,rotmat)
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
            hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,maxrot,rotmat)
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

end subroutine cova3


subroutine block_covariance(xdb,ydb,zdb, ndb, &
                            nst,it,c0,cc,aa,aa1,aa2,ang1,ang2,ang3, &
                            unbias,cbb)
    !-----------------------------------------------------------------------

    !                    Block Covariance 
    !                    *****************************

    ! This subroutine calculate the block covariance associated with a variogram
    ! model specified by a nugget effect and nested varigoram structures. 
    ! The anisotropy definition can be different for each nested structure.
    ! The block size is defined with input discretization points with 
    ! arbitrary locations (for example regular discretization, random 
    ! discretization or discretization points in an irregular polygon)



    ! INPUT VARIABLES:

    !   xdb,ydb,zdb      coordinates of discretization points
    !   ndb              number of discretization points
    !   nst(ivarg)       number of nested structures (maximum of 4)
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
    !   aa,aa1,aa2       parameters "a" of each nested structure.
    !   aa,aa1,aa2       rotation angles for each nested structure
    ! 
    ! OUTPUT VARIABLES:

    !   unbias           unbias variable, need internally for kriging 
    !   cbb              block covariance



    ! Notes 
    
    ! This functions was created from code lines in kt3d program
    ! 
    !
    ! Adrian Martinez  2015
    !-----------------------------------------------------------------------
    
    IMPLICIT NONE 

    ! in 
    ! variogram
    integer, intent(in)                 :: nst
    real*8, intent(in),   dimension(1)  :: c0
    integer, intent(in), dimension(nst) :: it 
    real*8, intent(in), dimension(nst)  :: cc,aa, aa1, aa2,ang1,ang2,ang3
    !discretization points of the last block
    integer, intent (in) :: ndb   ! number of discretization points
    real*8, intent (in), dimension(ndb) :: xdb,ydb,zdb      !coordinates of discretization points

    !out 
    real*8, intent(out) :: unbias, cbb


    ! internal 
    ! variogram 
    real*8 :: cmax,covmax, cov
    real*8, dimension(nst) :: anis1, anis2
    real*8, dimension(nst,3,3) :: rotmat
    real*8 ::  EPSLON=1.e-20, PMX=999.
    integer :: i, is, j
    
    
    
    do i=1,nst
        anis1(i) = aa1(i) / max(aa(i),EPSLON)
        anis2(i) = aa2(i) / max(aa(i),EPSLON)
        
        print *, aa(i), aa1(i), aa2(i), anis1(i), anis2(i)
        
        if(it(i).eq.4) then
              if(aa(i).lt.0.0) stop ' INVALID power variogram'
              if(aa(i).gt.2.0) stop ' INVALID power variogram'
        end if
    end do
    
    
    
    ! get the rotation matrix
    do is=1,nst
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,nst,rotmat)
        if(it(is) == 4) then
            covmax = covmax + PMX
        else
            covmax = covmax + cc(is)
        endif
    end do

    ! Calculate Block Covariance. Check for point kriging.
    call cova3(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst, &
                c0,it,cc,aa,1,nst,rotmat,cmax,cov)

    ! Set the 'unbias' variable so that the matrix solution is more stable

    unbias = cov
    cbb    = dble(cov)
    if(ndb > 1) then
        cbb = 0.0
        do i=1,ndb
            do j=1,ndb
                call cova3(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j), &
                           1,nst,c0,it,cc,aa,1,nst,rotmat,cmax,cov)
                if(i == j) cov = cov - c0 (1)
                cbb = cbb + dble(cov)
            end do
        end do
        cbb = cbb/dble(ndb*ndb)
    end if

    return
    
end subroutine block_covariance



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


subroutine sortem(ib,ie,a,iperm,b,c,d,e,f,g,h)
    !-----------------------------------------------------------------------

    !                      Quickersort Subroutine
    !                      **********************

    ! This is a subroutine for sorting a real array in ascending order. This
    ! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
    ! in collected algorithms of the ACM.

    ! The method used is that of continually splitting the array into parts
    ! such that all elements of one part are less than all elements of the
    ! other, with a third part in the middle consisting of one element.  An
    ! element with value t is chosen arbitrarily (here we choose the middle
    ! element). i and j give the lower and upper limits of the segment being
    ! split.  After the split a value q will have been found such that
    ! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
    ! performs operations on the two segments (i,q-1) and (q+1,j) as follows
    ! The smaller segment is split and the position of the larger segment is
    ! stored in the lt and ut arrays.  If the segment to be split contains
    ! two or fewer elements, it is sorted and another segment is obtained
    ! from the lt and ut arrays.  When no more segments remain, the array
    ! is completely sorted.


    ! INPUT PARAMETERS:

    !   ib,ie        start and end index of the array to be sorteda
    !   a            array, a portion of which has to be sorted.
    !   iperm        0 no other array is permuted.
    !                1 array b is permuted according to array a
    !                2 arrays b,c are permuted.
    !                3 arrays b,c,d are permuted.
    !                4 arrays b,c,d,e are permuted.
    !                5 arrays b,c,d,e,f are permuted.
    !                6 arrays b,c,d,e,f,g are permuted.
    !                7 arrays b,c,d,e,f,g,h are permuted.
    !               >7 no other array is permuted.

    !   b,c,d,e,f,g,h  arrays to be permuted according to array a.

    ! OUTPUT PARAMETERS:

    !    a      = the array, a portion of which has been sorted.

    !    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)

    ! NO EXTERNAL ROUTINES REQUIRED:

    !-----------------------------------------------------------------------
    real*8 :: a,b,c,d,e,f,g,h
    dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)

    real*8 :: ta,tb,tc,td,te,tf,tg,th
    real*8 :: xa,xb,xc,xd,xe,xf,xg,xh

    ! The dimensions for lt and ut have to be at least log (base 2) n

    integer ::   lt(64),ut(64),i,j,k,m,p,q
 

    ! Initialize:

    j     = ie
    m     = 1
    i     = ib
    iring = iperm+1
    if (iperm > 7) iring=1

    ! If this segment has more than two elements  we split it

    10 if (j-i-1) 100,90,15

    ! p is the position of an arbitrary element in the segment we choose the  
    ! middle element. Under certain circumstances it may be advantageous
    ! to choose p at random.

    15 p    = (j+i)/2
    ta   = a(p)
    a(p) = a(i)
    go to (21,19,18,17,16,161,162,163),iring
    163 th   = h(p)
    h(p) = h(i)
    162 tg   = g(p)
    g(p) = g(i)
    161 tf   = f(p)
    f(p) = f(i)
    16 te   = e(p)
    e(p) = e(i)
    17 td   = d(p)
    d(p) = d(i)
    18 tc   = c(p)
    c(p) = c(i)
    19 tb   = b(p)
    b(p) = b(i)
    21 continue

    ! Start at the beginning of the segment, search for k such that a(k)>t

    q = j
    k = i
    20 k = k+1
    if(k > q)     go to 60
    if(a(k) <= ta) go to 20

    ! Such an element has now been found now search for a q such that a(q)<t
    ! starting at the end of the segment.

    30 continue
    if(a(q) < ta) go to 40
    q = q-1
    if(q > k)     go to 30
    go to 50

    ! a(q) has now been found. we interchange a(q) and a(k)

    40 xa   = a(k)
    a(k) = a(q)
    a(q) = xa
    go to (45,44,43,42,41,411,412,413),iring
    413 xh   = h(k)
    h(k) = h(q)
    h(q) = xh
    412 xg   = g(k)
    g(k) = g(q)
    g(q) = xg
    411 xf   = f(k)
    f(k) = f(q)
    f(q) = xf
    41 xe   = e(k)
    e(k) = e(q)
    e(q) = xe
    42 xd   = d(k)
    d(k) = d(q)
    d(q) = xd
    43 xc   = c(k)
    c(k) = c(q)
    c(q) = xc
    44 xb   = b(k)
    b(k) = b(q)
    b(q) = xb
    45 continue

    ! Update q and search for another pair to interchange:

    q = q-1
    go to 20
    50 q = k-1
    60 continue

    ! The upwards search has now met the downwards search:

    a(i)=a(q)
    a(q)=ta
    go to (65,64,63,62,61,611,612,613),iring
    613 h(i) = h(q)
    h(q) = th
    612 g(i) = g(q)
    g(q) = tg
    611 f(i) = f(q)
    f(q) = tf
    61 e(i) = e(q)
    e(q) = te
    62 d(i) = d(q)
    d(q) = td
    63 c(i) = c(q)
    c(q) = tc
    64 b(i) = b(q)
    b(q) = tb
    65 continue

    ! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
    ! store the position of the largest segment in lt and ut

    if (2*q <= i+j) go to 70
    lt(m) = i
    ut(m) = q-1
    i = q+1
    go to 80
    70 lt(m) = q+1
    ut(m) = j
    j = q-1

    ! Update m and split the new smaller segment

    80 m = m+1
    go to 10

    ! We arrive here if the segment has  two elements we test to see if
    ! the segment is properly ordered if not, we perform an interchange

    90 continue
    if (a(i) <= a(j)) go to 100
    xa=a(i)
    a(i)=a(j)
    a(j)=xa
    go to (95,94,93,92,91,911,912,913),iring
    913 xh   = h(i)
    h(i) = h(j)
    h(j) = xh
    912 xg   = g(i)
    g(i) = g(j)
    g(j) = xg
    911 xf   = f(i)
    f(i) = f(j)
    f(j) = xf
    91 xe   = e(i)
    e(i) = e(j)
    e(j) = xe
    92 xd   = d(i)
    d(i) = d(j)
    d(j) = xd
    93 xc   = c(i)
    c(i) = c(j)
    c(j) = xc
    94 xb   = b(i)
    b(i) = b(j)
    b(j) = xb
    95 continue

    ! If lt and ut contain more segments to be sorted repeat process:

    100 m = m-1
    if (m <= 0) go to 110
    i = lt(m)
    j = ut(m)
    go to 10
    110 continue
    return
end subroutine sortem



!*********************************************************************************
!     Subroutines for GSLIB datafile 
!*********************************************************************************
subroutine read_header(datafl, comment_line, varnames, nvar, error, maxnvar)
    ! read gslib data file header
    ! returns comment line, number of variables and variable names
    ! use this with the functions read_ndata (to get the number of lines in the file)
    ! and read_data to get an array with data
    ! warning: all values are converted to real and character variables are not admissible
    

    ! parameters:
    !--------------------------------
    ! input parameters:  
    !    datafl :  file path 
    !    maxnvar:  maximum number of variables (use a large number, ex. 500 or 500)
    !
    ! output: 
    !    comment_line : string with comments in the datafile
    !    varnames     : string array with varnames    
    !    nvar         : integer with number of variables in the file
    !    error        : integer with error code
    

    implicit none

    integer, intent(in) :: maxnvar    
    character(len=500), intent(in)  :: datafl
    character(len=500), intent(out) :: comment_line
    character, intent(out) , dimension (maxnvar,80) ::  varnames   
            
    integer,   intent(out)  ::  nvar
    
    character(len=80) str

    !error variable
    integer,   intent(out) ::   error
    
    integer i, j, lin

    error=0
    lin=1
    
    ! open the file and read in the header 
    open(unit=lin,file=datafl(1:500), status='old')

    read(lin,'(a500)',err=99) comment_line(1:500)
    read(lin,*,err=99)       nvar
    
    if (nvar > maxnvar) then
        goto 97
    end if
    
    do i=1,nvar
        read(lin,'(a80)',err=99) str
        do j=1,80
            varnames(i,j)= str(j:j)
        end do
    end do

    close(lin)
    
    return
    
    97 error=97 !'nvar > maxnvar'
    close(lin)
    return
    
    99 error=99 !'error in data file!'
    close(lin)
    return
    
end subroutine read_header


subroutine read_ndata(datafl, nvar, maxdat, error)
    ! find out number of rows in a gslib file
    ! use this with the functions read_header (to get nvar and varnames)
    ! and read_data to get an array with data
    ! warning: all values are converted to real and character variables are not admissible

    ! parameters:
    !--------------------------------
    ! input parameters:  
    !    datafl :  file path 
    !    nvar   :  integer with number of variables in the file
    !
    ! output: 
    !    maxdat : integer with number of rows in the datafile. 
    !    error  : integer with error code


    implicit none
  
    character(len=250), intent(in)  :: datafl
    integer,  intent(in)  :: nvar
    integer, intent(out) ::  maxdat
    real*8, dimension (nvar) ::  var       
    
    !error variable
    integer,   intent(out) ::   error
    
    integer i, j, lin

    error=0
    lin=1
    
    ! Open the file and read number of data points 
    open(lin,file=datafl,status='OLD')
    read(lin,*,err=99) !comment line 
    read(lin,*,err=99) !nvar 
    do i=1,nvar        !varnames 
        read(lin,*)
    end do
    
    !initialize number of data to zero
    maxdat = 0
    2 read(lin,*,end=4,err=99) (var(j),j=1,nvar)
        maxdat = maxdat + 1
        go to 2
    4 continue
    
    close(lin)
    
    return
       
    99 error=99 !'ERROR in data file!'
    close(lin)
    return
    
end subroutine read_ndata


subroutine read_data(datafl, nvar, maxdat, table, error)
    ! read data in gslib file
    ! returns array with data (REAL format)
    ! use this with the functions read_header (to get nvar and varnames)
    ! and read_ndata to get the number of rows
    ! warning: all values are converted to real and character variables are not admissible


    ! parameters:
    !--------------------------------
    ! input parameters:  
    !    datafl :  file path 
    !    nvar   :  integer with number of variables in the file
    !    maxdat :  integer with number of rows in the datafile. 
    !
    ! output: 
    !    table    :  real array (maxdat,nvar) with data
    !    varnames :  string array with varnames    
    !    error        : integer with error code


    IMPLICIT NONE
  
    character(len=250), intent(in)  :: datafl
    integer,  intent(in)  :: nvar, maxdat
    real*8, intent(out), dimension (maxdat,nvar) ::  table
    real*8, dimension (nvar) ::  var       
    
    !error variable
    integer,   intent(out) ::   error
    
    integer i, j, lin

    error=0
    lin=1
    
    ! Open the file and read number of data points 
    open(lin,file=datafl,status='OLD')
    read(lin,*,err=99) !comment line 
    read(lin,*,err=99) !nvar 
    do i=1,nvar        !varnames 
        read(lin,*)
    end do


    ! Now, read the data
    do i=1, maxdat
        read(lin,*,end=94,err=99) (var(j),j=1,nvar)
        table(i,:)=var
    end do
    

    close(lin)
    
    return

    94  error=94 !'unexpected end of line'
    close(lin)
    return
    
    99 error=99 !'ERROR in data file!'
    close(lin)
    return
    
end subroutine read_data


!*********************************************************************************
!     Subroutines auxiliary program
!*********************************************************************************
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






!*********************************************************************************
! 
! 
!     Subroutines for GSLIB gam/gamv (variograms)
! 
! 
!*********************************************************************************

!---------------------------------------------------------------------------
!     Subroutine writeout for gamv
!---------------------------------------------------------------------------
subroutine writeout(nvarg, ndir, nlag, np, dis, gam, hm, tm, hv, tv, &
                    pdis, pgam, phm, ptm, phv, ptv, pnump)  

    ! This function is required to transform the output of gam and gamv (1D array)
    ! into a readable 3D array with variogram results
    !
    !                  Write Out the Results of GAMV3
    !                  ******************************
    !
    ! An output will be prepared in variables which contains each directional
    ! variogram. 
    ! the output is: 
    !   pnump()             Number of pairs
    !   pdis()            Distance of pairs falling into this lag
    !   pgam()            Semivariogram, covariance, correlogram,... value
    !   phm()             Mean of the tail data
    !   ptm()             Mean of the head data
    !   phv()             Variance of the tail data
    !   ptv()             Variance of the head data
    !  
    !  the output is in a 3D matrix with dimensions (nvarg, ndir, nlag+2)
    !-----------------------------------------------------------------------
    
    real*8, intent(in), dimension(ndir*(nlag+2)*nvarg)  :: np, dis, gam, hm, tm, hv, tv
    integer, intent(in) :: nvarg, ndir, nlag
    real*8,    intent(out), dimension(nvarg, ndir, nlag+2) :: pdis,pgam, phm,ptm,phv,ptv
    integer, intent(out), dimension(nvarg, ndir, nlag+2) :: pnump
    
    integer :: iv, id, il


    ! Loop over all the variograms that have been computed:

    do iv=1,nvarg
        do id=1,ndir
            do il=1,nlag+2
            
                i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                
                pdis(iv,id,il)=dis(i)
                pgam(iv,id,il)=gam(i)
                phm(iv,id,il)=hm(i)
                ptm(iv,id,il)=tm(i)
                phv(iv,id,il)=hv(i)
                ptv(iv,id,il)=tv(i)
                pnump(iv,id,il)=int(np(i))
                
            end do
        end do
    end do

    return

end subroutine writeout


!---------------------------------------------------------------------------
!     Subroutine gam
!---------------------------------------------------------------------------
subroutine gamv(nd, x, y, z, bhid, nv, vr, &                           ! data array and coordinares
                tmin, tmax, nlag, xlag, xltol,  &                ! lag definition
                ndir, azm, atol, bandwh, dip, dtol, bandwd,  &   ! directions and parameters
                isill, sills, nvarg, ivtail, ivhead,ivtype,  &   ! variograms types and parameters
                np, dis, gam, hm, tm, hv, tv, &                  ! output
                cldi, cldj, cldg, cldh, l, maxclp )              ! extra parameters for cloud veriogram

    !----------------------------------------------------------------------
    ! This code was modified from original f77 GSLIB code (v.2)
    ! Mayor changes
    ! a) Fixed lenght arrays are now externally defined (example in python)
    ! b) The fixed parameters were replaced by new variables as follow
    !   MAXDAT => nd     :  number of data points
    !   MAXVAR => nv     :  number of variables
    !   MAXDIR => ndir   :  number of directions possible at one time
    !   MAXLAG => nlag*2 :  number of lags at one time
    !   MXVARG => nvarg  :  number of variograms possible at one time

    !   MXDLV  =>  ndir*(nlag*2)*nvarg 
    
    ! c) Support for downhole variogram was added (see bhid(nd)). To
    !    ignore downhole option use bhid(:)= constant
    !
    !     Comment: you can use this to avoid mix of data (example variograms)
    !              at both sides of a fault. 
    !
    ! d) Variogram cloud was added for the first direction and first variogram. 
    !    The number of variogram cloud required may be calculated externally 
    !    and may be equal to sum(np) in first direction of the first variogram.
    !    For this function the following variables were added
    !       maxclp : max number of points (this is sum(np(varg1,dir1, :)))
    !       cldi : point i ID, array('i') with bounds (maxclp)
    !       cldj : point j ID, array('i') with bounds (maxclp)
    !       cldg : pair var/covar value, array('f') with bounds (maxclp)
    !       cldd : distance i-j, array('f') with bounds (maxclp)
    !       l    : number of variogram cloud points calculated,  int
    !  
    ! Warning: The automatic implementation of indicator is not implemented here
    !          you may calculate the indicators variables externally  
    !
    ! TODO: -Add extra output parameters to plot box and whiskers per lag.
    !        Eventually this can be calculated with the varigram cloud output. 
    !----------------------------------------------------------------------
    

    !----------------------------------------------------------------------
    
    !              Variogram of 3-D Irregularly Spaced Data
    !              ****************************************
    ! This subroutine computes a variety of spatial continuity measures of a
    ! set for irregularly spaced data.  The user can specify any combination
    ! of direct and cross variograms using any of eight "variogram" measures
    !
    !
    ! INPUT VARIABLES:
    !
    !   nd               Number of data (no missing values)
    !   x(nd)            X coordinates of the data
    !   y(nd)            Y coordinates of the data
    !   z(nd)            Z coordinates of the data
    !   bhid(nd)         bhid ID (integer) to force downhole variogram
    !   nv               The number of variables
    !   vr(nd,nv)        Data values
    !   tmin,tmax        Trimming limits
    !   nlag             Number of lags to calculate
    !   xlag             Length of the unit lag
    !   xltol            Distance tolerance (if <0 then set to xlag/2)
    !   ndir             Number of directions to consider
    !   azm(ndir)        Azimuth angle of direction (measured positive
    !                      degrees clockwise from NS).
    !   atol(ndir)       Azimuth (half window) tolerances
    !   bandwh(ndir)     Maximum Horizontal bandwidth (i.e., the deviation
    !                      perpendicular to the defined azimuth).
    !   dip(ndir)        Dip angle of direction (measured in negative
    !                      degrees down from horizontal).
    !   dtol(ndir)       Dip (half window) tolerances
    !   bandwd(ndir)     Maximum "Vertical" bandwidth (i.e., the deviation
    !                      perpendicular to the defined dip).
    !   isill            1=attempt to standardize, 0=do not
    !   sills            the sills (variances) to standardize with
    !   nvarg            Number of variograms to compute
    !   ivtail(nvarg)    Variable for the tail of each variogram
    !   ivhead(nvarg)    Variable for the head of each variogram
    !   ivtype(nvarg)    Type of variogram to compute:
    !                      1. semivariogram
    !                      2. cross-semivariogram
    !                      3. covariance
    !                      4. correlogram
    !                      5. general relative semivariogram
    !                      6. pairwise relative semivariogram
    !                      7. semivariogram of logarithms
    !                      8. rodogram
    !                      9. indicator semivariogram (continuous)
    !                     10. indicator semivariogram (categorical)
    !
    !
    !
    ! OUTPUT VARIABLES:  The following arrays are ordered by direction,
    !                    then variogram, and finally lag, i.e.,
    !                      iloc = (id-1)*nvarg*MAXLG+(iv-1)*MAXLG+il
    !
    !   np()             Number of pairs
    !   dis()            Distance of pairs falling into this lag
    !   gam()            Semivariogram, covariance, correlogram,... value
    !   hm()             Mean of the tail data
    !   tm()             Mean of the head data
    !   hv()             Variance of the tail data
    !   tv()             Variance of the head data
    !
    !
    ! PROGRAM NOTES:
    !
    !   1. The file "gamv.inc" contains the dimensioning parameters.
    !      These may have to be changed depending on the amount of data
    !      and the requirements to compute different variograms.
    ! Original:  A.G. Journel                                           1978
    ! Revisions: K. Guertin                                             1980
    !-----------------------------------------------------------------------

    !for safety reason we don't want undeclared variables
    IMPLICIT NONE    

    integer, intent(in)                   ::nd, nv, ndir, nvarg, nlag, isill
    integer, intent(in), dimension(nvarg) :: ivtail, ivhead,ivtype
    real*8, intent(in), dimension(nd)       :: x, y, z
    
    real*8, intent(in), dimension(nd,nv)    :: vr
    real*8, intent(in), dimension(nv)       :: sills
    real*8, intent(in)                      :: tmin, tmax, xlag 
    real*8, intent(in)                      :: xltol
    real*8, intent(in), dimension(ndir)     :: azm, atol, bandwh, dip, dtol, bandwd
    real*8, intent(out), dimension(ndir*(nlag+2)*nvarg)  :: np, dis, gam, hm, tm, hv, tv

    !new variables
    integer, intent(in), dimension(nd)         :: bhid
    integer, intent(in)                        :: maxclp
    integer, intent(out), dimension(maxclp)    :: cldi,cldj
    real*8, intent(out), dimension(maxclp)     :: cldg, cldh
    integer, intent(out)                       :: l
               
    ! some general declarations
    real*8 :: PI, EPSLON
    parameter(PI=3.14159265, EPSLON=1.0e-20)
    
    real*8, dimension(ndir) :: uvxazm, uvyazm, uvzdec, uvhdec, csatol, csdtol
    logical               :: omni


    !Extra Variables not declared in the original library
    real*8  :: azmuth, declin, band, dcazm, dismxs, dx, dxs, dxy, dy, &
             dys, dz, dzs, dcdec, gamma, h, hs, htave, variance, & 
             vrh, vrhpr, vrt,  vrtpr, xltoll
    integer :: i, id, ii, iii, il, ilag, it, iv, j, lagbeg, lagend, nsiz, rnum
    
    
    !initialize counter for variogram cloud
    l=0
    

    
    !here we rename the variable xltol to xtoll to avoid reference conflict in (real, intent(in) :: xltol)
    xltoll=xltol

    ! Define the distance tolerance if it isn't already:

    if(xltoll <= 0.0) xltoll = 0.5 * xlag

    ! Define the angles and tolerance for each direction:

    do id=1,ndir
    
    ! The mathematical azimuth is measured counterclockwise from EW and
    ! not clockwise from NS as the conventional azimuth is:
    
        azmuth     = (90.0-azm(id))*PI/180.0
        uvxazm(id) = cos(azmuth)
        uvyazm(id) = sin(azmuth)
        if(atol(id) <= 0.0) then
            csatol(id) = cos(45.0*PI/180.0)
        else
            csatol(id) = cos(atol(id)*PI/180.0)
        endif
    
    ! The declination is measured positive down from vertical (up) rather
    ! than negative down from horizontal:
    
        declin     = (90.0-dip(id))*PI/180.0
        uvzdec(id) = cos(declin)
        uvhdec(id) = sin(declin)
        if(dtol(id) <= 0.0) then
            csdtol(id) = cos(45.0*PI/180.0)
        else
            csdtol(id) = cos(dtol(id)*PI/180.0)
        endif
    end do

! Initialize the arrays for each direction, variogram, and lag:

    nsiz = ndir*nvarg*(nlag+2)
    do i=1,nsiz
        np(i)  = 0.
        dis(i) = 0.0
        gam(i) = 0.0
        hm(i)  = 0.0
        tm(i)  = 0.0
        hv(i)  = 0.0
        tv(i)  = 0.0
    end do
    dismxs = ((dble(nlag) + 0.5 - EPSLON) * xlag) ** 2  
    
! MAIN LOOP OVER ALL PAIRS:

    do 3 i=1,nd
        do 4 j=i,nd
        
        !implement the downhole variogram  (change added on may 2015)
            if(bhid(j).NE.bhid(i)) go to 4
        
        ! Definition of the lag corresponding to the current pair:
        
            dx  = x(j) - x(i)
            dy  = y(j) - y(i)
            dz  = z(j) - z(i)
            dxs = dx*dx
            dys = dy*dy
            dzs = dz*dz
            hs  = dxs + dys + dzs
            if(hs > dismxs) go to 4
            if(hs < 0.0) hs = 0.0
            h   = sqrt(hs)
        
        ! Determine which lag this is and skip if outside the defined distance
        ! tolerance:
        
            if(h <= EPSLON) then
                lagbeg = 1
                lagend = 1
            else
                lagbeg = -1
                lagend = -1
                do ilag=2,nlag+2
                    if(h >= (xlag*dble(ilag-2)-xltoll) .AND. &
                    h <= (xlag*dble(ilag-2)+xltoll)) then
                        if(lagbeg < 0) lagbeg = ilag
                        lagend = ilag
                    end if
                end do
                if(lagend < 0) go to 4
            endif
        
        ! Definition of the direction corresponding to the current pair. All
        ! directions are considered (overlapping of direction tolerance cones
        ! is allowed):
        
            do 5 id=1,ndir
            
            ! Check for an acceptable azimuth angle:
            
                dxy = sqrt(max((dxs+dys),0.0))
                if(dxy < EPSLON) then
                    dcazm = 1.0
                else
                    dcazm = (dx*uvxazm(id)+dy*uvyazm(id))/dxy
                endif
                if(abs(dcazm) < csatol(id)) go to 5
            
            ! Check the horizontal bandwidth criteria (maximum deviation
            ! perpendicular to the specified direction azimuth):
            
                band = uvxazm(id)*dy - uvyazm(id)*dx
                if(abs(band) > bandwh(id)) go to 5
            
            ! Check for an acceptable dip angle:
            
                if(dcazm < 0.0) dxy = -dxy
                if(lagbeg == 1) then
                    dcdec = 0.0
                else
                    dcdec = (dxy*uvhdec(id)+dz*uvzdec(id))/h
                    if(abs(dcdec) < csdtol(id)) go to 5
                endif
            
            ! Check the vertical bandwidth criteria (maximum deviation perpendicular
            ! to the specified dip direction):
            
                band = uvhdec(id)*dz - uvzdec(id)*dxy
                if(abs(band) > bandwd(id)) go to 5
            
            ! Check whether or not an omni-directional variogram is being computed:
            
                omni = .FALSE. 
                if(atol(id) >= 90.0) omni = .TRUE. 
            
            ! This direction is acceptable - go ahead and compute all variograms:
            
                do 6 iv=1,nvarg
                
                ! For this variogram, sort out which is the tail and the head value:
                
                    it = ivtype(iv)
                    if(dcazm >= 0.0 .AND. dcdec >= 0.0) then
                        ii = ivtail(iv)
                        vrh   = vr(i,ii)
                        ii = ivhead(iv)
                        vrt   = vr(j,ii)
                        if(omni .OR. it == 2) then
                            ii    = ivhead(iv)
                            vrtpr = vr(i,ii)
                            ii    = ivtail(iv)
                            vrhpr = vr(j,ii)
                        endif
                    else
                        ii = ivtail(iv)
                        vrh   = vr(j,ii)
                        ii = ivhead(iv)
                        vrt   = vr(i,ii)
                        if(omni .OR. it == 2) then
                            ii    = ivhead(iv)
                            vrtpr = vr(j,ii)
                            ii    = ivtail(iv)
                            vrhpr = vr(i,ii)
                        endif
                    endif
                
                ! Reject this pair on the basis of missing values:
                
                    if(vrt < tmin .OR. vrh < tmin .OR. &
                    vrt > tmax .OR. vrh > tmax) go to 6
                    if(it == 2 .AND. (vrtpr < tmin .OR. vrhpr < tmin .OR. &
                    vrtpr > tmax .OR. vrhpr > tmax)) &
                    go to 6
                
                !             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
                
                
                ! The Semivariogram:
                
                    if(it == 1 .OR. it == 5 .OR. it >= 9) then
                        do il=lagbeg,lagend
                            ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                            np(ii)  = np(ii)  + 1.
                            dis(ii) = dis(ii) + dble(h)
                            tm(ii)  = tm(ii)  + dble(vrt)
                            hm(ii)  = hm(ii)  + dble(vrh)
                            gam(ii) = gam(ii) + dble((vrh-vrt)*(vrh-vrt))
                            if(omni) then
                                if(vrtpr >= tmin .AND. vrhpr >= tmin .AND. &
                                vrtpr < tmax .AND. vrhpr < tmax) then
                                    np(ii)  = np(ii)  + 1.
                                    dis(ii) = dis(ii) + dble(h)
                                    tm(ii)  = tm(ii)  + dble(vrtpr)
                                    hm(ii)  = hm(ii)  + dble(vrhpr)
                                    gam(ii) = gam(ii) + dble((vrhpr-vrtpr)* &
                                    (vrhpr-vrtpr))
                                endif
                            endif
                        end do
                    
                    ! The Traditional Cross Semivariogram:
                    
                    else if(it == 2) then
                        do il=lagbeg,lagend
                            ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                            np(ii)  = np(ii)  + 1.
                            dis(ii) = dis(ii) + dble(h)
                            tm(ii)  = tm(ii)  + dble(0.5*(vrt+vrtpr))
                            hm(ii)  = hm(ii)  + dble(0.5*(vrh+vrhpr))
                            gam(ii) = gam(ii) + dble((vrhpr-vrh)*(vrt-vrtpr))
                        end do
                    
                    ! The Covariance:
                    
                    else if(abs(it) == 3) then
                        do il=lagbeg,lagend
                            ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                            np(ii)  = np(ii)  + 1.
                            dis(ii) = dis(ii) + dble(h)
                            tm(ii)  = tm(ii)  + dble(vrt)
                            hm(ii)  = hm(ii)  + dble(vrh)
                            gam(ii) = gam(ii) + dble(vrh*vrt)
                            if(omni) then
                                if(vrtpr >= tmin .AND. vrhpr >= tmin .AND. &
                                vrtpr < tmax .AND. vrhpr < tmax) then
                                    np(ii)  = np(ii)  + 1.
                                    dis(ii) = dis(ii) + dble(h)
                                    tm(ii)  = tm(ii)  + dble(vrtpr)
                                    hm(ii)  = hm(ii)  + dble(vrhpr)
                                    gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
                                endif
                            endif
                        end do
                    
                    ! The Correlogram:
                    
                    else if(it == 4) then
                        do il=lagbeg,lagend
                            ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                            np(ii)  = np(ii)  + 1.
                            dis(ii) = dis(ii) + dble(h)
                            tm(ii)  = tm(ii)  + dble(vrt)
                            hm(ii)  = hm(ii)  + dble(vrh)
                            hv(ii)  = hv(ii)  + dble(vrh*vrh)
                            tv(ii)  = tv(ii)  + dble(vrt*vrt)
                            gam(ii) = gam(ii) + dble(vrh*vrt)
                            if(omni) then
                                if(vrtpr >= tmin .AND. vrhpr >= tmin .AND. &
                                vrtpr < tmax .AND. vrhpr < tmax) then
                                    np(ii)  = np(ii)  + 1.
                                    dis(ii) = dis(ii) + dble(h)
                                    tm(ii)  = tm(ii)  + dble(vrtpr)
                                    hm(ii)  = hm(ii)  + dble(vrhpr)
                                    hv(ii)  = hv(ii)  + dble(vrhpr*vrhpr)
                                    tv(ii)  = tv(ii)  + dble(vrtpr*vrtpr)
                                    gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
                                endif
                            endif
                        end do
                    
                    ! The Pairwise Relative:
                    
                    else if(it == 6) then
                        do il=lagbeg,lagend
                            ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                            if(abs(vrt+vrh) > EPSLON) then
                                np(ii)  = np(ii)  + 1.
                                dis(ii) = dis(ii) + dble(h)
                                tm(ii)  = tm(ii)  + dble(vrt)
                                hm(ii)  = hm(ii)  + dble(vrh)
                                gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                                gam(ii) = gam(ii) + dble(gamma*gamma)
                            endif
                            if(omni) then
                                if(vrtpr >= tmin .AND. vrhpr >= tmin .AND. &
                                vrtpr < tmax .AND. vrhpr < tmax) then
                                    if(abs(vrtpr+vrhpr) > EPSLON) then
                                        np(ii)  = np(ii)  + 1.
                                        dis(ii) = dis(ii) + dble(h)
                                        tm(ii)  = tm(ii)  + dble(vrtpr)
                                        hm(ii)  = hm(ii)  + dble(vrhpr)
                                        gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                                        gam(ii) = gam(ii) + dble(gamma*gamma)
                                    endif
                                endif
                            endif
                        enddo
                    
                    ! Variogram of Logarithms:
                    
                    else if(it == 7) then
                        do il=lagbeg,lagend
                            ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                            if(vrt > EPSLON .AND. vrh > EPSLON) then
                                np(ii)  = np(ii)  + 1.
                                dis(ii) = dis(ii) + dble(h)
                                tm(ii)  = tm(ii)  + dble(vrt)
                                hm(ii)  = hm(ii)  + dble(vrh)
                                gamma   = dlog(vrt)-dlog(vrh)
                                gam(ii) = gam(ii) + dble(gamma*gamma)
                            endif
                            if(omni) then
                                if(vrtpr >= tmin .AND. vrhpr >= tmin .AND. &
                                vrtpr < tmax .AND. vrhpr < tmax) then
                                    if(vrtpr > EPSLON .AND. vrhpr > EPSLON) then
                                        np(ii)  = np(ii)  + 1.
                                        dis(ii) = dis(ii) + dble(h)
                                        tm(ii)  = tm(ii)  + dble(vrtpr)
                                        hm(ii)  = hm(ii)  + dble(vrhpr)
                                        gamma   = dlog(vrt)-dlog(vrh)
                                        gam(ii) = gam(ii) + dble(gamma*gamma)
                                    endif
                                endif
                            endif
                        end do
                    
                    ! Madogram:
                    
                    else if(it == 8) then
                        do il=lagbeg,lagend
                            ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                            np(ii)  = np(ii)  + 1.
                            dis(ii) = dis(ii) + dble(h)
                            tm(ii)  = tm(ii)  + dble(vrt)
                            hm(ii)  = hm(ii)  + dble(vrh)
                            gam(ii) = gam(ii) + dble(abs(vrh-vrt))
                            if(omni) then
                                if(vrtpr >= tmin .AND. vrhpr >= tmin .AND. &
                                vrtpr < tmax .AND. vrhpr < tmax) then
                                    np(ii)  = np(ii)  + 1.
                                    dis(ii) = dis(ii) + dble(h)
                                    tm(ii)  = tm(ii)  + dble(vrtpr)
                                    hm(ii)  = hm(ii)  + dble(vrhpr)
                                    gam(ii) = gam(ii) + dble(abs(vrhpr-vrtpr))
                                endif
                            endif
                        end do
                    endif
                
                ! Finish loops over variograms, directions, and the double data loops:
                    
                ! Calculate variogram cloud for the first maxclp points   (change added on may 2015)

                    ! this is calculated only in the first variogram/first direction 
                    if (iv==1 .AND. id==1) then
                        if (l<=maxclp) then
                            !  report the semivariogram rather than variogram
                            if(it == 1 .OR. it == 2) then
                                l=l+1
                                cldi(l)=i
                                cldj(l)=j
                                cldh(l)=dble(h)
                                cldg(l) = 0.5 * gam(ii)

                            !these are pairwise relative, semivariogram of logarithms
                            !semimadogram and indicators. 
                            else if(it >= 6) then
                                l=l+1
                                cldi(l)=i
                                cldj(l)=j
                                cldh(l)=dble(h)
                                cldg(i) = 0.5 * gam(ii)
                            endif
                            !TODO implement it for other variogram types (3,4,5)?
                        endif
                    endif
                    
                    
                6 END DO
            5 END DO
        4 END DO
    3 END DO

! Get average values for gam, hm, tm, hv, and tv, then compute
! the correct "variogram" measure:

    do 7 id=1,ndir
        do 7 iv=1,nvarg
            do 7 il=1,nlag+2
                i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                if(np(i) <= 0.) go to 7
                rnum   = np(i)
                dis(i) = dis(i) / dble(rnum)
                gam(i) = gam(i) / dble(rnum)
                hm(i)  = hm(i)  / dble(rnum)
                tm(i)  = tm(i)  / dble(rnum)
                hv(i)  = hv(i)  / dble(rnum)
                tv(i)  = tv(i)  / dble(rnum)
                it = ivtype(iv)
            
            ! Attempt to standardize:
            
                if(isill == 1) then
                    if(ivtail(iv) == ivhead(iv)) then
                        iii = ivtail(iv)
                        if((it == 1 .OR. it >= 9) .AND. sills(iii) > 0.0) &
                        gam(i) = gam(i) / sills(iii)
                    end if
                end if
            
            !  1. report the semivariogram rather than variogram
            !  2. report the cross-semivariogram rather than variogram
            !  3. the covariance requires "centering"
            !  4. the correlogram requires centering and normalizing
            !  5. general relative requires division by lag mean
            !  6. report the semi(pairwise relative variogram)
            !  7. report the semi(log variogram)
            !  8. report the semi(madogram)
            
                if(it == 1 .OR. it == 2) then
                    gam(i) = 0.5 * gam(i)
                else if(abs(it) == 3) then
                    gam(i) = gam(i) - hm(i)*tm(i)
                    if(it < 0) then
                        if(sills(ivtail(iv)) < 0.0 .OR. &
                        sills(ivhead(iv)) < 0.0) then
                            gam(i) = -999.0
                        else
                            variance = ( sqrt(sills(ivtail(iv))) &
                            *   sqrt(sills(ivhead(iv))) )
                            gam(i) = variance - gam(i)
                        end if
                    end if
                else if(it == 4) then
                    hv(i)  = hv(i)-hm(i)*hm(i)
                    if(hv(i) < 0.0) hv(i) = 0.0
                    hv(i)  = dsqrt(hv(i))
                    tv(i)  = tv(i)-tm(i)*tm(i)
                    if(tv(i) < 0.0) tv(i) = 0.0
                    tv(i)  = dsqrt(tv(i))
                    if((hv(i)*tv(i)) < EPSLON) then
                        gam(i) = 0.0
                    else
                        gam(i) =(gam(i)-hm(i)*tm(i))/(hv(i)*tv(i))
                    endif
                
                ! Square "hv" and "tv" so that we return the variance:
                
                    hv(i)  = hv(i)*hv(i)
                    tv(i)  = tv(i)*tv(i)
                else if(it == 5) then
                    htave  = 0.5*(hm(i)+tm(i))
                    htave  = htave   *   htave
                    if(htave < EPSLON) then
                        gam(i) = 0.0
                    else
                        gam(i) = gam(i)/dble(htave)
                    endif
                else if(it >= 6) then
                    gam(i) = 0.5 * gam(i)
                endif
    7 END DO
    return


 end subroutine gamv    

!---------------------------------------------------------------------------
!     Subroutine gamma (this is gam program, variogram for data in grid)
!---------------------------------------------------------------------------

subroutine gamma(nx, ny, nz, bhid, nv, vr, &                   ! data array and coordinares
                tmin, tmax, nlag,   &                            ! lag definition
                ndir, ixd, iyd, izd,&                              ! directions and parameters
                isill, sills, nvarg, ivtail, ivhead,ivtype,  &   ! variograms types and parameters
                xsiz,ysiz,zsiz, &                                ! WARNING TODO dummy variable not used here, remove? 
                np, gam, hm, tm, hv, tv)                           ! output


    !-----------------------------------------------------------------------

    !                Variogram of Data on a Regular Grid
    !                ***********************************

    ! This subroutine computes any of eight different measures of spatial
    ! continuity for regular spaced 3-D data.  Missing values are allowed
    ! and the grid need not be cubic.



    ! INPUT VARIABLES:

    !   nlag             Maximum number of lags to be calculated
    !   nx               Number of units in x (number of columns)
    !   ny               Number of units in y (number of lines)
    !   nz               Number of units in z (number of levels)
    !   ndir             Number of directions to consider
    !   ixd(ndir)        X (column) indicator of direction - number of grid
    !                      columns that must be shifted to move from a node
    !                      on the grid to the next nearest node on the grid
    !                      which lies on the directional vector
    !   iyd(ndir)        Y (line) indicator of direction - similar to ixd,
    !                      number of grid lines that must be shifted to
    !                      nearest node which lies on the directional vector
    !   izd(ndir)        Z (level) indicator of direction - similar to ixd,
    !                      number of grid levels that must be shifted to
    !                      nearest node of directional vector
    !   nv               The number of variables
    !   vr(nx*ny*nz*nv)  Array of data
    !   tmin,tmax        Trimming limits
    !   isill            1=attempt to standardize, 0=do not
    !   sills            the sills (variances) to standardize with
    !   nvarg            Number of variograms to compute
    !   ivtail(nvarg)    Variable for the tail of the variogram
    !   ivhead(nvarg)    Variable for the head of the variogram
    !   ivtype(nvarg)    Type of variogram to compute:
    !                      1. semivariogram
    !                      2. cross-semivariogram
    !                      3. covariance
    !                      4. correlogram
    !                      5. general relative semivariogram
    !                      6. pairwise relative semivariogram
    !                      7. semivariogram of logarithms
    !                      8. madogram
    !                      9. indicator semivariogram: an indicator variable
    !                         is constructed in the main program.

    ! OUTPUT VARIABLES:  The following arrays are ordered by direction,
    !                    then variogram, and finally lag, i.e.,
    !                      iloc = (id-1)*nvarg*nlag+(iv-1)*nlag+il

    !   np()             Number of pairs
    !   gam()            Semivariogram, covariance, correlogram,... value
    !   hm()             Mean of the tail data
    !   tm()             Mean of the head data
    !   hv()             Variance of the tail data
    !   tv()             Variance of the head data



    ! Original:  A.G. Journel                                           1978
    ! Revisions: B.E. Buxton                                       Apr. 1982
    !-----------------------------------------------------------------------


    !for safety reason we don't want undeclared variables
    IMPLICIT NONE    

    integer, intent(in)                   :: nx, ny, nz,  nv, ndir, nvarg, nlag, isill
    integer, intent(in), dimension(nvarg) :: ivtail, ivhead,ivtype
    integer, intent(in), dimension(ndir)  :: ixd, iyd, izd
    
    real*8, intent(in), dimension(nx*ny*nz*nv) :: vr
    real*8, intent(in), dimension(nv)       :: sills
    real*8, intent(in)                      :: tmin, tmax, xsiz,ysiz,zsiz    
    real*8, intent(out), dimension(ndir*(nlag+2)*nvarg)  :: np, gam, hm, tm, hv, tv

    ! TODO: np is real here, see if we may declare this as integer

    !new variables
    integer, intent(in), dimension(nx*ny*nz)         :: bhid             ! not implemented (use similar to gamv)

               
    ! some general declarations
    real*8 :: PI, EPSLON
    parameter(PI=3.14159265, EPSLON=1.0e-20)
    

    !Extra Variables not declared in the original library
    real*8  ::  htave, variance, vrh, vrhpr, vrt,  vrtpr, tempvar
    integer :: i, id, iii, il, it, iv, nsiz, rnum
    integer :: ixinc, iyinc, izinc
    integer :: ix1, iy1, iz1, ix, iy, iz, index, nxy,nxyz


!moved from line 154, original file gam.f 
    nxy  = nx * ny
    nxyz = nx * ny * nz

               
! Initialize the summation arrays for each direction, variogram, and lag

    nsiz = ndir*nvarg*nlag
    do i=1,nsiz
        np(i)  = 0.
        gam(i) = 0.0
        hm(i)  = 0.0
        tm(i)  = 0.0
        hv(i)  = 0.0
        tv(i)  = 0.0
    end do

! First fix the location of a seed point on the grid (ix,iy,iz):

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz
            
            ! For the fixed seed point, loop through all directions:
            
                do id=1,ndir
                    ixinc = ixd(id)
                    iyinc = iyd(id)
                    izinc = izd(id)
                    ix1   = ix
                    iy1   = iy
                    iz1   = iz
                
                ! For this direction, loop over all the lags:
                
                    do il=1,nlag
                    
                    ! Check to be sure that the point being considered is still in the
                    ! grid - if not, then finished with this direction:
                    
                        ix1 = ix1 + ixinc
                        if(ix1 < 1 .OR. ix1 > nx) go to 3
                        iy1 = iy1 + iyinc
                        if(iy1 < 1 .OR. iy1 > ny) go to 3
                        iz1 = iz1 + izinc
                        if(iz1 < 1 .OR. iz1 > nz) go to 3
                    
                    ! For this direction and lag, loop over all variograms:
                    
                        do iv=1,nvarg
                            it = ivtype(iv)
                        
                        ! Get the head value, skip this value if missing:
                        
                            i     = ivhead(iv)
                            index = ix+(iy-1)*nx+(iz-1)*nxy+(i-1)*nxyz
                            vrt   = vr(index)
                            if(vrt < tmin .OR. vrt >= tmax) go to 5
                        
                        ! Get the tail value, skip this value if missing:
                        
                            i     = ivtail(iv)
                            index = ix1+(iy1-1)*nx+(iz1-1)*nxy+(i-1)*nxyz
                            vrh   = vr(index)
                            if(vrh < tmin .OR. vrh >= tmax) go to 5
                        
                        ! Need increment for the cross semivariogram only:
                        
                            if(it == 2) then
                                i     = ivtail(iv)
                                index = ix+(iy-1)*nx+(iz-1)*nxy+(i-1)*nxyz
                                vrhpr = vr(index)
                                if(vrhpr < tmin .OR. vrhpr >= tmax) go to 5
                                i     = ivhead(iv)
                                index = ix1+(iy1-1)*nx+(iz1-1)*nxy+(i-1)*nxyz
                                vrtpr = vr(index)
                                if(vrtpr < tmin .OR. vrtpr >= tmax) go to 5
                            endif
                        
                        ! We have an acceptable pair, therefore accumulate all the statistics
                        ! that are required for the variogram:
                        
                            i      = (id-1)*nvarg*nlag+(iv-1)*nlag+il
                            np(i)  = np(i) + 1.
                            tm(i)  = tm(i) + dble(vrt)
                            hm(i)  = hm(i) + dble(vrh)
                        
                        ! Choose the correct variogram type and keep relevant sums:
                        
                            if(it == 1 .OR. it >= 9) then
                                gam(i) = gam(i) + dble((vrh-vrt)*(vrh-vrt))
                            else if(it == 2) then
                                gam(i) = gam(i) + dble((vrhpr-vrh)*(vrt-vrtpr))
                            else if(abs(it) == 3) then
                                gam(i) = gam(i) +  dble(vrh*vrt)
                            else if(it == 4) then
                                gam(i) = gam(i) +  dble(vrh*vrt)
                                hv(i)  = hv(i)  +  dble(vrh*vrh)
                                tv(i)  = tv(i)  +  dble(vrt*vrt)
                            else if(it == 5) then
                                gam(i) = gam(i) + dble((vrh-vrt)*(vrh-vrt))
                            else if(it == 6) then
                                if((vrt+vrh) < EPSLON) then
                                    np(i)  = np(i) - 1.
                                    tm(i)  = tm(i) - dble(vrt)
                                    hm(i)  = hm(i) - dble(vrh)
                                else
                                    tempvar= 2.0*(vrt-vrh)/(vrt+vrh)
                                    gam(i) = gam(i) + dble(tempvar*tempvar)
                                endif
                            else if(it == 7) then
                                if(vrt < EPSLON .OR. vrh < EPSLON) then
                                    np(i)  = np(i) - 1.
                                    tm(i)  = tm(i) - dble(vrt)
                                    hm(i)  = hm(i) - dble(vrh)
                                else
                                    tempvar= dlog(vrt)-dlog(vrh)
                                    gam(i) = gam(i) + dble(tempvar*tempvar)
                                endif
                            else if(it == 8) then
                                gam(i) = gam(i) + dble(abs(vrt-vrh))
                            endif
                            5 continue
                        end do
                        4 continue
                    end do
                    3 continue
                end do
            end do
        end do
    end do

! Get average values for gam, hm, tm, hv, and tv, then compute
! the correct "variogram" measure:

    do id=1,ndir
        do iv=1,nvarg
            do il=1,nlag
                i = (id-1)*nvarg*nlag+(iv-1)*nlag+il
                if(np(i) == 0.) go to 6
                rnum   = int(np(i))
                gam(i) = gam(i) / dble(rnum)
                hm(i)  = hm(i)  / dble(rnum)
                tm(i)  = tm(i)  / dble(rnum)
                hv(i)  = hv(i)  / dble(rnum)
                tv(i)  = tv(i)  / dble(rnum)
                it     = ivtype(iv)
            
            ! Attempt to standardize:
            
                if(isill == 1) then
                    if(ivtail(iv) == ivhead(iv)) then
                        iii = ivtail(iv)
                        if((it == 1 .OR. it >= 9) .AND. sills(iii) > 0.0) &
                        gam(i) = gam(i) / sills(iii)
                    end if
                end if
            
            ! 1. report the semivariogram rather than variogram
            ! 2. report the cross-semivariogram rather than variogram
            ! 3. the covariance requires "centering"
            ! 4. the correlogram requires centering and normalizing
            ! 5. general relative requires division by lag mean
            ! 6. report the semi(pairwise relative variogram)
            ! 7. report the semi(log variogram)
            ! 8. report the semi(madogram)
            
                if(it == 1 .OR. it == 2) then
                    gam(i) = 0.5 * gam(i)
                else if(abs(it) == 3) then
                    gam(i) = gam(i) - hm(i)*tm(i)
                    if(it < 0) then
                        if(sills(ivtail(iv)) < 0.0 .OR. &
                        sills(ivhead(iv)) < 0.0) then
                            gam(i) = -999.0
                        else
                            variance = ( sqrt(sills(ivtail(iv))) &
                            *   sqrt(sills(ivhead(iv))) )
                            gam(i) = variance - gam(i)
                        end if
                    end if
                else if(it == 4) then
                    hv(i)  = hv(i)-hm(i)*hm(i)
                    if(hv(i) <= 0.0) hv(i) = 0.0
                    hv(i)  = sqrt(hv(i))
                    tv(i)  = tv(i)-tm(i)*tm(i)
                    if(tv(i) <= 0.0) tv(i) = 0.0
                    tv(i)  = sqrt(tv(i))
                    if((hv(i)*tv(i)) < EPSLON) then
                        gam(i) = 0.0
                    else
                        gam(i) =(gam(i)-hm(i)*tm(i))/(hv(i)*tv(i))
                    endif
                
                ! Square "hv" and "tv" so that we return the variance:
                
                    hv(i)  = hv(i)*hv(i)
                    tv(i)  = tv(i)*tv(i)
                else if(it == 5) then
                    htave  = 0.5*(hm(i)+tm(i))
                    htave  = htave   *   htave
                    if(htave < EPSLON) then
                        gam(i) = 0.0
                    else
                        gam(i) = gam(i)/dble(htave)
                    endif
                else if(it >= 6) then
                    gam(i) = 0.5 * gam(i)
                endif
                6 continue
            end do
        end do
    end do
    return
end subroutine gamma


!-----------------------------------------------------------------------
!     Subroutine write out for gam
!-----------------------------------------------------------------------
subroutine writeout_gam (nvarg, ndir, nlag, ixd, xsiz, iyd, ysiz, &
                         izd, zsiz, np, gam, hm, tm, hv, tv, & 
                         pdis, pgam, phm, ptm, phv, ptv, pnump)

    !-------------------------------------------------------------------

    !                  Write Out the Results of GAM
    !                  ****************************

    ! An output file will be written which contains each directional
    ! variogram ordered by direction and then variogram (the directions
    ! cycle fastest then the variogram number). For each variogram there
    ! will be a one line description and then "nlag" lines with:

    !        a) lag number (increasing from 1 to nlag)
    !        b) separation distance
    !        c) the "variogram" value
    !        d) the number of pairs for the lag
    !        e) the mean of the data contributing to the tail
    !        f) the mean of the data contributing to the head
    !        g) IF the correlogram - variance of tail values
    !        h) IF the correlogram - variance of head values

    !-----------------------------------------------------------------------


    real*8, intent(in), dimension(ndir*(nlag+2)*nvarg)  :: np,  gam, hm, tm, hv, tv
    real*8, intent(in)  :: xsiz, ysiz, zsiz
    integer, intent(in) :: nvarg, ndir, nlag 
    integer, intent(in), dimension(ndir) :: ixd, iyd, izd
    real*8,    intent(out), dimension(nvarg, ndir, nlag+2) :: pdis,pgam, phm,ptm,phv,ptv
    integer, intent(out), dimension(nvarg, ndir, nlag+2) :: pnump
    

    integer :: iv, id, il, i
    real*8 :: dis

    ! Loop over all the variograms that have been computed:

    do iv=1,nvarg
       
        ! Loop over all the directions (note the direction in the title):
    
        do id=1,ndir
        
            ! Compute the unit lag distance along the directional vector:
        
            dis = sqrt( max(((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2 + &
            (izd(id)*zsiz)**2),0.0) )
        
            ! Write out all the lags:
        
            do il=1,nlag
                i = (id-1)*nvarg*nlag+(iv-1)*nlag+il
                pdis(iv,id,il) = dble(il)*dis
                pgam(iv,id,il)=gam(i)
                phm(iv,id,il)=hm(i)
                ptm(iv,id,il)=tm(i)
                phv(iv,id,il)=hv(i)
                ptv(iv,id,il)=tv(i)
                pnump(iv,id,il)=int(np(i))

            end do
        end do
    end do

    return

 end subroutine writeout_gam


!***********************************************************************
! 
! 
!     Subroutines for interpolation GSLIB kt3d
! 
! 
!***********************************************************************

subroutine ktsol(n  ,ns  ,nv,a,b,x,ktilt,maxeq)
    !-----------------------------------------------------------------------

    ! Solution of a system of linear equations by gaussian elimination with
    ! partial pivoting.  Several right hand side matrices and several
    ! variables are allowed.


    !         NOTE: All input matrices must be in double precision


    ! INPUT/OUTPUT VARIABLES:

    !   n                Number of equations
    !   ns               Number of right hand side matrices
    !   nv               Number of variables.
    !   a(n*n*nv)        left hand side matrices versus columnwise.
    !   b(n*ns*nv)       input right hand side matrices.
    !   x(n*ns*nv)       solution matrices.
    !   ktilt            indicator of singularity
    !                      =  0  everything is ok.
    !                      = -1 n.le.1
    !                      =  k  a null pivot appeared at the kth iteration.
    !   tol              used in test for null pivot. depends on machine
    !                      precision and can also be set for the tolerance
    !                      of an ill-defined kriging system.


!-----------------------------------------------------------------------
    
    
    implicit real*8 (a-h,o-z)
    real*8 :: x(maxeq),a(maxeq*maxeq),b(maxeq)

    ! Make sure there are equations to solve:

    if(n <= 1) then
        ktilt = -1
        return
    endif

    ! Initialization:

    tol   = 0.1e-10
    ktilt = 0
    ntn   = n*n
    nm1   = n-1

    ! Triangulation is done variable by variable:

    do iv=1,nv
    
        ! Indices of location in vectors a and b:
    
        nva = ntn*(iv-1)
        nvb = n*ns*(iv-1)
    
        ! Gaussian elimination with partial pivoting:
    
        do k=1,nm1
            kp1 = k+1
        
            ! Indice of the diagonal element in the kth row:
        
            kdiag = nva+(k-1)*n+k
        
            ! Find the pivot - interchange diagonal element/pivot:
        
            npiv = kdiag
            ipiv = k
            i1   = kdiag
            do i=kp1,n
                i1 = i1+1
                if(abs(a(i1)) > abs(a(npiv))) then
                    npiv = i1
                    ipiv = i
                endif
            end do
            t        = a(npiv)
            a(npiv)  = a(kdiag)
            a(kdiag) = t
        
            ! Test for singularity:
        
            if(abs(a(kdiag)) < tol) then
                ktilt=k
                return
            endif
        
            ! Compute multipliers:
        
            i1 = kdiag
            do i=kp1,n
                i1    = i1+1
                a(i1) = -a(i1)/a(kdiag)
            end do
        
            ! Interchange and eliminate column per column:
        
            j1 = kdiag
            j2 = npiv
            do j=kp1,n
                j1    = j1+n
                j2    = j2+n
                t     = a(j2)
                a(j2) = a(j1)
                a(j1) = t
                i1    = j1
                i2    = kdiag
                do i=kp1,n
                    i1    = i1+1
                    i2    = i2+1
                    a(i1) = a(i1)+a(i2)*a(j1)
                end do
            end do
        
            ! Interchange and modify the ns right hand matrices:
        
            i1 = nvb+ipiv
            i2 = nvb+k
            do i=1,ns
                t     = b(i1)
                b(i1) = b(i2)
                b(i2) = t
                j1    = i2
                j2    = kdiag
                do j=kp1,n
                    j1    = j1+1
                    j2    = j2+1
                    b(j1) = b(j1)+b(i2)*a(j2)
                end do
                i1 = i1+n
                i2 = i2+n
            end do
        end do
    
        ! Test for singularity for the last pivot:
    
        kdiag = ntn*iv
        if(abs(a(kdiag)) < tol) then
            ktilt = n
            return
        endif
    end do

    ! End of triangulation. Now, solve back variable per variable:

    do iv=1,nv
    
    ! Indices of location in vectors a and b:
    
        nva  = ntn*iv
        nvb1 = n*ns*(iv-1)+1
        nvb2 = n*ns*iv
    
        ! Back substitution with the ns right hand matrices:
    
        do il=1,ns
            do k=1,nm1
                nmk = n-k
            
                ! Indice of the diagonal element of the (n-k+1)th row and of
                ! the (n-k+1)th element of the left hand side.
            
                kdiag = nva-(n+1)*(k-1)
                kb    = nvb2-(il-1)*n-k+1
                b(kb) = b(kb)/a(kdiag)
                t     = -b(kb)
                i1    = kb
                i2    = kdiag
                do i=1,nmk
                    i1    = i1-1
                    i2    = i2-1
                    b(i1) = b(i1)+a(i2)*t
                end do
            end do
            kdiag = kdiag-n-1
            kb    = kb-1
            b(kb) = b(kb)/a(kdiag)
        end do
    
        ! End of back substitution:
    
    end do

    ! Restitution of the solution:

    itot = n*ns*nv
    do i=1,itot
        x(i) = b(i)
    end do

    ! Finished:

    return
end subroutine ktsol


subroutine kt3d( &   
                 na,xa,ya,za,vra, vea,                     &            ! input parameters (Data within search ellipse)
                 ndb, xdb, ydb, zdb, extest, cbb,          &            ! input grid block and block covariance (use discretization points for block estimate or center point for point estimate)
                 radius,                                   &            ! input search parameters (this is a dummy parameter for rescaling, TODO: find the way to remove it)
                 nst,c0,it,cc,aa,aa1,aa2,ang1,ang2,ang3,   &            ! input variogram
                 ktype, skmean,  UNEST,                    &            ! input kriging parameters 
                 idrift,                                   &            ! input drift parameters                         
                 kneq,                                     &            ! number of kriging equations (can be calculated with kt3d_getmatrix_size )
                 est, estv, estt, estvt,                   &            ! output: estimate, estimation variance on a single block (including estimate of the trend)
                 w, wt,  error,                            &            ! output: weight (and trend wight) and error (if error error > 0)
                 kmatrix, kvector, ksolution)                           ! output: system of kriging equations, only returned if kneq> 0

    
    !-----------------------------------------------------------------------
    !                Krige a 3-D Grid of Rectangular Blocks
    !                **************************************
    !
    ! This subroutine estimates point or block values of one variable by
    ! simple, ordinary, or kriging with a trend model.  It is also possible
    ! to estimate the trend directly.

    ! Original:  A.G. Journel and C. Lemmer                             1981
    ! Revisions: A.G. Journel and C. Kostov                             1984
    !
    !
    ! 2015 changes
    ! this is only to calculate in a single block
    ! the block is defined by discretization points (this functions may estimate polygos... )
    ! all data input is used
    ! This is a function to be used from python... 
    ! TODO: add comments 
    !
    !
    ! PARAMETERS:
    ! Input: 
    !  *Data points for estimation
    !     na                    - integer:  number of rows in the data
    !     xa(na),ya(na),za(na)  - [double]: coordinates 
    !     vra(nd),vea(nd)       - [double]: variable and external drift
    !  *Grid parameters (no rotation allowed)
    !     ndb                   - integer: number of discretization points
    !     xdb(ndb)              - [double]: x coordinates of discretization points
    !     ydb(ndb)              - [double]: y coordinates of discretization points
    !     zdb(ndb)              - [double]: z coordinates of discretization points 
    !     extest                - double:  external drift on block/point
    !     cbb                   - double:  block covariance
    !  *Search parameters                      
    !     radius                - double: ellipse radius (this is dummy, TODO: remove it)
    !  *Variogram 
    !     nst(1)                   - [integer]: number of structures 
    !     c0(1)                    - [double]: nugget effect    
    !     it(nst)                  - [integer]: structure type 
    !                                  1. spherical (range a)
    !                                  2. exponential (p'range 3a)
    !                                  3. gaussian (p'range a*sqrt(3))
    !                                  4. power (0<a<2), if linear model, a=1,c=slope.
    !                                  5. hole effect
    !     cc(nst)                  - [double]: structure variance
    !     aa, aa1, aa2(nst)        - [double]: structure ranges
    !                                 ** aa is used to calculate anisotroy
    !     ang1,ang2,ang3(nst)      - [double]: structure rotation angles
    !                                 ** variable angles per structure allowed    
    !  *Kriging parameters
    !     ktype                    - integer: 0=SK,1=OK,2=non-st SK,3=exdrift
    !     UNEST                    - double: non estimated values (ex. numpy.nan)
    !     skmean                   - double: constant used for simple kriging mean
    !                                      *warning this is an inout parameter*
    !     kneq                     - integer: number of kriging equations 
    !                                         if 0 no kriging equations reported
    !                                         if equal to the actual number of k equations k equations reported
    !                                         if > 0 and wrong number report error
    !                                         Note: use kt3d_getmatrix_size to calculate kneq
    !  *Drift
    !     idrift(9)                 - [integer]: if true will use or ignore
    !                                           the following drift terms: 
    !                                           x,y,z,xx,yy,zz,xy,xz,zy 
    ! Output: 
    !  *Estimate
    !     est, estv,               - double: estimate and kriging variance
    !     estt, estvt              - double: drift estimate and kriging variance (drift)
    !  *weight
    !    w(na), wt(na)             - [double] : estimate wight and trend weight
    !  *Kriging equations
    !    kmatrix(kneq,kneq)       - [double] : kriging matriz 
    !    kvector(kneq)            - [double] : kriging vector
    !    ksolution(kneq)          - [double] : kriging solution vector
    !-----------------------------------------------------------------------
    

	implicit none

	! input variables
	! target
	integer, intent(in) :: ndb  ! total number of discretization points 
	real*8,  intent(in), dimension (ndb) :: xdb, ydb, zdb
	real*8,  intent(in) :: extest   ! this is the external drift at target location
	real*8,  intent(in) :: radius   ! this is for internal rescal factor, TODO: find way to remove it
	real*8,  intent(in) :: UNEST

	! kriging parameters
	integer, intent(in) :: ktype
	real*8, intent(in) :: cbb
	real*8, intent(inout) :: skmean

	! drift 
	integer, intent(in), dimension (9) :: idrift

	! variogram
	integer, intent(in) :: nst
	real*8,  intent(in), dimension (1) :: c0
	real*8,  intent(in), dimension (nst) :: cc, aa, aa1, aa2, ang1, ang2, ang3
	integer, intent(in), dimension (nst) :: it 
	
	! data
	integer, intent(in) :: na
	real*8,  intent(in), dimension (na) :: xa, ya, za, vra, vea

	! ksystem of equations
	integer, intent(in) :: kneq


	! output variables 
	real*8,  intent(out):: est, estv ! estimate and estimate variance
	real*8,  intent(out):: estt, estvt ! estimate and estimate variance (trend)
	integer, intent(out) :: error ! 1=> allocation error
	real*8,  intent(out), dimension (na) :: w, wt  ! kriging weight (variable and trend) asigned to each variable
	! ksystem of equations
	real*8,  intent(out), dimension (kneq,kneq) :: kmatrix
	real*8,  intent(out), dimension (1,kneq) :: kvector, ksolution


	! internal variables
	real*8 :: PMX, covmax, EPSLON, radsqd , resc, xk, &
              vk, xkmae, xkmse, cb, cmax, cb1, cov, dx, dy, dz, &
			  unbias, resce
	real*8, dimension (nst,3,3) :: rotmat
	real*8, dimension (nst) :: anis1,anis2
	real*8, dimension (9) :: bv
    integer, dimension (9) :: idrif
	integer :: is, ie, i, j, k , maxrot, mdt, nk, neq, im, test, &
			   kadim, ksdim, nrhs, nv, ising

	! internal and dinamic
	real*8, allocatable,  dimension (:) :: a, at, att, r, rr, rt, rrt, s, st
	! TODO: find how to provide a,r,rt,s,st as outpot




	error = 0
	PMX    = 999.0
	EPSLON = 0.00000001


	! ********** this is the actual function ***************

	!calculate anisotropy factor 
    do i=1,nst
        anis1(i) = aa1 (i) / max(aa(i),EPSLON)
        anis2(i) = aa2 (i) / max(aa(i),EPSLON)
    end do
    
    ! put drift in a new variable to avoid inout 
    do i=1, 9
        idrif(i) = idrift(i)
    end do

	! Set up the rotation/anisotropy matrices that are needed for the
	! variogram.  Also compute the maximum covariance for
	! the rescaling factor:
    
	maxrot = nst
    covmax = c0(1)	
    do is=1,nst
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,MAXROT,rotmat)
        if(it(is) == 4) then
            covmax = covmax + PMX
        else
            covmax = covmax + cc(is)
        endif
    end do


	! Finish computing the rescaling factor and stop if unacceptable:
	radsqd = radius * radius
    if(radsqd < 1.0) then
        resc = 2.0 * radius / max(covmax,0.0001)
    else
        resc =(4.0 * radsqd)/ max(covmax,0.0001)
    endif
    if(resc <= 0.0) then
        write(*,*) 'ERROR KT3D: The rescaling value is wrong ',resc
        write(*,*) '            Maximum covariance: ',covmax
        write(*,*) '            search radius:      ',radius
        error = 700   ! rescale factor error
        return
    endif
    resc = 1.0 / resc


	! Compute the number of drift terms, if an external drift is being
	! considered then it is one more drift term, if SK is being considered
	! then we will set all the drift terms off and mdt to 0):

    mdt = 1
    do i=1,9
        if(ktype == 0 .OR. ktype == 2) idrif(i) = 0
        if(idrif(i) < 0 .OR. idrif(i) > 1) then
            write(*,*) 'ERROR KT3D: invalid drift term',idrif(i)
            write(*,*) 'Using idrif(i) = 0 on drift i= ', i
            idrif(i) = 0
        endif
        mdt = mdt + idrif(i)
    end do
    if(ktype == 3) mdt = mdt + 1
    if(ktype == 0) mdt = 0
    if(ktype == 2) mdt = 0
	
	! print *, 'mdt : ', mdt

	! Set up the discretization points per block.  Figure out how many
	! are needed, the spacing, and fill the xdb,ydb, and zdb arrays with
	! the offsets relative to the block center (this only gets done once):

	! In all cases the offsets are relative to the lower left corner.
	! This is done for rescaling the drift terms in the kriging matrix.

	! in this version we input arbitrary xdb,ydb, zdb values... 


	! Initialize accumulators:

    nk    = 0
    xk    = 0.0
    vk    = 0.0
    xkmae = 0.0
    xkmse = 0.0

	!
	! Calculate point Covariance. The block covariance is calculated externally
	!
	call cova3(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst, &
		       c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
	!
	! Set the ``unbias'' variable so that the matrix solution is more stable
	!
    
	unbias = cov


	! Mean values of the drift functions:

    do i=1,9
        bv(i) = 0.0
    end do
    do i=1,ndb
        bv(1) = bv(1) + xdb(i)
        bv(2) = bv(2) + ydb(i)
        bv(3) = bv(3) + zdb(i)
        bv(4) = bv(4) + xdb(i)*xdb(i)
        bv(5) = bv(5) + ydb(i)*ydb(i)
        bv(6) = bv(6) + zdb(i)*zdb(i)
        bv(7) = bv(7) + xdb(i)*ydb(i)
        bv(8) = bv(8) + xdb(i)*zdb(i)
        bv(9) = bv(9) + ydb(i)*zdb(i)
    end do
    do i=1,9
        bv(i) = (bv(i) / real(ndb)) * resc
    end do

	! Test if there are enough samples to estimate all drift terms:

    if(na >= 1 .AND. na <= mdt) then
        est  = UNEST
        estv = UNEST
		error = 100        ! no enough samples error
        return 
    end if

	! There are enough samples - proceed with estimation.

    if(na <= 1) then
    
    	! Handle the situation of only one sample:
    
        call cova3(xa(1),ya(1),za(1),xa(1),ya(1),za(1),1,nst, &
        c0,it,cc,aa,1,maxrot,rotmat,cmax,cb1)
    
    	! Establish Right Hand Side Covariance:
    
        if(ndb <= 1) then
            call cova3(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),1, &
            nst,c0,it,cc,aa,1,maxrot,rotmat,cmax,cb)
        else
            cb  = 0.0
            do i=1,ndb
                call cova3(xa(1),ya(1),za(1),xdb(i),ydb(i), &
                zdb(i),1,nst,c0,it,cc,aa,1, &
                MAXROT,rotmat,cmax,cov)
                cb = cb + cov
                dx = xa(1) - xdb(i)
                dy = ya(1) - ydb(i)
                dz = za(1) - zdb(i)
                if((dx*dx+dy*dy+dz*dz) < EPSLON) cb=cb-c0(1)
            end do
            cb = cb / real(ndb)
        end if
        est  = vra(1)
        estv = real(cbb) - 2.0*cb + cb1
        nk   = nk + 1
        xk   = xk + vra(1)
        vk   = vk + vra(1)*vra(1)
		error = 900000               ! warning, estimate with one sample
        return 
    end if

	! Go ahead and set up the OK portion of the kriging matrix:

    neq = mdt+na

	! Initialize the main kriging matrix:

    allocate( a(neq*neq), att(neq*neq), at(neq*neq), r(neq), rr(neq), rt(neq), rrt(neq), s(neq), st(neq),  stat = test)
	if(test.ne.0)then
    	error = 1   ! allocation error
		return        
    end if


    do i=1,neq*neq
        a(i) = 0.0
    end do
	

    do i=1,neq
        r(i) = 0.0
		rr(i) = 0.0
		rt(i) = 0.0
		rrt(i) = 0.0
		s(i) = 0.0
		st(i) = 0.0
    end do

	! Fill in the kriging matrix:

    do i=1,na
        do j=i,na
            call cova3(xa(i),ya(i),za(i),xa(j),ya(j),za(j),1,nst, &
            c0,it,cc,aa,1,maxrot,rotmat,cmax,cov)
            a(neq*(i-1)+j) = dble(cov)
            a(neq*(j-1)+i) = dble(cov)
        end do
    end do

	! Fill in the OK unbiasedness portion of the matrix (if not doing SK):

    if(neq > na) then
        do i=1,na
            a(neq*(i-1)+na+1) = dble(unbias)
            a(neq*na+i)       = dble(unbias)
        end do
    endif


	! Set up the right hand side:

    do i=1,na
        if(ndb <= 1) then  ! point kriging
			! print *, 'doing point kriging'
            call cova3(xa(i),ya(i),za(i),xdb(1),ydb(1),zdb(1),1, &
            nst,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
        else
			! print *, 'doing block kriging'
            cb  = 0.0
            do j=1,ndb
                call cova3(xa(i),ya(i),za(i),xdb(j),ydb(j), &
                zdb(j),1,nst,c0,it,cc,aa,1, &
                MAXROT,rotmat,cmax,cov)
                cb = cb + cov
                dx = xa(i) - xdb(j)
                dy = ya(i) - ydb(j)
                dz = za(i) - zdb(j)
                if((dx*dx+dy*dy+dz*dz) < EPSLON) cb=cb-c0(1)
            end do
            cb = cb / real(ndb)
        end if
        r(i) = dble(cb)
    end do
    if(neq > na) r(na+1) = dble(unbias)

	! Add the additional unbiasedness constraints:

	im = na + 1

	! First drift term (linear in "x"):

    if(idrif(1) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(xa(k)*resc)
            a(neq*(k-1)+im) = dble(xa(k)*resc)
        end do
        r(im) = dble(bv(1))
    endif

	! Second drift term (linear in "y"):

    if(idrif(2) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(ya(k)*resc)
            a(neq*(k-1)+im) = dble(ya(k)*resc)
        end do
        r(im) = dble(bv(2))
    endif

	! Third drift term (linear in "z"):

    if(idrif(3) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(za(k)*resc)
            a(neq*(k-1)+im) = dble(za(k)*resc)
        end do
        r(im) = dble(bv(3))
    endif

	! Fourth drift term (quadratic in "x"):

    if(idrif(4) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(xa(k)*xa(k)*resc)
            a(neq*(k-1)+im) = dble(xa(k)*xa(k)*resc)
        end do
        r(im) = dble(bv(4))
    endif

	! Fifth drift term (quadratic in "y"):

    if(idrif(5) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(ya(k)*ya(k)*resc)
            a(neq*(k-1)+im) = dble(ya(k)*ya(k)*resc)
        end do
        r(im) = dble(bv(5))
    endif

	! Sixth drift term (quadratic in "z"):

    if(idrif(6) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(za(k)*za(k)*resc)
            a(neq*(k-1)+im) = dble(za(k)*za(k)*resc)
        end do
        r(im) = dble(bv(6))
    endif

	! Seventh drift term (quadratic in "xy"):

    if(idrif(7) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(xa(k)*ya(k)*resc)
            a(neq*(k-1)+im) = dble(xa(k)*ya(k)*resc)
        end do
        r(im) = dble(bv(7))
    endif

	! Eighth drift term (quadratic in "xz"):

    if(idrif(8) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(xa(k)*za(k)*resc)
            a(neq*(k-1)+im) = dble(xa(k)*za(k)*resc)
        end do
        r(im) = dble(bv(8))
    endif

	! Ninth drift term (quadratic in "yz"):

    if(idrif(9) == 1) then
        im=im+1
        do k=1,na
            a(neq*(im-1)+k) = dble(ya(k)*za(k)*resc)
            a(neq*(k-1)+im) = dble(ya(k)*za(k)*resc)
        end do
        r(im) = dble(bv(9))
    endif

	! External drift term (specified by external variable):

    if(ktype == 3) then
        im=im+1
		resce  = covmax / max(extest,0.0001)
        do k=1,na
            a(neq*(im-1)+k) = dble(vea(k)*resce)
            a(neq*(k-1)+im) = dble(vea(k)*resce)
        end do
        r(im) = dble(extest*resce)
		! print *, 'r(im)', r(im), im
		! print *, 'covmax, extest, resce ', covmax , extest, resce
    endif


 	! Copy the right hand side to compute the kriging variance later:
	! this is because ksolve may change r values... 
    do k=1,neq
        rr(k) = r(k)
		rt(k) = r(k)
		rrt(k)= r(k) 
    end do
	! doing the same with a
    do k=1,neq*neq
        at(k) = a(k)
		att(k) = a(k) 
    end do
    kadim = neq * neq
    ksdim = neq
    nrhs  = 1
    nv    = 1

   
	! To estimate the trend we reset all the right hand side terms=0.0:
    do i=1,na
          rt(i)  = 0.0
          rrt(i) = 0.0
    end do


	! Solve the kriging system for data estimate
    call ktsol(neq,nrhs,nv,a,r,s,ising,neq)

	! Solve the kriging system for trend estimate
	call ktsol(neq,nrhs,nv,at,rt,st,ising,neq)


	! Compute the solution:
    if(ising /= 0) then
        est  = UNEST
        estv = UNEST
        estt  = UNEST
        estvt = UNEST
        error = 20    ! singular matrix
        deallocate( a, r, rr, rt, rrt, s, st,  stat = test)
        return
    else
        est  = 0.0
        estv = real(cbb)
        estt  = 0.0
        estvt = real(cbb)
        if(ktype == 2) skmean = extest
        do j=1,neq
            estv = estv - real(s(j))*rr(j)
			estvt = estvt - real(st(j))*rrt(j)
            if(j <= na) then
                if(ktype == 0 .OR. ktype == 2) then
                    est = est + real(s(j))*(vra(j)-skmean)
					estt = estt + real(st(j))*(vra(j)-skmean)
                else
                    est = est + real(s(j))*vra(j)
					estt = estt + real(st(j))*vra(j)
                endif
                w(j) = s(j)
                wt(j) = st(j)
            endif
        end do
        if(ktype == 0 .OR. ktype == 2) then
			est = est + skmean
			estt = estt + skmean			
		end if
 
    end if


	!
	! Write out the kriging equations
	!

	! The matrix has the right size?
	if (kneq>0 .and. kneq/=neq) then
		error = 10  ! wrong size for the kriging matrix, use kt3d_getmatrix_size to calculate right size
		deallocate( a, r, rr, rt, rrt, s, st,  stat = test)
		return
	end if	

	! then we populate the external matrix/vectors with values
    if(kneq>0) then
		is = 1 - neq

		do i=1, neq
			is = 1 + (i-1)*neq
			ie = is + neq - 1
			kvector(1,i)   = rr(i)
			ksolution(1,i) =  s(i)
			k=0
		    do j=is,ie
				k=k+1
		        kmatrix (i,k) = att(j) 
		        kmatrix (k,i) = att(j)
		    end do
			!kmatrix (neq,neq) = a(neq*neq)
		end do
	endif
           
	!dealocate arrays
    deallocate( a, r, rr, rt, rrt, s, st,  stat = test)
	if(test.ne.0)then
    	error = 2    ! deallocation error
		return           
    end if

	return

end subroutine kt3d



subroutine kt3d_getmatrix_size ( &
			  ktype, idrift , na, & 
			  mdt, kneq, error)
    !-----------------------------------------------------------------------
    !                Gets the size of the kriging equations from parameters
    !                **************************************
    !
    !
    ! PARAMETERS:
    ! Input: 
    !  *Data points for estimation
    !     na                    - integer:  number of rows in the data
    !  *Kriging parameters
    !     ktype                    - integer: 0=SK,1=OK,2=non-st SK,3=exdrift

    !  *Drift
    !     idrift(9)                 - [logical]: if true will use or ignore
    !                                           the following drift terms: 
    !                                           x,y,z,xx,yy,zz,xy,xz,zy 
    ! Output: 
    !     mdt                      - integer: size of the unbias terms
    !     kneq                     - integer: number of kriging equations (na + mdt)
    !-----------------------------------------------------------------------


	implicit none

	! input variables
	! target

	! kriging parameters
	integer, intent(in) :: ktype

	! drift 
	integer, intent(in), dimension (9) :: idrift
	
	! data
	integer, intent(in) :: na

	! output variables 
	integer, intent(out) :: error ! 1=> allocation error
	integer, intent(out) ::	mdt, kneq

	! internal variables
	integer :: i
    integer,  dimension (9) :: idrif

	error = 0

    ! put drift in a new variable to avoid inout 
    do i=1, 9
        idrif(i) = idrift(i)
    end do


	! Compute the number of drift terms, if an external drift is being
	! considered then it is one more drift term, if SK is being considered
	! then we will set all the drift terms off and mdt to 0):

    mdt = 1
    do i=1,9
        if(ktype == 0 .OR. ktype == 2) idrif(i) = 0
        if(idrif(i) < 0 .OR. idrif(i) > 1) then
            write(*,*) 'ERROR KT3D: invalid drift term',idrif(i)
            write(*,*) 'Using idrif(i) = 0 on drift i= ', i
            error = 1
            idrif(i) = 0
        endif
        mdt = mdt + idrif(i)
    end do
    if(ktype == 3) mdt = mdt + 1
    if(ktype == 0) mdt = 0
    if(ktype == 2) mdt = 0
	

	! Test if there are enough samples to estimate all drift terms:

    if(na >= 1 .AND. na <= mdt) then
		error = 100        ! no enough samples error
        return 
    end if

	! Go ahead and set up the OK portion of the kriging matrix:

    kneq = mdt+na

	return

end subroutine kt3d_getmatrix_size 


!***********************************************************************
! 
! 
!     Subroutines for data statistics and transformation
! 
! 
!***********************************************************************

!---------------------------------------------------------------------------
!     Subroutine declus (this is the declus program for declustering)
!---------------------------------------------------------------------------
subroutine declus( &
				 x,y,z,vr, nd, anisy,anisz, minmax, ncell, cmin, cmax, noff, MAXCEL,  & ! input
				 wtopt, vrop,wtmin,wtmax, error, &         ! out
				 xinc, yinc, zinc, rxcs, rycs, rzcs,rvrcr)
	! this is modified from opengeostat 2015, MIT lisence. .
	! this is developing test... 

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

	!         DECLUS: a three dimensional cell declustering program
	!         *****************************************************

	! See paper in Computers and Geosciences Vol 15 No 3 (1989) pp 325-332

	! INPUT Parameters:

	!   x,y,z,vr        columns for X, Y, Z, and variable
	!   not using this !         tmin,tmax       trimming limits
	!   anisy,anisz     Y and Z cell anisotropy (Ysize=size*Yanis)
	!   minmax          0=look for minimum declustered mean (1=max)
	!   ncell,cmin,cmax number of cell sizes, min size, max size
	!   noff            number of origin offsets


	! OUTPUT 
	!   wtopt                    declustering wight 
	!   rxcs, rycs, rzcs,rvrcr   cell size and declustered mean


	! PROGRAM NOTES:
	!   3. This program requires knowledge of whether the samples are
	!      clustered in high or low values or, alternately, knowledge
	!      of a ``natural'' cell size, e.g., an underlying regular data
	!      spacing.



	! The following Parameters control static dimensioning:
	!   MAXCEL    maximum number of cells.  The number of cells is a
	!             function of the cell size and the size of the area of
	!             interest.  In many cases a larger minimum cell size will
	!             remove the need to increase MAXCEL.




	!-----------------------------------------------------------------------
    

	implicit none

	! Input
	real*8, intent(in) :: anisy,anisz, cmin
	integer, intent(in) :: nd, minmax, ncell, noff
    real*8, intent(in), dimension(nd)  ::  x(nd),y(nd),z(nd),vr(nd)
	
	! Inout
	real*8, intent(in) :: cmax 
	integer, intent(in) :: MAXCEL

	! Output
	real*8, intent(out), dimension(nd)  :: wtopt
	real*8, intent(out), dimension(ncell+1)  :: rxcs, rycs, rzcs,rvrcr 
    

	real*8, intent(out) ::vrop,wtmin,wtmax , &
						  xinc, yinc, zinc
	integer, intent(out) :: error

	! Internal 
	real*8 :: vrmin, vrmax, roff, vrav, best, facto, &
			  xo1, yo1, zo1, xo, yo, zo, sumw, sumwg, xfac, yfac, zfac
    logical :: min
	real*8 :: wt(nd)
    integer :: index(nd), lp, i, kp,  & 
			   icell, icellx, icelly, icellz, &
			   ipoint, ncellt, ncellx, ncelly, ncellz, &
			   test
    real*8, allocatable :: cellwt(:)
    real*8 :: xmin,ymin,zmin,xmax,ymax,zmax, xcs, ycs, zcs,vrcr
	
	real*8 :: tcmax 
	integer:: tMAXCEL


    xmin = 1.0e21
	ymin = 1.0e21
	zmin = 1.0e21
    xmax =-1.0e21
	ymax =-1.0e21
	zmax =-1.0e21
	
	error = 0

    tcmax = cmax
    tMAXCEL = MAXCEL


	! Doing only one/final declustering calculation? 
    if(ncell == 1) tcmax = cmin

	! Some Initialization:

    min  = .TRUE. 
    if(minmax == 1) min = .FALSE. 
    roff = real(noff)

	! compute min, max, and average:

    vrav = 0.0
    do i=1,nd
        wtopt(i) = 1.0
        vrav  = vrav + vr(i)
        if(x(i) < xmin) xmin=x(i)
        if(x(i) > xmax) xmax=x(i)
        if(y(i) < ymin) ymin=y(i)
        if(y(i) > ymax) ymax=y(i)
        if(z(i) < zmin) zmin=z(i)
        if(z(i) > zmax) zmax=z(i)
    end do
    vrav = vrav / real(nd)

	! initialize the "best" weight values:

    vrop = vrav
    best = 0.0

	! define a "lower" origin to use for the cell sizes:

    xo1 = xmin - 0.01
    yo1 = ymin - 0.01
    zo1 = zmin - 0.01

	! define the increment for the cell size:

    xinc = (tcmax-cmin) / real(ncell)
    yinc = anisy * xinc
    zinc = anisz * xinc

	! loop over "ncell+1" cell sizes in the grid network:

    xcs =  cmin        - xinc
    ycs = (cmin*anisy) - yinc
    zcs = (cmin*anisz) - zinc

	! MAIN LOOP over cell sizes:

    do lp=1,ncell+1
        xcs = xcs + xinc
        ycs = ycs + yinc
        zcs = zcs + zinc
    
    	! initialize the weights to zero:
    
        do i=1,nd
            wt(i) = 0.0
        end do
    
    	! determine the maximum number of grid cells in the network:
    
        ncellx = int((xmax-(xo1-xcs))/xcs)+1
        ncelly = int((ymax-(yo1-ycs))/ycs)+1
        ncellz = int((zmax-(zo1-zcs))/zcs)+1
        ncellt = real(ncellx*ncelly*ncellz)
    
    	! check the array MAXCEL dimensions:  (warning: if MAXCEL<1 we don't check this)
    
        if(ncellt > real(tMAXCEL) .and. tMAXCEL>1 ) then
            error = 10   ! ncellt > MAXCEL ' check for outliers - increase cmin and/or MAXCEL'
            return 
        end if
    
    	allocate( cellwt(ncellt),  stat = test)
		if(test.ne.0)then
    		error = 1    ! allocation error
			return        
		end if

    	! loop over all the origin offsets selected:
    
        xfac = amin1((xcs/roff),(0.5*(xmax-xmin)))
        yfac = amin1((ycs/roff),(0.5*(ymax-ymin)))
        zfac = amin1((zcs/roff),(0.5*(zmax-zmin)))
        do kp=1,noff
            xo = xo1 - (real(kp)-1.0)*xfac
            yo = yo1 - (real(kp)-1.0)*yfac
            zo = zo1 - (real(kp)-1.0)*zfac
        
        	! initialize the cumulative weight indicators:
        
            do i=1,ncellt
                cellwt(i) = 0.0
            end do
        
        	! determine which cell each datum is in:
        
            do i=1,nd
                icellx = int((x(i) - xo)/xcs) + 1
                icelly = int((y(i) - yo)/ycs) + 1
                icellz = int((z(i) - zo)/zcs) + 1
                icell  = icellx + (icelly-1)*ncellx &
                + (icellz-1)*ncelly*ncellx
                index(i)      = icell
                cellwt(icell) = cellwt(icell) + 1.0
            end do
        
        	! The weight assigned to each datum is inversely proportional to the
        	! number of data in the cell.  We first need to get the sum of weights
        	! so that we can normalize the weights to sum to one:
        
            sumw = 0.0
            do i=1,nd
                ipoint = index(i)
                sumw   = sumw + (1.0 / cellwt(ipoint))
            end do
            sumw = 1.0 / sumw
        
        	! Accumulate the array of weights (that now sum to one):
        
            do i=1,nd
                ipoint = index(i)
                wt(i) = wt(i) + (1.0/cellwt(ipoint))*sumw
            end do
        
        	! End loop over all offsets:
        
        end do

		deallocate( cellwt,  stat = test)
		if(test.ne.0)then
			error = 2    ! deallocation error
			return 
		end if  
    
    	! compute the weighted average for this cell size:
    
        sumw  = 0.0
        sumwg = 0.0
        do i=1,nd
            sumw  = sumw  + wt(i)
            sumwg = sumwg + wt(i)*vr(i)
        end do
        vrcr  = sumwg / sumw
		
		! this is the report for each cell size
        rxcs(lp) = xcs
		rycs(lp) = ycs
		rzcs(lp) = xcs
		rvrcr(lp) = vrcr                                                    
    
    	! see if this weighting is optimal:
    
        if((min .AND. vrcr < vrop) .OR. ( .NOT. min .AND. vrcr > vrop) .OR. &
        (ncell == 1)) then
            best = xcs
            vrop = vrcr
            do i=1,nd
                wtopt(i) = wt(i)
            end do
        end if
    
    	! END MAIN LOOP over all cell sizes:
    
    end do
 

	! Get the optimal weights:

    sumw = 0.0
    do i=1,nd
        sumw = sumw + wtopt(i)
    end do
    wtmin = 99999.
    wtmax =-99999.
    facto = real(nd) / sumw
    do i = 1,nd
        wtopt(i) = wtopt(i) * facto
        if(wtopt(i) < wtmin) wtmin = wtopt(i)
        if(wtopt(i) > wtmax) wtmax = wtopt(i)
    end do



    return

end subroutine declus
