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


!*********************************************************************************
!     Subroutines in GSLIB (auxiliary functions)
!*********************************************************************************
include 'setrot.f90'
include 'cova3.f90'
include 'sqdist.f90'
! include 'getindx.f90'
! include 'locate.f90'
! this includes sortem and dsortem
! include 'sortem.f90'
! include 'acorni.f90'
! include 'gauinv.f90'
! include 'powint.f90'
! include 'gcum.f90'


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
    integer, intent(out) ::    mdt, kneq

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
