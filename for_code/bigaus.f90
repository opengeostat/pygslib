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


real function f(x, zc)
    !-----------------------------------------------------------------------

    ! This function calculates the values of the function
    !       f( zc, xi)= exp ( -zc**2/ (1+sin(x))) for zc and x

    !-----------------------------------------------------------------------
    real ::       x, zc

    f = exp ( -zc**2 / (1+sin(x)))

    return

end function f


real function g(x, zc)
    real :: x, zc
    g= exp ( - x**2/2.0)
    return
end function g

real function simpson(alim,blim,tol,nitn,f,zc)
    !-----------------------------------------------------------------------

    !       simpson performs numerical integration over the interval
    !       a,b of the function f. the method employed is modified
    !       simpson's rule. the procedure was adapted from algorithm
    !       182 of "collected algorithms from acm" vol. 1.
    !       the algorithm has been tested and found particularly useful
    !       for integration of a strongly peaked function. the function
    !       used for testing was not able to be suitably evaluated using
    !       gauss quadrature or romberg integration.

    ! PARAMETERS:

    !       a    -    lower limit of integration      (real)
    !       b    -    upper limit of integration      (real)
    !       eps  -    required tolerance            (real)
    !       nitn -    maximum level of subdivision (integer)


    !-----------------------------------------------------------------------
    
    ! external, intent(hide) :: f
    
     
    double precision :: absarea
    parameter (num1 = 100)
    dimension dx(num1) , x2(num1) , x3(num1) , f2(num1) , &
    f3(num1) , epsp(num1) , f4(num1) , fmp(num1) , &
    fbp(num1) , est2(num1) , est3(num1) , &
    pval(num1,3) , irtrn(num1)

    !       initialise the return level array.

    do 30 i = 1,num1
        irtrn(i) = 0
    30 END DO
    !if(nitn > 50) then
    !    error = 1000 ! 'Sorry, I only do 50 iterations!'
    !end if    

    a = alim
    b = blim
    eps = tol
    lvl = 0
    est = 1.0
    absarea = 1.0
    da = b - a
    fm = 4.0 * f((a+b)/2.0, zc)
    fa = f(a,zc)
    fb = f(b,zc)

    10 continue
    lvl = lvl + 1
    dx(lvl) = da / 3.0
    sx = dx(lvl) / 6.0
    f1 = 4.0 * f(a+dx(lvl)/2.0,zc)
    x2(lvl) = a + dx(lvl)
    f2(lvl) = f(x2(lvl),zc)
    x3(lvl) = x2(lvl) + dx(lvl)
    f3(lvl) = f(x3(lvl),zc)
    epsp(lvl) = eps
    f4(lvl) = 4.0 * f(x3(lvl) + dx(lvl)/2.0,zc)
    fmp(lvl) = fm
    est1 = ( fa + f1 + f2(lvl)) * sx
    fbp(lvl) = fb
    est2(lvl) = (f2(lvl) + f3(lvl) + fm) * sx
    est3(lvl) = (f3(lvl) + f4(lvl) + fb) * sx
    sum = est1 + est2(lvl) + est3(lvl)
    absarea = absarea - abs(est) + abs(est1) + &
    abs(est2(lvl)) + abs(est3(lvl))

    if ((abs(est-sum) <= epsp(lvl) * absarea &
     .AND. est /= 1.0) &
     .OR. lvl >= nitn) then
        lvl = lvl - 1
        pval(lvl,irtrn(lvl)) = sum
        ipoint = irtrn(lvl) + 1
        goto (1,2,3,4) ipoint
    endif

    1 da = dx(lvl)
    fm = f1
    fb = f2(lvl)
    est = est1
    eps = epsp(lvl) / 1.7
    irtrn(lvl) = 1

    goto 10

    2 da = dx(lvl)
    fa = f2(lvl)
    fm = fmp(lvl)
    fb = f3(lvl)
    eps = epsp(lvl) / 1.7
    est = est2(lvl)
    a = x2(lvl)
    irtrn(lvl) = 2

    goto 10

    3 da = dx(lvl)
    fa = f3(lvl)
    fm = f4(lvl)
    fb = fbp(lvl)
    eps = epsp(lvl) / 1.7
    est = est3(lvl)
    a = x3(lvl)
    irtrn(lvl) = 3

    goto 10

    4 sum = pval(lvl,1) + pval(lvl,2) + pval(lvl,3)

    if ( lvl > 1 ) then
        lvl = lvl - 1
        pval(lvl,irtrn(lvl)) = sum
        ipoint = irtrn(lvl) + 1
        goto (1,2,3,4) ipoint
    endif

    simpson = sum
    return
 
end function simpson


subroutine bigaus(ncut, zcc, ndir,nlag, nst, azm,dip,xlag, &
                it, c0,cc,ang1,ang2,ang3, aa,aa1,aa2, &
                zci, lag, dir, hh, ri, ci, rop, error, rotmat)
    !          Indicator Variograms for BiGaussian Distribution
    !          ************************************************

    ! This program returns the values of the theoretical indicator
    ! semivariograms corresponding to an input normal scores semivariogram

    ! INPUT PARAMETERS:
    !   ncut                   number of cutoffs
    !   cut()                  the cutoffs (in cdf units)
    !   ndir,nlag              number of directions, number of lags
    !   xoff,yoff,zoff         for each direction: the specification
    !   nst,c0                 Normal scores variogram: nst, nugget
    !   it,aa,cc               type, a parameter, c parameter
    !   ang1,ang2,...,anis2    anisotropy definition


    ! EXTERNAL REFERENCES:

    !       f        to calculate exp(-zc**2/(1+x))/sqrt(1-x**2)
    !       g        to calculate exp(-x**2/2)
    !       simpson  to do the numerical calculation
    !       gauinv   to get standard normal deviate

    ! Original:  H. Xiao                                                1975
    !-----------------------------------------------------------------------
    
    implicit none
    
    external   f
    external   g

    ! input 
    integer, intent(in) :: ncut, ndir,nlag, nst
    real, intent(in) :: c0
    real, intent(in), dimension(ncut) :: zcc
    real, intent(in), dimension(ndir) :: azm,dip,xlag
    real, intent(in), dimension(nst) :: it,cc,ang1,ang2,ang3, aa,aa1,aa2

    ! output (similar to writeout)
    real, intent(out), dimension(ncut, ndir, nlag) ::  zci, lag, dir, hh, ri, ci, rop
    real*8, intent(out) ::     rotmat(nst,3,3)
    
    integer, intent(out) ::  error

    ! internal
    real, dimension(ncut) :: zc
    real, dimension(nlag*ndir) :: gam, ro, h
    integer :: i, ii, is, icut, il, l, j, id, ierr
    real, dimension(ndir) :: xoff, yoff, zoff
     real, dimension(nst) :: anis1,anis2
    real :: PI, DEG2RAD, EPSLON
    real ::       b, p,maxcov, ci0, cmax, cov, simpson, zz, xx, xp, yy
    


    PI=3.14159265
    DEG2RAD=PI/180.0
    EPSLON = 1.0e-20
    !    MLAG=50, MDIR=4, MAXNST=4, MAXROT=4


    if(ncut < 0) then
        error = 10    !ncut > 0?
        return
    endif

    do i=1,ncut
        zc(i)=zcc(i)
    end do

    do i=1,ndir
        xoff(i) = sin(DEG2RAD*azm(i))*xlag(i)*cos(DEG2RAD*dip(i))
        yoff(i) = cos(DEG2RAD*azm(i))*xlag(i)*cos(DEG2RAD*dip(i))
        zoff(i) = sin(DEG2RAD*dip(i))*xlag(i)
    end do

    if(nst <= 0) then
        error = 100 ! ' nst must be at least 1'
        return 
    endif


    do i=1,nst
        anis1(i) = aa1(i) / max(aa(i),EPSLON)
        anis2(i) = aa2(i) / max(aa(i),EPSLON)
        if(it(i) == 4) then
            if(aa(i) < 0.0) then
                error = 200 ! ' INVALID power variogram'
            end if 
            if(aa(i) > 2.0) then
                error=210 !' INVALID power variogram'
            end if
        end if
    end do


    ! Set up the rotation matrices:

    rotmat(:, :, :) =0

    do is=1,nst
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,nst,rotmat)
    end do

    ! Convert cutoffs cumulative probabilities to standard normal cutoffs:

    do icut=1,ncut
        call gauinv(dble(zc(icut)),xp,ierr)
        zc(icut) = xp
    end do

    ! Get ready to go:

    call cova3(0.0,0.0,0.0,0.0,0.0,0.0,1,nst,c0,it,cc,aa, &
                1,nst,rotmat,cmax,maxcov)
    ! Set the variogram data, direction by direction up to ndir directions:

    i = 0
    do id=1,ndir
        xx  = -xoff(id)
        yy  = -yoff(id)
        zz  = -zoff(id)
        do il=1,nlag
            xx  = xx + xoff(id)
            yy  = yy + yoff(id)
            zz  = zz + zoff(id)
            call cova3(0.0,0.0,0.0,xx,yy,zz,1,nst,c0,it, &
            cc,aa,1,nst,rotmat,cmax,cov)
            i      = i + 1
            gam(i) = maxcov - cov
            ro(i)  = cov/maxcov
            h(i)   = sqrt(xx*xx+yy*yy+zz*zz)
            hh (1,id,il) =  h(i)
        end do
    end do

    do i=2,ncut
        do l=1,ndir
            do j=1,nlag
                hh (i,l,j) = hh (1,l,j)
            end do
        end do
    end do

    ! Now, loop over all the cutoffs:

    do i=1,ncut
        ci0 = simpson(0.0, pi/2.0, 1.e-6, 40, f, zc(i))
        p   = simpson(0.0,zc(i),1.e-6,40,g, zc(i))
        p   = 0.5 + p/sqrt(2*PI)
        do l=1,ndir
            do j=1,nlag
                ii = (l-1)*nlag + j
                b = ro(ii)
                b = asin(b)
                ci(i,l,j) = simpson (0.0, b, 1.e-6, 40, f, zc(i))
                ri(i,l,j) = ( ci0-ci(i,l,j)) / (2*pi)
                rop(i,l,j) = ri(i,l,j)/(p*(1-p))
                lag(i,l,j) = j
                dir(i,l,j) = l
                zci(i,l,j) = zc(i)
                
            end do
        end do
    end do

end subroutine bigaus
