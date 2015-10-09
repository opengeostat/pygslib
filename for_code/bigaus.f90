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
include 'gauinv.f90'
! include 'powint.f90'
! include 'gcum.f90'


subroutine bigaus(c0,it,cc,aa,aa1,aa2,ang1,ang2,ang3,nst, &
                  ncut,zc,ndir,nlag,azm,dip,xlag, &
                  out_j,out_l,out_h,out_ri,out_ci,out_rop,out_p,zcc,error)
    !-----------------------------------------------------------------------

    !          Indicator Variograms for BiGaussian Distribution
    !          ************************************************


    ! This program returns the values of the theoretical indicator
    ! semivariograms corresponding to an input normal scores semivariogram

    ! INPUT/OUTPUT PARAMETERS:

    !   outfl                  the output file for vargplt
    !   ncut                   number of cutoffs
    !   cut()                  the cutoffs (in cdf units)
    !   ndir,nlag              number of directions, number of lags
    !   xoff,yoff,zoff         for each direction: the specification
    !   nst,c0                 Normal scores variogram: nst, nugget
    !   it,aa,cc               type, a parameter, c parameter
    !   ang1,ang2,...,anis2    anisotropy definition



    ! PROGRAM NOTES:

    !   1. Setting the number of cutoffs to a -1 times the number of
    !      cutoffs will cause the variograms to be standardized to a sill
    !      of one.



    ! EXTERNAL REFERENCES:

    !       f        to calculate exp(-zc**2/(1+x))/sqrt(1-x**2)
    !       g        to calculate exp(-x**2/2)
    !       simpson  to do the numerical calculation
    !       gauinv   to get standard normal deviate



    ! Original:  H. Xiao                                                1975
    !-----------------------------------------------------------------------


    implicit none    


    !input 
    integer, intent(in) :: ncut,ndir,nlag,nst
    
    real, intent(in), dimension(nst) :: aa,aa1,aa2,cc,ang1,ang2,ang3
    integer, intent(in), dimension(nst) :: it
    real, intent(in), dimension(1) :: c0
    real, intent(in), dimension(ncut) :: zcc
    real, intent(in), dimension(ndir) :: xlag,azm,dip
    

    !ouput
    integer, intent(out), dimension(nlag,ndir,ncut) :: out_j, out_l
    real, intent(out), dimension(nlag,ndir,ncut) :: out_h,out_ri,out_ci,out_rop
    real, intent(out), dimension(ncut) :: out_p, zc
    
    integer, intent(out) :: error

    !internal
    real ::  b, p, maxcov
    real, dimension(ndir) :: xoff,yoff,zoff
    real, dimension(nlag*ndir) :: h,ri,ci,ro,gam,rop
    real*8, dimension(nst,3,3) :: rotmat
    integer :: i, id, j, l, is, icut, ii, il, ierr
    real :: xp, cov, zz, yy, xx, cmax, ci0
    real :: PI, DEG2RAD, EPSLON
    real, dimension(nst) :: anis1,anis2
    real ::  simpsong, simpsonf

    PI=3.14159265
    DEG2RAD=PI/180.0
    EPSLON = 1.0e-20


    if(ncut < 0) then
        error = -10 ! ncut may be positive 
        return 
    endif

    do i=1,ndir
        xoff(i) = sin(DEG2RAD*azm(i))*xlag(i)*cos(DEG2RAD*dip(i))
        yoff(i) = cos(DEG2RAD*azm(i))*xlag(i)*cos(DEG2RAD*dip(i))
        zoff(i) = sin(DEG2RAD*dip(i))*xlag(i)
        write(*,*) ' x,y,z offsets = ',xoff(i),yoff(i),zoff(i)
    end do

    do i=1,nst
        anis1(i) = aa1(i) / max(aa(i),EPSLON)
        anis2(i) = aa2(i) / max(aa(i),EPSLON)
    end do

    ! Set up the rotation matrices:

    do is=1,nst
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,nst,rotmat)
    end do


    ! Convert cutoffs cumulative probabilities to standard normal cutoffs:

    do icut=1,ncut
        zc(icut)=zcc(icut)
        call gauinv(dble(zc(icut)),xp,ierr)
        zc(icut) = xp
    end do

    ! Get ready to go:

    call cova3(0.0,0.0,0.0,0.0,0.0,0.0,1,nst,c0,it,cc,aa, &
    &            1,nst,rotmat,cmax,maxcov)


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
            print *, '@@@', maxcov , cov
        end do
    end do

    ! Now, loop over all the cutoffs:

    do i=1,ncut
        ci0 = simpsonf(0.0, pi/2.0, 1.e-6, 40, zc(i))
        p   = simpsong(0.0,zc(i),1.e-6,40)
        p   = 0.5 + p/sqrt(2*PI)
        out_p(i) = p
        do l=1,ndir
            do j=1,nlag
                ii = (l-1)*nlag + j
                b = ro(ii)
                b = asin(b)
                ci(ii) = simpsonf (0.0, b, 1.e-6, 40, zc(i))
                ri(ii) = ( ci0-ci(ii)) / (2*pi)
                rop(ii) = ri(ii)/(p*(1-p))

                out_j(j,l,i) = j
                out_h(j,l,i) = h(ii)
                out_ri(j,l,i) = ri(ii)
                out_l(j,l,i) = l
                out_ci(j,l,i) = ci(ii)
                out_rop(j,l,i) = rop(ii)

            end do
        end do
    end do

    ! Finished:

end subroutine bigaus



real function f(x, zc)
    !-----------------------------------------------------------------------

    ! This function calculates the values of the function
    !       f( zc, xi)= exp ( -zc**2/ (1+sin(x))) for zc and x

    !-----------------------------------------------------------------------
    real ::       x, zc
    f = exp ( -zc**2 / (1+sin(x)))
    return
end function f


real function g(x)
    real :: x
    g= exp ( - x**2/2.0)
    return
end function g



function simpsonf(alim,blim,tol,nitn,zc)
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
    if(nitn > 50) stop 'Sorry, I only do 50 iterations!'

    a = alim
    b = blim
    eps = tol
    lvl = 0
    est = 1.0
    absarea = 1.0
    da = b - a
    fm = 4.0 * f((a+b)/2.0,zc)
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

    simpsonf = sum
    return
end function simpsonf

function simpsong(alim,blim,tol,nitn)
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
    if(nitn > 50) stop 'Sorry, I only do 50 iterations!'

    a = alim
    b = blim
    eps = tol
    lvl = 0
    est = 1.0
    absarea = 1.0
    da = b - a
    fm = 4.0 * g((a+b)/2.0)
    fa = g(a)
    fb = g(b)

    10 continue
    lvl = lvl + 1
    dx(lvl) = da / 3.0
    sx = dx(lvl) / 6.0
    f1 = 4.0 * g(a+dx(lvl)/2.0)
    x2(lvl) = a + dx(lvl)
    f2(lvl) = g(x2(lvl))
    x3(lvl) = x2(lvl) + dx(lvl)
    f3(lvl) = g(x3(lvl))
    epsp(lvl) = eps
    f4(lvl) = 4.0 * g(x3(lvl) + dx(lvl)/2.0)
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

    simpsong = sum
    return
end function simpsong
