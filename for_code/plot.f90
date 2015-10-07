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
! include 'setrot.f90'
! include 'cova3.f90'
! include 'getindx.f90'
include 'locate.f90'
! this includes sortem and dsortem
include 'sortem.f90'
! include 'acorni.f90'
! include 'gauinv.f90'
! include 'powint.f90'
! include 'gcum.f90'


!***********************************************************************
! 
! 
!     Subroutines for plotting data using weight (qq-pp, cdf and histogram)
! 
! 
!***********************************************************************

subroutine  histplt(hmin,hmax, ncl, iwt, ilog, icum, nd, va, wt,  &
                    binval, nincls, cl, &
                    xpt025, xlqt, xmed, xuqt, xpt975, xmin, &
                    xmax, xcvr, xmen, xvar, xfrmx, dcl, &
                    error)
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

    !                           Histogram Plot
    !                           **************

    ! This program generates a PostScript file with a histogram and summary
    ! statistics.

    ! INPUT Parameters:

    !   hmin,hmax   plotting limts (will choose automatically if hmax<hmin)
    !   ncl         the number of classes
    !   ilog        1=log scale, 0=arithmetic
    !   icum        1=cumulative histogram, 0=frequency
    !   iwt         use weight variable?

    ! OUTPUT:
    !  xpt025, xlqt, xmed, xuqt, xpt975, xmin, xmax, xcvr, xmean, xvarm xfrmx
    !  ar3,ar4

    ! PROGRAM NOTES:



    ! The following Parameters control static dimensioning:

    !   MAXDAT    maximum number of data
    !   MAXCLS    maximum number of histogram classes



    !-----------------------------------------------------------------------
    
    implicit none    

    ! input 
    real*8, intent(in), dimension(nd) :: va, wt
    integer, intent(in) ::  ncl, nd,  ilog, icum, iwt
    real*8, intent(in) :: hmin, hmax

    ! output
    integer, intent(out) :: error
    ! the statistics
    real*8, intent(out) :: xpt025, xlqt, xmed, xuqt, xpt975, xmin, xmax, xcvr, xmen, xvar, xfrmx, dcl
    ! the actual histogram
    real*8, intent(out), dimension(ncl) :: binval, nincls, cl

    ! internal 
    logical ::  reghist,cumhist,connum
    real*8 :: EPSLON, vrmin, vrmax, art, xtwti, wt1, oldcp,cp, xtwt, thmin, thmax
    real*8, dimension (1) :: h, g, f, e, d, c

    real*8::       ar1(nd),ar2(nd)
    integer ::     j, i

    ! ar1 ==> data, ar2 ==> wight, binval 

    EPSLON=1.0e-20
    vrmin = 1.0e21
    vrmax =-1.0e21
    reghist = .TRUE. 
    cumhist = .FALSE. 
    connum  = .FALSE. 
    error = 0
    xtwt = 0.0
    xmen = 0.0
    thmin = hmin
    thmax = hmax

    if(icum /= 0) reghist = .FALSE. 
    if(icum == 1) cumhist = .TRUE.  
    if( .NOT. cumhist .AND. .NOT. connum) reghist = .TRUE. 


    if(nd <= 1) then 
        error = 1  ! 'ERROR: there is less than one datum ',nd
        return
    endif
    
    if(ncl < 1) then
        error = 100           ! ncl >=1 required 
        return
    endif

    ! Assign the defaults if necessary:

    if(hmax <= hmin) then
        thmin = vrmin
        thmax = vrmax
    endif

    do j=1, nd 
    
        ar1(j) = va(j)
        vrmin = min(ar1(nd),vrmin)
        vrmax = max(ar1(nd),vrmax)

        if(iwt >= 1) then
            ! Invalid wight?
            if(wt(j) <= EPSLON) then
                error = 10           ! weight too low, filter low wight before
                return
            endif
            ar2(j) = wt(j)
        else
            ar2(j) = 1.0
        endif

        xmen = xmen + ar1(j)*ar2(j)
        xtwt = xtwt + ar2(j)
        
    end do

    ! Get mean and total weight:
    xtwti = 1.0  / xtwt
    xmen  = xmen * xtwti

    ! Get the variance:
    xvar = 0.0
    do i=1,nd
        xvar = xvar + (ar1(i)-xmen) * (ar1(i)-xmen) * ar2(i)
    end do
    xvar  = xvar * xtwti


    ! Infer some of the histogram parameters:
    dcl = (thmax-thmin)/real(ncl+2)
    if(ilog == 1) then
        if(thmin <= 0.0) thmin = vrmin
        if(thmax <= 0.0) thmax = vrmax
        thmin = real(int(dlog10(max(thmin,EPSLON))-0.9))
        thmax = real(int(dlog10(max(thmax,EPSLON))+0.9))
    endif
    dcl  = (thmax-thmin)/real(ncl)
    

    ! Determine the histogram class structure:
    
    do i=1,ncl
        binval(i) = 0.0
        nincls(i) = 0
        cl(i) = thmin + i*dcl
        if(ilog == 1) cl(i) = 10**cl(i) 
    end do
    
    
    ! Here is where we build the histogram
    do i=1,nd
        if(ilog == 1) then
            art = dlog10(max(ar1(i),EPSLON))
        else
            art = ar1(i)
        endif
        wt1 = ar2(i) * xtwti
        j = ((art-thmin)/dcl)
        if(j < 1)   j = 1
        if(j > ncl) j = ncl
        binval(j)    =binval(j)    + wt1
        nincls(j) =nincls(j) + 1
    end do

    ! Sort the Data in Ascending Order:

    call dsortem(1,nd,ar1,1,ar2,c,d,e,f,g,h)

    ! Turn the weights into a cdf:

    oldcp = 0.0
    cp    = 0.0
    do i=1,nd
        cp     = cp + ar2(i) * xtwti
        ar2(i) =(cp + oldcp) * 0.5
        oldcp  = cp
    end do

    ! Obtain the quantiles:

    call dlocate(ar2,nd,1,nd,0.025,i)
    if(i == 0) then
        xpt025 = ar1(1)
    else if(i == nd) then
        xpt025 = ar1(nd)
    else
        xpt025 = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.025-ar2(i))/(ar2(i+1)-ar2(i))
    endif
    call dlocate(ar2,nd,1,nd,0.25,i)
    if(i == 0) then
        xlqt = ar1(1)
    else if(i == nd) then
        xlqt = ar1(nd)
    else
        xlqt = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.25-ar2(i))/(ar2(i+1)-ar2(i))
    endif
    call dlocate(ar2,nd,1,nd,0.50,i)
    if(i == 0) then
        xmed = ar1(1)
    else if(i == nd) then
        xmed = ar1(nd)
    else
        xmed = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.50-ar2(i))/(ar2(i+1)-ar2(i))
    endif
    call dlocate(ar2,nd,1,nd,0.75,i)
    if(i == 0) then
        xuqt = ar1(1)
    else if(i == nd) then
        xuqt = ar1(nd)
    else
        xuqt = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.75-ar2(i))/(ar2(i+1)-ar2(i))
    endif
    call dlocate(ar2,nd,1,nd,0.975,i)
    if(i == 0) then
        xpt975 = ar1(1)
    else if(i == nd) then
        xpt975 = ar1(nd)
    else
        xpt975 = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.975-ar2(i))/(ar2(i+1)-ar2(i))
    endif

    ! The coeficient of variation:

    xmin = ar1(1)
    xmax = ar1(nd)
    if(xmin < 0.0 .OR. xmen <= EPSLON) then
        xcvr = -1.0
    else
        xcvr = sqrt(max(xvar,0.0))/xmen
    endif

    ! Find the Maximum Class Frequency:

    xfrmx  = binval(1)
    do i=2,ncl
        xfrmx  = max(xfrmx,binval(i))
    end do


    ! Change things if we are considering a cumulative histogram:

    if(cumhist) then
        xfrmx = 1.0
        do i=2,ncl
            binval(i) = binval(i-1) + binval(i)
        end do
    end if


end subroutine  histplt



subroutine probplt( iwt, nd, va, wt,  &
                    binval, cl, &
                    xpt025, xlqt, xmed, xuqt, xpt975, xmin, &
                    xmax, xcvr, xmen, xvar, &
                    error)
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

    !                  Normal/Lognormal Probability Plot
    !                  *********************************

    ! Displays a set of data values on a probability plot with either an
    ! arithmetic or logarithmic scaling.

    ! INPUT/OUTPUT Parameters:
    !  iwt, nd, va, wt         


    ! PROGRAM NOTES:

    ! 1. The program is executed with no command line arguments.  The user
    !    will be prompted for the name of a parameter file.  The parameter
    !    file is described in the documentation (see example scatplt.par)

    ! 2. The logarithmic scaling is base 10

    ! 3. All acceptable values outside the plotting limits (including
    !    zero or negative values when a logarithmic scale is used) are
    !    used to establish the cdf but are not shown.



    ! The following Parameters control static dimensioning:

    !   MAXDAT    maximum number of data



    !-----------------------------------------------------------------------

    implicit none    

    ! input 
    real*8, intent(in), dimension(nd) :: va, wt
    integer, intent(in) ::  nd, iwt


    ! output
    integer, intent(out) :: error
    ! the statistics
    real*8, intent(out) :: xmen, xvar, xmed, xcvr, xmin, xmax, xpt025, xlqt, xuqt, xpt975
    ! the actual histogram
    real*8, intent(out), dimension(nd) :: binval, cl

    ! internal 
    real*8 :: EPSLON, xtwti, oldcp,cp, xtwt
    real*8, dimension (1) :: h, g, f, e, d, c
    real*8::       ar1(nd),ar2(nd)
    integer ::     j, i
 

    ! ar1 ==> data, ar2 ==> wight

    EPSLON=1.0e-20
    error = 0
    xtwt = 0.0
    xmen = 0.0


   if(nd <= 1) then 
        error = 1  ! 'ERROR: there is less than one datum ',nd
        return
    endif
    

    xtwt = 0.0
    xmen = 0.0
    do j=1, nd    
        ar1(j)=va(j)
        if(iwt >= 1) then
            ! Invalid wight?
            if(wt(j) <= EPSLON) then
                error = 10           ! weight too low, filter low wight before
                return
            endif
            ar2(j) = wt(j)
        else
            ar2(j) = 1.0
        endif
        ! Get mean and total weight:
        xmen = xmen + ar1(j)*ar2(j)
        xtwt = xtwt + ar2(j)
    end do
  
    if(xtwt < EPSLON) then
        error = 10         ! 'Cumulative Probability too LOW'
        return
    end if
    xtwti = 1.0  / xtwt
    xmen  = xmen * xtwti

    ! Get the variance:

    xvar = 0.0
    do i=1,nd
        xvar = xvar + (ar1(i)-xmen) * (ar1(i)-xmen) * ar2(i)
    end do
    xvar  = xvar * xtwti

   
    ! Sort the Data in Ascending Order:

    call dsortem(1,nd,ar1,1,ar2,c,d,e,f,g,h)


    ! Get cumulative probability and normalize:

    oldcp = 0.0
    cp    = 0.0
    do i=1,nd
        cp     = cp + ar2(i)*xtwti
        ar2(i) = 0.5*(cp+oldcp)
        oldcp  = cp
        binval(i) = ar2(i)
        cl(i) =      ar1(i)
    end do
    
    call dlocate(ar2,nd,1,nd,0.50,i)
    if(i == 0) then
        xmed = ar1(1)
    else if(i == nd) then
        xmed = ar1(nd)
    else
        xmed = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.50-ar2(i))/(ar2(i+1)-ar2(i))
    endif

    ! Obtain the quantiles:

    call dlocate(ar2,nd,1,nd,0.025,i)
    if(i == 0) then
        xpt025 = ar1(1)
    else if(i == nd) then
        xpt025 = ar1(nd)
    else
        xpt025 = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.025-ar2(i))/(ar2(i+1)-ar2(i))
    endif
    call dlocate(ar2,nd,1,nd,0.25,i)
    if(i == 0) then
        xlqt = ar1(1)
    else if(i == nd) then
        xlqt = ar1(nd)
    else
        xlqt = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.25-ar2(i))/(ar2(i+1)-ar2(i))
    endif
    call dlocate(ar2,nd,1,nd,0.75,i)
    if(i == 0) then
        xuqt = ar1(1)
    else if(i == nd) then
        xuqt = ar1(nd)
    else
        xuqt = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.75-ar2(i))/(ar2(i+1)-ar2(i))
    endif
    call dlocate(ar2,nd,1,nd,0.975,i)
    if(i == 0) then
        xpt975 = ar1(1)
    else if(i == nd) then
        xpt975 = ar1(nd)
    else
        xpt975 = ar1(i) +      (ar1(i+1)-ar1(i)) * &
        (0.975-ar2(i))/(ar2(i+1)-ar2(i))
    endif

    ! The coeficient of variation:
    xmin = ar1(1)
    xmax = ar1(nd)
    if(xmin < 0.0 .OR. xmen <= EPSLON) then
        xcvr = -1.0
    else
        xcvr = sqrt(max(xvar,0.0))/xmen
    endif

end subroutine probplt

subroutine qpplt(qqorpp, npts, &
                 n1, n2, va1,va2,wt1, wt2, &
                 vr1,vr2, error)
    !-----------------------------------------------------------------------

    !                        Q-Q, or P-P Plots
    !                        *****************

    ! This program generates values for Q-Q or P-P plot

    ! INPUT Parameters:

    !   va1,wt1     variable and the weight
    !   va2,wt2     variable and the weight
    !   qqorpp      0 = Q-Q plot, 1 = P-P plot
    !   npts        number of points to label

    ! Output 

    !   vr1, vr2   (npts), arrays with Q & Q or P & P, depending on qqorpp option

    !-----------------------------------------------------------------------

    implicit none


    ! input 
    integer, intent(in) :: qqorpp, npts, n1, n2
    real*8, intent(in), dimension(n1) :: va1, wt1
    real*8, intent(in), dimension(n2) :: va2, wt2

    ! output
    integer, intent(out) :: error
    real*8, intent(out), dimension(npts) :: vr1, vr2


    ! internal
    real*8 :: EPSLON, ccdf, cp1, cp2, cpinc, & 
             oldcp, cp, zz1,zz2, zz, zmin, zmax, zinc
    real*8,  dimension(n1) :: z1, p1
    real*8,  dimension(n2) :: z2, p2
    integer :: i, j1, j2, nd, npoints
    real*8, dimension (1) :: h, g, f, e, d, c  ! this are not used but required in sortem


    EPSLON=1.0e-20

    npoints = npts

    npoints = min(npoints,n1,n2)
    
    if (npoints /= npts) then 
        error = 100  ! the number of points may be equal or lowr than the number of points in the smaller dataset
        return
    end if

    ! using internal coppy of the input data to avoid inout variables
    do i=1, n1
        z1(i) = va1(i)
        p1(i) = wt1(i)
    end do
    do i=1, n2
        z2(i) = va2(i)
        p2(i) = wt2(i)
    end do


! Create Cumulative Probabilities out of Variables:


! CDF for first data set:

    call dsortem(1,n1,z1,1,p1,c,d,e,f,g,h)
    ccdf = 0
    do i=1,n1
        ccdf = ccdf + p1(i)
    end do
    if(ccdf < EPSLON) then
        error = 10  ! 'Cumulative Probability too LOW'
        return  
    end if
    oldcp = 0.0
    cp    = 0.0
    do i=1,n1
        cp    = cp + p1(i)/ccdf
        p1(i) = 0.5*(cp+oldcp)
        oldcp = cp
    end do

    ! CDF for second data set:

    call dsortem(1,n2,z2,1,p2,c,d,e,f,g,h)
    ccdf = 0
    do i=1,n2
        ccdf = ccdf + p2(i)
    end do
    if(ccdf < EPSLON) stop 'Cumulative Probability too LOW'
    oldcp = 0.0
    cp    = 0.0
    do i=1,n2
        cp    = cp + p2(i)/ccdf
        p2(i) = 0.5*(cp+oldcp)
        oldcp = cp
    end do

    ! Set up either a Q-Q or a P-P plot:

    nd = 0
    if(qqorpp == 0) then
    
        ! Q-Q plot:
    
        cpinc = 1.0 / real(npoints)
        cp    = -0.5*cpinc
        do i=1,npoints
            cp = cp + cpinc
            call dlocate(p2,n2,1,n2,cp,j2)
            call dlocate(p1,n1,1,n1,cp,j1)
            j2 = min((n2-1),max(1,j2))
            j1 = min((n1-1),max(1,j1))
            zz1 = z1(j1)+(z1(j1+1)-z1(j1))*(cp-p1(j1)) &
            / (p1(j1+1)-p1(j1))
            zz2 = z2(j2)+(z2(j2+1)-z2(j2))*(cp-p2(j2)) &
            / (p2(j2+1)-p2(j2))
            nd = nd + 1
            vr1(nd) = zz1
            vr2(nd) = zz2
        end do
    
    else
        ! P-P plot
        zmin = max(z1(1), z2(1) )
        zmax = min(z1(n1),z2(n2))
        zinc = (zmax-zmin) / real(npoints+1)
        zz   = zmin - 0.5*zinc
        do i=1,npoints
            zz = zz + zinc
            call dlocate(z1,n1,1,n1,zz,j1)
            call dlocate(z2,n2,1,n2,zz,j2)
            j2 = min((n2-1),max(1,j2))
            j1 = min((n1-1),max(1,j1))
            cp1 = p1(j1)+(p1(j1+1)-p1(j1))*(zz-z1(j1)) &
            / (z1(j1+1)-z1(j1))
            cp2 = p2(j2)+(p2(j2+1)-p2(j2))*(zz-z2(j2)) &
            / (z2(j2+1)-z2(j2))
            nd = nd + 1
            vr1(nd) = cp1
            vr2(nd) = cp2
        end do
    end if

end subroutine qpplt
