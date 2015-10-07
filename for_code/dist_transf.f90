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
include 'acorni.f90'
include 'gauinv.f90'
include 'powint.f90'
include 'gcum.f90'

!***********************************************************************
! 
! 
!     Subroutines for data transformation (nscore and bactransform)
! 
! 
!***********************************************************************

!---------------------------------------------------------------------------
!     Subroutine ns_ttable (this one part of the nscore program)
!---------------------------------------------------------------------------
subroutine ns_ttable(va,wt, nd, transin, transout, error)
    ! 
    ! This code is based on GSLIB nscore  
    !
    !-----------------------------------------------------------------------

    !                Compute Normal Scores transformation table
    !                ***********************************

    ! PROGRAM NOTES:

    !   1. Create a transformation table from data or reference distribution

    ! Version from Adrian martinez Vargas, 2015

    !-----------------------------------------------------------------------

    implicit none    

    ! input 
    real*8, intent(in), dimension (nd) :: va, wt
    integer, intent(in) :: nd

    ! output
    real*8, intent (out), dimension (nd)  :: transin, transout
    integer, intent(out) :: error

    ! internal
    real*8 :: EPSLON
    real*8 :: ixv (12)
    real*8, dimension (nd) :: vr, wt_ns
    real*8 ::  twt,wtfac,w,cp,oldcp,vrg,acorni, p
    real :: vrrg
    real*8, dimension(1) :: c,d,e,f,g,h    ! these are dummies for dsortem
    integer :: KORDEI, MAXOP1, i, istart, iend, ierr, j
   

    EPSLON=1.0e-6

    ! ACORN parameters:
    KORDEI=12 ! this is hard coded... ixv (KORDEI=12)
    MAXOP1 = KORDEI+1


    do i=1,MAXOP1
        ixv(i) = 0.0
    end do
    ixv(1) = 69069
    do i=1,1000
        p = real(acorni(ixv, KORDEI))
    end do


    ! calculat total wight
    twt=0
    do i=1, nd
        if(wt(nd) <= 1.0e-10) then
               error = 10         ! Wight too small review your data and try again  
            return 
        end if
        wt_ns(i) = wt(i)
        vr(i) = va(i) + acorni(ixv, KORDEI)*dble(EPSLON)
        twt = twt + wt_ns(i)
    end do

    if(nd <= 1 .OR. real(twt) <= EPSLON) then
        error = 100       !'ERROR: too few data'
        return
    endif

    ! Sort data by value:

    istart = 1
    iend   = nd
    call dsortem(istart,iend,vr,1,wt_ns,c,d,e,f,g,h)

    ! Compute the cumulative probabilities and write transformation table

    wtfac = 1.0/twt
    oldcp = 0.0
    cp    = 0.0
    do j=istart,iend
        w     =  wtfac*wt_ns(j)
        cp    =  cp + w
        wt_ns(j) = (cp + oldcp)/2.0
        call gauinv(wt_ns(j),vrrg,ierr)
        vrg = dble(vrrg)
        transin(j) = vr(j)
        transout(j) = vrg
        oldcp =  cp
        ! Now, reset the weight to the normal scores value:
        ! wt_ns(j) = vrg   
    end do

end subroutine ns_ttable

subroutine nscore(va, nd, transin, transout, nt, getrank , nsc) 

    !-----------------------------------------------------------------------

    !                Compute Normal Scores of a Data Set
    !                ***********************************

    ! PROGRAM NOTES:

    !   2. Random Despiking (THIS IS WHY WE NEED A RANDOM NUMBER GENERATOR)



    !-----------------------------------------------------------------------

    implicit none

    ! input 
    integer, intent(in) :: nd, nt
    logical, intent(in):: getrank
    real*8, intent (in), dimension (nt)  :: transin, transout
    real*8, intent(in), dimension (nd) :: va

    ! output
    real*8, intent(out), dimension(nd):: nsc
    

    ! internal 
    real*8 :: doubone, dpowint
    real*8 :: EPSLON
    real*8 :: ixv (12)
    integer :: KORDEI, MAXOP1, i, j
    real*8 :: vrg,acorni, p, yy, vrr
    real ::  pp, gcum
    
    doubone=1.0
    EPSLON=1.0e-6

    ! ACORN parameters:
    KORDEI=12 ! this is hard coded... ixv (KORDEI=12)
    MAXOP1 = KORDEI+1

    ! TODO: review this part... ixv was suppose to be in shared memory space (use inout instad?)
    do i=1,MAXOP1
        ixv(i) = 0.0
    end do
    ixv(1) = 69069
    do i=1,1000
        p = real(acorni(ixv, KORDEI))
    end do

    ! Normal Scores Transform:     
    do i=1, nd 
    
        vrr = dble(va(i))+acorni(ixv, KORDEI)*dble(EPSLON)
    
        ! Now, get the normal scores value for "vrr"
    
        call dlocate(transin,nt,1,nt,vrr,j)
        j   = min(max(1,j),(nd-1))
        vrg = dpowint(transin(j),transin(j+1),transout(j),transout(j+1),vrr,doubone)

        ! Write the rank order instead of the rank?

        if(getrank) then
            pp  = real(vrg)
            yy  = gcum(pp)
            vrg = dble(yy)
        end if
    
        nsc(i) = vrg

    end do

end subroutine nscore

real*8 function backtr(vrgs,nt,vr,vrg,zmin,zmax,ltail,ltpar, utail,utpar)
    !-----------------------------------------------------------------------

    !           Back Transform Univariate Data from Normal Scores
    !           *************************************************

    ! This subroutine backtransforms a standard normal deviate from a
    ! specified back transform table and option for the tails of the
    ! distribution.  Call once with "first" set to true then set to false
    ! unless one of the options for the tail changes.



    ! INPUT VARIABLES:

    !   vrgs             normal score value to be back transformed
    !   nt               number of values in the back transform tbale
    !   vr(nt)           original data values that were transformed
    !   vrg(nt)          the corresponding transformed values
    !   zmin,zmax        limits possibly used for linear or power model
    !   ltail            option to handle values less than vrg(1):
    !   ltpar            parameter required for option ltail
    !   utail            option to handle values greater than vrg(nt):
    !   utpar            parameter required for option utail



    !-----------------------------------------------------------------------

    real*8, dimension(nt) :: vr, vrg
    real*8 ::    ltpar,utpar,lambda, EPSLON, vrgs,zmin,zmax, dpowint, &
                 cdflo, cdfbt, cpow, cdfhi
    real:: gcum
    integer ::   ltail,utail, nt

    EPSLON=1.0e-20

    ! Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):

    if(vrgs <= vrg(1)) then
        backtr = vr(1)
        cdflo  = gcum(vrg(1))
        cdfbt  = gcum(vrgs)
        if(ltail == 1) then
            backtr = dpowint(dble(0.0),cdflo,zmin,vr(1),cdfbt,dble(1.0))
        else if(ltail == 2) then
            cpow   = 1.0 / ltpar
            backtr = dpowint(dble(0.0),cdflo,zmin,vr(1),cdfbt,cpow)
        endif
    
    ! Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
    
    else if(vrgs >= vrg(nt)) then
        backtr = vr(nt)
        cdfhi  = gcum(vrg(nt))
        cdfbt  = gcum(vrgs)
        if(utail == 1) then
            backtr = dpowint(cdfhi,dble(1.0),vr(nt),zmax,cdfbt,dble(1.0))
        else if(utail == 2) then
            cpow   = 1.0 / utpar
            backtr = dpowint(cdfhi,dble(1.0),vr(nt),zmax,cdfbt,cpow)
        else if(utail == 4) then
            lambda = (vr(nt)**utpar)*(1.0-gcum(vrg(nt)))
            backtr = (lambda/(1.0-gcum(vrgs)))**(1.0/utpar)
        endif
    else
    
    ! Value within the transformation table:
    
        call dlocate(vrg,nt,1,nt,vrgs,j)
        j = max(min((nt-1),j),1)
        backtr = dpowint(vrg(j),vrg(j+1),vr(j),vr(j+1),vrgs,dble(1.0))
    endif

    return 

end function backtr


!---------------------------------------------------------------------------
!     Subroutine bacnscore (this is the bactr program for bac transformation of normal score variables)
!---------------------------------------------------------------------------
subroutine bacnscore(vnsc, nd, transin, transout, nt, ltail,utail, ltpar,utpar, zmin,zmax, getrank, va, error)
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

    !                     Gaussian Back Transformation
    !                     ****************************

    ! PROGRAM NOTES:

    !  1. ltail, utail options: 1=linear interpolation, 2=power model
    !     interpolation, and 4=hyperbolic model interpolation (only for
    !     upper tail)



    ! EXTERNAL REFERENCES:

    !   gcum     Inverse of Gaussian cumulative distribution function
    !   locate   locate a position in an array
    !   dpowint   power law interpolation (real*8)



    !-----------------------------------------------------------------------



    implicit none

    ! input 
    integer, intent(in) :: nd, nt
    logical, intent(in) ::  getrank    
    integer, intent(in) :: ltail,utail
    real*8,  intent(in) :: ltpar,utpar, zmin,zmax
    real*8, intent(in),   dimension(nd) :: vnsc
    real*8, intent (in), dimension (nt) :: transin, transout


    ! output
    real*8, intent(out), dimension(nd):: va
    integer, intent(out) :: error
    
    ! internal 
    real*8 :: backtr, p, bac
    real*8 :: EPSLON
    integer :: i, ierr

    

    real*8, dimension(nd) :: vr,vrg


    EPSLON=0.00001


    ! Check for error situation:

    if(ltail /= 1 .AND. ltail /= 2) then
        error = 10 ! 'ERROR invalid lower tail option only allow 1 or 2 - see manual '
        return
    endif
    if(utail /= 1 .AND. utail /= 2 .AND. utail /= 4) then
        error = 11 ! 'ERROR invalid upper tail option, only allow 1,2 or 4 - see manual '
        return
    endif
    if(utail == 4 .AND. utpar < 1.0) then
        error = 12 ! 'ERROR invalid power for hyperbolic tail, must be greater than 1.0!'
        return
    endif

    ! Read in the transformation table:

    do i=1, nt
        vr(i)  = transin(i)
        vrg(i) = transout(i)
        if(i > 1) then
            if(vr(i) < vr(i-1) .OR. vrg(i) < vrg(i-1)) then
                error = 100 ! 'ERROR transformation table must be monotonic increasing! '
                return 
            endif
        endif
    end do

    ! Check for error situation:

    if(utail == 4 .AND. vr(nt) <= 0.0) then
        error = 200   ! 'ERROR can not use hyperbolic tail with negative values! - see manual '
        return
    endif
    if(zmin > vr(1)) then
        error = 210   ! 'ERROR zmin should be no larger than the first entry in the transformation table '
        return
    endif
    if(zmax < vr(nt)) then
        error = 220   ! 'ERROR zmax should be no less than the last entry in the transformation table '
        return
    endif

    ! Now do the transformation using backtr function

    do i=1, nd

        if(getrank) then
            p = vnsc(i)
            call gauinv(p,vnsc(i),ierr)
        end if
        bac = backtr(vnsc(i),nt,vr,vrg,zmin,zmax,ltail,ltpar,utail,utpar)
        if(bac < zmin) bac = zmin
        if(bac > zmax) bac = zmax
        
        va(i)=bac

    end do
    
    return

end subroutine bacnscore
