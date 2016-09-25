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
! include 'sqdist.f90'
! include 'getindx.f90'
include 'locate.f90'
! this includes sortem and dsortem
! include 'sortem.f90'
include 'acorni.f90'
! include 'gauinv.f90'
! include 'powint.f90'
! include 'gcum.f90'
! include 'getz.f90'

subroutine draw(vr,wt,nd,nv,rseed,ndraw,vo,sumwts,error)
    !-------------------------------------------------------------------

    !       Simple Stochastic Simulation (Random Drawing) Program
    !       *****************************************************

    ! Monte Carlo simulation with replacement of observations in an input
    ! data file - could be easily modified to include a "bootstrap"
    ! statistic



    !-------------------------------------------------------------------

    implicit none

    ! input 
    integer, intent (in):: nd, nv,rseed, ndraw
    real, intent (in), dimension(nd) :: wt
    real, intent (in), dimension(nd,nv) :: vr


    ! output 
    real, intent (out):: sumwts, error
    real, intent (out), dimension(ndraw,nv) :: vo
    

    ! internal 
    real*8 ::  acorni, ixv(12)
    real :: zz, EPSLON, prob(nd), cdf, cd, oldcp, cp
    integer :: irepo, i, j, KORDEI, MAXOP1, ivar


    EPSLON=1.0e-10

    ! ACORN parameters:
    KORDEI=12 ! this is hard coded... ixv (KORDEI=12)
    MAXOP1 = KORDEI+1


    do i=1,MAXOP1
        ixv(i) = 0.0
    end do
    ixv(1) = rseed
    do i=1,10000
        zz = real(acorni(ixv, KORDEI))
    end do   
  
    ! Read through the input file:


    sumwts = 0.0

    do i=1, nd
        sumwts   = sumwts +  wt(i)
        prob(i) = wt(i)
    end do


    if(nd < 2 .OR. sumwts < EPSLON) then
        error = 10 ! ' too few data or too low sum of weights '
        stop
    end if

    ! Turn the data distribution into a CDF:

    oldcp   = 0.0
    cp      = 0.0
    sumwts  = 1.0 / sumwts
    do i=1,nd
        cp      = cp + prob(i) * sumwts
        prob(i) =(cp + oldcp)  * 0.5
        oldcp   = cp
    end do


    ! Make the drawings and write them to the output file:

    irepo = max(1,min((ndraw/10),10000))
    do i=1,ndraw

        cdf = real(acorni(ixv, KORDEI))
        call locate(prob,nd,1,nd,cdf,j)
        j = max(min(j,nd),1)
        do ivar = 1, nv
            vo(i,ivar) =vr(j,ivar)    ! try to improve this? the best practice for slicing is vr(:,l)...
        end do 
    end do
end subroutine draw


