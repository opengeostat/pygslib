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
! include 'locate.f90'
! this includes sortem and dsortem
! include 'sortem.f90'
! include 'acorni.f90'
! include 'gauinv.f90'
! include 'powint.f90'
! include 'gcum.f90'


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
    !   noff            number of origin offsets (regular spacing)


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
    real*8 :: roff, vrav, best, facto, &
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
