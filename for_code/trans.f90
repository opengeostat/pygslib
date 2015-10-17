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
include 'sortem.f90'
include 'acorni.f90'
! include 'gauinv.f90'
include 'powint.f90'
! include 'gcum.f90'
include 'getz.f90'

subroutine trans(ivtype,nr,vr,wt,no,vo,wo,nx,ny,nz,wx,wy,wz, &
                 nxyza,zmin,zmax,ltpar,utpar,ltail,utail,ldata,kv,ef,rseed, &
                 gmedian,rvr,rcdf,ncut,zval,error)

    !-----------------------------------------------------------------------

    !                      Univariate Transformation
    !                      *************************

    ! Reads in a reference distribution and a number of other distributions
    ! and then transforms the values in each of the second distributions
    ! such that their histograms match that of the reference distribution.



    ! INPUT/OUTPUT Parameters:

    !   ivtype      variable type (1=continuous, 0=categorical)
    !   nx, ny, nz  size of categorical variable realizations to transform
    !   wx, wy, wz  window size for breaking ties
    !   nxyz        size to of continuous variable data set to transform
    !   zmin,zmax   minimum and maximum data values
    !   ltail,ltpar lower tail: option, parameter
    !   utail,utpar upper tail: option, parameter
    !   ldata       honor local data (1=yes, 0=no)
    !   wtfac       control parameter


    ! This version use only one ralization

    !-----------------------------------------------------------------------

    implicit none 

    ! input 
    integer, intent(in) :: ivtype              ! type of distribution 1=continuous, 0=categorical
    integer, intent(in) :: nr
    real, intent(in), dimension(nr) :: vr, wt  ! reference distribution, variable and weight
    integer, intent(in) :: no
    real, intent(in), dimension(no) :: vo, wo  ! original distributions variable and weight                            
    integer, intent(in) :: nx, ny, nz          ! categorical: nx, ny, nz: size of 3-D model
    integer, intent(in) :: wx, wy, wz             !     window size for tie-breaking
    integer, intent(in) :: nxyza                ! continuous: number to transform per "set"
    real, intent(in) :: zmin,zmax              !     minimum and maximum values
    real, intent(in) :: ltpar,utpar            !     lower, upper tail parameter
    integer, intent(in) :: ltail,utail         !     lower, upper tail function type 
    integer, intent(in) :: ldata               ! honor local data? (1=yes, 0=no)
    real, intent(in), dimension(nr) :: kv      ! estimation variance
    real, intent(in) ::    ef               ! control parameter or scaling factor ( 0.33 < w < 3.0 )
    integer, intent(in) ::     rseed           ! random number seed (conditioning cat.)


    ! output 
    integer, intent(out) :: error
    real, intent(out) :: gmedian                  ! median (continous)
    real, intent(out), dimension(nr) :: rvr,rcdf  ! category and cumulative proportion (categorical)
    integer, intent(out) ::                ncut      ! num of categories (use rvr[1:ncut] or rcdf[1:ncut] )
    real, intent(out), dimension(no) :: zval      ! the transformed variable


    ! internal 
    real ::  dcdf(no),dvr(no), &
             indx(no),fuzzcat(no),catcdf(nr)
    integer ::   category(nr), KORDEI, MAXOP1, i, j, &
                 icut,iix,iiy,iiz,isim,iwt,iwtd,ix,iy,iz, &
                 jx,jy,jz, nloc, num, nxy, nxyz
    real :: EPSLON, wtfac, wtt, tcdf, evmax, cp, oldcp ,rn,tmax,powint,getz
    real*8 ::  acorni, ixv(12)
    real, dimension(1) :: c,d,e,f,g,h

    wtfac = ef  
    EPSLON=1.0e-12
    tmax = 1.0e21
    nxyz = nxyza

    ! ACORN parameters:
    KORDEI=12 ! this is hard coded... ixv (KORDEI=12)
    MAXOP1 = KORDEI+1
    

    if(ldata == 1) then    
        ! Scale from 0 to 1  --> 0.33 to 3.0
        wtfac = 0.33 + wtfac*(3.0-0.33)

        do i=1,MAXOP1
            ixv(i) = 0.0
        end do
        ixv(1) = rseed
        do i=1,10000
            rn = real(acorni(ixv, KORDEI))
        end do    
    end if

    ! Check for error situation:

    if(ltail /= 1 .AND. ltail /= 2) then
        error =10 ! 'ERROR invalid lower tail option, only allow 1 or 2 - see manual '
        stop
    endif
    if(utail /= 1 .AND. utail /= 2 .AND. utail /= 4) then
        error =11 ! 'ERROR invalid upper tail option, only allow 1,2 or 4 - see manual '
        stop
    endif
    if(utail == 4 .AND. utpar < 1.0) then
        error =12 ! 'ERROR invalid power for hyperbolic tail must be greater than 1.0!'
        stop
    endif


    ! Read as much data for reference distribution as possible:

    ncut = 0
    tcdf = 0
    do j=1, nr

        wtt = 1.0
        if(iwt >= 1) wtt = wt(j)

        if(ivtype == 1) then
            ncut = ncut + 1
            rvr(ncut)  = vr(j)
            rcdf(ncut) = wtt
        else
            do i=1,ncut
                if(int(vr(j)+0.5) == int(rvr(i))) then
                    icut = i
                    go to 4
                end if
            end do
            ncut = ncut + 1
            rvr(ncut)  = vr(j)
            rcdf(ncut) = 0.0
            icut       = ncut
            4 continue
            rcdf(icut) = rcdf(icut) + wtt
        end if
        tcdf = tcdf + wtt

    end do

    ! Sort the Reference Distribution and Check for error situation:

    call sortem(1,ncut,rvr,1,rcdf,c,d,e,f,g,h)
    if(ncut <= 1 .OR. tcdf <= EPSLON) then
        error =100 ! 'ERROR: too few data or too low weight'
        return
    endif
    if(ivtype == 1 .AND. utail == 4 .AND. rvr(ncut) <= 0.0) then
        error = 200 !'ERROR can not use hyperbolic tail with negative values! - see manual '
        return 
    endif

    ! Turn the (possibly weighted) distribution into a cdf that is useful:

    tcdf  = 1.0 / tcdf
    if(ivtype == 1) then
        oldcp = 0.0
        cp    = 0.0
        do i=1,ncut
            cp     = cp + rcdf(i) * tcdf
            rcdf(i) =(cp + oldcp) * 0.5
            oldcp  = cp
        end do
    else
        do i=1,ncut
            rcdf(i) = rcdf(i) * tcdf
        end do
    end if

    ! Write Some of the Statistics to the screen:

    if(ivtype == 1) then
        call locate(rcdf,ncut,1,ncut,0.5,j)
        gmedian = powint(rcdf(j),rcdf(j+1),rvr(j),rvr(j+1),0.5,1.0)
    else
        do i=1,ncut
            catcdf(i) = 0.0
            if(i > 1) rcdf(i) = rcdf(i) + rcdf(i-1)
        end do
    end if


    ! MAIN Loop over all the increments to transform:

    if(ivtype == 0) then
        nxyz = nx*ny*nz
        nxy  = nx*ny
    end if

    if(nxyz/=nxyza) then
        error =-1 ! 'Warning:  nx*ny*nz/=nxyza'
    endif
    
    if(nxyz/=nr) then
        error =-2 ! 'Warning:  nx*ny*nz/=nr'
    endif
    
    
    ! Read in the data values:

    tcdf = 0.0
    num  = 0
    do i=1,no
        num = num + 1
        if (no<num) then 
            error = 1000  ! bad parameters, exceding array length no<num 
            return 
        end if 
        dvr(num)  = vo(num)
        indx(num) = real(num)
        wtt = 1.0
        if(iwtd >= 1) wtt = vo(num)
        dcdf(num) = wtt
        tcdf      = tcdf + wtt

    
        ! Keep track of the proportions if working with a categorical variable:
    
        if(ivtype == 0) then
            do j=1,ncut
                if(int(dvr(num)+0.5) == int(rvr(j))) then
                    icut = j
                    go to 6
                end if
            end do
            error=2000 ! 'Found a code not found in reference'   write(*,*) dvr(num)
            return 
            6 continue
            catcdf(icut)  = catcdf(icut) + wtt
        end if
    end do

    if(tcdf <= EPSLON) then
        error=3000 ! 'ERROR: no data'
        return
    endif
    if(ivtype == 0 .AND. num /= nxyz) then
        error=400 ! 'ERROR: you must have nxyz for a categorical transformation'
        return
    endif

    ! For categorical transformation we need a fuzzy category which is the
    ! local average category - this is for breaking ties in the
    ! transformation procedure.

    if(ivtype == 0) then
        do i=1,nxyz
            fuzzcat(i) = 0.0
            nloc       = 0
            iz = int((i-1)/nxy) + 1
            iy = int((i-(iz-1)*nxy-1)/nx) + 1
            ix = i- (iz-1)*nxy - (iy-1)*nx
            do iix=-wx,wx
                do iiy=-wy,wy
                    do iiz=-wz,wz
                        jx = ix + iix
                        jy = iy + iiy
                        jz = iz + iiz
                        if(jx >= 1 .AND. jx <= nx .AND. &
                        jy >= 1 .AND. jy <= ny .AND. &
                        jz >= 1 .AND. jz <= nz) then
                            j = jx + (jy-1)*nx + (jz-1)*nxy
                            nloc = nloc + 1
                            fuzzcat(i) = fuzzcat(i) + dvr(j)
                        end if
                    end do
                end do
            end do
            if(nloc > 0) then
                fuzzcat(i) = fuzzcat(i) / real(nloc)
            else
                fuzzcat(i) = 2*tmax
            end if
        end do
    end if

    ! Turn the (possibly weighted) data distribution into a useful cdf:

    if(ivtype == 1) then
        call sortem(1,num,dvr,2,dcdf,indx,d,e,f,g,h)
    else
        call sortem(1,num,fuzzcat,3,dcdf,dvr,indx,e,f,g,h)
    end if
    oldcp = 0.0
    cp    = 0.0
    tcdf  = 1.0 / tcdf
    do i=1,num
        cp     = cp + dcdf(i) * tcdf
        dcdf(i) =(cp + oldcp) * 0.5
        oldcp  = cp
    end do

    ! Now, get the right order back:

    if(ivtype == 1) then
        call sortem(1,num,indx,2,dcdf,dvr,d,e,f,g,h)
    else
        call sortem(1,num,indx,3,dcdf,dvr,fuzzcat,e,f,g,h)
    end if

    ! Read in the kriging variance to array "indx" if we have to honor
    ! local data:

    if(ldata == 1) then
        evmax = -1.0e21
        do i=1,num
            indx(i) = kv(i)
            if(indx(i) >= 0.0) then
                indx(i) = sqrt(max(indx(i),0.0))
                if(indx(i) > evmax) evmax = indx(i)
            end if
        end do
    end if

    ! Go through all the data back transforming them to the reference CDF:

    do i=1,num

        if(ivtype == 1) then
            zval(i) = getz(dcdf(i),ncut,rvr,rcdf,zmin, &
            zmax,ltail,ltpar,utail,utpar)
        else
            zval(i) = rvr(1)
            do j=2,ncut
                if(dcdf(i) > rcdf(j-1)) &
                zval(i) = rvr(j)
            end do
        end if
    
        ! Now, do we have to honor local data?
    
        if(ldata == 1) then
            if(indx(i) >= 0) then
                wtt = (indx(i)/evmax)**wtfac
                if(ivtype == 1) then
                    zval(i) = dvr(i)+wtt*(zval(i)-dvr(i))
                else
                    rn = real(acorni(ixv, KORDEI))
                    if(rn > wtt) zval(i) = dvr(i)
                end if
            end if
        end if

    end do
    

    


end subroutine trans
