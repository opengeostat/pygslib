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
                rnum   = int(np(i))
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



!---------------------------------------------------------------------------
!     Subroutine gamv
!---------------------------------------------------------------------------
subroutine gamv3D(nd, x, y, z, bhid, nv, vr, &                   ! data array and coordinares
                tmin, tmax, nlag, xlag, &                        ! lag definition
                ndir, ndip, orgdir, orgdip, &                     ! directions and parameters
                isill, sills, nvarg, ivtail, ivhead,ivtype,  &   ! variograms types and parameters
                np, dis, gam, hm, tm, hv, tv)                    ! output

    !----------------------------------------------------------------------
    !
    ! Calculate variograms in all posible direction and place it in a
    ! 3D spheric regular grid. The result can be ploted later in a
    ! VTK structured grid
    !
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! This code was modified from original f77 GSLIB code (v.2) gamv function
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
    !   ndir, ndip       Number of directions, inclinations to consider
    !                    angles of direction  are measured positive degrees clockwise from NS.
    !   orgdir, orgdip   clockwise rotation added to dir and
    !                    down dip positive angle added to dip
    !                    dir rotation is applied and the dip rotation
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
    ! OUTPUT VARIABLES:  The following arrays are stored in (nlag,ndir,ndip,nvarg) arrays

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

    integer, intent(in)                     :: nlag, ndir, ndip, nvarg
    integer, intent(in)                     :: nd, nv, isill
    integer, intent(in), dimension(nvarg)   :: ivtail, ivhead,ivtype
    real*8, intent(in), dimension(nd)       :: x, y, z

    real*8, intent(in), dimension(nd,nv)    :: vr
    real*8, intent(in), dimension(nv)       :: sills
    real*8, intent(in)                      :: tmin, tmax, xlag, orgdir, orgdip

    real*8, intent(out), dimension(nlag,ndir,ndip,nvarg)  :: np, dis, gam, hm, tm, hv, tv



    !new variables
    integer, intent(in), dimension(nd)         :: bhid


    ! some general declarations
    real*8 :: PI, EPSLON
    parameter(PI=3.14159265, EPSLON=1.0e-20)

    real*8                :: uvxazm, uvyazm, uvzdec, uvhdec, csatol, csdtol
    logical               :: omni


    !Extra Variables not declared in the original library
    real*8  :: azmuth, declin, band, dcazm, dismxs, dx, dxs, dxy, dy, &
             dys, dz, dzs, dcdec, gamma, h, hs, htave, variance, &
             vrh, vrhpr, vrt,  vrtpr, xltoll
    integer :: i, id, ii, iii, il, ilag, it, iv, j, lagbeg, lagend, nsiz, rnum

    real*8 ::  dazm, ddip, dp, az, azimdg
    integer :: idp, jdp, k, kdip, jdir, l

    write (*,*) 'Nlag: ', nlag
    write (*,*) 'Ndir: ', ndir
    write (*,*) 'Ndir: ', ndip
    write (*,*) 'Nvarg:', nvarg

    ! the azimuth separation
    if (ndir<1) then
        dazm = 180
    else
        dazm = 180./ndir
    end if

  ! the dip separation
    if (ndip<1) then
        ddip = 180.
    else
        ddip = 180./(ndip)
    end if

    ! Initialize the arrays for each direction, variogram, and lag:

    do l=1,nvarg
        do k=1,ndip
            do i=1,nlag
                do j=1,ndir
                    np(i,j,k,l)  = 0
                    dis(i,j,k,l) = 0.0
                    gam(i,j,k,l) = 0.0
                    hm(i,j,k,l)  = 0.0
                    tm(i,j,k,l)  = 0.0
                    hv(i,j,k,l)  = 0.0
                    tv(i,j,k,l)  = 0.0
                end do
            end do
        end do
    end do
    dismxs = ((dble(nlag) + 0.5 - EPSLON) * xlag) ** 2



! MAIN LOOP OVER ALL PAIRS:

    write (*,*) 'i', 'j', 'x(i)','y(i)','z(i)', 'x(j)','y(j)','z(j)', 'h', 'ilag', 'az', 'jdir', 'dp','kdip'

    do 3 i=1,nd
        do 4 j=i+1,nd

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
            if(hs > dismxs) goto 4
            if(hs < 0.0) hs = 0.0
            h   = sqrt(hs)

            ! now we determine the cell of the 3D variogram

            ! a) the ilag position = int (h/xlag + 1)

            if(h <= EPSLON) then
                goto 4 !no zero distance allowed
            else
                ilag = int(h/xlag + 1)
            endif

            ! b) the jdir direction (math angle)

            az = azimdg(dx,dy)+dazm/2.+orgdir     ! we add tolerance angle and rotation of direction origing to compare as firs lag from zero to alpha+beta
            az = modulo(az,360.)                  ! fix angle to interval [0,360[


            jdir = int(az/dazm + 1)                          ! get cell array location. int(az/dazm + 1) -> direction in 360
            jdir = modulo(real(jdir),real(ndir)+EPSLON)      ! this is to correct directions to 180 degrees.

            ! c) get kdip direction

            dxy = sqrt(max((dxs+dys),0.0))
            dp = modulo(azimdg(dxy,dz) + ddip/2. + orgdip, 180.)

            kdip = int(dp/ddip+1)                ! get cell array location properly corrected to first cell starting at [0,ddip[


            !write (*,*) i, j, x(i),y(i),z(i), x(j),y(j),z(j), h, ilag, az, jdir, dp,kdip




            do 6 iv=1,nvarg

                ! For this variogram, sort out which is the tail and the head value:

                it = ivtype(iv)
                ii = ivtail(iv)
                vrh   = vr(i,ii)
                ii = ivhead(iv)   ! TODO: Review this, using same variable for tail and head
                vrt   = vr(j,ii)
                if(it == 2) then
                    ii    = ivhead(iv)  ! TODO: Review this, using same variable for tail and head
                    vrtpr = vr(i,ii)
                    ii    = ivtail(iv)
                    vrhpr = vr(j,ii)
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

                    np(ilag,jdir,kdip,iv)  =  np(ilag,jdir,kdip,iv) + 1
                    dis(ilag,jdir,kdip,iv) = dis(ilag,jdir,kdip,iv) + dble(h)
                    tm(ilag,jdir,kdip,iv)  =  tm(ilag,jdir,kdip,iv) + dble(vrt)
                    hm(ilag,jdir,kdip,iv)  =  hm(ilag,jdir,kdip,iv) + dble(vrh)
                    gam(ilag,jdir,kdip,iv) = gam(ilag,jdir,kdip,iv) + dble((vrh-vrt)*(vrh-vrt))

                ! The Traditional Cross Semivariogram:

                else if(it == 2) then

                    np(ilag,jdir,kdip,iv)  = np(ilag,jdir,kdip,iv)  + 1
                    dis(ilag,jdir,kdip,iv) = dis(ilag,jdir,kdip,iv) + dble(h)
                    tm(ilag,jdir,kdip,iv)  = tm(ilag,jdir,kdip,iv)  + dble(0.5*(vrt+vrtpr))
                    hm(ilag,jdir,kdip,iv)  = hm(ilag,jdir,kdip,iv)  + dble(0.5*(vrh+vrhpr))
                    gam(ilag,jdir,kdip,iv) = gam(ilag,jdir,kdip,iv) + dble((vrhpr-vrh)*(vrt-vrtpr))


                ! The Covariance:

                else if(it == 3) then

                    np(ilag,jdir,kdip,iv)  = np(ilag,jdir,kdip,iv)  + 1
                    dis(ilag,jdir,kdip,iv) = dis(ilag,jdir,kdip,iv) + dble(h)
                    tm(ilag,jdir,kdip,iv)  = tm(ilag,jdir,kdip,iv)  + dble(vrt)
                    hm(ilag,jdir,kdip,iv)  = hm(ilag,jdir,kdip,iv)  + dble(vrh)
                    gam(ilag,jdir,kdip,iv) = gam(ilag,jdir,kdip,iv) + dble(vrh*vrt)



                ! The Correlogram:

                else if(it == 4) then

                    np(ilag,jdir,kdip,iv)  = np(ilag,jdir,kdip,iv)  + 1
                    dis(ilag,jdir,kdip,iv) = dis(ilag,jdir,kdip,iv) + dble(h)
                    tm(ilag,jdir,kdip,iv)  = tm(ilag,jdir,kdip,iv)  + dble(vrt)
                    hm(ilag,jdir,kdip,iv)  = hm(ilag,jdir,kdip,iv)  + dble(vrh)
                    hv(ilag,jdir,kdip,iv)  = hv(ilag,jdir,kdip,iv)  + dble(vrh*vrh)
                    tv(ilag,jdir,kdip,iv)  = tv(ilag,jdir,kdip,iv)  + dble(vrt*vrt)
                    gam(ilag,jdir,kdip,iv) = gam(ilag,jdir,kdip,iv) + dble(vrh*vrt)


                ! The Pairwise Relative:

                else if(it == 6) then

                    if(abs(vrt+vrh) > EPSLON) then
                        np(ilag,jdir,kdip,iv)  = np(ilag,jdir,kdip,iv)  + 1
                        dis(ilag,jdir,kdip,iv) = dis(ilag,jdir,kdip,iv) + dble(h)
                        tm(ilag,jdir,kdip,iv)  = tm(ilag,jdir,kdip,iv)  + dble(vrt)
                        hm(ilag,jdir,kdip,iv)  = hm(ilag,jdir,kdip,iv)  + dble(vrh)
                        gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                        gam(ilag,jdir,kdip,iv) = gam(ilag,jdir,kdip,iv) + dble(gamma*gamma)
                    endif


                ! Variogram of Logarithms:

                else if(it == 7) then
                    if(vrt > EPSLON .AND. vrh > EPSLON) then
                        np(ilag,jdir,kdip,iv)  = np(ilag,jdir,kdip,iv)  + 1
                        dis(ilag,jdir,kdip,iv) = dis(ilag,jdir,kdip,iv) + dble(h)
                        tm(ilag,jdir,kdip,iv)  = tm(ilag,jdir,kdip,iv)  + dble(vrt)
                        hm(ilag,jdir,kdip,iv)  = hm(ilag,jdir,kdip,iv)  + dble(vrh)
                        gamma   = dlog(vrt)-dlog(vrh)
                        gam(ilag,jdir,kdip,iv) = gam(ilag,jdir,kdip,iv) + dble(gamma*gamma)
                    endif


                ! Madogram:

                else if(it == 8) then

                        np(ilag,jdir,kdip,iv)  = np(ilag,jdir,kdip,iv)  + 1
                        dis(ilag,jdir,kdip,iv) = dis(ilag,jdir,kdip,iv) + dble(h)
                        tm(ilag,jdir,kdip,iv)  = tm(ilag,jdir,kdip,iv)  + dble(vrt)
                        hm(ilag,jdir,kdip,iv)  = hm(ilag,jdir,kdip,iv)  + dble(vrh)
                        gam(ilag,jdir,kdip,iv) = gam(ilag,jdir,kdip,iv) + dble(abs(vrh-vrt))

                endif

            ! Finish loops over variograms, directions, and the double data loops:

            6 END DO
            ! 5 END DO
        4 END DO
    3 END DO

   ! Get average values for gam, hm, tm, hv, and tv, then compute
   ! the correct "variogram" measure:


    do 7 id=1,ndir
        do 7 idp=1,ndip
            do 7 iv=1,nvarg
                do 7 il=1,nlag

                    ! write(*,*) il,id,idp,iv, np(il,id,idp,iv), dis(il,id,idp,iv), gam(il,id,idp,iv)

                    if(np(il,id,idp,iv) <= 0.) go to 7

                    rnum   = np(il,id,idp,iv)
                    dis(il,id,idp,iv) = dis(il,id,idp,iv) / dble(rnum)
                    gam(il,id,idp,iv) = gam(il,id,idp,iv) / dble(rnum)
                    hm(il,id,idp,iv)  = hm(il,id,idp,iv)  / dble(rnum)
                    tm(il,id,idp,iv)  = tm(il,id,idp,iv)  / dble(rnum)
                    hv(il,id,idp,iv)  = hv(il,id,idp,iv)  / dble(rnum)
                    tv(il,id,idp,iv)  = tv(il,id,idp,iv)  / dble(rnum)
                    it = ivtype(iv)

                    ! Attempt to standardize:

                    if(isill == 1) then
                        if(ivtail(iv) == ivhead(iv)) then
                            iii = ivtail(iv)
                            if((it == 1 .OR. it >= 9) .AND. sills(iii) > 0.0) &
                            gam(il,id,idp,iv) = gam(il,id,idp,iv) / sills(iii)
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
                        gam(il,id,idp,iv) = 0.5 * gam(il,id,idp,iv)
                    else if(it == 3) then
                        gam(il,id,idp,iv) = gam(il,id,idp,iv) - hm(il,id,idp,iv)*tm(il,id,idp,iv)
                    else if(it == 4) then
                        hv(il,id,idp,iv)  = hv(il,id,idp,iv)-hm(il,id,idp,iv)*hm(il,id,idp,iv)
                        if(hv(il,id,idp,iv) < 0.0) hv(il,id,idp,iv) = 0.0
                        hv(il,id,idp,iv)  = dsqrt(hv(il,id,idp,iv))
                        tv(il,id,idp,iv)  = tv(il,id,idp,iv)-tm(il,id,idp,iv)*tm(il,id,idp,iv)
                        if(tv(il,id,idp,iv) < 0.0) tv(il,id,idp,iv) = 0.0
                        tv(il,id,idp,iv)  = dsqrt(tv(il,id,idp,iv))
                        if((hv(il,id,idp,iv)*tv(il,id,idp,iv)) < EPSLON) then
                            gam(il,id,idp,iv) = 0.0
                        else
                            gam(il,id,idp,iv) =(gam(il,id,idp,iv)-hm(il,id,idp,iv)*tm(il,id,idp,iv))/ &
                                               (hv(il,id,idp,iv)*tv(il,id,idp,iv))
                        endif

                        ! Square "hv" and "tv" so that we return the variance:

                        hv(il,id,idp,iv)  = hv(il,id,idp,iv)*hv(il,id,idp,iv)
                        tv(il,id,idp,iv)  = tv(il,id,idp,iv)*tv(il,id,idp,iv)
                    else if(it == 5) then
                        htave  = 0.5*(hm(il,id,idp,iv)+tm(il,id,idp,iv))
                        htave  = htave   *   htave
                        if(htave < EPSLON) then
                            gam(il,id,idp,iv) = 0.0
                        else
                            gam(il,id,idp,iv) = gam(il,id,idp,iv)/dble(htave)
                        endif
                    else if(it >= 6) then
                        gam(il,id,idp,iv) = 0.5 * gam(il,id,idp,iv)
                    endif



    7 END DO


    return


 end subroutine gamv3D


 real*8 function azimdg(dx,dy)

    real*8, intent(in) :: dx, dy
    real*8 :: PI
    parameter(PI=3.14159265)


    azimdg=180./PI*atan2(dy,dx)

    ! here we fix the angle to be 0<=alpha<360
    azimdg=modulo(azimdg,360.)

    return

end function azimdg
