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
subroutine bicalib(ndp,vval,nd,u,v,wt, ncutu, cutu, ncutv, cutv, &
                  ssqu, avgu, umin, umax, ssqv, avgv, vmin, vmax, &
                  pdfrep, fract, yx, em, vm, nm, b,&
                  lcdf, error)
    !-----------------------------------------------------------------------

    !               Calibration for Markov/Bayes Simulation
    !               ***************************************

    ! This program calibrates a set of primary and secondary data for input
    ! to the ``mbsim'' program.  Collocated primary (u) and secondary (v)
    ! samples are input data and the output is the local prior cdfs for the
    ! primary variable (u) given that the secondary variable (v) belongs to
    ! specific classes.  Other calibration information is also written
    ! to the output file.


    ! The program is executed with no command line arguments.  The user
    ! will be prompted for the name of a parameter file.  The parameter
    ! file is described in the documentation (see the example bicalib.par)



    ! DIMENSIONING PARAMETERS:

    !      nDAT       maximum number of calibration data
    !      nUCUT       maximum number of cutoffs on u
    !      nVCUT       maximum number of cutoffs on v



    ! Original: H. Zhu                                   Date:     July 1990
    !-----------------------------------------------------------------------

    implicit none 

    ! input 
    integer, intent(in) :: ndp,nd,ncutu,ncutv
    real, intent(in), dimension(ndp) :: vval        !data
    real, intent(in), dimension(nd) :: u,v,wt       !auxiliar demse file
    real, intent(in), dimension(ncutu) :: cutu
    real, intent(in), dimension(ncutv) :: cutv


    ! output
    !   primary
    real, intent (out) :: ssqu, avgu, umin, umax !variance, average, minimum, maximum
    !   secundary
    real, intent (out) :: ssqv, avgv, vmin, vmax !variance, average, minimum, maximum
    !  error 
    integer, intent(out) :: error
    !  report output table
    real, intent(out), dimension(ncutu):: fract, b 
    real, intent(out), dimension(ncutv+1,ncutu+1) :: pdfrep, yx
    real, intent(out), dimension(ncutu,0:1) :: em, vm
    integer, intent(out), dimension(ncutu,0:1) :: nm
    real, intent(out), dimension(ndp,ncutu):: lcdf    

    ! internal

    real :: pdf(ncutv+1,ncutu+1), soft(ncutu,0:1,nd), softw(ncutu,0:1,nd)
    integer :: i, j, icls, k, n
    real :: cum, dev, sum, xd
    real, dimension(0:ncutu+1) :: ucut
    real, dimension(0:ncutv+1) :: vcut

    ! Initialize for some statistics:

    avgu = 0.0
    avgv = 0.0
    ssqu = 0.0
    ssqv = 0.0
    umin = 1.0e10
    vmin = 1.0e10
    umax =-1.0e10
    vmax =-1.0e10
    sum  = 0.0

    do i=1, ncutu
        ucut (i) = cutu(i)
    end do
    do i=1, ncutv
        vcut (i) = cutv(i)
    end do


    ! Read in as much data as the allocated storage will allow:

    do i=1, nd 
        sum  = sum  + wt(i)
        avgu = avgu + u(i)*wt(i)
        avgv = avgv + v(i)*wt(i)
        ssqv = ssqv + v(i)*v(i)*wt(i)
        ssqu = ssqu + u(i)*u(i)*wt(i)
        if(u(i) < umin) umin = u(i)
        if(v(i) < vmin) vmin = v(i)
        if(u(i) > umax) umax = u(i)
        if(v(i) > vmax) vmax = v(i)
    end do

    ! There has to be at least two data to go ahead with calibration:

    if(nd <= 1 .OR. sum <= 0.001) then
        error = 10 ! 'TOO FEW DATA to go ahead  or  sum of weights too small
        stop
    endif

    ! write(*,*) 'Calculating the calibration parameters'

    xd = 1.0 / real(sum)

    ! Compute the averages and variances as an error check for the user:

    avgu = avgu * xd
    avgv = avgv * xd
    ssqu = ssqu * xd - avgu * avgu
    ssqv = ssqv * xd - avgv * avgv


    ! Establish lower and upper bounds:

    ucut(0)       = umin - 1. - abs(umax)
    ucut(ncutu+1) = umax + abs(umax+1.)
    vcut(0)       = vmin - 1. - abs(vmax)
    vcut(ncutv+1) = vmax + abs(vmax+1.)
    do i=1,ncutu
        fract(i)  = 0.0
    end do

    ! Calculate the indicator mean for each cutoff:

    do i=1,nd
        do j=1,ncutu
            if((u(i) > ucut(j-1)) .AND. (u(i) <= ucut(j)))then
                fract(j)=fract(j)+wt(i)
                go to 2
            endif
        end do
        2 continue
    end do

   ! Turn fract() into a cdf:

    i  = 1
    fract(i) = fract(i) * xd
    do i=2,ncutu
        fract(i) = fract(i-1) + fract(i) * xd
    end do
    

   ! Compute the pdf table:

    do k=1,ncutv+1
        do j=1,ncutu+1
            pdf(k,j) = 0.0
            pdfrep(k,j) = 0.0
        end do
    end do
    do i=1,nd
        do k=1,ncutv+1
            if((v(i) > vcut(k-1)) .AND. (v(i) <= vcut(k)))then
                do j=1,ncutu+1
                    if((u(i) > ucut(j-1)) .AND. (u(i) <= ucut(j)))then
                        pdf(k,j) = pdf(k,j) + wt(i)
                        go to 6
                    endif
                end do
            endif
        end do
        6 continue
    end do

    ! Turn the bivariate pdf into conditional cdf, i.e. the local prior cdf

    do i=1,ncutv+1
        do j=1,ncutu+1
            pdf(i,j) = pdf(i,j) * xd
            pdfrep(i,j) =pdf(i,j)
        end do
    end do


    ! Loop over all the v cutoffs:

    do i=1,ncutv+1
        cum = 0.0
        do j=1,ncutu+1
            cum = cum + pdf(i,j)
        end do
        if(cum > 0.0) then
            do j=1,ncutu+1
                pdf(i,j) = pdf(i,j) / cum
            end do
        endif
        do j=1,ncutu+1
            yx(i,j) = 0.0
            do k=1,j
                yx(i,j) = yx(i,j) + pdf(i,k)
            end do
        end do
    end do

    ! Calculate the calibration parameters from the local prior cdf table
    ! and the input data:  First, initialize counters:

    do j=1,ncutu
        do k=0,1
            nm(j,k)   = 0
            em(j,k)   = 0.0
            vm(j,k)   = 0.0
            do n=1,nd
                soft(j,k,n)  = 0.0
                softw(j,k,n) = 0.0
            end do
        end do
    end do

    ! MAIN Loop over all the u cutoffs to sort out the classification of
    ! all the data:

    do j=1,ncutu
    
    ! Now, loop over all the data:
    
        do i=1,nd
            if(u(i) <= ucut(j))then
                do k=1,ncutv+1
                    if((v(i) > vcut(k-1)) .AND. &
                    (v(i) <= vcut(k)))  then
                        nm(j,1)            = nm(j,1) + 1
                        soft(j,1,nm(j,1))  = yx(k,j)
                        softw(j,1,nm(j,1)) = wt(i)
                        go to 11
                    endif
                end do
            else
                do k=1,ncutv+1
                    if((v(i) > vcut(k-1)) .AND. &
                    (v(i) <= vcut(k)))  then
                        nm(j,0)            = nm(j,0) + 1
                        soft(j,0,nm(j,0))  = yx(k,j)
                        softw(j,0,nm(j,0)) = wt(i)
                        go to 11
                    endif
                end do
            endif
            11 continue
        end do
    end do

    ! Now, get some statistics on the calibration.  Loop over all the u
    ! cutoffs and then the two classication categories:

    do j=1,ncutu
        do k=0,1
            if(nm(j,k) >= 1) then
                em(j,k) = 0.0
                sum     = 0.0
                do i=1,nm(j,k)
                    em(j,k) = em(j,k) + soft(j,k,i)* &
                    softw(j,k,i)
                    sum     = sum     + softw(j,k,i)
                end do
                em(j,k) = em(j,k) / max(real(sum),0.000001)
                vm(j,k) = 0.0
                sum     = 0.0
                do i=1,nm(j,k)
                    dev     = soft(j,k,i) - em(j,k)
                    vm(j,k) = vm(j,k) + dev*dev* &
                    softw(j,k,i)
                    sum     = sum     + softw(j,k,i)
                end do
                vm(j,k) = vm(j,k) / max(real(sum),0.000001)
            endif
        end do
    end do

    ! Hardness coefficients:

    do j=1,ncutu
        b(j)   = em(j,1) - em(j,0)
        if(nm(j,1) < 1 .OR. nm(j,0) < 1) b(j) = -9.0
    end do

    ! Loop through all observations in the input file:

    do j=1, ndp
        do i=1,ncutv+1
            if((vval(j) > vcut(i-1)) .AND. (vval(j) <= vcut(i))) then
                icls = i
                go to 42
            endif
        end do
        42 continue
        do i=1,ncutu
            lcdf(j,i) = yx(icls,i)
        end do
    end do

end subroutine bicalib

