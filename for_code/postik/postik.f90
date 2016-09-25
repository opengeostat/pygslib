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


Module Commons


    real :: UNEST    ! asign this externally at init , otherwise will be zero by default


End Module Commons  


subroutine set_unest(unest_)

    use Commons

    real, intent(in) :: unest_

    UNEST = unest_

end subroutine set_unest  


subroutine get_unest(unest_)

    use Commons
    
    real, intent(out) :: unest_

    unest_ = UNEST 
    
    if (isnan(UNEST)) write(*,*) 'UNEST is a NaN'
    
    if (UNEST==anan) write(*,*) 'UNEST == NaN'
    
    return

end subroutine get_unest  


subroutine postik( &
                  iout,outpar, &  ! output option, output parameter
                  nccut, ccut1, & ! number of thresholds, the thresholds
                  ivol,ivtyp,varred, & ! volume support?, type, varred
                  zmin,zmax, & ! minimum and maximum Z value
                  ltail,ltpar, & ! lower tail: option, parameter
                  middle,mpar, & ! middle    : option, parameter
                  utail,utpar, & ! upper tail: option, parameter
                  maxdis, & ! maximum discretization
                  vr, wt, nc, & ! global distribution
                  p, na, &  ! input data 
                  out1,out2,out3, &  ! output data 
                  error) ! to control the error
                  
    !-----------------------------------------------------------------------

    !                   Post Process IK Distributions
    !                   *****************************

    ! Reads IK-generated distributions and post processes them
    ! according to user specifications.

    !      - affine or indirect lognormal volume support correction
    !      - E-type mean
    !      - probability of exceeding specified threshold and the mean
    !           value above threshold.
    !      - compute a z-threshold value for a specified CDF value.

    ! See the text for detailed discussion of parameters to extrapolate the
    ! discretely coded cdf values.


    ! INPUT/OUTPUT Parameters:

    !   distin           the input distributions (output from IK3D)
    !   outfl            the output file for E-type,....
    !   iout,outpar         =1 E-type,
    !                       =2 prob and grade > and <= outpar;
    !                       =3 Z percentile corresponding to outpar,
    !                       =4 conditional variance.
    !   nccut            the number of cutoffs
    !   ccut(i)          cutoffs
    !   ivol             =1, the consider volume support corrections
    !   ivtyp            =1, affine; =2 indirect lognormal
    !   varred           variance adjustment factor, between 0 and 1
    !   datafl           the global data distribution
    !   ivr,iwt,tmin     parameters to read in global data distribution
    !   varred           variance reduction parameter
    !   zmin,zmax        minimum and maximum Z values
    !   ltail,ltpar      option and parameter to handle values in lower tail
    !   middle,mpar      option and parameter to handle values in the middle
    !   utail,utpar      option and parameter to handle values in upper tail
    !   maxdis           discretization to compute E-type and mean > cutoff



    !-----------------------------------------------------------------------


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    iout == 1 out1 -> 'mean', out2 -> man, out3 -> man
    !    iout == 2 out1 -> 'prob > cutoff', out2 -> 'mean > cutoff', out3 -> 'mean < cutoff'
    !    iout == 3 out1 -> 'Z value corresponding to CDF = value', out2 -> man, out3 -> man
    !    iout == 4 out1 -> 'Conditional Variance', out2 -> man, out3 -> man
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    

    use commons


    ! external variables (in)
    real, intent(in) :: outpar, varred, zmin,zmax,ltpar,mpar,utpar
    integer, intent(in) :: iout, nccut, ivol,ivtyp, ltail,middle,utail,maxdis,nc, na
    real, intent(in), dimension(nccut) :: ccut1
    real, intent(in), dimension(nc) :: vr, wt
    real, intent(in), dimension(na,nccut) :: p

    ! external variables (out)
    real, intent(out), dimension(na) :: out1,out2,out3
    integer, intent(out) :: error
    
    ! local variables     
    real :: EPSLON, outpar_, varred_
    real, dimension(nc) :: cdf,cut
    real, dimension(nccut) :: ccdf,ccdf1,ccdf2,ccut
    logical ::   testfl

    EPSLON=1.0e-6
    error = 0
    outpar_ = outpar
    varred_ = varred

    write(*,*) ' output option & par = ',iout,outpar_
    if(iout /= 1 .AND. iout /= 2 .AND. iout /= 3 .AND. iout /= 4) then
        write(*,*) ' ERROR: invalid output option ',iout
        error = 1
        return
    end if
    if(iout == 3) then
        if(outpar_ < 0.0) then
            write (*,*) 'Invalid p-value for iout=3'
            error = 2
            return ! stop 'Invalid p-value for iout=3'
        end if
        if(outpar_ > 1.0) then
            write (*,*) 'Invalid p-value for iout=3'
            error = 3
            return ! 'Invalid p-value for iout=3'
        end if
        outpar_ = min(max(outpar_,EPSLON),(1.0-EPSLON))
    end if

    write(*,*) ' number of cutoffs = ',nccut
    write(*,*) ' cutoffs = ',(ccut1(i),i=1,nccut)
    write(*,*) ' volume variance = ',ivol,ivtyp,varred_
    if(varred_ < 0.0 .OR. varred_ > 1.0) then
        write(*,*) ' ERROR: invalid variance reduction ',varred_
        write(*,*) '        must be between 0 and 1'
        error = 3
        return
    end if

    write(*,*) ' minimum and maximum = ',zmin,zmax
    if(iout == 2 .AND. outpar_ < zmin) then
        write(*,*) 'Invalid z-value for iout=2'
        error = 4
    end if
    if(iout == 2 .AND. outpar_ > zmax) then
        write(*,*) 'Invalid z-value for iout=2'
        error = 5
    end if

    write(*,*) ' ltail, ltpar = ',ltail,ltpar
    write(*,*) ' middle, mpar = ',middle,mpar
    write(*,*) ' utail, utpar = ',utail, utpar
    write(*,*) ' discretization = ',maxdis


    ! Do we have a global distribution and do we need it anyway?

    ncut = 0
    if(ltail /= 3 .AND. middle /= 3 .AND. utail /= 3) testfl = .FALSE. 

    ! Read in the global cdf ("cut" and "cdf" arrays):

    if(testfl) then
        do l=1, nc         
            tcdf = 0.0
            gmean = 0.0
            ncut  = 0
            if(wt(l) <= 0.0 .or. wt(l) == UNEST .or. isnan(wt(l))) then
                write(*,*) ' ERROR: negative weights at index ', l
                error = 6
                return 
            endif
            
            if(vr(l) == UNEST .or. isnan(vr(l)) ) then
                write(*,*) ' ERROR: nans in globa CDF at index ', l
                error = 6
                return 
            endif
            
            ncut = ncut + 1
            cut(ncut) = vr(l)
            gmean     = gmean + vr(l)
            cdf(ncut) = wt(l)
            tcdf = tcdf + cdf(ncut)
                        
            if(tcdf <= 0) then
                write (*,*) ' total global CDF <= 0'
                return
            end if
            
            gmean = gmean / tcdf
        
            ! Turn the (possibly weighted) distribution into a cdf that is useful:
        
            call sortem(1,ncut,cut,1,cdf,c,d,e,f,g,h)
            oldcp = 0.0
            cp    = 0.0
            tcdf  = 1.0 / tcdf
            do i=1,ncut
                cp     = cp + cdf(i) * tcdf
                cdf(i) =(cp + oldcp) * 0.5
                oldcp  = cp
            end do
        
            ! Get median and write some info to the screen:
        
            call locate(cdf,ncut,1,ncut,0.5,j)
            gmedian = powint(cdf(j),cdf(j+1),cut(j),cut(j+1),0.5,1.)
            write(*,*) 'Global cdf from file: ',datafl
            write(*,*) '   number of data: ',ncut
            write(*,*) '   global mean:    ',gmean
            write(*,*) '   global median:  ',gmedian
            write(*,*)
        end do
    endif


    ! If we are doing the affine correction then compute the square root
    ! of the variance ``adjustment'' factor:

    if(ivol == 1 .AND. ivtyp == 1) varred_ = sqrt(max(varred_,0.0))

    ! BIG LOOP reading in each distribution in turn:

    nproc = 0
    procm = 0.0

    do l=1, na
        do i=1,nccut
            ccdf(i) = p(l,i)
        end do 

        ! initialize output
        out1(l) = UNEST
        out2(l) = UNEST
        out3(l) = UNEST

        ! Check for missing values:

        if(ccdf(nccut) < -0.1) cycle  ! skip this iteration and do next iteration 
        if(ccdf(nccut) == UNEST) cycle  ! skip this iteration and do next iteration 
        if(isnan(ccdf(nccut))) cycle

        ! Reinstate "ccut1" because volume support may have changed "ccut":

        do i=1,nccut
            ccut(i) = ccut1(i)
        end do

        ! Correct Order Relations (the distributions coming from IK3D have
        ! already been corrected):

        do i=1,nccut
            ccdf(i)=max(min(1.0,ccdf(i)),0.0)
        end do
        ccdf1(1) = ccdf(1)
        do i=2,nccut
            ccdf1(i) = max(ccdf1(i-1),ccdf(i))
        end do
        ccdf2(nccut) = ccdf(nccut)
        do i=nccut-1,1,-1
            ccdf2(i) = min(ccdf2(i+1),ccdf(i))
        end do
        do i=1,nccut
            ccdf(i) = 0.5*(ccdf1(i)+ccdf2(i))
        end do

        ! Volume support correction with parameters (ivol,ivtyp,varred_):

        if(ivol == 1) then
            dis    = 1.0 / real(maxdis)
            cdfval = -0.5*dis
            etype  = 0.0
            ecv    = 0.0
            do i=1,maxdis
                cdfval  = cdfval + dis
                zval    = -1.0
                call beyond(1,nccut,ccut,ccdf,ncut,cut,cdf,zmin, &
                zmax,ltail,ltpar,middle,mpar,utail, &
                utpar,zval,cdfval,ierr)
                if(ierr /= 0) write(*,*) 'ERROR: ',ierr,' continuing1'
                etype = etype + zval
                ecv   = ecv   + zval*zval
            end do
            etype = etype / real(maxdis)
            ecv   = sqrt(max((ecv/real(maxdis)-etype*etype),0.0)) &
            / max(etype,EPSLON)
            if(etype == 0.0) then
                write(*,*) 'NO support correction with 0 mean'
                etype = EPSLON
            endif
        
            ! Affine Correction:
        
            if(ivtyp == 1) then
                do i=1,nccut
                    ccut(i) = etype+varred_*(ccut(i)-etype)
                end do
            else
            
                ! Indirect lognormal:  1. Compute parameters "a" and "b"
                !                      2. Correct quantiles
                !                      3. recompute etype
                !                      4. Recorrect quantiles to correct mean
                
                
                !            Parameters "a" and "b":
            
                b = sqrt(max((alog(varred_*ecv*ecv+1) &
                /alog(ecv*ecv+1)),0.0))
                a = (etype/sqrt(max((varred_*ecv*ecv+1),0.0))) * &
                (sqrt(max((ecv*ecv+1),0.0))/etype)**b
            
                !            Correct quantiles:
            
                do i=1,nccut
                    ccut(i) = a*ccut(i)**b
                end do
            
                !            New etype:
            
                cdfval = -0.5*dis
                enew  = 0.0
                do i=1,maxdis
                    cdfval  = cdfval + dis
                    zval    = -1.0
                    call beyond(1,nccut,ccut,ccdf,ncut,cut,cdf,zmin, &
                    zmax,ltail,ltpar,middle,mpar,utail, &
                    utpar,zval,cdfval,ierr)
                    if(ierr /= 0) then
                        write(*,*) 'ERROR: ',ierr,' continuing2'
                    endif
                    enew = enew + zval
                end do
                enew = enew / real(maxdis)
            
                !            Recorrect "quantiles":
            
                do i=1,nccut
                    ccut(i) = (etype/enew)*ccut(i)
                end do
            endif
        endif

        ! Compute mean of local distribution if E-type (iout=1), or need the
        ! mean above (below) a cutoff, or if performing volume support:

        if(iout <= 2 .OR. iout == 4) then
            dis    = 1.0 / real(maxdis)
            cdfval = -0.5*dis
            eabove = 0.0
            nabove = 0
            ebelow = 0.0
            nbelow = 0
            etype  = 0.0
            ecv    = 0.0
            
            do i=1,maxdis
                cdfval  = cdfval + dis
                zval    = -1.0
                call beyond(1,nccut,ccut,ccdf,ncut,cut,cdf,zmin, &
                zmax,ltail,ltpar,middle,mpar,utail, &
                utpar,zval,cdfval,ierr)
                if(ierr /= 0) then
                    write(*,*) 'ERROR: ',ierr,' continuing3'
                endif
                etype = etype + zval
                ecv   = ecv   + zval*zval
                if(zval <= outpar_) then
                    nbelow = nbelow + 1
                    ebelow = ebelow + zval
                else
                    nabove = nabove + 1
                    eabove = eabove + zval
                end if
            end do
        
            ! e-type and conditional variance:
        
            etype = etype / real(maxdis)
            ecv   = ecv/real(maxdis)-etype*etype
        endif

        ! Write out E-type?

        if(iout == 1) out1(l) = etype

        ! Do we need probability and mean value above a threshold?

        if(iout == 2) then
        
        !      Get probability:
        
            cdfval  = -1.0
            call beyond(1,nccut,ccut,ccdf,ncut,cut,cdf,zmin, &
            zmax,ltail,ltpar,middle,mpar,utail, &
            utpar,outpar_,cdfval,ierr)
            if(ierr /= 0) then
                write(*,*) 'ERROR: ',ierr,' continuing4'
            endif
            prob   = 1.0 - cdfval
            eabove = eabove / real(max(1,nabove))
            ebelow = ebelow / real(max(1,nbelow))
            out1(l) = prob
            out1(l) = eabove
            out1(l) = ebelow
        endif

        ! Do we need the "Z" value corresponding to a particular CDF?

        if(iout == 3) then
            zval    = -1.0
            call beyond(1,nccut,ccut,ccdf,ncut,cut,cdf,zmin, &
            zmax,ltail,ltpar,middle,mpar,utail, &
            utpar,zval,outpar_,ierr)
            if(ierr /= 0) then
                write(*,*) 'ERROR: ',ierr,' continuing5'
            endif
            out1(l) = zval
        endif

        ! Write out conditional variance:

        if(iout == 4) out1(l) = ecv

        ! Return for another:

        nproc = nproc + 1
        procm = procm + etype

    end do

    ! Finished:

    procm = procm / max(real(nproc),1.0)
    write(*,*)
    write(*,*) 'Number of distributions ',nproc
    if(procm /= 0.0) write(*,*) 'Overall mean:           ',procm
    if(iout == 1) write(*,*) 'Local Means (E-type):   ',outfl
    if(iout == 2) write(*,*) 'Prob and mean > cutoff: ',outfl
    if(iout == 3) write(*,*) 'Z values for outpar in: ',outfl
    write(*,*)

    ! Finished:

    write (*,*) ' POSTIK Finished'

    return 

end subroutine postik
