    subroutine beyond(ivtype,nccut,ccut,ccdf,ncut,cut,cdf,zmin,zmax, &
    ltail,ltpar,middle,mpar,utail,utpar,zval,cdfval,ierr)
!-----------------------------------------------------------------------

!                     Go Beyond a Discrete CDF
!                     ************************

! This subroutine is a general purpose subroutine to interpolate within
! and extrapolate beyond discrete points on a conditional CDF.  If the
! Z value "zval" is specified then the corresponding CDF value "cdfval"
! will be computed, if the CDF value "cdfval" is specified the
! corresponding Z value "zval" will be computed.



! INPUT/OUTPUT VARIABLES:

!   ivtype           variable type (1=continuous, 0=categorical)
!   nccut            number of cutoffs defining the conditional CDF
!   ccut()           real array of the nccut cutoffs
!   ccdf()           real array of the conditional cdf values
!   ncut             number of cutoffs defining the global CDF
!   cut()            real array of the ncut cutoffs
!   cdf()            real array of the global cdf values

!   zmin,zmax        minimum and maximum allowable data values
!   ltail            option to handle values in lower tail
!   ltpar            parameter required for option ltail
!   middle           option to handle values in the middle
!   mpar             parameter required for option middle
!   utail            option to handle values in upper tail
!   utpar            parameter required for option utail

!   zval             interesting cutoff (if -1 then it is calculated)
!   cdfval           interesting CDF (if -1 then it is calculated)


!-----------------------------------------------------------------------
    parameter(EPSLON=1.0e-20,UNEST=-1.0)
    dimension ccut(nccut),ccdf(nccut),cut(1),cdf(1)
    real ::      utpar,mpar,ltpar,lambda
    integer ::   ltail,utail,middle,cclow,cchigh

! Check for both "zval" and "cdfval" defined or undefined:

    ierr  = 1
    if(zval > UNEST .AND. cdfval > UNEST) return
    if(zval <= UNEST .AND. cdfval <= UNEST) return

! Handle the case of a categorical variable:

    if(ivtype == 0) then
        cum = 0
        do i=1,nccut
            cum = cum + ccdf(i)
            if(cdfval <= cum) then
                zval = ccut(i)
                return
            endif
        end do
        return
    end if

! Figure out what part of distribution: ipart = 0 - lower tail
!                                       ipart = 1 - middle
!                                       ipart = 2 - upper tail
    ierr  = 0
    ipart = 1
    if(zval > UNEST) then
        if(zval <= ccut(1))       ipart = 0
        if(zval >= ccut(nccut))   ipart = 2
    else
        if(cdfval <= ccdf(1))     ipart = 0
        if(cdfval >= ccdf(nccut)) ipart = 2
    endif

! ARE WE IN THE LOWER TAIL?

    if(ipart == 0) then
        if(ltail == 1) then
        
        ! Straight Linear Interpolation:
        
            powr = 1.0
            if(zval > UNEST) then
                cdfval = powint(zmin,ccut(1),0.0,ccdf(1), &
                zval,powr)
            else
                zval = powint(0.0,ccdf(1),zmin,ccut(1), &
                cdfval,powr)
            endif
        else if(ltail == 2) then
        
        ! Power Model interpolation to lower limit "zmin"?
        
            if(zval > UNEST) then
                cdfval = powint(zmin,ccut(1),0.0,ccdf(1), &
                zval,ltpar)
            else
                powr = 1.0 / ltpar
                zval = powint(0.0,ccdf(1),zmin,ccut(1), &
                cdfval,powr)
            endif
        
        ! Linear interpolation between the rescaled global cdf?
        
        else if(ltail == 3) then
            if(zval > UNEST) then
            
            ! Computing the cdf value. Locate the point and the class bound:
            
                call locate(cut,ncut,1,ncut,zval,idat)
                call locate(cut,ncut,1,ncut,ccut(1),iupp)
            
            ! Straight linear interpolation if no data; otherwise, linear:
            
                if(idat <= 0 .OR. idat >= ncut .OR. &
                iupp <= 0 .OR. iupp >= ncut) then
                    cdfval = powint(zmin,cut(1),0.0,cdf(1), &
                    zval,1.)
                else
                    temp   = powint(cut(idat),cut(idat+1), &
                    cdf(idat),cdf(idat+1),zval,1.)
                    cdfval = temp*ccdf(1)/cdf(iupp)
                endif
            else
            
            ! Computing Z value: Are there any data out in the tail?
            
                call locate(cut,ncut,1,ncut,ccut(1),iupp)
            
            ! Straight linear interpolation if no data; otherwise, local linear
            ! interpolation:
            
                if(iupp <= 0 .OR. iupp >= ncut) then
                    zval = powint(0.0,cdf(1),zmin,cut(1), &
                    cdfval,1.)
                else
                    temp = cdfval*cdf(iupp)/ccdf(1)
                    call locate(cdf,ncut,1,ncut,temp,idat)
                    if(idat <= 0 .OR. idat >= ncut) then
                        zval = powint(0.0,cdf(1),zmin, &
                        cut(1),cdfval,1.)
                    else
                        zval = powint(cdf(idat),cdf(idat+1), &
                        cut(idat),cut(idat+1),temp,1.)
                    end if
                endif
            endif
        else
        
        ! Error situation - unacceptable option:
        
            ierr = 2
            return
        endif
    endif

! FINISHED THE LOWER TAIL,  ARE WE IN THE MIDDLE?

    if(ipart == 1) then
    
    ! Establish the lower and upper limits:
    
        if(zval > UNEST) then
            call locate(ccut,nccut,1,nccut,zval,cclow)
        else
            call locate(ccdf,nccut,1,nccut,cdfval,cclow)
        endif
        cchigh = cclow + 1
        if(middle == 1) then
        
        ! Straight Linear Interpolation:
        
            powr = 1.0
            if(zval > UNEST) then
                cdfval = powint(ccut(cclow),ccut(cchigh), &
                ccdf(cclow),ccdf(cchigh),zval,powr)
            else
                zval = powint(ccdf(cclow),ccdf(cchigh), &
                ccut(cclow),ccut(cchigh),cdfval,powr)
            endif
        
        ! Power interpolation between class bounds?
        
        else if(middle == 2) then
            if(zval > UNEST) then
                cdfval = powint(ccut(cclow),ccut(cchigh), &
                ccdf(cclow),ccdf(cchigh),zval,mpar)
            else
                powr = 1.0 / mpar
                zval = powint(ccdf(cclow),ccdf(cchigh), &
                ccut(cclow),ccut(cchigh),cdfval,powr)
            endif
        
        ! Linear interpolation between the rescaled global cdf?
        
        else if(middle == 3) then
            call locate(cut,ncut,1,ncut,ccut(cclow),ilow)
            call locate(cut,ncut,1,ncut,ccut(cchigh),iupp)
            if(cut(ilow) < ccut(cclow))  ilow = ilow + 1
            if(cut(iupp) > ccut(cchigh)) iupp = iupp - 1
            if(zval > UNEST) then
                call locate(cut,ncut,1,ncut,zval,idat)
            
            ! Straight linear interpolation if no data; otherwise, local linear
            ! interpolation:
            
                if(idat <= 0 .OR. idat >= ncut .OR. &
                ilow <= 0 .OR. ilow >= ncut .OR. &
                iupp <= 0 .OR. iupp >= ncut .OR. &
                iupp <= ilow) then
                    cdfval=powint(ccut(cclow),ccut(cchigh), &
                    ccdf(cclow),ccdf(cchigh),zval,1.)
                else
                    temp = powint(cut(idat),cut(idat+1), &
                    cdf(idat),cdf(idat+1),zval,1.)
                    cdfval=powint(cdf(ilow),cdf(iupp), &
                    ccdf(cclow),ccdf(cchigh),temp,1.)
                endif
            else
            
            ! Straight linear interpolation if no data; otherwise, local linear
            ! interpolation:
            
                if(ilow <= 0 .OR. ilow >= ncut .OR. &
                iupp <= 0 .OR. iupp >= ncut .OR. &
                iupp <= ilow) then
                    zval=powint(ccdf(cclow),ccdf(cchigh), &
                    ccut(cclow),ccut(cchigh),cdfval,1.)
                else
                    temp=powint(ccdf(cclow),ccdf(cchigh), &
                    cdf(ilow),cdf(iupp),cdfval,1.)
                    call locate(cdf,ncut,1,ncut,temp,idat)
                    if(cut(idat) < ccut(cclow)) idat=idat+1
                    if(idat <= 0 .OR. idat >= ncut .OR. &
                    cut(idat+1) > ccut(cchigh)) then
                        zval = powint(ccdf(cclow), &
                        ccdf(cchigh),ccut(cclow), &
                        ccut(cchigh),cdfval,1.)
                    else
                        zval = powint(cdf(idat),cdf(idat+1), &
                        cut(idat),cut(idat+1),temp,1.)
                    end if
                    zval = powint(cdf(idat),cdf(idat+1), &
                    cut(idat),cut(idat+1),temp,1.)
                endif
            endif
        else
        
        ! Error situation - unacceptable option:
        
            ierr = 2
            return
        endif
    endif

! FINISHED THE MIDDLE,  ARE WE IN THE UPPER TAIL?

    if(ipart == 2) then
        if(utail == 1) then
            powr = 1.0
            if(zval > UNEST) then
                cdfval = powint(ccut(nccut),zmax,ccdf(nccut), &
                &                                   1.0,zval,powr)
            else
                zval   = powint(ccdf(nccut),1.0,ccut(nccut), &
                zmax,cdfval,powr)
            endif

        else if(utail == 2) then
        
        ! Power interpolation to upper limit "utpar"?
        
            if(zval > UNEST) then
                cdfval = powint(ccut(nccut),zmax,ccdf(nccut), &
                &                                   1.0,zval,utpar)
            else
                powr = 1.0 / utpar
                zval   = powint(ccdf(nccut),1.0,ccut(nccut), &
                zmax,cdfval,powr)
            endif
        
        ! Linear interpolation between the rescaled global cdf?
        
        else if(utail == 3) then
            if(zval > UNEST) then
            
            ! Approximately Locate the point and the class bound:
            
                call locate(cut,ncut,1,ncut,zval,idat)
                call locate(cut,ncut,1,ncut,ccut(nccut),ilow)
                if(cut(idat) < zval)        idat = idat + 1
                if(cut(ilow) < ccut(nccut)) ilow = ilow + 1
            
            ! Straight linear interpolation if no data; otherwise, local linear
            ! interpolation:
            
                if(idat <= 0 .OR. idat >= ncut .OR. &
                ilow <= 0 .OR. ilow >= ncut) then
                    cdfval = powint(ccut(nccut),zmax, &
                    ccdf(nccut),1.0,zval,1.)
                else
                    temp   = powint(cut(idat),cut(idat+1), &
                    cdf(idat),cdf(idat+1),zval,1.)
                    cdfval = powint(cdf(ilow),1.0, &
                    ccdf(nccut),1.0,temp,1.)
                endif
            else
            
            ! Computing Z value: Are there any data out in the tail?
            
                call locate(cut,ncut,1,ncut,ccut(nccut),ilow)
                if(cut(ilow) < ccut(nccut)) ilow = ilow + 1
            
            ! Straight linear interpolation if no data; otherwise, local linear
            ! interpolation:
            
                if(ilow <= 0 .OR. ilow >= ncut) then
                    zval   = powint(ccdf(nccut),1.0, &
                    ccut(nccut),zmax,cdfval,1.)
                else
                    temp = powint(ccdf(nccut),1.0, &
                    cdf(ilow),1.0,cdfval,1.)
                    call locate(cdf,ncut,1,ncut,temp,idat)
                    if(cut(idat) < ccut(nccut)) idat=idat+1
                    if(idat >= ncut) then
                        zval   = powint(ccdf(nccut),1.0, &
                        ccut(nccut),zmax,cdfval,1.)
                    else
                        zval = powint(cdf(idat),cdf(idat+1), &
                        cut(idat),cut(idat+1),temp,1.)
                    endif
                endif
            endif
        
        ! Fit a Hyperbolic Distribution?
        
        else if(utail == 4) then
        
        ! Figure out "lambda" and required info:
        
            lambda = (ccut(nccut)**utpar)*(1.0-ccdf(nccut))
            if(zval > UNEST) then
                cdfval = 1.0 - (lambda/(zval**utpar))
            else
                zval = (lambda/(1.0-cdfval))**(1.0/utpar)
            endif
        else
        
        ! Error situation - unacceptable option:
        
            ierr = 2
            return
        endif
    endif
    if(zval < zmin) zval = zmin
    if(zval > zmax) zval = zmax

! All finished - return:

    return
    end subroutine beyond
