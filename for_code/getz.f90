    real function getz(pval,nt,vr,cdf,zmin,zmax,ltail,ltpar, &
    utail,utpar)
!-----------------------------------------------------------------------

!           Back Transform Univariate Data from Normal Scores
!           *************************************************

! This subroutine backtransforms a standard normal deviate from a
! specified back transform table and option for the tails of the
! distribution.  Call once with "first" set to true then set to false
! unless one of the options for the tail changes.



! INPUT VARIABLES:

!   pval             probability value to use
!   nt               number of values in the back transform tbale
!   vr(nt)           original data values that were transformed
!   cdf(nt)          the corresponding transformed values
!   zmin,zmax        limits possibly used for linear or power model
!   ltail            option to handle values less than cdf(1)
!   ltpar            parameter required for option ltail
!   utail            option to handle values greater than cdf(nt)
!   utpar            parameter required for option utail



!-----------------------------------------------------------------------
    parameter(EPSLON=1.0e-20)
    dimension vr(nt),cdf(nt)
    real ::      ltpar,utpar,lambda
    integer ::   ltail,utail

! Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):

    if(pval <= cdf(1)) then
        getz = vr(1)
        if(ltail == 1) then
            getz = powint(0.0,cdf(1),zmin,vr(1),pval,1.0)
        else if(ltail == 2) then
            cpow = 1.0 / ltpar
            getz = powint(0.0,cdf(1),zmin,vr(1),pval,cpow)
        endif
    
    ! Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
    
    else if(pval >= cdf(nt)) then
        cdfhi  = cdf(nt)
        getz   = vr(nt)
        if(utail == 1) then
            getz   = powint(cdfhi,1.0,vr(nt),zmax,pval,1.0)
        else if(utail == 2) then
            cpow   = 1.0 / utpar
            getz   = powint(cdfhi,1.0,vr(nt),zmax,pval,cpow)
        else if(utail == 4) then
            lambda = (vr(nt)**utpar)*(1.0-cdf(nt))
            getz   = (lambda/(1.0-pval))**(1.0/utpar)
        endif
    else
    
    ! Value within the transformation table:
    
        call locate(cdf,nt,1,nt,pval,j)
        j    = max(min((nt-1),j),1)
        getz = powint(cdf(j),cdf(j+1),vr(j),vr(j+1),pval,1.0)
    endif
    if(getz < zmin) getz = zmin
    if(getz > zmax) getz = zmax
    return
    end function getz
