    real function powint(xlow,xhigh,ylow,yhigh,xval,pow)
!-----------------------------------------------------------------------

! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
!                 for a value of x and a power pow.

!-----------------------------------------------------------------------
    parameter(EPSLON=1.0e-20)

    if((xhigh-xlow) < EPSLON) then
        powint = (yhigh+ylow)/2.0
    else
        powint = ylow + (yhigh-ylow)* &
        (((xval-xlow)/(xhigh-xlow))**pow)
    end if

    return
    end function powint
