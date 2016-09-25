    subroutine locate(xx,n,is,ie,x,j)
!-----------------------------------------------------------------------

! Given an array "xx" of length "n", and given a value "x", this routine
! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
! returned to indicate that x is out of range.

! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
!-----------------------------------------------------------------------
    dimension xx(n)

! Initialize lower and upper methods:

    if(is <= 0) is = 1
    jl = is-1
    ju = ie
    if(xx(n) <= x) then
        j = ie
        return
    end if

! If we are not done then compute a midpoint:

    10 if(ju-jl > 1) then
        jm = (ju+jl)/2
    
    ! Replace the lower or upper limit with the midpoint:
    
        if((xx(ie) > xx(is)).eqv.(x > xx(jm))) then
            jl = jm
        else
            ju = jm
        endif
        go to 10
    endif

! Return with the array index:

    j = jl
    return
    end subroutine locate
