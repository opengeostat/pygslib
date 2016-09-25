    subroutine getindx(n,min,siz,loc,index,inflag)
!-----------------------------------------------------------------------

!     Gets the coordinate index location of a point within a grid
!     ***********************************************************


! n       number of "nodes" or "cells" in this coordinate direction
! min     origin at the center of the first cell
! siz     size of the cells
! loc     location of the point being considered
! index   output index within [1,n]
! inflag  true if the location is actually in the grid (false otherwise
!         e.g., if the location is outside then index will be set to
!         nearest boundary



!-----------------------------------------------------------------------
    integer ::   n,index
    real ::      min,siz,loc
    logical ::   inflag

! Compute the index of "loc":

    index = int( (loc-min)/siz + 1.5 )

! Check to see if in or out:

    if(index < 1) then
        index  = 1
        inflag = .FALSE. 
    else if(index > n) then
        index  = n
        inflag = .FALSE. 
    else
        inflag = .TRUE. 
    end if

! Return to calling program:

    return
    end subroutine getindx
