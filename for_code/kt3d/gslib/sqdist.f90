
real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
!-----------------------------------------------------------------------

!    Squared Anisotropic Distance Calculation Given Matrix Indicator
!    ***************************************************************

! This routine calculates the anisotropic distance between two points
!  given the coordinates of each point and a definition of the
!  anisotropy.


! INPUT VARIABLES:

!   x1,y1,z1         Coordinates of first point
!   x2,y2,z2         Coordinates of second point
!   ind              The rotation matrix to use
!   MAXROT           The maximum number of rotation matrices dimensioned
!   rotmat           The rotation matrices



! OUTPUT VARIABLES:

!   sqdist           The squared distance accounting for the anisotropy
!                      and the rotation of coordinates (if any).


! NO EXTERNAL REFERENCES


!-----------------------------------------------------------------------
    real*8 rotmat(MAXROT,3,3),cont,dx,dy,dz

! Compute component distance vectors and the squared distance:

    dx = dble(x1 - x2)
    dy = dble(y1 - y2)
    dz = dble(z1 - z2)
    sqdist = 0.0
    do i=1,3
        cont   = rotmat(ind,i,1) * dx &
        + rotmat(ind,i,2) * dy &
        + rotmat(ind,i,3) * dz
        sqdist = sqdist + cont * cont
    end do
    return

end function
