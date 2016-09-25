    subroutine picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup, &
    irot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr, &
    iysbtosr,izsbtosr)
!-----------------------------------------------------------------------

!             Establish Which Super Blocks to Search
!             **************************************

! This subroutine establishes which super blocks must be searched given
! that a point being estimated/simulated falls within a super block
! centered at 0,0,0.



! INPUT VARIABLES:

!   nxsup,xsizsup    Definition of the X super block grid
!   nysup,ysizsup    Definition of the Y super block grid
!   nzsup,zsizsup    Definition of the Z super block grid
!   irot             index of the rotation matrix for searching
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!   radsqd           squared search radius



! OUTPUT VARIABLES:

!   nsbtosr          Number of super blocks to search
!   ixsbtosr         X offsets for super blocks to search
!   iysbtosr         Y offsets for super blocks to search
!   izsbtosr         Z offsets for super blocks to search



! EXTERNAL REFERENCES:

!   sqdist           Computes anisotropic squared distance



!-----------------------------------------------------------------------
    real*8 ::  rotmat(MAXROT,3,3),hsqd,sqdist,shortest
    integer :: ixsbtosr(*),iysbtosr(*),izsbtosr(*)

! MAIN Loop over all possible super blocks:

    nsbtosr = 0
    do i=-(nxsup-1),(nxsup-1)
        do j=-(nysup-1),(nysup-1)
            do k=-(nzsup-1),(nzsup-1)
                xo = real(i)*xsizsup
                yo = real(j)*ysizsup
                zo = real(k)*zsizsup
            
            ! Find the closest distance between the corners of the super blocks:
            
                shortest = 1.0e21
                do i1=-1,1
                    do j1=-1,1
                        do k1=-1,1
                            do i2=-1,1
                                do j2=-1,1
                                    do k2=-1,1
                                        if(i1 /= 0 .AND. j1 /= 0 .AND. k1 /= 0 .AND. &
                                        i2 /= 0 .AND. j2 /= 0 .AND. k2 /= 0) then
                                            xdis = real(i1-i2)*0.5*xsizsup + xo
                                            ydis = real(j1-j2)*0.5*ysizsup + yo
                                            zdis = real(k1-k2)*0.5*zsizsup + zo
                                            hsqd = sqdist(0.0,0.0,0.0,xdis,ydis,zdis, &
                                            irot,MAXROT,rotmat)
                                            if(hsqd < shortest) shortest = hsqd
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            
            ! Keep this super block if it is close enoutgh:
            
                if(real(shortest) <= radsqd) then
                    nsbtosr = nsbtosr + 1
                    ixsbtosr(nsbtosr) = i
                    iysbtosr(nsbtosr) = j
                    izsbtosr(nsbtosr) = k
                end if
            end do
        end do
    end do

! Finished:

    return
    end subroutine picksup
