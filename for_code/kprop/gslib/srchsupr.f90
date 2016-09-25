    subroutine srchsupr(xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat, &
    nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd, &
    x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup, &
    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup, &
    nclose,close,infoct)
!-----------------------------------------------------------------------

!              Search Within Super Block Search Limits
!              ***************************************


! This subroutine searches through all the data that have been tagged in
! the super block subroutine.  The close data are passed back in the
! index array "close".  An octant search is allowed.



! INPUT VARIABLES:

!   xloc,yloc,zloc   location of point being estimated/simulated
!   radsqd           squared search radius
!   irot             index of the rotation matrix for searching
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!   nsbtosr          Number of super blocks to search
!   ixsbtosr         X offsets for super blocks to search
!   iysbtosr         Y offsets for super blocks to search
!   izsbtosr         Z offsets for super blocks to search
!   noct             If >0 then data will be partitioned into octants
!   nd               Number of data
!   x(nd)            X coordinates of the data
!   y(nd)            Y coordinates of the data
!   z(nd)            Z coordinates of the data
!   tmp(nd)          Temporary storage to keep track of the squared
!                      distance associated with each data
!   nisb()                Array with cumulative number of data in each
!                           super block.
!   nxsup,xmnsup,xsizsup  Definition of the X super block grid
!   nysup,ymnsup,ysizsup  Definition of the X super block grid
!   nzsup,zmnsup,zsizsup  Definition of the X super block grid



! OUTPUT VARIABLES:

!   nclose           Number of close data
!   close()          Index of close data
!   infoct           Number of informed octants (only computes if
!                      performing an octant search)



! EXTERNAL REFERENCES:

!   sqdist           Computes anisotropic squared distance
!   sortem           Sorts multiple arrays in ascending order



!-----------------------------------------------------------------------
    real ::    x(*),y(*),z(*),tmp(*),close(*)
    real*8 ::  rotmat(MAXROT,3,3),hsqd,sqdist
    integer :: nisb(*),inoct(8)
    integer :: ixsbtosr(*),iysbtosr(*),izsbtosr(*)
    logical :: inflag

! Determine the super block location of point being estimated:

    call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
    call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
    call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)

! Loop over all the possible Super Blocks:

    nclose = 0
    do 1 isup=1,nsbtosr
    
    ! Is this super block within the grid system:
    
        ixsup = ix + ixsbtosr(isup)
        iysup = iy + iysbtosr(isup)
        izsup = iz + izsbtosr(isup)
        if(ixsup <= 0 .OR. ixsup > nxsup .OR. &
        iysup <= 0 .OR. iysup > nysup .OR. &
        izsup <= 0 .OR. izsup > nzsup) go to 1
    
    ! Figure out how many samples in this super block:
    
        ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
        if(ii == 1) then
            nums = nisb(ii)
            i    = 0
        else
            nums = nisb(ii) - nisb(ii-1)
            i    = nisb(ii-1)
        endif
    
    ! Loop over all the data in this super block:
    
        do 2 ii=1,nums
            i = i + 1
        
        ! Check squared distance:
        
            hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot, &
            MAXROT,rotmat)
            if(real(hsqd) > radsqd) go to 2
        
        ! Accept this sample:
        
            nclose = nclose + 1
            close(nclose) = real(i)
            tmp(nclose)  = real(hsqd)
        2 END DO
    1 END DO

! Sort the nearby samples by distance to point being estimated:

    call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)

! If we aren't doing an octant search then just return:

    if(noct <= 0) return

! PARTITION THE DATA INTO OCTANTS:

    do i=1,8
        inoct(i) = 0
    end do

! Now pick up the closest samples in each octant:

    nt = 8*noct
    na = 0
    do j=1,nclose
        i  = int(close(j))
        h  = tmp(j)
        dx = x(i) - xloc
        dy = y(i) - yloc
        dz = z(i) - zloc
        if(dz < 0.) go to 5
        iq=4
        if(dx <= 0.0 .AND. dy > 0.0) iq=1
        if(dx > 0.0 .AND. dy >= 0.0) iq=2
        if(dx < 0.0 .AND. dy <= 0.0) iq=3
        go to 6
        5 iq=8
        if(dx <= 0.0 .AND. dy > 0.0) iq=5
        if(dx > 0.0 .AND. dy >= 0.0) iq=6
        if(dx < 0.0 .AND. dy <= 0.0) iq=7
        6 continue
        inoct(iq) = inoct(iq) + 1
    
    ! Keep this sample if the maximum has not been exceeded:
    
        if(inoct(iq) <= noct) then
            na = na + 1
            close(na) = i
            tmp(na)   = h
            if(na == nt) go to 7
        endif
    end do

! End of data selection. Compute number of informed octants and return:

    7 nclose = na
    infoct = 0
    do i=1,8
        if(inoct(i) > 0) infoct = infoct + 1
    end do

! Finished:

    return
    end subroutine srchsupr
