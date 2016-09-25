    subroutine setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z, &
    vr,tmp,nsec,sec1,sec2,sec3,MAXSBX,MAXSBY, &
    MAXSBZ,nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup, &
    ysizsup,nzsup,zmnsup,zsizsup)
!-----------------------------------------------------------------------

!           Establish Super Block Search Limits and Sort Data
!           *************************************************

! This subroutine sets up a 3-D "super block" model and orders the data
! by super block number.  The limits of the super block is set to the
! minimum and maximum limits of the grid; data outside are assigned to
! the nearest edge block.

! The idea is to establish a 3-D block network that contains all the
! relevant data.  The data are then sorted by their index location in
! the search network, i.e., the index location is given after knowing
! the block index in each coordinate direction (ix,iy,iz):
!          ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix
! An array, the same size as the number of super blocks, is constructed
! that contains the cumulative number of data in the model.  With this
! array it is easy to quickly check what data are located near any given
! location.



! INPUT VARIABLES:

!   nx,xmn,xsiz      Definition of the X grid being considered
!   ny,ymn,ysiz      Definition of the Y grid being considered
!   nz,zmn,zsiz      Definition of the Z grid being considered
!   nd               Number of data
!   x(nd)            X coordinates of the data
!   y(nd)            Y coordinates of the data
!   z(nd)            Z coordinates of the data
!   vr(nd)           Variable at each location.
!   tmp(nd)          Temporary storage to keep track of the super block
!                      index associated to each data (uses the same
!                      storage already allocated for the simulation)
!   nsec             Number of secondary variables to carry with vr
!   sec1(nd)         First secondary variable (if nsec >= 1)
!   sec2(nd)         Second secondary variable (if nsec >= 2)
!   sec3(nd)         Third secondary variable (if nsec = 3)
!   MAXSB[X,Y,Z]     Maximum size of super block network



! OUTPUT VARIABLES:

!   nisb()                Array with cumulative number of data in each
!                           super block.
!   nxsup,xmnsup,xsizsup  Definition of the X super block grid
!   nysup,ymnsup,ysizsup  Definition of the Y super block grid
!   nzsup,zmnsup,zsizsup  Definition of the Z super block grid



! EXTERNAL REFERENCES:

!   sortem           Sorting routine to sort the data



!-----------------------------------------------------------------------
    real ::    x(*),y(*),z(*),vr(*),tmp(*),sec1(*),sec2(*),sec3(*)
    integer :: nisb(*)
    logical :: inflag

! Establish the number and size of the super blocks:

    nxsup   = min(nx,MAXSBX)
    nysup   = min(ny,MAXSBY)
    nzsup   = min(nz,MAXSBZ)
    xsizsup = real(nx)*xsiz/real(nxsup)
    ysizsup = real(ny)*ysiz/real(nysup)
    zsizsup = real(nz)*zsiz/real(nzsup)
    xmnsup  = (xmn-0.5*xsiz)+0.5*xsizsup
    ymnsup  = (ymn-0.5*ysiz)+0.5*ysizsup
    zmnsup  = (zmn-0.5*zsiz)+0.5*zsizsup

! Initialize the extra super block array to zeros:

    do i=1,nxsup*nysup*nzsup
        nisb(i) = 0
    end do

! Loop over all the data assigning the data to a super block and
! accumulating how many data are in each super block:

    do i=1,nd
        call getindx(nxsup,xmnsup,xsizsup,x(i),ix,inflag)
        call getindx(nysup,ymnsup,ysizsup,y(i),iy,inflag)
        call getindx(nzsup,zmnsup,zsizsup,z(i),iz,inflag)
        ii = ix + (iy-1)*nxsup + (iz-1)*nxsup*nysup
        tmp(i)   = ii
        nisb(ii) = nisb(ii) + 1
    end do

! Sort the data by ascending super block number:

    nsort = 4 + nsec
    call sortem(1,nd,tmp,nsort,x,y,z,vr,sec1,sec2,sec3)

! Set up array nisb with the starting address of the block data:

    do i=1,(nxsup*nysup*nzsup-1)
        nisb(i+1) = nisb(i) + nisb(i+1)
    end do

! Finished:

    return
    end subroutine setsupr
