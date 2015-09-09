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
    implicit none


	!external references
	real*8,external :: sqdist

    !input
    integer, intent(in) :: nxsup, nysup, nzsup
    integer, intent(in) :: irot, MAXROT, nsbtosr, noct, nd
    real*8, intent(in), dimension (MAXROT,3,3) :: rotmat
    real*8, intent(in) ::  xloc,yloc,zloc,radsqd 
    integer, intent(in), dimension(nsbtosr) ::ixsbtosr, iysbtosr, izsbtosr
    real*8, intent(in), dimension(nd) :: x, y, z
    integer, intent(in), dimension (nxsup*nysup*nzsup) :: nisb
    real*8, intent(in) :: xmnsup,xsizsup,ymnsup, ysizsup, zmnsup, zsizsup

	!inout
    real*8, intent(inout), dimension(nd) :: tmp

    !out
    real*8, intent(out), dimension(nd) :: close
    integer, intent(out) :: infoct, nclose


    !internal
    real*8 ::  hsqd, h, dx, dy, dz
    integer :: inoct(8)
    logical :: inflag
    integer :: isup,ixsup, iysup, izsup, i, ii, nt, na, nums, j, ix, iy, iz, iq

    ! TODO: this is for sortem, see if this can be modified and avoid declaring non used arrays (e, f, g, hh)
    real*8, dimension(nd) :: c, d, e, f, g, hh


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
            if(dble(hsqd) > radsqd) go to 2
        
            ! Accept this sample:
        
            nclose = nclose + 1
            close(nclose) = dble(i)
            tmp(nclose)  = dble(hsqd)
        2 END DO
    1 END DO

    ! Sort the nearby samples by distance to point being estimated:

    call sortem(1,nclose,tmp,1,close,c,d,e,f,g,hh)

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
    !                         super block.
    !   nxsup,xmnsup,xsizsup  Definition of the X super block grid
    !   nysup,ymnsup,ysizsup  Definition of the Y super block grid
    !   nzsup,zmnsup,zsizsup  Definition of the Z super block grid



    ! EXTERNAL REFERENCES:

    !   sortem           Sorting routine to sort the data



    !-----------------------------------------------------------------------

    implicit none
    !input
    real*8, intent(in), dimension (nd) ::  x,y,z,vr,sec1,sec2,sec3
    integer , intent(in) :: nx, ny, nz, nd, nsec, MAXSBX,MAXSBY,MAXSBZ
    real*8 , intent(in) :: xmn,xsiz,ymn,ysiz,zmn,zsiz

    !inout (warning...)
    real*8, intent(inout), dimension (nd) ::  tmp    

    !output
    !nisb(nxsup*nysup*nzsup) not possible, assigning large value or MAXSBX*MAXSBY*MAXSBZ
    integer, intent(out), dimension (MAXSBX*MAXSBY*MAXSBZ) :: nisb
    integer, intent(out) :: nxsup, nysup, nzsup
    real*8, intent(out) :: xmnsup,xsizsup,ymnsup, ysizsup,zmnsup,zsizsup


    !local
    logical :: inflag
    integer :: i, ii, nsort, ix, iy, iz

    ! Establish the number and size of the super blocks:

    nxsup   = min(nx,MAXSBX)
    nysup   = min(ny,MAXSBY)
    nzsup   = min(nz,MAXSBZ)
    xsizsup = dble(nx)*xsiz/dble(nxsup)
    ysizsup = dble(ny)*ysiz/dble(nysup)
    zsizsup = dble(nz)*zsiz/dble(nzsup)
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
    implicit none


    !external references
    real*8,external :: sqdist
    

    !input
    real*8, intent(in), dimension (MAXROT,3,3) :: rotmat
    integer, intent(in) :: nxsup, nysup, nzsup, irot, MAXROT
    real*8, intent(in) ::  xsizsup, ysizsup, zsizsup, radsqd


    !output
    integer, intent(out) :: nsbtosr
    integer, intent(out) :: ixsbtosr(*),iysbtosr(*),izsbtosr(*)
    

    !local
    integer :: i, j, k, i1, j1, k1, i2, j2, k2
    real*8 :: xo, yo,zo, hsqd, shortest, xdis, ydis, zdis



    ! MAIN Loop over all possible super blocks:

    nsbtosr = 0
    do i=-(nxsup-1),(nxsup-1)
        do j=-(nysup-1),(nysup-1)
            do k=-(nzsup-1),(nzsup-1)
                xo = dble(i)*xsizsup
                yo = dble(j)*ysizsup
                zo = dble(k)*zsizsup
            
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
                                            xdis = dble(i1-i2)*0.5*xsizsup + xo
                                            ydis = dble(j1-j2)*0.5*ysizsup + yo
                                            zdis = dble(k1-k2)*0.5*zsizsup + zo
                                            hsqd = sqdist(dble(0.0),dble(0.0),dble(0.0),xdis,ydis,zdis, &
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
            
                if(dble(shortest) <= radsqd) then
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
    implicit none
    

    !input
    integer, intent(in) :: n
    real*8, intent(in) ::  min, siz, loc

    !output
    integer, intent(out) :: index 
    logical, intent(out) :: inflag



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


subroutine sortem(ib,ie,a,iperm,b,c,d,e,f,g,h)
    !-----------------------------------------------------------------------

    !                      Quickersort Subroutine
    !                      **********************

    ! This is a subroutine for sorting a real array in ascending order. This
    ! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
    ! in collected algorithms of the ACM.

    ! The method used is that of continually splitting the array into parts
    ! such that all elements of one part are less than all elements of the
    ! other, with a third part in the middle consisting of one element.  An
    ! element with value t is chosen arbitrarily (here we choose the middle
    ! element). i and j give the lower and upper limits of the segment being
    ! split.  After the split a value q will have been found such that
    ! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
    ! performs operations on the two segments (i,q-1) and (q+1,j) as follows
    ! The smaller segment is split and the position of the larger segment is
    ! stored in the lt and ut arrays.  If the segment to be split contains
    ! two or fewer elements, it is sorted and another segment is obtained
    ! from the lt and ut arrays.  When no more segments remain, the array
    ! is completely sorted.


    ! INPUT PARAMETERS:

    !   ib,ie        start and end index of the array to be sorteda
    !   a            array, a portion of which has to be sorted.
    !   iperm        0 no other array is permuted.
    !                1 array b is permuted according to array a
    !                2 arrays b,c are permuted.
    !                3 arrays b,c,d are permuted.
    !                4 arrays b,c,d,e are permuted.
    !                5 arrays b,c,d,e,f are permuted.
    !                6 arrays b,c,d,e,f,g are permuted.
    !                7 arrays b,c,d,e,f,g,h are permuted.
    !               >7 no other array is permuted.

    !   b,c,d,e,f,g,h  arrays to be permuted according to array a.

    ! OUTPUT PARAMETERS:

    !    a      = the array, a portion of which has been sorted.

    !    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)

    ! NO EXTERNAL ROUTINES REQUIRED:

    !-----------------------------------------------------------------------
    real*8 :: a,b,c,d,e,f,g,h
    dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)

    ! The dimensions for lt and ut have to be at least log (base 2) n

    integer ::   lt(64),ut(64),i,j,k,m,p,q

    ! Initialize:

    j     = ie
    m     = 1
    i     = ib
    iring = iperm+1
    if (iperm > 7) iring=1

    ! If this segment has more than two elements  we split it

    10 if (j-i-1) 100,90,15

    ! p is the position of an arbitrary element in the segment we choose the  
    ! middle element. Under certain circumstances it may be advantageous
    ! to choose p at random.

    15 p    = (j+i)/2
    ta   = a(p)
    a(p) = a(i)
    go to (21,19,18,17,16,161,162,163),iring
    163 th   = h(p)
    h(p) = h(i)
    162 tg   = g(p)
    g(p) = g(i)
    161 tf   = f(p)
    f(p) = f(i)
    16 te   = e(p)
    e(p) = e(i)
    17 td   = d(p)
    d(p) = d(i)
    18 tc   = c(p)
    c(p) = c(i)
    19 tb   = b(p)
    b(p) = b(i)
    21 continue

    ! Start at the beginning of the segment, search for k such that a(k)>t

    q = j
    k = i
    20 k = k+1
    if(k > q)     go to 60
    if(a(k) <= ta) go to 20

    ! Such an element has now been found now search for a q such that a(q)<t
    ! starting at the end of the segment.

    30 continue
    if(a(q) < ta) go to 40
    q = q-1
    if(q > k)     go to 30
    go to 50

    ! a(q) has now been found. we interchange a(q) and a(k)

    40 xa   = a(k)
    a(k) = a(q)
    a(q) = xa
    go to (45,44,43,42,41,411,412,413),iring
    413 xh   = h(k)
    h(k) = h(q)
    h(q) = xh
    412 xg   = g(k)
    g(k) = g(q)
    g(q) = xg
    411 xf   = f(k)
    f(k) = f(q)
    f(q) = xf
    41 xe   = e(k)
    e(k) = e(q)
    e(q) = xe
    42 xd   = d(k)
    d(k) = d(q)
    d(q) = xd
    43 xc   = c(k)
    c(k) = c(q)
    c(q) = xc
    44 xb   = b(k)
    b(k) = b(q)
    b(q) = xb
    45 continue

    ! Update q and search for another pair to interchange:

    q = q-1
    go to 20
    50 q = k-1
    60 continue

    ! The upwards search has now met the downwards search:

    a(i)=a(q)
    a(q)=ta
    go to (65,64,63,62,61,611,612,613),iring
    613 h(i) = h(q)
    h(q) = th
    612 g(i) = g(q)
    g(q) = tg
    611 f(i) = f(q)
    f(q) = tf
    61 e(i) = e(q)
    e(q) = te
    62 d(i) = d(q)
    d(q) = td
    63 c(i) = c(q)
    c(q) = tc
    64 b(i) = b(q)
    b(q) = tb
    65 continue

    ! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
    ! store the position of the largest segment in lt and ut

    if (2*q <= i+j) go to 70
    lt(m) = i
    ut(m) = q-1
    i = q+1
    go to 80
    70 lt(m) = q+1
    ut(m) = j
    j = q-1

    ! Update m and split the new smaller segment

    80 m = m+1
    go to 10

    ! We arrive here if the segment has  two elements we test to see if
    ! the segment is properly ordered if not, we perform an interchange

    90 continue
    if (a(i) <= a(j)) go to 100
    xa=a(i)
    a(i)=a(j)
    a(j)=xa
    go to (95,94,93,92,91,911,912,913),iring
    913 xh   = h(i)
    h(i) = h(j)
    h(j) = xh
    912 xg   = g(i)
    g(i) = g(j)
    g(j) = xg
    911 xf   = f(i)
    f(i) = f(j)
    f(j) = xf
    91 xe   = e(i)
    e(i) = e(j)
    e(j) = xe
    92 xd   = d(i)
    d(i) = d(j)
    d(j) = xd
    93 xc   = c(i)
    c(i) = c(j)
    c(j) = xc
    94 xb   = b(i)
    b(i) = b(j)
    b(j) = xb
    95 continue

    ! If lt and ut contain more segments to be sorted repeat process:

    100 m = m-1
    if (m <= 0) go to 110
    i = lt(m)
    j = ut(m)
    go to 10
    110 continue
    return
end subroutine sortem

subroutine cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,c0,it,cc,aa, &
    irot,maxrot,rotmat,cmax,cova)
    !-----------------------------------------------------------------------

    !                    Covariance Between Two Points
    !                    *****************************

    ! This subroutine calculated the covariance associated with a variogram
    ! model specified by a nugget effect and nested varigoram structures.
    ! The anisotropy definition can be different for each nested structure.



    ! INPUT VARIABLES:

    !   x1,y1,z1         coordinates of first point
    !   x2,y2,z2         coordinates of second point
    !   nst(ivarg)       number of nested structures (maximum of 4)
    !   ivarg            variogram number (set to 1 unless doing cokriging
    !                       or indicator kriging and in this case use same number of structures in all variograms)
    !   c0(ivarg)        isotropic nugget constant
    !   it(i)            type of each nested structure:
    !                      1. spherical model of range a;
    !                      2. exponential model of parameter a;
    !                           i.e. practical range is 3a
    !                      3. gaussian model of parameter a;
    !                           i.e. practical range is a*sqrt(3)
    !                      4. power model of power a (a must be gt. 0  and
    !                           lt. 2).  if linear model, a=1,c=slope.
    !                      5. hole effect model
    !   cc(i)            multiplicative factor of each nested structure.
    !                      (sill-c0) for spherical, exponential,and gaussian
    !                      slope for linear model.
    !   aa(i)            parameter "a" of each nested structure.
    !   irot             index of the rotation matrix for the first nested
    !                    structure (the second nested structure will use
    !                    irot+1, the third irot+2, and so on)
    !   maxrot           size of rotation matrix arrays
    !   rotmat           rotation matrices
    ! 
    !  Note that aa, it and cc are 1D arrays with size (MXVARG*MAXNST).     
    !  MAXNST was removed from this code and recalculated as MAXROT=MAXNST+1
    !  MXVARG is equal to ivarg
    ! 
    ! OUTPUT VARIABLES:

    !   cmax             maximum covariance
    !   cova             covariance between (x1,y1,z1) and (x2,y2,z2)



    ! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
    !                      rotmat    computes rotation matrix for distance
    !-----------------------------------------------------------------------

    implicit none
    !external references
    real*8,external :: sqdist


    ! input
    real*8, intent(in) :: x1,y1,z1,x2,y2,z2
    integer, intent(in) :: ivarg, irot,maxrot
    integer, intent(in), dimension(ivarg) :: nst, it
    real*8, intent(in), dimension(ivarg) :: c0
    real*8, intent(in), dimension((maxrot-1)*ivarg) :: cc, aa
    real*8, intent(in), dimension(maxrot,3,3) :: rotmat

    ! output
    real*8, intent(out) :: cmax, cova

    ! internal variables
    real*8 ::    hsqd, maxnst, h, hr
    integer :: ir, is, ist, istart
    
    !parameters
    real*8 :: DEG2RAD,EPSLON,PI,PMX
 
    DEG2RAD=3.141592654/180.0
    EPSLON=1.e-20
    PI=3.14159265
    PMX=999.
    EPSLON=1.e-10


 
    ! Calculate the maximum covariance value (used for zero distances and
    ! for power model covariance):

    maxnst=maxrot-1
    istart = 1 + (ivarg-1)*maxnst
    cmax   = c0(ivarg)
    do is=1,nst(ivarg)
        ist = istart + is - 1
        if(it(ist) == 4) then
            cmax = cmax + PMX
        else
            cmax = cmax + cc(ist)
        endif
    end do

    ! Check for "zero" distance, return with cmax if so:

    hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,maxrot,rotmat)
    if(dble(hsqd) < EPSLON) then
        cova = cmax
        return
    endif

    ! Loop over all the structures:

    cova = 0.0
    do is=1,nst(ivarg)
        ist = istart + is - 1
    
    ! Compute the appropriate distance:
    
        if(ist /= 1) then
            ir = min((irot+is-1),maxrot)
            hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,maxrot,rotmat)
        end if
        h = dble(dsqrt(hsqd))
    
    ! Spherical Variogram Model?
    
        if(it(ist) == 1) then
            hr = h/aa(ist)
            if(hr < 1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
        
        ! Exponential Variogram Model?
        
        else if(it(ist) == 2) then
            cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
        
        ! Gaussian Variogram Model?
        
        else if(it(ist) == 3) then
            cova = cova + cc(ist)*exp(-(3.0*h/aa(ist)) &
            *(3.0*h/aa(ist)))
        
        ! Power Variogram Model?
        
        else if(it(ist) == 4) then
            cova = cova + cmax - cc(ist)*(h**aa(ist))
        
        ! Hole Effect Model?
        
        else if(it(ist) == 5) then
        !                 d = 10.0 * aa(ist)
        !                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
            cova = cova + cc(ist)*cos(h/aa(ist)*PI)
        endif
    end do

    ! Finished:

    return

end subroutine cova3

!*********************************************************************************
!     Subroutines in GSLIB (auxiliary functions)
!*********************************************************************************
subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,maxrot,rotmat)
    !-----------------------------------------------------------------------

    !              Sets up an Anisotropic Rotation Matrix
    !              **************************************

    ! Sets up the matrix to transform Cartesian coordinates to coordinates
    ! accounting for angles and anisotropy (see manual for a detailed
    ! definition):


    ! INPUT PARAMETERS:

    !   ang1             Azimuth angle for principal direction
    !   ang2             Dip angle for principal direction
    !   ang3             Third rotation angle
    !   anis1            First anisotropy ratio
    !   anis2            Second anisotropy ratio
    !   ind              matrix indicator to initialize
    !   maxrot           maximum number of rotation matrices dimension 
    !                    for example maxrot = number of structures in  variogram + 1
    !   rotmat           rotation matrices


    ! This code was modified from original f77 GSLIB code (v.2)
    ! Mayor changes
    ! rotmat is dynamically defined
    ! DEG2RAD,EPSLON redefined as variables
    ! maxrot is now variable defined externally... 


    !-----------------------------------------------------------------------

    implicit none

    
    
    ! input
    real*8, intent(in) :: ang1,ang2,ang3,anis1,anis2
    integer, intent(in) :: ind, maxrot

    ! output
    real*8, intent(out), dimension(maxrot,3,3) :: rotmat

    ! internal variables
    real*8 ::  afac1,afac2,sina,sinb,sint, cosa,cosb,cost, alpha, beta, theta
    
    !parameters
    real*8 :: DEG2RAD,EPSLON

    DEG2RAD=3.141592654/180.0
    EPSLON=1.e-20

    ! Converts the input angles to three angles which make more
    !  mathematical sense:

    !         alpha   angle between the major axis of anisotropy and the
    !                 E-W axis. Note: Counter clockwise is positive.
    !         beta    angle between major axis and the horizontal plane.
    !                 (The dip of the ellipsoid measured positive down)
    !         theta   Angle of rotation of minor axis about the major axis
    !                 of the ellipsoid.

    if(ang1 >= 0.0 .AND. ang1 < 270.0) then
        alpha = (90.0   - ang1) * DEG2RAD
    else
        alpha = (450.0  - ang1) * DEG2RAD
    endif
    beta  = -1.0 * ang2 * DEG2RAD
    theta =        ang3 * DEG2RAD

    ! Get the required sines and cosines:

    sina  = dble(sin(alpha))
    sinb  = dble(sin(beta))
    sint  = dble(sin(theta))
    cosa  = dble(cos(alpha))
    cosb  = dble(cos(beta))
    cost  = dble(cos(theta))

    ! Construct the rotation matrix in the required memory:

    afac1 = 1.0 / dble(max(anis1,EPSLON))
    afac2 = 1.0 / dble(max(anis2,EPSLON))
    rotmat(ind,1,1) =       (cosb * cosa)
    rotmat(ind,1,2) =       (cosb * sina)
    rotmat(ind,1,3) =       (-sinb)
    rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
    rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
    rotmat(ind,2,3) = afac1*( sint * cosb)
    rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
    rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
    rotmat(ind,3,3) = afac2*(cost * cosb)

    ! Return to calling program:

    return
end subroutine setrot



real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,maxrot,rotmat)
    !-----------------------------------------------------------------------
    !
    !    Squared Anisotropic Distance Calculation Given Matrix Indicator
    !    ***************************************************************
    !
    ! This routine calculates the anisotropic distance between two points
    !  given the coordinates of each point and a definition of the
    !  anisotropy.
    !
    !
    ! INPUT VARIABLES:
    !
    !   x1,y1,z1         Coordinates of first point
    !   x2,y2,z2         Coordinates of second point
    !   ind              The rotation matrix to use
    !   maxrot           The maximum number of rotation matrices dimensioned
    !   rotmat           The rotation matrices
    !
    !
    !
    ! OUTPUT VARIABLES:
    !
    !   sqdis           The squared distance accounting for the anisotropy
    !                      and the rotation of coordinates (if any).
    !
    !
    ! NO EXTERNAL REFERENCES
    !
    !
    !-----------------------------------------------------------------------
     
    implicit none
    
    ! input
    integer, intent(in) ::  maxrot, ind
    real*8, intent(in):: x1,y1,z1,x2,y2,z2
    real*8, intent(in), dimension(maxrot,3,3) ::  rotmat

    
    ! output
    ! real*8 :: sqdis

    ! Internal 
    real*8 :: cont,dx,dy,dz
    integer  :: i

    !
    ! Compute component distance vectors and the squared distance:
    !
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
end function sqdist


!
! #####################################################################
! #####################################################################
!

subroutine kt3d(radius,radius1,radius2,sang1,sang2,sang3, &  
                nst,c0,it,cc,aa,aa1,aa2,ang1,ang2,ang3,   &
                nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, extve,  &
                nd,x,y,z,vr, ve, &
                nclose, close, infoct, noct, &
                MAXSBX, MAXSBY, MAXSBZ, tmin, tmax,  &
                ktype, skmean,  koption, itrend, nxdis,nydis,nzdis,xdb,ydb,zdb, &
                bv,ndmax, ndmin, xa,ya,za,vra, vea, na, &
                ndjak, xlj, ylj, zlj, vrlj, extvj, &
                UNEST)
    !-----------------------------------------------------------------------

    !                Krige a 3-D Grid of Rectangular Blocks
    !                **************************************

    ! This subroutine estimates point or block values of one variable by
    ! simple, ordinary, or kriging with a trend model.  It is also possible
    ! to estimate the trend directly.



    ! PROGRAM NOTES:

    !   1. The data and parameters are passed in common blocks defined
    !      in kt3d.inc.  Local storage is allocated in the subroutine
    !      for kriging matrices, i.e.,
    !         - xa,ya,za,vra   arrays for data within search neighborhood
    !         - a,r,rr,s       kriging arrays
    !         - xdb,ydb,zdb    relative position of discretization points
    !         - cbb            block covariance
    !   2. The kriged value and the kriging variance is written to Fortran
    !      unit number "lout".




    ! Original:  A.G. Journel and C. Lemmer                             1981
    ! Revisions: A.G. Journel and C. Kostov                             1984
    
    ! 2015 changes
    ! This is a funtion to be used from python... 
    ! TODO: add comments 
    
    ! input variables
    !   - AXSAM    maximum number of data points to use in one kriging system
    !   - MAXDT     maximum number of drift terms
    !   - MAXEQ    maximum numb of equations
    !   - MAXDIS maximum number of discretization points per block
    !   - MAXNST (using nst)  maximum number of nested structures
    !   - MAXDAT    maximum number of data points
    !   - MAXSBX    maximum super block nodes in X direction
    !   - MAXSBY    maximum super block nodes in Y direction    
    !   - MAXSBZ    maximum super block nodes in Z direction
    !   - 
            
  
    
    
    ! output 
    !   - cbb block covariance
    !   - xa,ya,za,vra (MAXSAM) array with data used last block estimated
    !   - a,r,rr,s              kriging arrays
    !   - xdb,ydb,zdb           relative position of discretization points  (last block?)
    !   - bv(9)                 mean drift value 
    
    !-----------------------------------------------------------------------


    !for safety reason we don't want undeclared variables
    IMPLICIT NONE 


    ! in 
    ! variogram
    integer, intent(in), dimension(1)     :: nst
    integer, intent(in), dimension(1)     :: it
    real*8, intent(in), dimension(1)      :: c0
    real*8, intent(in), dimension(nst(1)) :: cc,aa, aa1, aa2,ang1,ang2,ang3
    ! search
    real*8, intent(in) :: radius, radius1, radius2, sang1, sang2, sang3
    integer, intent(in) :: MAXSBX, MAXSBY, MAXSBZ
    integer, intent(in) :: ndmax, ndmin
    ! grid 
    integer, intent(in) :: nx,ny,nz
    real*8 , intent(in) :: xmn,ymn,zmn,xsiz,ysiz,zsiz
    real*8 , intent(in), dimension(nx*ny*nz) :: extve
    !Jacknife
    integer, intent(in)               :: ndjak
    real*8, intent(in), dimension(nd) :: xlj, ylj, zlj, vrlj, extvj
    ! data 
    integer, intent(in)               :: nd
    real*8, intent(in), dimension(nd) :: x,y,z,vr, ve                   !vr->variable, ve->external drift
    real*8, intent(in) :: tmin, tmax
    !kriging parameters     
    integer, intent(in) :: ktype, koption, noct, itrend
    integer, intent(inout) :: nxdis,nydis,nzdis                      
    real*8, intent(in) :: UNEST
    real*8, intent(inout) :: skmean
    

    ! out
    ! kriging parameters
    ! discretization points of the last block
    real*8, intent (out), dimension(nxdis*nydis*nzdis) :: xdb,ydb,zdb      ! MAXDIS=nxdis*nydis*nzdis number of discretization points (variable)
    ! drift
    real*8, intent(out) ::  bv(9)
    ! the samples around the block (last block) we may use this to test estimate in a single block
    real*8, intent(out), dimension(nd) :: close
    integer, intent(out) :: nclose, infoct, na
    real*8, intent(out), dimension(ndmax)  ::  xa,ya,za,vra, vea
    
  
    ! internal 
    ! constants
    real*8 :: EPSLON=0.000001,PMX = 999.0
    integer :: MAXROT
    ! search 
    real*8, dimension(nst(1)+1,3,3) :: rotmat
    real*8 :: radsqd,resc,sanis1, sanis2
    integer, dimension (MAXSBX*MAXSBY*MAXSBZ) :: nisb                   ! Array with cumulative number of data in each super block.
    integer            :: nxsup, nysup, nzsup                           ! Definition of the X super block grid
    real*8             :: xmnsup,xsizsup,ymnsup, ysizsup,zmnsup,zsizsup ! Definition of the X super block grid
    integer, dimension (8*MAXSBX*MAXSBY*MAXSBZ) ::ixsbtosr,iysbtosr,izsbtosr  ! X offsets for super blocks to search
    integer :: nsbtosr                                                        ! Number of super blocks to search
    ! variogram 
    real*8 :: covmax, unbias, cbb
    real*8, dimension(nst(1)) :: anis1 , anis2
    integer :: isrot 
    ! data
    real*8, dimension(nd) :: tmp, sec2, sec3                            !required in supper-block functions
    ! drift 
    integer :: mdt
    real*8, dimension(9) :: idrif                                       !MAXDT=9     maximum number of drift terms
    ! kriging parameters
    real*8 :: xdis, ydis, zdis, xloc, yloc, zloc
    real*8 :: true, secj, extest
    real*8,  dimension((ndmax+9+2)*(ndmax+9+2))  ::  a                  ! this is MAXDT=9 and MAXEQ=MAXSAM+MAXDT+2
    real*8,  dimension(ndmax+9)  ::  r, rr, s
    integer :: ndb, iktype = 0                                          ! I'm removing this part from the code IK with KT3D
    ! other
    integer :: i, j, is, nsec, ix, iy, iz, index, ind, neq, im, k, &
               kadim, ksdim, nrhs, nv, maxeq, ising
    integer :: nxy, nxyz, nloop, irepo, nk, xkmae, xkmse        !some variables to iterate in main 
    real*8 :: est, estv, resce, cb, cb1, cmax, cov, dx, dy, dz, err, xk, vk
    logical :: accept
 
    MAXROT = nst(1) + 1
    maxeq = ndmax+9+2
    
    ! calculate anisotropy
    do i=1,nst(1)

            anis1(i) = aa1(i) / max(aa(i),EPSLON)
            anis2(i) = aa2(i) / max(aa(i),EPSLON)
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i), ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1(i),aa2(i)
            write(*,*) ' anis1 anis2: ',anis1(i),anis2(i)
            
            if(it(i).eq.4) then
                  if(aa(i).lt.0.0) stop ' INVALID power variogram'
                  if(aa(i).gt.2.0) stop ' INVALID power variogram'
            end if
    end do
   
    ! Set up the rotation/anisotropy matrices that are needed for the
    ! variogram and search.  Also compute the maximum covariance for
    ! the rescaling factor:

    write(*,*) 'Setting up rotation matrices for variogram and search'
    
    radsqd = radius * radius ! the search first radius
    if(radius.lt.EPSLON) stop 'radius must be greater than zero'
    radsqd = radius  * radius
    sanis1 = radius1 / radius
    sanis2 = radius2 / radius
    
 
    covmax = c0(1)
    do is=1,nst(1)
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,MAXROT,rotmat)
        if(it(is) == 4) then
            covmax = covmax + PMX
        else
            covmax = covmax + cc(is)
        endif
    end do
    isrot = nst(1) + 1
    
    
    call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)     ! this is the rotation matrix for the search ellipse
    

    ! Finish computing the rescaling factor and stop if unacceptable:
    
    if(radsqd < 1.0) then
        resc = 2.0 * radius / max(covmax,0.0001)
    else
        resc =(4.0 * radsqd)/ max(covmax,0.0001)
    endif
    if(resc <= 0.0) then
        write(*,*) 'ERROR KT3D: The rescaling value is wrong ',resc
        write(*,*) '            Maximum covariance: ',covmax
        write(*,*) '            search radius:      ',radius
        stop                                                            ! remove this guy
    endif
    resc = 1.0 / resc
   
   
    ! Set up for super block searching:

    write(*,*) 'Setting up super block search strategy'
    nsec = 1
    call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z, &
                 vr,tmp,nsec,ve,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,nisb, &
                 nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup, &
                 zmnsup,zsizsup)
            
                
    call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup, &
                isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr, &
                iysbtosr,izsbtosr)
    
 
    ! Compute the number of drift terms, if an external drift is being
    ! considered then it is one more drift term, if SK is being considered
    ! then we will set all the drift terms off and mdt to 0):

    mdt = 1
    do i=1,9
        if(ktype == 0 .OR. ktype == 2) idrif(i) = 0
        if(idrif(i) < 0 .OR. idrif(i) > 1) then
            write(*,*) 'ERROR KT3D: invalid drift term',idrif(i)
            stop
        endif
        mdt = mdt + idrif(i)
    end do
    if(ktype == 3) mdt = mdt + 1
    if(ktype == 0) mdt = 0
    if(ktype == 2) mdt = 0

    ! Set up the discretization points per block.  Figure out how many
    ! are needed, the spacing, and fill the xdb,ydb, and zdb arrays with
    ! the offsets relative to the block center (this only gets done once):

    ! In all cases the offsets are relative to the lower left corner.
    ! This is done for rescaling the drift terms in the kriging matrix.

    if(nxdis < 1) nxdis = 1
    if(nydis < 1) nydis = 1
    if(nzdis < 1) nzdis = 1
    ndb = nxdis * nydis * nzdis
    xdis = xsiz  / max(dble(nxdis),1.0)
    ydis = ysiz  / max(dble(nydis),1.0)
    zdis = zsiz  / max(dble(nzdis),1.0)    
    i    = 0
    xloc = -0.5*(xsiz+xdis)
    do ix =1,nxdis
        xloc = xloc + xdis
        yloc = -0.5*(ysiz+ydis)
        do iy=1,nydis
            yloc = yloc + ydis
            zloc = -0.5*(zsiz+zdis)
            do iz=1,nzdis
                zloc = zloc + zdis
                i = i+1
                xdb(i) = xloc + 0.5*xsiz
                ydb(i) = yloc + 0.5*ysiz
                zdb(i) = zloc + 0.5*zsiz
            end do
        end do
    end do
  
  
    ! Calculate Block Covariance. Check for point kriging
    call block_covariance(xdb,ydb,zdb, ndb, &
                            nst,it,c0,cc,aa, aa1, aa2,ang1,ang2,ang3, &
                            unbias, cbb)

    do i=1,9
        bv(i) = 0.0
    end do
    do i=1,ndb
        bv(1) = bv(1) + xdb(i)
        bv(2) = bv(2) + ydb(i)
        bv(3) = bv(3) + zdb(i)
        bv(4) = bv(4) + xdb(i)*xdb(i)
        bv(5) = bv(5) + ydb(i)*ydb(i)
        bv(6) = bv(6) + zdb(i)*zdb(i)
        bv(7) = bv(7) + xdb(i)*ydb(i)
        bv(8) = bv(8) + xdb(i)*zdb(i)
        bv(9) = bv(9) + ydb(i)*zdb(i)
    end do
    do i=1,9
        bv(i) = (bv(i) / dble(ndb)) * resc
    end do

    ! Report on progress from time to time:

    if(koption == 0) then
        nxy   = nx*ny
        nxyz  = nx*ny*nz
        nloop = nxyz
        irepo = max(1,min((nxyz/10),10000))
    else
        nloop = ndjak
        irepo = max(1,min((ndjak/10),10000))
    end if
    write(*,*)
    write(*,*) 'Working on the kriging '


    ! Initialize accumulators:

    nk    = 0
    xk    = 0.0
    vk    = 0.0
    xkmae = 0.0
    xkmse = 0.0

! #########################################################

    do index=1,nloop
        if((int(index/irepo)*irepo) == index) write(*,103) index
        103 format('   currently on estimate ',i9)
    
        ! Where are we making an estimate?
    
        if(koption == 0) then
            iz   = int((index-1)/nxy) + 1
            iy   = int((index-(iz-1)*nxy-1)/nx) + 1
            ix   = index - (iz-1)*nxy - (iy-1)*nx
            xloc = xmn + dble(ix-1)*xsiz
            yloc = ymn + dble(iy-1)*ysiz
            zloc = zmn + dble(iz-1)*zsiz
            print *, ' estimating in grid', xloc, yloc, zloc
        else
            xloc = xmn
            yloc = ymn
            zloc = zmn
            true = UNEST
            secj = UNEST
            xloc   = xlj(index)
            yloc   = ylj(index)
            zloc   = zlj(index)
            true   = vrlj(index)
            extest = extvj(index)
            print *, ' estimating in point', xloc, yloc, zloc
        end if

        ! Read in the external drift variable for this grid node if needed:
    
        if(ktype == 2 .OR. ktype == 3) then
            if(koption == 0) then
                extest = extve(index)
            end if
            if(extest < tmin .OR. extest >= tmax) then
                est  = UNEST
                estv = UNEST
                go to 1
            end if
            resce  = covmax / max(extest,0.0001)
        endif
        
        ! Find the nearest samples:
    
        call srchsupr(xloc,yloc,zloc,radsqd,isrot,MAXROT,rotmat,nsbtosr, &
        ixsbtosr,iysbtosr,izsbtosr,noct,nd,x,y,z,tmp, &
        nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup, &
        nzsup,zmnsup,zsizsup,nclose,close,infoct)

        ! Load the nearest data in xa,ya,za,vra,vea:
        ! TODO: Improve this with drillhole column!!!!!!!
        ! TODO: put this search part in a different function? 
        na = 0
        do i=1,nclose
            ind    = int(close(i)+0.5)
            accept = .TRUE. 
            
            if(koption /= 0 .AND. &
            (abs(x(ind)-xloc)+abs(y(ind)-yloc)+ abs(z(ind)-zloc)) &
             < EPSLON) accept = .FALSE. 
            if(accept) then
                if(na < ndmax) then
                    na = na + 1
                    xa(na)  = x(ind) - xloc + 0.5*xsiz                   !this shift the data with centre at block location
                    ya(na)  = y(ind) - yloc + 0.5*ysiz
                    za(na)  = z(ind) - zloc + 0.5*zsiz
                    vra(na) = vr(ind)
                    vea(na) = ve(ind)
                end if
            end if
        end do

        ! Test number of samples found:
    
        if(na < ndmin) then
            est  = UNEST
            estv = UNEST
            go to 1
        end if

        ! Test if there are enough samples to estimate all drift terms:
    
        if(na >= 1 .AND. na <= mdt) then
            ! pot warning/error 999 here
            est  = UNEST
            estv = UNEST
            go to 1
        end if
        ! 999 format(' Encountered a location where there were too few data ',/, &
        ! ' to estimate all of the drift terms but there would be',/, &
        ! ' enough data for OK or SK.   KT3D currently leaves ',/, &
        ! ' these locations unestimated.',/, &
        ! ' This message is only written once - the first time.',/)
    
        ! There are enough samples - proceed with estimation.
        
        if(na <= 1) then
        
            ! Handle the situation of only one sample:
        
            call cova3(xa(1),ya(1),za(1),xa(1),ya(1),za(1),1,nst, &
            c0,it,cc,aa,1,maxrot,rotmat,cmax,cb1)
       
            
            ! Establish Right Hand Side Covariance:
        
            if(ndb <= 1) then
                call cova3(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),1, &
                nst,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
            else
                cb  = 0.0
                do i=1,ndb
                    call cova3(xa(1),ya(1),za(1),xdb(i),ydb(i), &
                    zdb(i),1,nst,c0,it,cc,aa,1, &
                    MAXROT,rotmat,cmax,cov)
                    cb = cb + cov
                    dx = xa(1) - xdb(i)
                    dy = ya(1) - ydb(i)
                    dz = za(1) - zdb(i)
                    if((dx*dx+dy*dy+dz*dz) < EPSLON) cb=cb-c0(1)
                end do
                cb = cb / dble(ndb)
            end if
            est  = vra(1)
            estv = dble(cbb) - 2.0*cb + cb1
            nk   = nk + 1
            xk   = xk + vra(1)
            vk   = vk + vra(1)*vra(1)
            go to 1
        end if
    
        ! Go ahead and set up the OK portion of the kriging matrix:

        neq = mdt+na
    
        ! Initialize the main kriging matrix:
        do i=1,neq*neq
            a(i) = 0.0
        end do
    
        ! Fill in the kriging matrix:
    
        do i=1,na
            do j=i,na
                call cova3(xa(i),ya(i),za(i),xa(j),ya(j),za(j),1,nst, &
                c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                a(neq*(i-1)+j) = dble(cov)
                a(neq*(j-1)+i) = dble(cov)
            end do
        end do
        
    
        ! Fill in the OK unbiasedness portion of the matrix (if not doing SK):
    
        if(neq > na) then
            do i=1,na
                a(neq*(i-1)+na+1) = dble(unbias)
                a(neq*na+i)       = dble(unbias)
            end do
        endif
        
        ! Set up the right hand side:
    
        do i=1,na
            if(ndb <= 1) then
                call cova3(xa(i),ya(i),za(i),xdb(1),ydb(1),zdb(1),1, &
                nst,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
            else
                cb  = 0.0
                do j=1,ndb
                    call cova3(xa(i),ya(i),za(i),xdb(j),ydb(j), &
                    zdb(j),1,nst,c0,it,cc,aa,1, &
                    MAXROT,rotmat,cmax,cov)
                    cb = cb + cov
                    dx = xa(i) - xdb(j)
                    dy = ya(i) - ydb(j)
                    dz = za(i) - zdb(j)
                    if((dx*dx+dy*dy+dz*dz) < EPSLON) cb=cb-c0(1)
                end do
                cb = cb / dble(ndb)
            end if
            r(i) = dble(cb)
        end do
        if(neq > na) r(na+1) = dble(unbias)
    
        ! Add the additional unbiasedness constraints:  
    
        im = na + 1
    
        ! First drift term (linear in "x"):
    
        if(idrif(1) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(xa(k)*resc)
                a(neq*(k-1)+im) = dble(xa(k)*resc)
            end do
            r(im) = dble(bv(1))
        endif
    
        ! Second drift term (linear in "y"):
    
        if(idrif(2) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(ya(k)*resc)
                a(neq*(k-1)+im) = dble(ya(k)*resc)
            end do
            r(im) = dble(bv(2))
        endif
    
        ! Third drift term (linear in "z"):
    
        if(idrif(3) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(za(k)*resc)
                a(neq*(k-1)+im) = dble(za(k)*resc)
            end do
            r(im) = dble(bv(3))
        endif
    
    ! Fourth drift term (quadratic in "x"):
    
        if(idrif(4) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(xa(k)*xa(k)*resc)
                a(neq*(k-1)+im) = dble(xa(k)*xa(k)*resc)
            end do
            r(im) = dble(bv(4))
        endif
    
        ! Fifth drift term (quadratic in "y"):
    
        if(idrif(5) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(ya(k)*ya(k)*resc)
                a(neq*(k-1)+im) = dble(ya(k)*ya(k)*resc)
            end do
            r(im) = dble(bv(5))
        endif
    
        ! Sixth drift term (quadratic in "z"):
    
        if(idrif(6) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(za(k)*za(k)*resc)
                a(neq*(k-1)+im) = dble(za(k)*za(k)*resc)
            end do
            r(im) = dble(bv(6))
        endif
    
        ! Seventh drift term (quadratic in "xy"):
    
        if(idrif(7) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(xa(k)*ya(k)*resc)
                a(neq*(k-1)+im) = dble(xa(k)*ya(k)*resc)
            end do
            r(im) = dble(bv(7))
        endif
    
        ! Eighth drift term (quadratic in "xz"):
    
        if(idrif(8) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(xa(k)*za(k)*resc)
                a(neq*(k-1)+im) = dble(xa(k)*za(k)*resc)
            end do
            r(im) = dble(bv(8))
        endif
    
        ! Ninth drift term (quadratic in "yz"):
    
        if(idrif(9) == 1) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(ya(k)*za(k)*resc)
                a(neq*(k-1)+im) = dble(ya(k)*za(k)*resc)
            end do
            r(im) = dble(bv(9))
        endif
    
        ! External drift term (specified by external variable):
    
        if(ktype == 3) then
            im=im+1
            do k=1,na
                a(neq*(im-1)+k) = dble(vea(k)*resce)
                a(neq*(k-1)+im) = dble(vea(k)*resce)
            end do
            r(im) = dble(extest*resce)
        endif
    
        ! Copy the right hand side to compute the kriging variance later:
    
        do k=1,neq
            rr(k) = r(k)
        end do
        kadim = neq * neq
        ksdim = neq
        nrhs  = 1
        nv    = 1
    
        ! If estimating the trend then reset all the right hand side terms=0.0:
    
        if(itrend >= 1) then
            do i=1,na
                r(i)  = 0.0
                rr(i) = 0.0
            end do
        endif
    
        ! Write out the kriging Matrix if Seriously Debugging:
    
!~         if(idbg == 3) then
!~             write(ldbg,*) 'Estimating node index : ',ix,iy,iz
!~             is = 1 - neq
!~             do i=1,neq
!~                 is = 1 + (i-1)*neq
!~                 ie = is + neq - 1
!~                 write(ldbg,100) i,r(i),(a(j),j=is,ie)
!~                 100 format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
!~             end do
!~         endif
    
        ! Solve the kriging system:
    
        call ktsol(neq,nrhs,nv,a,r,s,ising,maxeq)
    
        ! Compute the solution:
    
        if(ising /= 0) then
!~             if(idbg >= 3) write(ldbg,*) ' Singular Matrix ',ix,iy,iz
            est  = UNEST
            estv = UNEST
        else
            est  = 0.0
            estv = dble(cbb)
            if(ktype == 2) skmean = extest
            do j=1,neq
                estv = estv - dble(s(j))*rr(j)
                if(j <= na) then
                    if(ktype == 0 .OR. ktype == 2) then
                        est = est + dble(s(j))*(vra(j)-skmean)
                    else
                        est = est + dble(s(j))*vra(j)
                    endif
                endif
            end do
            if(ktype == 0 .OR. ktype == 2) est = est + skmean
            nk   = nk + 1
            xk   = xk + est
            vk   = vk + est*est
        
            ! Write the kriging weights and data if debugging level is above 2:
        
!~             if(idbg >= 2) then
!~                 write(ldbg,*) '       '
!~                 write(ldbg,*) 'BLOCK: ',ix,iy,iz,' at ',xloc,yloc,zloc
!~                 write(ldbg,*) '       '
!~                 if(ktype /= 0) &
!~                 write(ldbg,*) '  Lagrange : ',s(na+1)*unbias
!~                 write(ldbg,*) '  BLOCK EST: x,y,z,vr,wt '
!~                 do i=1,na
!~                     xa(i) = xa(i) + xloc - 0.5*xsiz
!~                     ya(i) = ya(i) + yloc - 0.5*ysiz
!~                     za(i) = za(i) + zloc - 0.5*zsiz
!~                     write(ldbg,'(5f12.3)') xa(i),ya(i),za(i), &
!~                     vra(i),s(i)
!~                 end do
!~                 write(ldbg,*) '  estimate, variance  ',est,estv
!~             endif
        endif
        
        ! END OF MAIN KRIGING LOOP:
    
        1 continue

        if(iktype == 0) then
            if(koption == 0) then
!~                 write(lout,'(f9.3,1x,f9.3)') est,estv
            else
                err = UNEST
                if(true /= UNEST .AND. est /= UNEST)err=est-true
!~                 write(lout,'(7(f12.3,1x))') xloc,yloc,zloc,true, &
!~                 est,estv,err
                xkmae = xkmae + abs(err)
                xkmse = xkmse + err*err
            end if
!~         else
!~         
!~         ! Work out the IK-type distribution implicit to this data configuration
!~         ! and kriging weights:
!~         
!~             do icut=1,ncut
!~                 cdf(icut) = -1.0
!~             end do
!~             wtmin = 1.0
!~             do i=1,na
!~                 if(s(i) < wtmin) wtmin = s(i)
!~             end do
!~             sumwt = 0.0
!~             do i=1,na
!~                 s(i)  = s(i) - wtmin
!~                 sumwt = sumwt + s(i)
!~             end do
!~             do i=1,na
!~                 s(i) = s(i) / max(0.00001,sumwt)
!~             end do
!~             if(na > 1 .AND. sumwt > 0.00001) then
!~                 do icut=1,ncut
!~                     cdf(icut) = 0.0
!~                     do i=1,na
!~                         if(vra(i) <= cut(icut)) &
!~                         cdf(icut)=cdf(icut)+s(i)
!~                     end do
!~                 end do
!~             end if
!~             if(koption == 0) then
!~                 write(lout,'(30(f8.4))') (cdf(i),i=1,ncut)
!~             else
!~                 write(lout,'(30(f8.4))') (cdf(i),i=1,ncut),true
!~             end if
        end if
    end do
    2 continue
!~     if(koption > 0) close(ljack)

    ! Write statistics of kriged values:

     
    if(nk > 0 ) then
        xk    = xk/dble(nk)
        vk    = vk/dble(nk) - xk*xk
        xkmae = xkmae/dble(nk)
        xkmse = xkmse/dble(nk)
        ! write(ldbg,105) nk,xk,vk
        print *, 'Estimated blocks ', nk
        print *, 'Average          ', xk
        print *, 'Variance         ', vk
        if(koption /= 0) then
            print *, 'Mean error', xkmae
            print *, 'Mean sqd e', xkmse
        end if
    endif

    ! All finished the kriging:
    
! ###########################################################   


  
    return
    
    ! 96 stop 'ERROR in jackknife file!'                                ! remove this guy
    
end subroutine kt3d

subroutine ktsol(n,ns,nv,a,b,x,ktilt,maxeq)
    !-----------------------------------------------------------------------

    ! Solution of a system of linear equations by gaussian elimination with
    ! partial pivoting.  Several right hand side matrices and several
    ! variables are allowed.


    !         NOTE: All input matrices must be in double precision


    ! INPUT/OUTPUT VARIABLES:

    !   n                Number of equations
    !   ns               Number of right hand side matrices
    !   nv               Number of variables.
    !   a(n*n*nv)        left hand side matrices versus columnwise.
    !   b(n*ns*nv)       input right hand side matrices.
    !   x(n*ns*nv)       solution matrices.
    !   ktilt            indicator of singularity
    !                      =  0  everything is ok.
    !                      = -1 n.le.1
    !                      =  k  a null pivot appeared at the kth iteration.
    !   tol              used in test for null pivot. depends on machine
    !                      precision and can also be set for the tolerance
    !                      of an ill-defined kriging system.


!-----------------------------------------------------------------------
    
    
    implicit real*8 (a-h,o-z)
    real*8 :: x(maxeq),a(maxeq*maxeq),b(maxeq)

    ! Make sure there are equations to solve:

    if(n <= 1) then
        ktilt = -1
        return
    endif

    ! Initialization:

    tol   = 0.1e-10
    ktilt = 0
    ntn   = n*n
    nm1   = n-1

    ! Triangulation is done variable by variable:

    do iv=1,nv
    
        ! Indices of location in vectors a and b:
    
        nva = ntn*(iv-1)
        nvb = n*ns*(iv-1)
    
        ! Gaussian elimination with partial pivoting:
    
        do k=1,nm1
            kp1 = k+1
        
            ! Indice of the diagonal element in the kth row:
        
            kdiag = nva+(k-1)*n+k
        
            ! Find the pivot - interchange diagonal element/pivot:
        
            npiv = kdiag
            ipiv = k
            i1   = kdiag
            do i=kp1,n
                i1 = i1+1
                if(abs(a(i1)) > abs(a(npiv))) then
                    npiv = i1
                    ipiv = i
                endif
            end do
            t        = a(npiv)
            a(npiv)  = a(kdiag)
            a(kdiag) = t
        
            ! Test for singularity:
        
            if(abs(a(kdiag)) < tol) then
                ktilt=k
                return
            endif
        
            ! Compute multipliers:
        
            i1 = kdiag
            do i=kp1,n
                i1    = i1+1
                a(i1) = -a(i1)/a(kdiag)
            end do
        
            ! Interchange and eliminate column per column:
        
            j1 = kdiag
            j2 = npiv
            do j=kp1,n
                j1    = j1+n
                j2    = j2+n
                t     = a(j2)
                a(j2) = a(j1)
                a(j1) = t
                i1    = j1
                i2    = kdiag
                do i=kp1,n
                    i1    = i1+1
                    i2    = i2+1
                    a(i1) = a(i1)+a(i2)*a(j1)
                end do
            end do
        
            ! Interchange and modify the ns right hand matrices:
        
            i1 = nvb+ipiv
            i2 = nvb+k
            do i=1,ns
                t     = b(i1)
                b(i1) = b(i2)
                b(i2) = t
                j1    = i2
                j2    = kdiag
                do j=kp1,n
                    j1    = j1+1
                    j2    = j2+1
                    b(j1) = b(j1)+b(i2)*a(j2)
                end do
                i1 = i1+n
                i2 = i2+n
            end do
        end do
    
        ! Test for singularity for the last pivot:
    
        kdiag = ntn*iv
        if(abs(a(kdiag)) < tol) then
            ktilt = n
            return
        endif
    end do

    ! End of triangulation. Now, solve back variable per variable:

    do iv=1,nv
    
    ! Indices of location in vectors a and b:
    
        nva  = ntn*iv
        nvb1 = n*ns*(iv-1)+1
        nvb2 = n*ns*iv
    
        ! Back substitution with the ns right hand matrices:
    
        do il=1,ns
            do k=1,nm1
                nmk = n-k
            
                ! Indice of the diagonal element of the (n-k+1)th row and of
                ! the (n-k+1)th element of the left hand side.
            
                kdiag = nva-(n+1)*(k-1)
                kb    = nvb2-(il-1)*n-k+1
                b(kb) = b(kb)/a(kdiag)
                t     = -b(kb)
                i1    = kb
                i2    = kdiag
                do i=1,nmk
                    i1    = i1-1
                    i2    = i2-1
                    b(i1) = b(i1)+a(i2)*t
                end do
            end do
            kdiag = kdiag-n-1
            kb    = kb-1
            b(kb) = b(kb)/a(kdiag)
        end do
    
        ! End of back substitution:
    
    end do

    ! Restitution of the solution:

    itot = n*ns*nv
    do i=1,itot
        x(i) = b(i)
    end do

    ! Finished:

    return
end subroutine ktsol


subroutine block_covariance(xdb,ydb,zdb, ndb, &
                            nst,it,c0,cc,aa,aa1,aa2,ang1,ang2,ang3, &
                            unbias,cbb)
    ! #################################################################
    ! modified from kt3d 
    ! 
    ! to calculate the block covariance
    ! it may work with random/any location of discretization points
    !
    ! Adrian Martinez  2015
    ! ################################################################# 
    
    IMPLICIT NONE 

    ! in 
    ! variogram
    integer, intent(in)                 :: nst(1)
    real*8, intent(in)                  :: c0(1)
    integer, intent(in), dimension(nst(1)) :: it 
    real*8, intent(in), dimension(nst(1))  :: cc,aa, aa1, aa2,ang1,ang2,ang3
    !discretization points of the last block
    integer, intent (in) :: ndb   ! number of discretization points
    real*8, intent (in), dimension(ndb) :: xdb,ydb,zdb      !coordinates of discretization points

    !out 
    real*8, intent(out) :: unbias, cbb


    ! internal 
    ! variogram 
    real*8 :: cmax,covmax, cov
    real*8, dimension(nst(1)) :: anis1, anis2
    real*8, dimension(nst(1),3,3) :: rotmat
    real*8 ::  EPSLON=1.e-20, PMX=999.
    integer :: i, is, j, MAXROT
    
    do i=1,nst(1)
        anis1(i) = aa1(i) / max(aa(i),EPSLON)
        anis2(i) = aa2(i) / max(aa(i),EPSLON)
        
        print *, aa(i), aa1(i), aa2(i), anis1(i), anis2(i)
        
        if(it(i).eq.4) then
              if(aa(i).lt.0.0) stop ' INVALID power variogram'
              if(aa(i).gt.2.0) stop ' INVALID power variogram'
        end if
    end do
    
    
    
    ! get the rotation matrix
    do is=1,nst(1)
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,nst(1),rotmat)
        if(it(is) == 4) then
            covmax = covmax + PMX
        else
            covmax = covmax + cc(is)
        endif
    end do

    ! Calculate Block Covariance. Check for point kriging.
    call cova3(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst, &
                c0,it,cc,aa,1,nst(1),rotmat,cmax,cov)

    ! Set the 'unbias' variable so that the matrix solution is more stable

    unbias = cov
    cbb    = dble(cov)
    if(ndb > 1) then
        cbb = 0.0
        do i=1,ndb
            do j=1,ndb
                call cova3(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j), &
                           1,nst,c0,it,cc,aa,1,nst(1),rotmat,cmax,cov)
                if(i == j) cov = cov - c0 (1)
                cbb = cbb + dble(cov)
            end do
        end do
        cbb = cbb/dble(ndb*ndb)
    end if

    return

end subroutine block_covariance




program main

    integer :: nst(1), it(2), nx, ny, nz, nd, koption, ndjak, noct, &
               nclose, ndmax, ndmin, na, itrend
    integer :: MAXSBX =50, MAXSBY=50, MAXSBZ=50  ! this are supperbloock parameters
    
    real*8 :: radius,radius1,radius2,sang1,sang2,sang3, &
              c0(1),cc(2),aa(2),aa1(2),aa2(2),ang1(2),ang2(2),ang3(2), &
              xmn,ymn,zmn, xsiz,ysiz,zsiz, extve(1), &
              x(4),y(4),z(4),vr(4), ve(4), xdb(9),ydb(9),zdb(9), &
              unbias, cbb, bv(9), xa(3),ya(3),za(3),vra(3), vea(3), &
              UNEST, tmin, tmax, skmean , &
              xlj(5), ylj(5), zlj(5), vrlj(5), extvj(5), close(4)
    

    !variogram
    nst(1) = 2  ! note that dimension = 2 is because nst =2
    c0 = 0.25
    cc(1:2) = (/ 0.5, 0.25 /)
    aa(1:2) = (/ 20, 200 /)
    aa1(1:2) = (/ 15, 150 /)
    aa2(1:2) = (/ 10, 100 /)
    ang1(:)= 45 
    ang2(:)= 20
    ang3(:)= 0
    it(:) = 1 
    !search
    radius = 20
    radius1 = 15
    radius2 = 10
    sang1 = 45
    sang2 = 20
    sang3 = 0
    ndmax = 3
    ndmin = 2
    !grid
    nx = 1
    ny = 1
    nz = 1
    xmn =0
    ymn =0
    zmn =0
    xsiz = 10
    ysiz = 10
    zsiz = 5
    extve(:)=1
    
    !jacknife (This is different, here we pass the array)
    ndjak = 5
    xlj = (/ 1 , 2 , 3 ,  1, 2 /)
    ylj = (/ 1 , 1 , 2 ,  2, 3 /)
    zlj = (/ 0 , 0 , 0 ,  0, 0 /)
    vrlj = (/ 0.2 , 2.1 , 3. ,  1.8, 2.5 /)
    extvj = (/ 1 , 1 , 1 ,  1, 1 /)
    
    !data
    nd = 4
    x  = (/ -5 , -5 , 5 ,  5 /)
    y  = (/ -5 ,  5 , -5 , 5 /)
    z(:)  = 0
    vr  = (/ 0.2 , 1.5 , 5.21 , 2.8 /)
    ve(:) = 1  ! this is the external drift, use 1 if no ext drift...
    close(:) = -99999 
    tmin =  99999 
    tmax = -99999 
    ! kriging parameters (see declaration of xdb,ydb,zdb )
    ktype =0
    koption= 0
    nxdis = 3
    nydis = 3
    nzdis = 1
    UNEST = -99999    
    noct = 10
    itrend = 0  ! estimate the drift ? 
    skmean =0

    
    call   kt3d(radius,radius1,radius2,sang1,sang2,sang3, &  
                nst,c0,it,cc,aa,aa1,aa2,ang1,ang2,ang3,   &
                nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, extve,  &
                nd,x,y,z,vr, ve, &
                nclose, close, infoct, noct, &
                MAXSBX, MAXSBY, MAXSBZ, tmin, tmax,  &
                ktype, skmean,  koption, itrend, nxdis,nydis,nzdis,xdb,ydb,zdb, &
                bv,ndmax, ndmin, xa,ya,za,vra, vea, na, &
                ndjak, xlj, ylj, zlj, vrlj, extvj, &
                UNEST)


    ! xdb,ydb,zdb are output. 
    do i=1,9
        print *, xdb(i),ydb(i),zdb(i)
    end do
    
    
    call block_covariance(xdb,ydb,zdb, 9, &
                            nst,it,c0,cc,aa, aa1, aa2,ang1,ang2,ang3, &
                            unbias, cbb)
                            
    print *, 'unbias', unbias
    print *, 'cbb', cbb
    
    print *, 'index of samples collected' , int(close(1:nclose))
    print *, 'samples retained X' , xa (1:na)
    print *, 'samples retained Y' , ya (1:na)
    print *, 'samples retained Z' , za (1:na)
    

    
end 
