!  kt3d.f90
!  
!  Copyright 2016 Adrian Martinez Vargas <adrian.martinez@opengeostat.com>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
!  
!  
! Version ### 5

! compile with 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  gfortran -W -g -fbacktrace -ffpe-trap=zero,overflow,underflow kt3d.f90 gslib/*.f90 -o kt3d
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




! ======================================================================
! Module with global variables
! Note it will affect all python objects in a python section
! Set variable values with set_varname functions
!
! Note:  linked automatically to python but do not write it directly,
!        use set_varname instead
!
! ====================================================================== 
Module Commons

    ! User Adjustable Parameters:

    integer :: MAXSBX, MAXSBY, MAXSBZ, MAXDT, MAXSAM, MAXCUT

    real :: UNEST                                                       ! asign this externally at init , otherwise will be zero by default

    parameter(MAXSBX =    21, MAXSBY =  21, MAXSBZ =  21, &
    MAXDT  =   9, MAXSAM =  64, MAXCUT=3)

    ! Fixed Parameters:

    parameter(MAXSB=MAXSBX*MAXSBY*MAXSBZ,MXSXY=4*MAXSBX*MAXSBY, &
    MXSX=2*MAXSBX,EPSLON=0.000001, &
    VERSION=2.000)


End Module Commons          

! ======================================================================
!
! Functions to set Global "System" variable values 
!
! Note:  link these functions manually in setup.py
!
! ====================================================================== 
subroutine set_unest(unest_)

    use Commons

    real, intent(in) :: unest_

    UNEST = unest_

end subroutine set_unest


! ======================================================================
!
! Python wrapper, this is what we link to python 
!
! Note:  link these functions manually in setup.py
!
! ====================================================================== 
subroutine pykt3d( nd, x,y,z,vr,ve, &                                   ! input data 
           nx,ny,nz,xmn,ymn,zmn, xsiz,ysiz,zsiz, nxdis,nydis,nzdis,&    ! block model definition (including  discretization)
           radius,radius1,radius2, ndmax,ndmin,noct,sang1,sang2,sang3, & ! search parameters
           idrif,&                                                      ! drift terms
           itrend,ktype,skmean,koption, iktype,ncut,cut, &              ! kriging options
           nst,c0,it,cc,aa,aa1,aa2,ang1,ang2,ang3,  &                   ! variogram parameters                                                       
           nout, outx, outy, outz, outextve, outest, outkvar, &         ! output variable  
           outcdf, &                                                    ! output in case of cdf
           idbg, cbb, neq, na, dbgxdat,dbgydat,dbgzdat, &                         ! debug level and data
           dbgvrdat,dbgwt,dbgxtg,dbgytg,dbgztg, dbgkvector, dbgkmatrix, & 
           errors, warns)                                               ! Error output for python




    use commons


    ! ================
    ! input variables
    ! ================       
    integer, intent (in) :: nd
    real, intent(in), dimension(nd) :: x                                ! only x is compulsory, if y,z,vr, ve not provided will be initialized to zero                      
    real, intent(in), dimension(nd), optional :: y,z,vr,ve 

    
    
    integer, intent(in) :: nx,ny,nz                                     ! block model definition \nx,xmn,xsiz ...
    real, intent(in) ::    xmn,ymn,zmn, xsiz,ysiz,zsiz 
    integer, intent(in) :: nxdis,nydis,nzdis                            ! block discretization

    real, intent(in) :: radius,radius1,radius2                          ! search parameters
    real, intent(in), optional  :: sang1,sang2,sang3       
    integer, intent(in) :: ndmax,ndmin
    integer, intent(in), optional ::noct

    integer, intent(in), dimension (MAXDT), optional  :: idrif          ! drift terms (external)
     
    integer, intent(in), optional :: itrend,ktype,koption               ! kriging options
    real, intent(in), optional :: skmean

    integer, intent(in), optional ::  iktype                       ! especial parameters for indicator krigimg  iktype =1       
    integer, intent(in) :: ncut
    real, intent(in), dimension(ncut) ::cut
    
    integer, intent(in) :: nst                                          ! Variogram parameters
    real, intent(in) :: c0
    integer, intent(in), dimension(nst):: it
    real, intent(in), dimension(nst):: cc,aa,aa1,aa2
    real, intent(in), dimension(nst), optional:: ang1,ang2,ang3
    
    integer, intent(in), optional :: idbg                                
    
    ! ==================
    ! output variables
    ! ==================
    integer, intent (in) :: nout
    real, intent(in), dimension(nout), optional :: outy, outz, outextve 
    real, intent(in), dimension(nout) :: outx    
    real, intent(out), dimension(nout) :: outest, outkvar                ! output variable with the estimate, kvar and an indicator of success for a given block 
    real, intent(out), dimension(nout,ncut) :: outcdf
    
    ! debug
    real*8, intent(out) :: cbb
    integer, intent(out) ::neq ,na
     real, intent(out), dimension(ndmax):: dbgxdat,dbgydat,dbgzdat,dbgvrdat
    real, intent(out), dimension(ndmax+MAXDT+2):: dbgwt,dbgkvector
    real, intent(out), dimension((ndmax+MAXDT+2),(ndmax+MAXDT+2)):: dbgkmatrix
 
    real, intent (out) :: dbgxtg,dbgytg,dbgztg
 
    character(LEN=250), intent(out)  :: errors, warns
                

               
    
    ! ==================
    ! internal variables
    ! ==================
    real, dimension(nst):: anis1,anis2     ! anisotropy in variogram structures
    real                   :: sanis1, sanis2  ! anisotropy in the search ellipse 
    integer,  dimension (MAXDT) :: idrif_                 ! drift terms (local)
    real :: skmean_


    ! create a local copy of the array
    
    idrif_ = idrif
    skmean_ = skmean

    ! some of the check before in read parameter
    
    if(nxdis < 1) stop 'Error in parameters nxdis,nxdis,nydis < 1'
    
    if(nst <= 0) then
        write(*,9997) nst
        9997 format(' nst must be at least 1, it has been set to ',i4,/, &
        ' The c or a values can be set to zero')
        errors = 'Error in parameters: nst < 1'
        return
    endif   


    ! compute anisotropy on search ellipse 
    
    if(radius < EPSLON) then 
        errors = 'Error in parameters: radius < 0'
        return
    end if 
    sanis1 = radius1 / radius
    sanis2 = radius2 / radius  

    

    !compute anis1 and anis2 and check the variogram 

    do i=1,nst
        anis1(i) = aa1(i) / max(aa(i),EPSLON)
        anis2(i) = aa2(i) / max(aa(i),EPSLON)
        if(it(i) == 4) then
            if(aa(i) < 0.0) then
                 errors = 'Error in parameters: aa(i) < 0.0 invalid power variogram'
                 return 
            end if
            if(aa(i) > 2.0) then
                 errors = 'Error in parameters: aa(i) > 2.0 invalid power variogram'
                 return 
            end if
        end if
    end do    



    ! Compute the averages and variances as an error check for the user:

    ! >>>>>>  TODO: put this optional, it may slowdown a sequential run

    av = 0.0
    ss = 0.0
    do i=1, nd
        av = av + vr(i)
        ss = ss + vr(i)*vr(i)
    end do

    av = av / max(real(nd),1.0)
    ss =(ss / max(real(nd),1.0)) - av * av
    write(*,*) '  Number   = ',nd
    write(*,*) '  Average  = ',av
    write(*,*) '  Variance = ',ss
    
    ! report the input

    write(*,*) ' kriging option = ',koption
    write(*,*) ' iktype option = ',iktype

    if(iktype == 1) then
        write(*,*) ' Doing Median Indicator Kriging Estimate (iktype == 1) '
        write(*,*) ' number of cutoffs = ',ncut
        write(*,*) ' cutoffs = ',(cut(i),i=1,ncut)
    end if

    write(*,*) ' debugging level = ',idbg
    write(*,*) ' nx, xmn, xsiz = ',nx,xmn,xsiz
    write(*,*) ' ny, ymn, ysiz = ',ny,ymn,ysiz
    write(*,*) ' nz, zmn, zsiz = ',nz,zmn,zsiz
    write(*,*) ' block discretization:',nxdis,nydis,nzdis
    write(*,*) ' ndmin,ndmax = ',ndmin,ndmax
    write(*,*) ' max per octant = ',noct
    write(*,*) ' search radii = ',radius,radius1,radius2
    write(*,*) ' search anisotropy = ', sanis1,sanis2, '** computed as radius secondary / radius primary **'
    write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3
    write(*,*) ' ktype, skmean =',ktype,skmean
    write(*,*) ' drift terms = ',(idrif(i),i=1,9)
    write(*,*) ' itrend = ',itrend
    write(*,*) ' nst, c0 = ',nst,c0

    do i=1,nst
        write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i), &
                                ang1(i),ang2(i),ang3(i)
        write(*,*) ' a1 a2 a3 anis1 anis2: ',aa(i),aa1(i),aa2(i), anis1(i), anis2(i)
    end do

    write (*,*) '       ** anis computed as a2/a1  **' 


    ! call kriging

    call   kt3d(   nd, x,y,z,vr,ve, &                                        ! input data 
                   nx,ny,nz,xmn,ymn,zmn, xsiz,ysiz,zsiz, nxdis,nydis,nzdis,& ! block model definition (including  discretization)
                   radius, ndmax,ndmin,noct,sang1,sang2,sang3,sanis1,sanis2, &              ! search parameters
                   idrif,&                                                   ! drift terms
                   itrend,ktype,skmean,koption, iktype,ncut,cut, &           ! kriging options
                   nst,it,c0,cc,aa,ang1,ang2,ang3, anis1,anis2, &            ! variogram parameters
                   idbg, &
                   nout, outx, outy, outz, outextve, outest, outkvar, &      ! output variable with the estimate, kvar and an indicator of success for a given block  
                   outcdf, & 
                   cbb,neq, na, dbgxdat,dbgydat,dbgzdat, &
                   dbgvrdat,dbgwt,dbgxtg,dbgytg,dbgztg,dbgkvector, dbgkmatrix, &                  ! debug
                   errors, warns) 



    return


end subroutine pykt3d



! ======================================================================
!
! Modified GSLIB Kriging functions   
!
! Note:  Do not link these to python, use wraper (pykt3d) instead
!
! ======================================================================
subroutine kt3d(   nd, x,y,z,vr,ve, &                                        ! input data 
                   nx,ny,nz,xmn,ymn,zmn, xsiz,ysiz,zsiz, nxdis,nydis,nzdis,& ! block model definition (including  discretization)
                   radius, ndmax,ndmin,noct,sang1,sang2,sang3,sanis1,sanis2, &              ! search parameters
                   idrif,&                                                   ! drift terms
                   itrend,ktype,skmean,koption, iktype,ncut,cut, &           ! kriging options
                   nst_,it,c0_,cc,aa,ang1,ang2,ang3, anis1,anis2, &            ! variogram parameters
                   idbg, &
                   nout, outx, outy, outz, outextve, outest, outkvar, &                                  ! output variable with the estimate, kvar and an indicator of success for a given block  
                   outcdf, & 
                   cbb,neq, na, dbgxdat,dbgydat,dbgzdat, &
                   dbgvrdat,dbgwt,dbgxtg,dbgytg,dbgztg,dbgkvector, dbgkmatrix, &                  ! debug
                   errors, warns) 

    !-----------------------------------------------------------------------

    !                Krige a 3-D Grid of Rectangular Blocks
    !                **************************************

    ! This subroutine estimates point or block values of one variable by
    ! simple, ordinary, or kriging with a trend model.  It is also possible
    ! to estimate the trend directly.

    ! This is a modified version, adapted to be linked to python
    ! all the output is trough variables (file output removed)
    
    !
    ! Modified by : Adrian Martinez Vargas                             2016
    !
    

    ! Original:  A.G. Journel and C. Lemmer                             1981
    ! Revisions: A.G. Journel and C. Kostov                             1984
    !-----------------------------------------------------------------------               
    
    use commons
    
           
    ! ==================
    ! input variables
    ! ==================              
           
    integer, intent (in) :: nd                               
    real, intent(in), dimension(nd) :: x,y,z,vr,ve                  ! input data
    
    integer, intent(in) :: nx,ny,nz                                 ! block model definition \nx,xmn,xsiz ...
    real, intent(in) ::    xmn,ymn,zmn, xsiz,ysiz,zsiz 
    integer, intent(in) :: nxdis,nydis,nzdis                     ! block discretization
    
    real, intent(in) :: radius, sang1,sang2,sang3,sanis1,sanis2                   ! search radius
    integer, intent(in) :: ndmax,ndmin,noct

    integer, intent(inout), dimension (MAXDT) :: idrif              ! drift terms 

    integer, intent(in) :: itrend,ktype,koption                     ! kriging options
    real, intent(inout) :: skmean
    
    integer, intent(in) ::  iktype,ncut                             ! especial parameters for indicator krigimg  iktype =1       
    real, intent(in), dimension(ncut)::cut


    integer, intent(in) :: nst_                                   ! Variogram parameters
    real, intent(in) :: c0_
    integer, intent(in), dimension(nst_):: it
    real, intent(in), dimension(nst_):: cc,aa,ang1,ang2,ang3, &
                                          anis1,anis2
    
    integer, intent(in) :: idbg  ! TODO: >>>>> remove this    

    
    ! ==================
    ! output variables
    ! ==================
    real, intent(in), dimension(nout) :: outx, outy, outz           ! centroid of the block where we want the estimation
                                                                     ! block size is parent  
                                                          
    real, intent(in), dimension(nout) :: outextve                   ! external drift                                  
                                                                    
    
    ! block covariance
    real*8, intent (out) ::     cbb                                                               
                                                                     
    
    integer, intent (in) :: nout
    real, intent(out), dimension(nout) :: outest, outkvar                                 ! output variable with the estimate, kvar and an indicator of success for a given block 
    real, intent(out), dimension(nout,ncut) :: outcdf
    
    ! debug
    integer, intent(out) ::neq , na
    real, intent(out), dimension(ndmax):: dbgxdat,dbgydat,dbgzdat,dbgvrdat
    real, intent(out), dimension(ndmax+MAXDT+2):: dbgwt,dbgkvector
    real, intent(out), dimension((ndmax+MAXDT+2),(ndmax+MAXDT+2)):: dbgkmatrix
    
    real, intent (out) :: dbgxtg,dbgytg,dbgztg
    
    character(LEN=250), intent(out)  :: errors, warns

    ! =================
    ! local variables
    ! =================
    
    integer :: nst(1)                                   ! Variogram parameters
    real :: c0(1)
    integer :: maxrot
    
    ! array of varianles to read from file >>>>>> remove 
    real ::       var(20)
    logical ::    first,fircon,accept
    data       fircon/ .TRUE. /
    
    ! set of array with the same size of data
    real, dimension(:), allocatable :: tmp, close, xa,ya,za,vra,vea
    integer :: AllocateStatus
    
    ! unbias and drift means
    real :: unbias, bv(MAXDT)
    
    ! cdf of the actual block estimate, if iktype = 1
    real :: cdf(ncut)
    
    ! kriging matrix and rotation matrix
    real*8 ::  rotmat(nst_+1,3,3)
    
    ! block discretization
    real :: xdb(nxdis*nydis*nzdis), &
            ydb(nxdis*nydis*nzdis), &
            zdb(nxdis*nydis*nzdis)

    ! superblock search (consider to use kdtree here)
    integer ::  nisb(MAXSB), &
    ixsbtosr(8*MAXSB),iysbtosr(8*MAXSB),izsbtosr(8*MAXSB)

    ! kriging system
    real*8, dimension(ndmax+MAXDT+2)::r,rr,s
    real*8, dimension((ndmax+MAXDT+2)*(ndmax+MAXDT+2)):: a

    integer :: nloop


    ! alocate local arrays
    
    allocate ( tmp(nd), close(nd),xa(nd),ya(nd),za(nd),vra(nd),vea(nd), STAT = AllocateStatus)
    if (AllocateStatus /= 0) then 
        errors = "Error Internal: There was a problem allocating arrays in memory"
        return 
    end if                                                                        

    c0(1)=c0_
    nst(1)=nst_
    maxrot = nst_+1

    ! =================
    ! computations start here
    ! =================

    ! we define the array where we estimate externally (like in jackknife) 
    nloop = nout


    ! Set up the rotation/anisotropy matrices that are needed for the
    ! variogram and search.  Also compute the maximum covariance for
    ! the rescaling factor:

    write(*,*) 'Setting up rotation matrices for variogram and search'
    radsqd = radius * radius
    PMX    = 999.0
    covmax = c0(1)
    do is=1,nst(1)
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,maxrot,rotmat)
        if(it(is) == 4) then
            covmax = covmax + PMX
        else
            covmax = covmax + cc(is)
        endif
    end do
    isrot = nst_ + 1
    call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,maxrot,rotmat)
    
    
    ! Finish computing the rescaling factor and stop if unacceptable:

    if(radsqd < 1.0) then
        resc = 2.0 * radius / max(covmax,0.0001)
    else
        resc =(4.0 * radsqd)/ max(covmax,0.0001)
    endif
! >>> handle error    
    if(resc <= 0.0) then
        write(*,*) 'ERROR KT3D: The rescaling value is wrong ',resc
        write(*,*) '            Maximum covariance: ',covmax
        write(*,*) '            search radius:      ',radius
        errors = "Error Internal: The rescaling value is wrong"
        return 
    endif
    resc = 1.0 / resc

    ! Set up for super block searching:                               ! here use kdtree 

    write(*,*) 'Setting up super block search strategy'
    nsec = 1
    call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z, &
    vr,tmp,nsec,ve,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,nisb, &
    nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup, &
    zmnsup,zsizsup)
    call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup, &
    isrot,maxrot,rotmat,radsqd,nsbtosr,ixsbtosr, &
    iysbtosr,izsbtosr)

    ! Compute the number of drift terms, if an external drift is being
    ! considered then it is one more drift term, if SK is being considered
    ! then we will set all the drift terms off and mdt to 0):

    mdt = 1
    do i=1,9
        if(ktype == 0 .OR. ktype == 2) idrif(i) = 0
        if(idrif(i) < 0 .OR. idrif(i) > 1) then
            write(*,*) 'ERROR KT3D: invalid drift term',idrif(i)
            errors = "Error Internal: Invalid drift term (idrif(i) < 0 .OR. idrif(i) > 1)"
            return 
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


    ndb = nxdis * nydis * nzdis
    xdis = xsiz  / max(real(nxdis),1.0)
    ydis = ysiz  / max(real(nydis),1.0)
    zdis = zsiz  / max(real(nzdis),1.0)
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

    ! output discretization 
    
    if (idbg > 0) then 
        write (*,*) ' for block with coordinate at ', 0,0,0
        write (*,*) '                and side size ', xsiz,xsiz,xsiz
        write (*,*) '  The discretization coordinates are : '
        write (*,*) ' xdb ', xdb
        write (*,*) ' ydb ', ydb
        write (*,*) ' zdb ', zdb
    end if


    ! Initialize accumulators:

    nk    = 0
    xk    = 0.0
    vk    = 0.0
    xkmae = 0.0
    xkmse = 0.0

    ! Calculate Block Covariance. Check for point kriging.

    call cova3(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst,nst_, &
    c0,it,cc,aa,1,maxrot,rotmat,cmax,cov)
    
    if (idbg>0)  write (*,*) ' The point covariance (cov) is ', cov 

    ! Set the ``unbias'' variable so that the matrix solution is more stable

    unbias = cov
    cbb    = dble(cov)
    if(ndb > 1) then
        cbb = 0.0
        do i=1,ndb
            do j=1,ndb
                call cova3(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j), &
                &                1,nst,nst_,c0,it,cc,aa,1,maxrot,rotmat,cmax,cov)
                if(i == j) cov = cov - c0(1)
                cbb = cbb + dble(cov)
            end do
        end do
        cbb = cbb/dble(real(ndb*ndb))
    end if
    
    if (idbg>0)  write (*,*) ' The block covariance (cbb) is ', cbb


    ! Mean values of the drift functions:

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
        bv(i) = (bv(i) / real(ndb)) * resc
    end do

    ! Report on progress from time to time:

    irepo = max(1,min((nloop/10),10000))

    write(*,*)
    write(*,*) 'Working on the kriging with nloop:  ', nloop

    ! MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:

    do index=1,nloop
        if((int(index/irepo)*irepo) == index) &
           write(*,*) '   currently on estimate [ind,x,y,z]', &
           index, outx(index), outy(index), outz(index)
            
        ! Where are we making an estimate?
             

        xloc = outx(index)
        yloc = outy(index)
        zloc = outz(index)
        true = UNEST
        secj = UNEST
        

    
        ! Read in the external drift variable for this grid node if needed:
    
        if(ktype == 2 .OR. ktype == 3) then
            extest = outextve(index)                                    ! warn if doing it in sparse data 
            resce  = covmax / max(extest,0.0001)
        endif
    
        ! Find the nearest samples:
    
        call srchsupr(xloc,yloc,zloc,radsqd,isrot,maxrot,rotmat,nsbtosr, &
        ixsbtosr,iysbtosr,izsbtosr,noct,nd,x,y,z,tmp, &
        nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup, &
        nzsup,zmnsup,zsizsup,nclose,close,infoct)
    
        ! Load the nearest data in xa,ya,za,vra,vea:
    
        na = 0
        do i=1,nclose
            ind    = int(close(i)+0.5)
            accept = .TRUE. 
            if(koption /= 0 .AND. &
            (abs(x(ind)-xloc)+abs(y(ind)-yloc)+ abs(z(ind)-zloc)) & ! this is what excludes the point in cross val and jackknife 
             < EPSLON) accept = .FALSE. 
            if(accept) then
                if(na < ndmax) then
                    na = na + 1
                    xa(na)  = x(ind) - xloc + 0.5*xsiz
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
            if(fircon) then
                warns = trim(warns) // ' Encountered a location where there were too few data to estimate all of the drift terms'
                fircon = .FALSE. 
            end if
            est  = UNEST
            estv = UNEST
            if (idbg>0)  write (*,*) ' The block ', index, 'not estimated due to lack of data, na:', na
            go to 1
        end if

    
        ! There are enough samples - proceed with estimation.
    
        if(na <= 1) then
        
            ! Handle the situation of only one sample:
        
            call cova3(xa(1),ya(1),za(1),xa(1),ya(1),za(1),1,nst,nst_, &
            c0,it,cc,aa,1,maxrot,rotmat,cmax,cb1)
        
            ! Establish Right Hand Side Covariance:
        
            if(ndb <= 1) then
                call cova3(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),1, &
                nst,MAXNST,c0,it,cc,aa,1,maxrot,rotmat,cmax,cb)
            else
                cb  = 0.0
                do i=1,ndb
                    call cova3(xa(1),ya(1),za(1),xdb(i),ydb(i), &
                    zdb(i),1,nst,nst_,c0,it,cc,aa,1, &
                    maxrot,rotmat,cmax,cov)
                    cb = cb + cov
                    dx = xa(1) - xdb(i)
                    dy = ya(1) - ydb(i)
                    dz = za(1) - zdb(i)
                    if((dx*dx+dy*dy+dz*dz) < EPSLON) cb=cb-c0(1)
                end do
                cb = cb / real(ndb)
            end if
            est  = vra(1)
            estv = real(cbb) - 2.0*cb + cb1
            nk   = nk + 1
            xk   = xk + vra(1)
            vk   = vk + vra(1)*vra(1)
            go to 1
        end if
    
        ! Go ahead and set up the OK portion of the kriging matrix:
    
        neq = mdt+na
    
        ! Initialize the main kriging matrix:
    
        first = .FALSE. 
        do i=1,neq*neq
            a(i) = 0.0
        end do
    
        ! Fill in the kriging matrix:
    
        do i=1,na
            do j=i,na
                call cova3(xa(i),ya(i),za(i),xa(j),ya(j),za(j),1,nst,nst_, &
                c0,it,cc,aa,1,maxrot,rotmat,cmax,cov)
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
                nst,nst_,c0,it,cc,aa,1,maxrot,rotmat,cmax,cb)
            else
                cb  = 0.0
                do j=1,ndb
                    call cova3(xa(i),ya(i),za(i),xdb(j),ydb(j), &
                    zdb(j),1,nst,nst_,c0,it,cc,aa,1, &
                    maxrot,rotmat,cmax,cov)
                    cb = cb + cov
                    dx = xa(i) - xdb(j)
                    dy = ya(i) - ydb(j)
                    dz = za(i) - zdb(j)
                    if((dx*dx+dy*dy+dz*dz) < EPSLON) cb=cb-c0(1)
                end do
                cb = cb / real(ndb)
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
    
        if(idbg>0) then
            write(*,*) 'Estimating node index : ', index
            is = 1 - neq
            do i=1,neq
                is = 1 + (i-1)*neq
                ie = is + neq - 1
                write(*,100) i,r(i),(a(j),j=is,ie)
                100 format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
                
                ! output kriging matrix of the las estimate
                dbgkvector(i) = r(i)
                dbgkmatrix(i,1:neq)= a(is:ie)
                
            end do
        endif
    
        ! Solve the kriging system:
    
        call ktsol(neq,nrhs,nv,a,r,s,ising,ndmax+MAXDT+2)
    
        ! Compute the solution:
    
        if(ising /= 0) then
            if(idbg >0) write(*,*) ' Singular Matrix ', index
            est  = UNEST
            estv = UNEST
        else
            est  = 0.0
            estv = real(cbb)
            if(ktype == 2) skmean = extest
            do j=1,neq
                estv = estv - real(s(j))*rr(j)
                if(j <= na) then
                    if(ktype == 0 .OR. ktype == 2) then
                        est = est + real(s(j))*(vra(j)-skmean)
                    else
                        est = est + real(s(j))*vra(j)
                    endif
                endif
            end do
            if(ktype == 0 .OR. ktype == 2) est = est + skmean
            nk   = nk + 1
            xk   = xk + est
            vk   = vk + est*est
        
            ! Write the kriging weights and data if debugging level is above 2:
        
            if(idbg > 0) then
                write(*,*) '       '
                write(*,*) 'BLOCK: ',index,' at ',xloc,yloc,zloc
                write(*,*) '       '
                if(ktype /= 0) &
                write(*,*) '  Lagrange : ',s(na+1)*unbias
                write(*,*) '  BLOCK EST: x,y,z,vr,wt '
                do i=1,na
                    xa(i) = xa(i) + xloc - 0.5*xsiz
                    ya(i) = ya(i) + yloc - 0.5*ysiz
                    za(i) = za(i) + zloc - 0.5*zsiz
                    
                    ! output sample location and wt of the last estimate
                    dbgxdat(i) = xa(i)
                    dbgydat(i) = ya(i)
                    dbgzdat(i) = za(i)
                    dbgvrdat(i) = vra(i)
                    dbgwt(i) = s(i)
                    
                    dbgxtg = xloc 
                    dbgytg = yloc
                    dbgztg = zloc
                    
                    write(*,'(5f12.3)') xa(i),ya(i),za(i), vra(i),s(i)
                end do
                do i=na+1,neq
                    dbgwt(i)=s(i)*unbias                                ! (note that here we report all lagranges, rescaled back to covariace scale)
                end do
                write(*,*) '  estimate, variance  ',est,estv
            endif
        endif
    
        ! END OF MAIN KRIGING LOOP:
    
        1 continue
        if(iktype == 0) then

                outest(index)=est
                outkvar(index)=estv

        else
        
            ! Work out the IK-type distribution implicit to this data configuration
            ! and kriging weights:
        
            do icut=1,ncut
                cdf(icut) = -1.0
            end do
            wtmin = 1.0
            do i=1,na
                if(s(i) < wtmin) wtmin = s(i)
            end do
            sumwt = 0.0
            do i=1,na
                s(i)  = s(i) - wtmin
                sumwt = sumwt + s(i)
            end do
            do i=1,na
                s(i) = s(i) / max(0.00001,sumwt)
            end do
            if(na > 1 .AND. sumwt > 0.00001) then
                do icut=1,ncut
                    cdf(icut) = 0.0
                    do i=1,na
                        if(vra(i) <= cut(icut)) &
                        cdf(icut)=cdf(icut)+s(i)
                    end do
                end do
            end if
            
            ! write the cdf in a 2D array
            do icut=1,ncut
                outcdf(index,icut)=cdf(icut)
            end do
            
        end if
    end do



    ! Write statistics of kriged values:

     
    if(nk > 0 .AND. idbg > 0) then
        xk    = xk/real(nk)
        vk    = vk/real(nk) - xk*xk
        xkmae = xkmae/real(nk)
        xkmse = xkmse/real(nk)
        ! write(ldbg,105) nk,xk,vk
        write(*,   105) nk,xk,vk
        105 format(/,'Estimated   ',i8,' blocks ',/, &
        '  average   ',f9.4,/,'  variance  ',f9.4,/)
    endif

    ! All finished the kriging:

    ! deallocate internal arrays 
    deallocate ( tmp, close,xa,ya,za,vra,vea, STAT = AllocateStatus)
    if (AllocateStatus /= 0) Then
        errors = "Error Internal: There was a problem deallocating arrays from memory"
        return 
    end if

    return

end subroutine kt3d



