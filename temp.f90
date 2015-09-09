

!***********************************************************************
! 
! 
!     Subroutines for GSLIB gam/gamv (variograms)
! 
! 
!***********************************************************************

subroutine kt3d(  &
MAXSAM, MAXDIS, MAXDT, MAXDAT, MAXSBX, MAXSBY, MAXSBZ, &
UNEST, &
nd,x,y,z,vr,ve,tmin,tmax,nx,ny,nz,xmn,ymn,zmn, &
xsiz,ysiz,zsiz,ndmax,ndmin,radius,noct,nxdis, & 
nydis,nzdis,idrif,itrend,ktype,skmean,koption, &
idbg, & !iout,lout,iext,lext,
iextve,ljack, &
ixlj,iylj,izlj,ivrlj,iextvj,nvarij, &
sang1,sang2,sang3,sanis1,sanis2,xa,ya, &
za,vra,vea,xdb,ydb,zdb,bv, &
iktype,ncut,cut,cdf, &
nst,it,c0,cc,aa,ang1,ang2,ang3, anis1,anis2, rotmat, &
cbb,a,r,rr,s)
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
 
    !input
    integer, intent(in) :: MAXSAM ,MAXDIS, MAXDT,&
                           MAXSBX, MAXSBY, MAXSBZ, MAXDAT, nst, nd, &
                           iktype,ncut, nx,ny,nz, ktype, itrend, &
                           nxdis, nydis, nzdis, noct, ndmax,ndmin, nvarij, &
                           koption, idbg
                           
    real, intent(in) :: radius , c0(1), UNEST, tmin,tmax, xmn,ymn,zmn, xsiz,ysiz,zsiz
    integer, intent(in), dimension(nst) :: it
    real, intent(in), dimension(nst) :: cc,aa,ang1,ang2,ang3, anis1,anis2
    real, intent(in), dimension(nd) :: x, y, z, vr, ve
    integer, intent(in), dimension(MAXDT) :: idrif
    integer, intent(in), dimension(MAXSBX*MAXSBY*MAXSBZ) :: nisb
    integer, intent(in), dimension(8*MAXSBX*MAXSBY*MAXSBZ) :: ixsbtosr, iysbtosr, izsbtosr
    
    
    !inout
    real*8, intent(inout), dimension(nst+1,3,3) :: rotmat
    real*8, intent(inout), dimension(nd) :: tmp, close
    real, intent(inout), dimension(ncut) :: cdf
        
    !output
    real*8, intent(out) ::  cbb, bv(9)
    real, intent(out), dimension(MAXSAM)  ::  xa,ya,za,vra, vea
    real*8, intent (out), dimension(MAXSAM+MAXDT+2) :: r,rr,s
    real, intent (out), dimension((MAXSAM+MAXDT+2)*(MAXSAM+MAXDT+2)) :: a
    real, intent (out), dimension(MAXDIS) :: xdb,ydb,zdb
    
    
    !internal
    integer :: MAXROT,  MXSXY, MXSX, MAXEQ, ndb, na, isrot
    real ::       var(20), PMX, covmax, EPSLON
    logical ::    first,fircon,accept
    real*8 :: unbias
    
    MAXROT=nst+1
    MXSXY=4*MAXSBX*MAXSBY
    MXSX=2*MAXSBX
    EPSLON=0.000001
    PMX    = 999.0
    fircon= .TRUE.
    MAXEQ=MAXSAM+MAXDT+2
 
    
    
    ! TODO: define this as optional parameters?
    ! MAXSBX =    21 , MAXSBY =  21, MAXSBZ =  21,
    ! MAXDAT = 10000, MAXDT  =   9, MAXSAM =  64,
    ! MAXDIS =    64, MAXNST =   4, 
    ! UNEST=-999.0, 

    ! remove this? (not using indicators here )
    ! MAXCUT =  50, cut(MAXCUT),cdf(MAXCUT)


    ! Set up the rotation/anisotropy matrices that are needed for the
    ! variogram and search.  Also compute the maximum covariance for
    ! the rescaling factor:

    write(*,*) 'Setting up rotation matrices for variogram and search'
    radsqd = radius * radius ! the search first radius
    PMX    = 999.0
    covmax = c0(1)
    do is=1,nst
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,MAXROT,rotmat)
        if(it(is) == 4) then
            covmax = covmax + PMX
        else
            covmax = covmax + cc(is)
        endif
    end do
    isrot = nst + 1
    call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)

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
        stop
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
    if(ndb > MAXDIS) then
        write(*,*) 'ERROR KT3D: Too many discretization points',ndb
        write(*,*) '            Increase MAXDIS or lower n[xyz]dis'
        stop
    endif
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

    ! Initialize accumulators:

    nk    = 0
    xk    = 0.0
    vk    = 0.0
    xkmae = 0.0
    xkmse = 0.0

    ! Calculate Block Covariance. Check for point kriging.

    call cova3(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst,nst, &
    c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)

    ! Set the 'unbias' variable so that the matrix solution is more stable

    unbias = cov
    cbb    = dble(cov)
    if(ndb > 1) then
        cbb = 0.0
        do i=1,ndb
            do j=1,ndb
                call cova3(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j), &
                &                1,nst,nst,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                if(i == j) cov = cov - c0(1)
                cbb = cbb + dble(cov)
            end do
        end do
        cbb = cbb/dble(real(ndb*ndb))
    end if
!~     if(idbg > 1) then
!~         write(ldbg,*) ' '
!~         write(ldbg,*) 'Block Covariance: ',cbb
!~         write(ldbg,*) ' '
!~     end if

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

    if(koption == 0) then
        nxy   = nx*ny
        nxyz  = nx*ny*nz
        nloop = nxyz
        irepo = max(1,min((nxyz/10),10000))
    else
        nloop = 10000000
        irepo = max(1,min((nd/10),10000))
    end if
    write(*,*)
    write(*,*) 'Working on the kriging '

    ! MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:

    do index=1,nloop
        if((int(index/irepo)*irepo) == index) write(*,103) index
        103 format('   currently on estimate ',i9)
    
        ! Where are we making an estimate?
    
        if(koption == 0) then
            iz   = int((index-1)/nxy) + 1
            iy   = int((index-(iz-1)*nxy-1)/nx) + 1
            ix   = index - (iz-1)*nxy - (iy-1)*nx
            xloc = xmn + real(ix-1)*xsiz
            yloc = ymn + real(iy-1)*ysiz
            zloc = zmn + real(iz-1)*zsiz
        else
            read(ljack,*,err=96,end=2) (var(i),i=1,nvarij)
            xloc = xmn
            yloc = ymn
            zloc = zmn
            true = UNEST
            secj = UNEST
            if(ixlj > 0)   xloc   = var(ixlj)
            if(iylj > 0)   yloc   = var(iylj)
            if(izlj > 0)   zloc   = var(izlj)
            if(ivrlj > 0)  true   = var(ivrlj)
            if(iextvj > 0) extest = var(iextvj)
        end if

    
        ! Read in the external drift variable for this grid node if needed:
    
        if(ktype == 2 .OR. ktype == 3) then
            if(koption == 0) then
!~                 read(lext,*) (var(i),i=1,iextve)
                extest = var(iextve)
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
!~                 write(ldbg,999)
                fircon = .FALSE. 
            end if
            est  = UNEST
            estv = UNEST
            go to 1
        end if
        999 format(' Encountered a location where there were too few data ',/, &
        ' to estimate all of the drift terms but there would be',/, &
        ' enough data for OK or SK.   KT3D currently leaves ',/, &
        ' these locations unestimated.',/, &
        ' This message is only written once - the first time.',/)
    
        ! There are enough samples - proceed with estimation.
    
        if(na <= 1) then
        
            ! Handle the situation of only one sample:
        
            call cova3(xa(1),ya(1),za(1),xa(1),ya(1),za(1),1,nst,nst, &
            c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb1)
        
            ! Establish Right Hand Side Covariance:
        
            if(ndb <= 1) then
                call cova3(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),1, &
                nst,nst,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
            else
                cb  = 0.0
                do i=1,ndb
                    call cova3(xa(1),ya(1),za(1),xdb(i),ydb(i), &
                    zdb(i),1,nst,nst,c0,it,cc,aa,1, &
                    MAXROT,rotmat,cmax,cov)
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
                call cova3(xa(i),ya(i),za(i),xa(j),ya(j),za(j),1,nst,nst, &
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
                nst,nst,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
            else
                cb  = 0.0
                do j=1,ndb
                    call cova3(xa(i),ya(i),za(i),xdb(j),ydb(j), &
                    zdb(j),1,nst,nst,c0,it,cc,aa,1, &
                    MAXROT,rotmat,cmax,cov)
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
    
        if(idbg == 3) then
!~             write(ldbg,*) 'Estimating node index : ',ix,iy,iz
            is = 1 - neq
            do i=1,neq
                is = 1 + (i-1)*neq
                ie = is + neq - 1
!~                 write(ldbg,100) i,r(i),(a(j),j=is,ie)
                100 format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
            end do
        endif
    
        ! Solve the kriging system:
    
        call ktsol(neq,nrhs,nv,a,r,s,ising,maxeq)
    
        ! Compute the solution:
    
        if(ising /= 0) then
!~             if(idbg >= 3) write(ldbg,*) ' Singular Matrix ',ix,iy,iz
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
        
            if(idbg >= 2) then
!~                 write(ldbg,*) '       '
!~                 write(ldbg,*) 'BLOCK: ',ix,iy,iz,' at ',xloc,yloc,zloc
!~                 write(ldbg,*) '       '
!~                 if(ktype /= 0) &
!~                 write(ldbg,*) '  Lagrange : ',s(na+1)*unbias
!~                 write(ldbg,*) '  BLOCK EST: x,y,z,vr,wt '
                do i=1,na
                    xa(i) = xa(i) + xloc - 0.5*xsiz
                    ya(i) = ya(i) + yloc - 0.5*ysiz
                    za(i) = za(i) + zloc - 0.5*zsiz
!~                     write(ldbg,'(5f12.3)') xa(i),ya(i),za(i), &
!~                     vra(i),s(i)
                end do
!~                 write(ldbg,*) '  estimate, variance  ',est,estv
            endif
        endif
    
        ! END OF MAIN KRIGING LOOP:
    
        1 continue
        if(iktype == 0) then
            if(koption == 0) then
!~                 write(lout,'(f9.3,1x,f9.3)') est,estv
            else
                err = UNEST
!~                 if(true /= UNEST .AND. est /= UNEST)err=est-true
!~                 write(lout,'(7(f12.3,1x))') xloc,yloc,zloc,true, &
!~                 est,estv,err
                xkmae = xkmae + abs(err)
                xkmse = xkmse + err*err
            end if
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
            if(koption == 0) then
!~                 write(lout,'(30(f8.4))') (cdf(i),i=1,ncut)
            else
!~                 write(lout,'(30(f8.4))') (cdf(i),i=1,ncut),true
            end if
        end if
    end do
    2 continue
    if(koption > 0) close(ljack)

    ! Write statistics of kriged values:

     
    if(nk > 0 .AND. idbg > 0) then
        xk    = xk/real(nk)
        vk    = vk/real(nk) - xk*xk
        xkmae = xkmae/real(nk)
        xkmse = xkmse/real(nk)
!~         write(ldbg,105) nk,xk,vk
        write(*,   105) nk,xk,vk
        105 format(/,'Estimated   ',i8,' blocks ',/, &
        '  average   ',f9.4,/,'  variance  ',f9.4,/)
        if(koption /= 0) then
            write(*,106) xkmae,xkmse
            106 format(/,'  mean error',f9.4,/,'  mean sqd e',f9.4)
        end if
    endif

    ! All finished the kriging:

    return
    96 stop 'ERROR in jackknife file!'
end subroutine kt3d
