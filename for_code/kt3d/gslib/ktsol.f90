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
