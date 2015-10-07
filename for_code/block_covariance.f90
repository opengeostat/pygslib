!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 2015 Adrian Martinez Vargas                            %
!                                                                      %
! This software may be modified and distributed under the terms        %
! of the MIT license.  See the LICENSE.txt file for details.           %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
! The subroutine below isbased on GSLIB code, 
! version 2.0 of the gslib code written in fortran 77
! 
!-----------------------------------------------------------------------


!*********************************************************************************
!     Subroutines in GSLIB (auxiliary functions)
!*********************************************************************************
include 'setrot.f90'
include 'cova3.f90'
include 'sqdist.f90'
! include 'getindx.f90'
! include 'locate.f90'
! this includes sortem and dsortem
! include 'sortem.f90'
! include 'acorni.f90'
! include 'gauinv.f90'
! include 'powint.f90'
! include 'gcum.f90'



subroutine block_covariance(xdb,ydb,zdb, ndb, &
                            nst,it,c0,cc,aa,aa1,aa2,ang1,ang2,ang3, &
                            unbias,cbb)
    !-----------------------------------------------------------------------

    !                    Block Covariance 
    !                    *****************************

    ! This subroutine calculate the block covariance associated with a variogram
    ! model specified by a nugget effect and nested varigoram structures. 
    ! The anisotropy definition can be different for each nested structure.
    ! The block size is defined with input discretization points with 
    ! arbitrary locations (for example regular discretization, random 
    ! discretization or discretization points in an irregular polygon)



    ! INPUT VARIABLES:

    !   xdb,ydb,zdb      coordinates of discretization points
    !   ndb              number of discretization points
    !   nst(ivarg)       number of nested structures (maximum of 4)
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
    !   aa,aa1,aa2       parameters "a" of each nested structure.
    !   aa,aa1,aa2       rotation angles for each nested structure
    ! 
    ! OUTPUT VARIABLES:

    !   unbias           unbias variable, need internally for kriging 
    !   cbb              block covariance



    ! Notes 
    
    ! This functions was created from code lines in kt3d program
    ! 
    !
    ! Adrian Martinez  2015
    !-----------------------------------------------------------------------
    
    IMPLICIT NONE 

    ! in 
    ! variogram
    integer, intent(in)                 :: nst
    real*8, intent(in),   dimension(1)  :: c0
    integer, intent(in), dimension(nst) :: it 
    real*8, intent(in), dimension(nst)  :: cc,aa, aa1, aa2,ang1,ang2,ang3
    !discretization points of the last block
    integer, intent (in) :: ndb   ! number of discretization points
    real*8, intent (in), dimension(ndb) :: xdb,ydb,zdb      !coordinates of discretization points

    !out 
    real*8, intent(out) :: unbias, cbb


    ! internal 
    ! variogram 
    real*8 :: cmax,covmax, cov
    real*8, dimension(nst) :: anis1, anis2
    real*8, dimension(nst,3,3) :: rotmat
    real*8 ::  EPSLON=1.e-20, PMX=999.
    integer :: i, is, j
    
    
    
    do i=1,nst
        anis1(i) = aa1(i) / max(aa(i),EPSLON)
        anis2(i) = aa2(i) / max(aa(i),EPSLON)
        
        ! print *, aa(i), aa1(i), aa2(i), anis1(i), anis2(i)
        
        if(it(i).eq.4) then
              if(aa(i).lt.0.0) stop ' INVALID power variogram'
              if(aa(i).gt.2.0) stop ' INVALID power variogram'
        end if
    end do
    
    
    
    ! get the rotation matrix
    do is=1,nst
        call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
        is,nst,rotmat)
        if(it(is) == 4) then
            covmax = covmax + PMX
        else
            covmax = covmax + cc(is)
        endif
    end do

    ! Calculate Block Covariance. Check for point kriging.
    call cova3(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst, &
                c0,it,cc,aa,1,nst,rotmat,cmax,cov)

    ! Set the 'unbias' variable so that the matrix solution is more stable

    unbias = cov
    cbb    = dble(cov)
    if(ndb > 1) then
        cbb = 0.0
        do i=1,ndb
            do j=1,ndb
                call cova3(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j), &
                           1,nst,c0,it,cc,aa,1,nst,rotmat,cmax,cov)
                if(i == j) cov = cov - c0 (1)
                cbb = cbb + dble(cov)
            end do
        end do
        cbb = cbb/dble(ndb*ndb)
    end if

    return
    
end subroutine block_covariance
