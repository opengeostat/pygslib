!-----------------------------------------------------------------------
!   PyGSLIB Desurvey, Module to calculate drillhole coordinates at 
!   interval tables and other drillhole relate process.  
! 
!   Copyright (C) 2015 Adrian Martinez Vargas 
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 3 of the License, or
!   any later version.
!    
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>
!-----------------------------------------------------------------------


! TODO: function to check relations between tables

integer function issort(nc,cid,ns,sid,sleng,na,aid,afrom,ato,tol,gaptol,inttol)

    implicit none

    ! input
    !  collar 
    integer, intent(in) :: nc
    integer, intent(in), dimension(nc) ::  cid
    !  survey 
    integer, intent(in) :: ns
    integer, intent(in), dimension(ns) ::  sid
    real, intent(in), dimension(ns) :: sleng
    !  assay
    integer, intent(in) :: na
    integer, intent(in), dimension(na) ::  aid
    real, intent(in), dimension(na) :: afrom, ato
    real, intent(in) :: tol,gaptol, inttol
    
    
    ! internal
    integer :: i
    real :: EPSLON
    
    EPSLON= 0.0001
    
    issort =0
    
    ! check collar
    if (nc>=2) then
        do i=2,nc
            ! id not sorted 
            if (cid(i)<cid(i-1)) then
                issort = 10 
                return
            end if
        end do
    end if

    ! check survey

    if (ns>=2) then
        do i=2,ns
            ! id not sorted 
            if (sid(i)<sid(i-1)) then
                issort = 20  
                return
            end if
            ! length not sorted 
            if (sid(i)==sid(i-1) .AND. sleng(i)+tol<=sleng(i-1)) then
                issort = 21   
                return
            end if
            ! data at length zero?
            if (sid(i)/=sid(i-1) .AND. sleng(i)>EPSLON) then
                issort = 23   
                return
            end if            
        end do
    end if
    ! data at length zero?
    if (sleng(1)>EPSLON) then
        issort = 23   
        return
    end if


    ! check assay
    if (na>=2) then
        do i=2,na
            ! id sorted ?
            if (aid(i)<aid(i-1)) then
                issort = 30   
                return
            end if
            ! from sorted ?
            if (aid(i)==aid(i-1) .AND. afrom(i)+tol<=afrom(i-1)) then
                issort = 31   
                return
            end if
            ! gaps?
            if (aid(i)==aid(i-1) .AND. ato(i-1)-afrom(i)>gaptol) then
                issort = 32   
                return
            end if
            ! overlap?
            if (aid(i)==aid(i-1) .AND. ato(i-1)-afrom(i)<-gaptol) then
                issort = 33  
                return
            end if
        end do
    end if

    !check assay zero interval or negative interval
    do i=1,na
        ! from > to?
        if (ato(i)-afrom(i)<-tol) then
            issort = 34   
            return
        end if
        ! from - to > inttol?
        if (abs(ato(i)-afrom(i))<inttol) then
            issort = 35   
            return
        end if
    end do

end function issort

subroutine ds_mincurb(nc,cid,cx,cy,cz,ns,sid,saz,sdip,sleng,na,aid,afrom,ato, &
                    ax,ay,az, error)
    ! Subroutine for drillhole desurvey with minimum curvature
    ! here we assume the data is sorted
    
    implicit none
    
    ! input
    !  collar 
    integer, intent(in) :: nc
    integer, intent(in), dimension(nc) ::  cid
    real, intent(in), dimension(nc) :: cx,cy,cz
    !  survey 
    integer, intent(in) :: ns
    integer, intent(in), dimension(ns) ::  sid
    real, intent(in), dimension(ns) :: saz,sdip,sleng !len is length from collar
    !  assay
    integer, intent(in) :: na
    integer, intent(in), dimension(na) ::  aid
    real, intent(in), dimension(na) :: afrom,ato
    
    ! output
    !  assay coordinates
    real, intent(out), dimension(na) :: ax,ay,az, sx,sy,sz   

    ! internal 
    integer :: i
    
    
    ! first desurvey at survey table (information at len 0 required)
    ! using formulas in http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf
    sx(:)=0
    sy(:)=0
    sz(:)=0
    


end subroutine ds_mincurb


subroutine ds_tangential(nc,cid,cx,cy,cz,ns,sid,saz,sdip,sleng,na,aid,afrom,ato, &
                    ax,ay,az)
    ! Subroutine for drillhole desurvey with minimum curvature
    ! here we assume the data is sorted
    
    implicit none
    
    ! input
    !  collar 
    integer, intent(in) :: nc
    integer, intent(in), dimension(nc) ::  cid
    real, intent(in), dimension(nc) :: cx,cy,cz
    !  survey 
    integer, intent(in) :: ns
    integer, intent(in), dimension(ns) ::  sid
    real, intent(in), dimension(ns) :: saz,sdip,sleng !len is length from collar
    !  assay
    integer, intent(in) :: na
    integer, intent(in), dimension(na) ::  aid
    real, intent(in), dimension(na) :: afrom,ato
    
    ! output
    !  assay coordinates
    real, intent(out), dimension(na) :: ax,ay,az, sx,sy,sz   
    integer, intent(out) :: error

    ! internal 
    integer :: i, nsurv, icoll
    real :: x,y,z
    
    
    
    if (ns<2) then
        error=20
        return
    end if 

    ! using formulas in http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf
    
    ! get collar
    isrv = 0
    do ic=1, nc
        nsurv=1
        sx(1)=0.
        sy(1)=0.
        sz(1)=0.
        do i=2, ns
            if (sid(i)==sid(i-1)) then 
                nsurv=nsurv+1
                sz(i)=
            end if 
        end do
    end do
    
Δ TVD = Δ MD cos I1
Δ N = Δ MD sin I1 cos A1
Δ E = Δ MD sin I1 Δ E Δ MD sin A1
    
end subroutine ds_tangential
