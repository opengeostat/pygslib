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

!-----------------------------------------------------------------------
! functions to interpolate angles at drillhole, ex, in assay from survey 
!----------------------------------------------------------------------- 

subroutine ang2cart(azm,dip,x,y,z)
    ! ------------------------------------------------------------------
    ! get azimuth and dip (downward positive) and convert it to x,y,z
    ! angles are in degrees
    ! x,y,z are vectors with origin 0,0,0
    ! for example: [0,x]
    ! ------------------------------------------------------------------
    
    implicit none
    
    ! inout 
    real, intent(in) :: azm,dip
    real, intent(out) :: x,y,z
    
    ! internal
    real :: razm,rdip, DEG2RAD
    
    DEG2RAD=3.141592654/180.0
    
    ! convert degree to rad and correct sign of dip
    razm = azm * DEG2RAD
    rdip = -dip * DEG2RAD

    ! do the conversion
    x = sin(razm) * cos(rdip)
    y = cos(razm) * cos(rdip)
    z = sin(rdip)
    
    return
    
end subroutine ang2cart

subroutine cart2ang(x,y,z,azm,dip)
    ! ------------------------------------------------------------------
    ! convert x,y,z to azimuth, dip (downward positive) 
    ! angles are in degrees
    ! x,y,z are assumed vectors with origin 0,0,0
    ! for example: [0,x]
    ! ------------------------------------------------------------------
    
    implicit none
    
    ! inout 
    real, intent(out) :: azm,dip
    real, intent(in) :: x,y,z
    
    ! internal
    real :: razm,rdip, RAD2DEG, pi
    
    RAD2DEG=180.0/3.141592654
    pi = 3.141592654

    if (x/=0. .and. y/= 0.) then 
        azm= atan2(y,x)
        if (azm<0) azm= azm + pi*2
        azm = azm * RAD2DEG
    else
        azm = 0
    end if 
    
    return 

end subroutine cart2ang


subroutine interp_ang1D(azm1,dip1,azm2,dip2,len12,d1, azm,dip)
    ! ------------------------------------------------------------------
    ! Interpolate the azimuth and dip angle over a line:
    !   given two points (p1, p2) over a line (1D problem);
    !   this subroutine calculate the average azimuth and dip of a point 
    !   between p1 and p2, located at a distance d1 from p1 one and a 
    !   distance len12-d1 from p2
    !   
    !   to do this we convert the (azimuth,dip) to (x,y,z), we 
    !   interpolate x,y,z and then we convert back to (azimuth,dip)
    ! 
    ! ------------------------------------------------------------------
    
    implicit none
    
    ! inout 
    real, intent(out) :: azm,dip
    real, intent(in) :: azm1,dip1,azm2,dip2,len12,d1
    
    ! internal
    real :: x1,y1,z1,x2,y2,z2,x,y,z
    
    
    ! convert angles to coordinates
    call ang2cart(azm1,dip1,x1,y1,z1)
    call ang2cart(azm2,dip2,x2,y2,z2)
    
    ! interpolate x,y,z
    x = x1*d1/len12 + x2*(len12-d1)/len12 
    y = y1*d1/len12 + y2*(len12-d1)/len12
    z = z1*d1/len12 + z2*(len12-d1)/len12
    
    ! get back the results as angles
    call cart2ang(x,y,z,azm,dip)
    
    return
    
end subroutine interp_ang1D

!-----------------------------------------------------------------------
! functions to put in assay x,y,z from collar
!----------------------------------------------------------------------- 


!-----------------------------------------------------------------------
! functions to put in assay az, dip from survey
!----------------------------------------------------------------------- 


!-----------------------------------------------------------------------
! functions to desurvey assay tables (assuming you have x,y,z,az,dip)
!----------------------------------------------------------------------- 

!subroutine ds_mincurb(nc,cid,cx,cy,cz,ns,sid,saz,sdip,sleng,na,aid,afrom,ato, &
!                    ax,ay,az, error)
!    ! using formulas in http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf
    
!    implicit none
    
!    ! input
!    !  collar 
!    integer, intent(in) :: nc
!    integer, intent(in), dimension(nc) ::  cid
!    real, intent(in), dimension(nc) :: cx,cy,cz
!    !  survey 
!    integer, intent(in) :: ns
!    integer, intent(in), dimension(ns) ::  sid
!    real, intent(in), dimension(ns) :: saz,sdip,sleng !len is length from collar
!    !  assay
!    integer, intent(in) :: na
!    integer, intent(in), dimension(na) ::  aid
!    real, intent(in), dimension(na) :: afrom,ato
    
!    ! output
!    !  assay coordinates
!    real, intent(out), dimension(na) :: ax,ay,az, sx,sy,sz   

!    ! internal 
!    integer :: i
    
    
!    ! first desurvey at survey table (information at len 0 required)
!    ! using formulas in http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf
!    sx(:)=0
!    sy(:)=0
!    sz(:)=0
    


!end subroutine ds_mincurb


!subroutine ds_tang(nc,cid,cx,cy,cz,ns,sid,saz,sdip,sleng,na,aid,afrom,ato, &
!                    ax,ay,az)
!    ! using formulas in http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf
    
!    implicit none
    
!    ! input
!    !  collar 
!    integer, intent(in) :: nc
!    integer, intent(in), dimension(nc) ::  cid
!    real, intent(in), dimension(nc) :: cx,cy,cz
!    !  survey 
!    integer, intent(in) :: ns
!    integer, intent(in), dimension(ns) ::  sid
!    real, intent(in), dimension(ns) :: saz,sdip,sleng !len is length from collar
!    !  assay
!    integer, intent(in) :: na
!    integer, intent(in), dimension(na) ::  aid
!    real, intent(in), dimension(na) :: afrom,ato
    
!    ! output
!    !  assay coordinates
!    real, intent(out), dimension(na) :: ax,ay,az, sx,sy,sz   
!    integer, intent(out) :: error

!    ! internal 
!    integer :: i, nsurv, icoll
!    real :: x,y,z
    
    
    
!    if (ns<2) then
!        error=20
!        return
!    end if 

!    ! using formulas in http://www.cgg.com/data//1/rec_docs/2269_MinimumCurvatureWellPaths.pdf
    
!    ! get collar
!    isrv = 0
!    do ic=1, nc
!        nsurv=1
!        sx(1)=0.
!        sy(1)=0.
!        sz(1)=0.
!        do i=2, ns
!            if (sid(i)==sid(i-1)) then 
!                nsurv=nsurv+1
!                sz(i)=
!            end if 
!        end do
!    end do
    
    
!end subroutine ds_tang
