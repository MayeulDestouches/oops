! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to ageom%ny jurisdiction.

module qg_convert_x_to_u_mod

use kinds
use qg_geom_mod

implicit none

private
public :: convert_x_to_u,convert_x_to_u_tl,convert_x_to_u_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Convert streafunction to zonal wind
subroutine convert_x_to_u(geom,x,x_north,x_south,u)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                            !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)    !< Streamfunction
real(kind_real),intent(in) :: x_north(geom%nz)              !< Streamfunction on northern wall
real(kind_real),intent(in) :: x_south(geom%nz)              !< Streamfunction on southern wall
real(kind_real),intent(inout) :: u(geom%nx,geom%ny,geom%nz) !< Zonal wind

! Local variables
integer :: iz

!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
  ! South bound
  u(:,1,iz) = (x_south(iz)-x(:,2,iz))/(2.0_kind_real*geom%deltay)

  ! Center
  u(:,2:geom%ny-1,iz) = (x(:,1:geom%ny-2,iz)-x(:,3:geom%ny,iz))/(2.0_kind_real*geom%deltay)

  ! North bound
  u(:,geom%ny,iz) = (x(:,geom%ny-1,iz)-x_north(iz))/(2.0_kind_real*geom%deltay)
enddo
!$omp end parallel do

end subroutine convert_x_to_u
! ------------------------------------------------------------------------------
!> Convert streafunction to zonal wind - tangent Linear
subroutine convert_x_to_u_tl(geom,x,u)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                          !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)  !< Streamfunction
real(kind_real),intent(out) :: u(geom%nx,geom%ny,geom%nz) !< Zonal wind

! Local variables
integer :: iz

!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
  ! South bound
  u(:,1,iz) = -x(:,2,iz)/(2.0_kind_real*geom%deltay)

  ! Center
  u(:,2:geom%ny-1,iz) = (x(:,1:geom%ny-2,iz)-x(:,3:geom%ny,iz))/(2.0_kind_real*geom%deltay)

  ! North bound
  u(:,geom%ny,iz) = x(:,geom%ny-1,iz)/(2.0_kind_real*geom%deltay)
enddo
!$omp end parallel do

end subroutine convert_x_to_u_tl
! ------------------------------------------------------------------------------
!> Convert streafunction to zonal wind - adjoint
subroutine convert_x_to_u_ad(geom,u,x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                             !< Geometry
real(kind_real),intent(in) :: u(geom%nx,geom%ny,geom%nz)     !< Zonal wind
real(kind_real),intent(inout)  :: x(geom%nx,geom%ny,geom%nz) !< Streamfunction

! Local variables
integer :: iz

do iz=1,geom%nz
  ! South bound
  x(:,2,iz) = x(:,2,iz)-u(:,1,iz)/(2.0_kind_real*geom%deltay)

  ! Center
  x(:,1:geom%ny-2,iz) = x(:,1:geom%ny-2,iz)+u(:,2:geom%ny-1,iz)/(2.0_kind_real*geom%deltay)
  x(:,3:geom%ny,iz) = x(:,3:geom%ny,iz)-u(:,2:geom%ny-1,iz)/(2.0_kind_real*geom%deltay)

  ! North bound
  x(:,geom%ny-1,iz) = x(:,geom%ny-1,iz)+u(:,geom%ny,iz)/(2.0_kind_real*geom%deltay)
enddo

end subroutine convert_x_to_u_ad
! ------------------------------------------------------------------------------
end module qg_convert_x_to_u_mod
