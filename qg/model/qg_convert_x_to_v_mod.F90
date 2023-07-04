! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to ageom%ny jurisdiction.

module qg_convert_x_to_v_mod

use kinds
use qg_geom_mod

implicit none

private
public :: convert_x_to_v,convert_x_to_v_tl,convert_x_to_v_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Convert streafunction to meridional wind
subroutine convert_x_to_v(geom,x,v)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                            !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)    !< Streamfunction
real(kind_real),intent(inout) :: v(geom%nx,geom%ny,geom%nz) !< Meridional wind

! Local variables
integer :: iz

!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
   ! West side (perdiodic)
   v(1,:,iz) = (x(2,:,iz)-x(geom%nx,:,iz))/(2.0_kind_real*geom%deltax)

   ! Center
   v(2:geom%nx-1,:,iz) = (x(3:geom%nx,:,iz)-x(1:geom%nx-2,:,iz))/(2.0_kind_real*geom%deltax)

   ! East side (periodic)
   v(geom%nx,:,iz) = (x(1,:,iz)-x(geom%nx-1,:,iz))/(2.0_kind_real*geom%deltax)
enddo
!$omp end parallel do

end subroutine convert_x_to_v
! ------------------------------------------------------------------------------
!> Convert streafunction to meridional wind - tangent Linear
subroutine convert_x_to_v_tl(geom,x,v)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                          !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)  !< Streamfunction
real(kind_real),intent(out) :: v(geom%nx,geom%ny,geom%nz) !< Meridional wind

! Local variables
integer :: iz

!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
   ! West side (perdiodic)
   v(1,:,iz) = (x(2,:,iz)-x(geom%nx,:,iz))/(2.0_kind_real*geom%deltax)

   ! Center
   v(2:geom%nx-1,:,iz) = (x(3:geom%nx,:,iz)-x(1:geom%nx-2,:,iz))/(2.0_kind_real*geom%deltax)

   ! East side (periodic)
   v(geom%nx,:,iz) = (x(1,:,iz)-x(geom%nx-1,:,iz))/(2.0_kind_real*geom%deltax)
enddo
!$omp end parallel do

end subroutine convert_x_to_v_tl
! ------------------------------------------------------------------------------
!> Convert streafunction to meridional wind - adjoint
subroutine convert_x_to_v_ad(geom,v,x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                             !< Geometry
real(kind_real),intent(in) :: v(geom%nx,geom%ny,geom%nz)     !< Meridional wind
real(kind_real),intent(inout)  :: x(geom%nx,geom%ny,geom%nz) !< Streamfunction

! Local variables
integer :: iz

do iz=1,geom%nz
   ! West side (perdiodic)
   x(2,:,iz) = x(2,:,iz)+v(1,:,iz)/(2.0_kind_real*geom%deltax)
   x(geom%nx,:,iz) = x(geom%nx,:,iz)-v(1,:,iz)/(2.0_kind_real*geom%deltax)

   ! Center
   x(3:geom%nx,:,iz) = x(3:geom%nx,:,iz)+v(2:geom%nx-1,:,iz)/(2.0_kind_real*geom%deltax)
   x(1:geom%nx-2,:,iz) = x(1:geom%nx-2,:,iz)-v(2:geom%nx-1,:,iz)/(2.0_kind_real*geom%deltax)

   ! East side (periodic)
   x(1,:,iz) = x(1,:,iz)+v(geom%nx,:,iz)/(2.0_kind_real*geom%deltax)
   x(geom%nx-1,:,iz) = x(geom%nx-1,:,iz)-v(geom%nx,:,iz)/(2.0_kind_real*geom%deltax)
enddo

end subroutine convert_x_to_v_ad
! ------------------------------------------------------------------------------
end module qg_convert_x_to_v_mod
