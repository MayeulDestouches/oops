! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_geom_mod

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use kinds
use iso_c_binding
use qg_constants_mod

implicit none

private
public :: qg_geom
public :: qg_geom_registry
public :: qg_geom_setup,qg_geom_model_setup,qg_geom_fill_atlas_fieldset,qg_geom_clone,qg_geom_delete,qg_geom_info

real(kind_real),parameter :: tol = 1.0e-2_kind_real !< x/y bounds tolerance

! ------------------------------------------------------------------------------
! TEMPORARY (should be removed once an updated version of ATLAS is used)
type :: projection
  type(fckit_configuration) :: conf
  real(kind_real) :: k_radius
  real(kind_real) :: inv_k_radius
contains
  procedure :: setup
  procedure :: copy
  procedure :: xy2lonlat
  procedure :: lonlat2xy
end type projection
! END TEMPORARY

type :: qg_geom
  ! Geometry parameters
  integer :: nx                                      !< Number of points in the zonal direction
  integer :: ny                                      !< Number of points in the meridional direction
  integer :: nz                                      !< Number of vertical levels
  real(kind_real) :: deltax                          !< Zonal cell size [m]
  real(kind_real) :: deltay                          !< Meridional cell size [m]
  real(kind_real),allocatable :: depths(:)           !< Depths of each layer [m]
  real(kind_real),allocatable :: x(:)                !< Zonal coordinate [m]
  real(kind_real),allocatable :: y(:)                !< Meridional coordinate [m]
  real(kind_real),allocatable :: z(:)                !< Altitude of each layer middle [m]
  real(kind_real),allocatable :: lat(:,:)            !< Latitude [rad]
  real(kind_real),allocatable :: lon(:,:)            !< Longitude [rad]
  real(kind_real),allocatable :: area(:,:)           !< Area [m^2]
  real(kind_real),allocatable :: ph_coeff            !< Perturbed Heating Coefficient
! TEMPORARY
  type(projection) :: aproj                          !< ATLAS projection
! ELSE
!  type(atlas_projection) :: aproj                    !< ATLAS projection TEMPORARY
! END TEMPORARY
  type(atlas_structuredgrid) :: agrid                !< ATLAS grid
  type(atlas_functionspace_structuredcolumns) :: afs !< ATLAS function space

  ! Interpolator
  character(len=1024) :: interpolator                !< Interpolator for changes of resolution

  ! Model-related geometry parameters
  real(kind_real),allocatable :: f(:,:)              !< Coefficients of PV operator
  real(kind_real),allocatable :: f_p(:,:)            !< Coefficients of PV operator, right eigenvectors
  real(kind_real),allocatable :: f_pinv(:,:)         !< Coefficients of PV operator, right eigenvectors inverse
  real(kind_real),allocatable :: f_d(:)              !< Coefficients of PV operator, eigenvalues
  real(kind_real),allocatable :: beta0y(:)           !< Beta plane linear term (beta0*y)
  real(kind_real),allocatable :: heat(:,:)           !< Heating term
end type qg_geom

#define LISTED_TYPE qg_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_geom_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! TEMPORARY
!> Projection setup
subroutine setup(proj)

implicit none

class(projection),intent(inout) :: proj

real(kind_real) :: lat1

! Get lat1
call proj%conf%get_or_die("latitude1", lat1)

! Projection radius
proj%k_radius = 6371229.0_kind_real*cos(lat1*deg_to_rad)

! Inverse projection radius
proj%inv_k_radius = 1.0_kind_real/proj%k_radius

end subroutine setup
! ------------------------------------------------------------------------------
!> Projection copy
subroutine copy(proj,other_proj)

implicit none

class(projection),intent(inout) :: proj
type(projection),intent(in) :: other_proj

! Copy
proj%conf = other_proj%conf
proj%k_radius = other_proj%k_radius
proj%inv_k_radius = other_proj%inv_k_radius

end subroutine copy
! ------------------------------------------------------------------------------
!> x/y to lon/lat
subroutine xy2lonlat(proj,x,y,lon,lat)

implicit none

! Passed variables
class(projection),intent(in) :: proj
real(kind_real),intent(in) :: x
real(kind_real),intent(in) :: y
real(kind_real),intent(out) :: lon
real(kind_real),intent(out) :: lat

lon = x*proj%inv_k_radius*rad_to_deg
lat = 90.0-2.0*atan(exp(-y*proj%inv_k_radius))*rad_to_deg

end subroutine xy2lonlat
! ------------------------------------------------------------------------------
!> lon/lat to x/y
subroutine lonlat2xy(proj,lon,lat,x,y)

implicit none

! Passed variables
class(projection),intent(in) :: proj
real(kind_real),intent(in) :: lon
real(kind_real),intent(in) :: lat
real(kind_real),intent(out) :: x
real(kind_real),intent(out) :: y

! Local variables
real(kind_real) :: sinlat

x = proj%k_radius*deg_to_rad*lon
sinlat = sin(lat*deg_to_rad)
y = proj%k_radius*0.5_kind_real*log((1.0_kind_real+sinlat)/(1.0_kind_real-sinlat))

end subroutine lonlat2xy
! END TEMPORARY
! ------------------------------------------------------------------------------
!> Setup geometry
subroutine qg_geom_setup(self,f_conf)

implicit none

! Passed variables
type(qg_geom),intent(inout) :: self            !< Geometry
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ix,iy,iz
real(kind_real) :: xy(2),lonlat(2),deltax,deltay
real(kind_real) :: lonmin,lonmax,latmin,latmax
real(kind_real),allocatable :: real_array(:)
real(kind_real) :: pcoeff
character(len=1024) :: record
character(len=:),allocatable :: str

! TEMPORARY
! Get projection parameters
call f_conf%get_or_die("projection", self%aproj%conf)
call self%aproj%setup
! END TEMPORARY

! Get geometry parameters
call f_conf%get_or_die("nx",self%nx)
call f_conf%get_or_die("ny",self%ny)
if (f_conf%has("depths")) then
  self%nz = f_conf%get_size("depths")
else
  self%nz = 1
endif

! Allocation
allocate(self%depths(self%nz))
allocate(self%x(self%nx))
allocate(self%y(self%ny))
allocate(self%z(self%nz))
allocate(self%lon(self%nx,self%ny))
allocate(self%lat(self%nx,self%ny))
allocate(self%area(self%nx,self%ny))
allocate(self%ph_coeff)

! Set perturbed heating coeff
call f_conf%get_or_die("perturbed heating",pcoeff)
if (pcoeff/=1.0) then
  call fckit_log%info('qg_geom_setup: Perturbed Heating ON')
else
  call fckit_log%info('qg_geom_setup: Perturbed Heating OFF')
end if
self%ph_coeff = pcoeff

! Get depths
if (f_conf%has("depths")) then
  call f_conf%get_or_die("depths",real_array)
  self%depths = real_array
  if (abs(sum(self%depths)-domain_depth)>0.0_kind_real) call abor1_ftn('sum of depths should add up to domain_depth')
else
  self%depths(1) = domain_depth
endif

! Define x/y grid properties
call f_conf%get_or_die("dx",self%deltax)
call f_conf%get_or_die("dy",self%deltay)
do ix=1,self%nx
  xy = self%agrid%xy(ix,1)
  self%x(ix) = xy(1)
enddo
do iy=1,self%ny
  xy = self%agrid%xy(1,iy)
  self%y(iy) = xy(2)
enddo
! Check x/y grid constants and effective bounds
if (abs(xmin-self%x(1))/abs(xmin)>tol) call abor1_ftn('wrong xmin')
if (abs(xmax-(self%x(self%nx)+1.0_kind_real*self%deltax))/abs(xmax)>tol) call abor1_ftn('wrong xmax')
if (abs(ymin-(self%y(1)-self%deltay))/abs(ymin)>tol) call abor1_ftn('wrong ymin')
if (abs(ymax-(self%y(self%ny)+self%deltay))/abs(ymax)>tol) call abor1_ftn('wrong ymax')

do iy=1,self%ny
  do ix=1,self%nx
    lonlat = self%agrid%lonlat(ix,iy)
    self%lon(ix,iy) = lonlat(1)
    self%lat(ix,iy) = lonlat(2)
  end do
end do
do iy=1,self%ny
  deltax = (self%lon(2,iy)-self%lon(1,iy))*deg_to_rad*req
  if (iy==1) then
     deltay = self%lat(1,2)-self%lat(1,1)
  elseif (iy==self%ny) then
     deltay = self%lat(1,self%ny)-self%lat(1,self%ny-1)
  else
     deltay = 0.5_kind_real*(self%lat(1,iy+1)-self%lat(1,iy-1))
  end if
  deltay = deltay*deg_to_rad*req
  self%area(:,iy) = deltax*deltay
end do

! Set heights
self%z(1) = 0.5_kind_real*self%depths(1)
do iz=2,self%nz
  self%z(iz) = sum(self%depths(1:iz-1))+0.5_kind_real*self%depths(iz)
end do

! Set interpolator for changes of variables
call f_conf%get_or_die("interpolator",str)
self%interpolator = str
call fckit_log%info('qg_geom_setup: interpolator: '//trim(self%interpolator))

! Print sizes and dx/dy
write(record,'(a,i3,a,i3,a,i3)') 'qg_geom_create: nx/ny/nz = ',self%nx,' /',self%ny,' /',self%nz
call fckit_log%info(record)
write(record,'(a,f7.2,a,f7.2,a)') '                deltax/deltay = ',self%deltax*1.0e-3_kind_real,' km / ', &
 & self%deltay*1.0e-3_kind_real,' km'
call fckit_log%info(record)

end subroutine qg_geom_setup
! ------------------------------------------------------------------------------
!> Setup model-related geometry
subroutine qg_geom_model_setup(self,f_conf)

implicit none

! Passed variables
type(qg_geom),intent(inout) :: self            !< Geometry
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ix,iy,iz,ix_c,iy_c,lwork,info
integer,allocatable :: ipiv(:),ipivsave(:)
real(kind_real) :: lonlat(2),distx,disty,f
real(kind_real),allocatable :: wi(:),vl(:,:),work(:)
real(kind_real),allocatable :: fsave(:,:),vrlu(:,:),vrlusave(:,:)
real(kind_real) :: norm
logical :: htype

! Allocation
allocate(self%f(self%nz,self%nz))
allocate(self%f_p(self%nz,self%nz))
allocate(self%f_pinv(self%nz,self%nz))
allocate(self%f_d(self%nz))
allocate(self%beta0y(self%ny))
allocate(self%heat(self%nx,self%ny))
allocate(wi(self%nz))
allocate(vl(self%nz,self%nz))
allocate(fsave(self%nz,self%nz))
allocate(vrlu(self%nz,self%nz))
allocate(vrlusave(self%nz,self%nz))
allocate(ipiv(self%nz))
allocate(ipivsave(self%nz))

! Define geophysical parameters
lonlat = 0.5_kind_real*(self%agrid%lonlat(1,1)+self%agrid%lonlat(self%nx,self%ny))*deg_to_rad

! Coefficients of PV operator
self%f = 0.0_kind_real
do iz=1,self%nz
  f = f0**2/(g*dlogtheta*self%depths(iz))
  if (iz>1) then
    self%f(iz,iz-1) = f
    self%f(iz,iz) = self%f(iz,iz)-f
  end if
  if (iz<self%nz) then
    self%f(iz,iz+1) = f
    self%f(iz,iz) = self%f(iz,iz)-f
  end if
enddo

! Compute eigendecomposition of ff
norm = maxval(abs(self%f))
fsave = self%f/norm
allocate(work(1))
call dgeev('V','V',self%nz,fsave,self%nz,self%f_d,wi,vl,self%nz,self%f_p,self%nz,work,-1,info)
if (info/=0) call abor1_ftn('error in dgeev, first pass')
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))
fsave = self%f/norm
call dgeev('V','V',self%nz,fsave,self%nz,self%f_d,wi,vl,self%nz,self%f_p,self%nz,work,lwork,info)
if (info/=0) call abor1_ftn('error in dgeev, second pass')
deallocate(work)

! Compute inverse of right eigenvectors of ff
vrlu = self%f_p
call dgetrf(self%nz,self%nz,vrlu,self%nz,ipiv,info)
if (info/=0) call abor1_ftn('error in dgetrf')
allocate(work(1))
vrlusave = vrlu
ipivsave = ipiv
call dgetri(self%nz,vrlusave,self%nz,ipivsave,work,-1,info)
if (info/=0) call abor1_ftn('error in dgetri, first pass')
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))
self%f_pinv = vrlu
ipivsave = ipiv
call dgetri(self%nz,self%f_pinv,self%nz,ipivsave,work,lwork,info)
if (info/=0) call abor1_ftn('error in dgetri, second pass')
deallocate(work)

! Re-apply normalization after eigendecomposition
self%f_d = self%f_d*norm
self%f_p = self%f_p*norm
self%f_pinv = self%f_pinv/norm

! Beta plane linear term (beta0*y)
do iy=1,self%ny
  self%beta0y(iy) = beta0*(self%lat(1,iy)*deg_to_rad-lat0)*req
enddo

! Set heating term
call f_conf%get_or_die("heating",htype)
if (.not. htype) then
  ! No heating term
  call fckit_log%info('qg_geom_setup: heating off')
  self%heat = 0.0_kind_real
else
  call fckit_log%info('qg_geom_setup: Gaussian heating on')
  ! Gaussian source
  ix_c = self%nx/4
  iy_c = 3*self%ny/4
  do iy=1,self%ny
    do ix=1,self%nx
      distx = abs(self%x(ix)-self%x(ix_c))
      if (distx>0.5_kind_real*domain_zonal) distx = domain_zonal-distx
      disty = abs(self%y(iy)-self%y(iy_c))
      self%heat(ix,iy) = heating_amplitude*exp(-(distx**2+disty**2)/heating_scale**2)
    enddo
  enddo
endif

end subroutine qg_geom_model_setup
! ------------------------------------------------------------------------------
!> Fill ATLAS fieldset
subroutine qg_geom_fill_atlas_fieldset(self,afieldset)

implicit none

! Passed variables
type(qg_geom),intent(inout) :: self             !< Geometry
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: ix,iy,iz,inode
real(kind_real),pointer :: real_ptr_1(:),real_ptr_2(:,:)
type(atlas_field) :: afield

! Add area
afield = self%afs%create_field(name='area',kind=atlas_real(kind_real),levels=0)
call afield%data(real_ptr_1)
inode = 0
do iy=1,self%ny
  do ix=1,self%nx
    inode = inode+1
    real_ptr_1(inode) = self%area(ix,iy)
  enddo
enddo
call afieldset%add(afield)
call afield%final()

! Add vertical unit
afield = self%afs%create_field(name='vunit',kind=atlas_real(kind_real),levels=self%nz)
call afield%data(real_ptr_2)
do iz=1,self%nz
  real_ptr_2(iz,1:self%nx*self%ny) = self%z(iz)
end do
call afieldset%add(afield)
call afield%final()

end subroutine qg_geom_fill_atlas_fieldset
! ------------------------------------------------------------------------------
!> Clone geometry
subroutine qg_geom_clone(self,other)

implicit none

! Passed variables
type(qg_geom),intent(inout) :: self !< Geometry
type(qg_geom),intent(in) :: other   !< Other geometry

! Copy dimensions
self%nx = other%nx
self%ny = other%ny
self%nz = other%nz

! Allocation
allocate(self%depths(self%nz))
allocate(self%x(self%nx))
allocate(self%y(self%ny))
allocate(self%z(self%nz))
allocate(self%lon(self%nx,self%ny))
allocate(self%lat(self%nx,self%ny))
allocate(self%area(self%nx,self%ny))
allocate(self%f(self%nz,self%nz))
allocate(self%f_p(self%nz,self%nz))
allocate(self%f_pinv(self%nz,self%nz))
allocate(self%f_d(self%nz))
allocate(self%beta0y(self%ny))
allocate(self%heat(self%nx,self%ny))
allocate(self%ph_coeff)

! Copy data
self%deltax = other%deltax
self%deltay = other%deltay
self%depths = other%depths
self%x = other%x
self%y = other%y
self%z = other%z
self%lon = other%lon
self%lat = other%lat
self%area = other%area
self%agrid = atlas_structuredgrid(other%agrid%c_ptr())
! TEMPORARY
call self%aproj%copy(other%aproj)
! ELSE
!self%aproj = atlas_projection(other%aproj%c_ptr())
! END TEMPORARY
self%afs = atlas_functionspace_structuredcolumns(other%afs%c_ptr())
self%interpolator = other%interpolator
self%f = other%f
self%f_p = other%f_p
self%f_pinv = other%f_pinv
self%f_d = other%f_d
self%beta0y = other%beta0y
self%heat = other%heat
self%ph_coeff = other%ph_coeff

end subroutine qg_geom_clone
! ------------------------------------------------------------------------------
!> Delete geometry
subroutine qg_geom_delete(self)

implicit none

! Passed variables
type(qg_geom),intent(inout) :: self !< Geometry

! Release memory
deallocate(self%depths)
deallocate(self%x)
deallocate(self%y)
deallocate(self%z)
deallocate(self%lon)
deallocate(self%lat)
deallocate(self%area)
call self%agrid%final()
! TEMPORARY
! ELSE
!call self%aproj%final()
! END TEMPORARY
call self%afs%final()
deallocate(self%f)
deallocate(self%f_p)
deallocate(self%f_pinv)
deallocate(self%f_d)
deallocate(self%beta0y)
deallocate(self%heat)
deallocate(self%ph_coeff)

end subroutine qg_geom_delete
! ------------------------------------------------------------------------------
!> Get geometry info
subroutine qg_geom_info(self,nx,ny,nz,deltax,deltay)

implicit none

! Passed variables
type(qg_geom),intent(in) :: self      !< Geometry
integer,intent(out) :: nx             !< Number of points in the zonal direction
integer,intent(out) :: ny             !< Number of points in the meridional direction
integer,intent(out) :: nz             !< Number of vertical levels
real(kind_real),intent(out) :: deltax !< Zonal cell size
real(kind_real),intent(out) :: deltay !< Meridional cell size

! Copy data
nx = self%nx
ny = self%ny
nz = self%nz
deltax = self%deltax
deltay = self%deltay

end subroutine qg_geom_info
! ------------------------------------------------------------------------------
end module qg_geom_mod
