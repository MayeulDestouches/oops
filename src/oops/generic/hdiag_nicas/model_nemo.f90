!----------------------------------------------------------------------
! Module: module_nemo.f90
!> Purpose: NEMO model routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module model_nemo

use netcdf
use tools_const, only: req,deg2rad,sphere_dist
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use tools_missing, only: msvalr,msr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_geom, only: geomtype,geom_alloc
use type_nam, only: namtype

implicit none

private
public :: model_nemo_coord,model_nemo_read,model_nemo_write

contains

!----------------------------------------------------------------------
! Subroutine: model_nemo_coord
!> Purpose: get NEMO coordinates
!----------------------------------------------------------------------
subroutine model_nemo_coord(nam,geom)

implicit none

! Passed variables
type(namtype),intent(in) :: nam      !< Namelist
type(geomtype),intent(inout) :: geom !< Geometry

! Local variables
integer :: il0
integer :: ncid,nlon_id,nlat_id,nlev_id,lon_id,lat_id,tmask_id,e1t_id,e2t_id
integer(kind=1),allocatable :: tmask(:,:,:)
real(kind=4),allocatable :: lon(:,:),lat(:,:),e1t(:,:,:),e2t(:,:,:)
character(len=1024) :: subr = 'model_nemo_coord'

! Open file and get dimensions
call ncerr(subr,nf90_open(trim(nam%datadir)//'/grid.nc',nf90_nowrite,ncid))
call ncerr(subr,nf90_inq_dimid(ncid,'x',nlon_id))
call ncerr(subr,nf90_inq_dimid(ncid,'y',nlat_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlon_id,len=geom%nlon))
call ncerr(subr,nf90_inquire_dimension(ncid,nlat_id,len=geom%nlat))
geom%nc0 = geom%nlon*geom%nlat
call ncerr(subr,nf90_inq_dimid(ncid,'z',nlev_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nlev_id,len=geom%nlev))

! Allocation
allocate(lon(geom%nlon,geom%nlat))
allocate(lat(geom%nlon,geom%nlat))
allocate(tmask(geom%nlon,geom%nlat,geom%nl0))
allocate(e1t(geom%nlon,geom%nlat,geom%nl0))
allocate(e2t(geom%nlon,geom%nlat,geom%nl0))

! Read data and close file
call ncerr(subr,nf90_inq_varid(ncid,'nav_lon',lon_id))
call ncerr(subr,nf90_inq_varid(ncid,'nav_lat',lat_id))
call ncerr(subr,nf90_inq_varid(ncid,'tmask',tmask_id))
call ncerr(subr,nf90_inq_varid(ncid,'e1t',e1t_id))
call ncerr(subr,nf90_inq_varid(ncid,'e2t',e2t_id))
call ncerr(subr,nf90_get_var(ncid,lon_id,lon,(/1,1/),(/geom%nlon,geom%nlat/)))
call ncerr(subr,nf90_get_var(ncid,lat_id,lat,(/1,1/),(/geom%nlon,geom%nlat/)))
do il0=1,geom%nl0
   call ncerr(subr,nf90_get_var(ncid,tmask_id,tmask(:,:,il0),(/1,1,nam%levs(il0),1/),(/geom%nlon,geom%nlat,1,1/)))
   call ncerr(subr,nf90_get_var(ncid,e1t_id,e1t(:,:,il0),(/1,1,1/),(/geom%nlon,geom%nlat,1/)))
   call ncerr(subr,nf90_get_var(ncid,e2t_id,e2t(:,:,il0),(/1,1,1/),(/geom%nlon,geom%nlat,1/)))
end do
call ncerr(subr,nf90_close(ncid))

! Convert to radian
lon = lon*real(deg2rad,kind=4)
lat = lat*real(deg2rad,kind=4)

! Pack
call geom_alloc(geom)
geom%lon = pack(real(lon,kind_real),mask=.true.)
geom%lat = pack(real(lat,kind_real),mask=.true.)
do il0=1,geom%nl0
   ! Land/sea mask
   geom%mask(:,il0) = pack(tmask(:,:,il0)>0,mask=.true.)
end do

! Compute normalized area
do il0=1,geom%nl0
   geom%area(il0) = sum(e1t(:,:,il0)*e2t(:,:,il0),mask=tmask(:,:,il0)>0.0)/req**2
end do

! Vertical unit
geom%vunit = float(nam%levs(1:geom%nl0))

! Release memory
deallocate(lon)
deallocate(lat)
deallocate(tmask)

end subroutine model_nemo_coord

!----------------------------------------------------------------------
! Subroutine: model_nemo_read
!> Purpose: read NEMO field
!----------------------------------------------------------------------
subroutine model_nemo_read(nam,geom,ncid,its,fld)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                              !< Namelist
type(geomtype),intent(in) :: geom                            !< Geometry
integer,intent(in) :: ncid                                   !< NetCDF file ID
integer,intent(in) :: its                                    !< Timeslot index
real(kind_real),intent(out) :: fld(geom%nc0,geom%nl0,nam%nv) !< Read field

! Local variables
integer :: iv,il0
integer :: fld_id
real(kind=8) :: fld_tmp(geom%nlon,geom%nlat),fld_loc(geom%nlon,geom%nlat)
character(len=1024) :: subr = 'model_nemo_read'

! Initialize field
call msr(fld)

do iv=1,nam%nv
   ! 3d variable

   ! Get variable id
   call ncerr(subr,nf90_inq_varid(ncid,trim(nam%varname(iv)),fld_id))

   do il0=1,nam%nl
      call ncerr(subr,nf90_get_var(ncid,fld_id,fld_tmp,(/1,1,nam%levs(il0),nam%timeslot(its)/),(/geom%nlon,geom%nlat,1,1/)))
      select case (trim(nam%varname(iv)))
      case ('un')
         fld_loc(:,1) = 0.5*(fld_tmp(geom%nlon,:)+fld_tmp(:,1))
         fld_loc(2:geom%nlon,:) = 0.5*(fld_tmp(1:geom%nlon-1,:)+fld_tmp(2:geom%nlon,:))
      case ('vn')
         fld_loc(:,1) = 0.5*(fld_tmp(:,geom%nlat)+fld_tmp(:,1))
         fld_loc(:,2:geom%nlat) = 0.5*(fld_tmp(:,1:geom%nlat-1)+fld_tmp(:,2:geom%nlat))
      case default
         fld_loc = fld_tmp
      end select
      fld(:,il0,iv) = pack(real(fld_loc,kind_real),mask=.true.)
   end do

   if (trim(nam%addvar2d(iv))/='') then
      ! 2d variable

      ! Get id
      call ncerr(subr,nf90_inq_varid(ncid,trim(nam%addvar2d(iv)),fld_id))

      ! Read data
      call ncerr(subr,nf90_get_var(ncid,fld_id,fld_loc,(/1,1,nam%timeslot(its)/),(/geom%nlon,geom%nlat,1/)))
      fld(:,geom%nl0,iv) = pack(real(fld_loc,kind_real),.true.)
   end if
end do

end subroutine model_nemo_read

!----------------------------------------------------------------------
! Subroutine: model_nemo_write
!> Purpose: write NEMO field
!----------------------------------------------------------------------
subroutine model_nemo_write(geom,ncid,varname,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                    !< Geometry
integer,intent(in) :: ncid                           !< NetCDF file ID
character(len=*),intent(in) :: varname               !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0,geom%nl0) !< Written field

! Local variables
integer :: il0,ic0,ierr
integer :: nlon_id,nlat_id,nlev_id,fld_id
real(kind_real) :: fld_tmp(geom%nc0),fld_loc(geom%nlon,geom%nlat)
logical :: mask_unpack(geom%nlon,geom%nlat)
character(len=1024) :: subr = 'model_nemo_write'

! Initialization
mask_unpack = .true.

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Define dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'x',nlon_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'x',geom%nlon,nlon_id))
   ierr = nf90_inq_dimid(ncid,'y',nlat_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'y',geom%nlat,nlat_id))
   ierr = nf90_inq_dimid(ncid,'z',nlev_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'z',geom%nl0,nlev_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
do il0=1,geom%nl0
   if (isanynotmsr(fld(:,il0))) then
      call msr(fld_tmp)
      do ic0=1,geom%nc0
         if (geom%mask(ic0,il0)) fld_tmp(ic0) = fld(ic0,il0)
      end do
      fld_loc = unpack(fld_tmp,mask=mask_unpack,field=fld_loc)
      call ncerr(subr,nf90_put_var(ncid,fld_id,fld_loc,(/1,1,il0/),(/geom%nlon,geom%nlat,1/)))
   end if
end do

end subroutine model_nemo_write

end module model_nemo
