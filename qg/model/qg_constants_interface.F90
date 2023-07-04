! (C) Copyright 2021 UCAR.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module qg_constants_interface

use kinds
use iso_c_binding
use qg_constants_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Get domain parameters
subroutine qg_domain_parameters_c(c_lat_min, c_lat_max, c_domain_zonal, c_domain_meridional, c_xmin, c_ymin, c_lat_proj,  &
 & c_domain_depth) bind(c,name='qg_domain_parameters_f90')

implicit none

! Passed variables
real(c_double),intent(out) :: c_lat_min           !< Minimum latitude [degrees]
real(c_double),intent(out) :: c_lat_max           !< Maximum latitude, chosen to get a 1/3 ratio between domain horizontal dimensions [degrees]
real(c_double),intent(out) :: c_domain_zonal      !< Zonal domain [m]
real(c_double),intent(out) :: c_domain_meridional !< Meridional domain [m]
real(c_double),intent(out) :: c_lat_proj          !< Projection latitude [degrees]
real(c_double),intent(out) :: c_xmin              !< Projection xmin [m]
real(c_double),intent(out) :: c_ymin              !< Projection ymin [m]
real(c_double),intent(out) :: c_domain_depth      !< Domain depth [m]

! Return constant values
c_lat_min = lat_min
c_lat_max = lat_max
c_domain_zonal = domain_zonal
c_domain_meridional = domain_meridional
c_lat_proj = lat_proj
c_xmin = xmin
c_ymin = ymin
c_domain_depth = domain_depth

end subroutine qg_domain_parameters_c
! ------------------------------------------------------------------------------
end module qg_constants_interface
