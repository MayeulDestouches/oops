! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_constants_mod

use kinds

implicit none

private
public :: pi,deg_to_rad,rad_to_deg
public :: req,omega,g
public :: lat0,f0,beta0
public :: dlogtheta,heating_amplitude,heating_scale
public :: lat_min,lat_max,domain_zonal,domain_meridional,lat_proj,xmin,xmax,ymin,ymax,domain_depth

! Mathematical parameters
real(kind_real),parameter :: pi = acos(-1.0_kind_real)                        !< Pi
real(kind_real),parameter :: deg_to_rad = pi/180.0_kind_real                  !< Degrees to radians
real(kind_real),parameter :: rad_to_deg = 180.0_kind_real/pi                  !< Radians to degrees

! Geophysical parameters
real(kind_real),parameter :: req = 6371229.0_kind_real                        !< Earth radius at equator [m]
real(kind_real),parameter :: omega = 7.2921e-5_kind_real                      !< Rotation rate of the Earth [rad.s^{-1}]
real(kind_real),parameter :: g = 9.81_kind_real                               !< Gravity [m^2.s^{-2}]

! Beta plane approximation
real(kind_real),parameter :: lat0 = 50.0_kind_real*deg_to_rad                 !< Beta plane center [radians]
real(kind_real),parameter :: f0 = 2.0_kind_real*omega*sin(lat0)               !< Coriolis parameter at the beta plane center [rad.s^{-1}]
real(kind_real),parameter :: beta0 = 2.0_kind_real*omega*cos(lat0)/req        !< Meridional gradient of f [s^{-1} m^{-1}]

! Model dynamics parameters
real(kind_real),parameter :: dlogtheta = 0.1_kind_real                        !< Difference in log(pot. T) across boundary
real(kind_real),parameter :: heating_amplitude = 5.0e-5_kind_real             !< Heating term amplitude [s^{-1}]
real(kind_real),parameter :: heating_scale = 1.0e6_kind_real                  !< Heating term scale [m]

! Independent domain parameters
real(kind_real),parameter :: lat_min = 20.0_kind_real                         !< Minimum latitude [degrees]
real(kind_real),parameter :: lat_max = 80.14351_kind_real                     !< Maximum latitude, chosen to get a 1/3 ratio between domain horizontal dimensions [degrees]
real(kind_real),parameter :: lat_proj = 43.0_kind_real                        !< Projection latitude [degrees]
real(kind_real),parameter :: domain_depth = 1.0e4_kind_real                   !< Domain depth [m]

! Deduced domain parameters
real(kind_real),parameter :: domain_zonal = 29277267.940857869_kind_real      !< Zonal domain [m]
real(kind_real),parameter :: domain_meridional = 9759089.3553319573_kind_real !< Meridional domain [m]
real(kind_real),parameter :: xmin = -14638633.970428934_kind_real             !< Projection xmin [m]
real(kind_real),parameter :: xmax = 14638633.970428934_kind_real              !< Projection xmax [m]
real(kind_real),parameter :: ymin = 1660589.0899409982_kind_real              !< Projection ymin [m]
real(kind_real),parameter :: ymax = 11419678.445272956_kind_real              !< Projection ymax [m]

end module qg_constants_mod
