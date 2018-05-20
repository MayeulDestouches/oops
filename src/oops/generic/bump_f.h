/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_BUMP_F_H_
#define OOPS_GENERIC_BUMP_F_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
extern "C" {
  void create_bump_f90(int &, const eckit::Configuration * const *, const int &, const int &,
                       const int &, const int &, const double *, const double *, const double *,
                       const double *, const int *, const int &, const double *);
  void delete_bump_f90(const int &);
  void bump_multiply_f90(const int &, const int &);
}
}  // namespace oops

#endif  // OOPS_GENERIC_BUMP_F_H_
