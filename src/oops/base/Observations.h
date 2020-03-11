/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVATIONS_H_
#define OOPS_BASE_OBSERVATIONS_H_

#include <cstddef>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/Departures.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

/// Observations Class.
/*!
 *  Contains observed values or their model equivalents
 */

// -----------------------------------------------------------------------------
template <typename MODEL> class Observations : public util::Printable {
  typedef Departures<MODEL>          Departures_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
/// \brief create Observations for all obs (read from ObsSpace if name is specified)
  explicit Observations(const ObsSpaces_ &, const std::string & name = "");
/// \brief create local Observations
  Observations(const ObsSpaces_ &, const Observations &);

/// destructor and copy/move constructor/assignments
  ~Observations() = default;
  Observations(const Observations &);
  Observations(Observations &&);
  Observations & operator=(const Observations &);
  Observations & operator=(Observations &&);

/// Access
  std::size_t size() const {return obs_.size();}
  ObsVector_ & operator[](const std::size_t ii) {return obs_.at(ii);}
  const ObsVector_ & operator[](const std::size_t ii) const {return obs_.at(ii);}

/// Interactions with Departures
  Departures_ operator-(const Observations & other) const;
  Observations & operator+=(const Departures_ &);

/// Save observations values
  void save(const std::string &) const;

/// Accumulator
  void zero();
  void accumul(const Observations &);
  Observations & operator*=(const double);

 private:
  void print(std::ostream &) const;

/// Data
  const ObsSpaces_ &      obsdb_;
  std::vector<ObsVector_> obs_;
};

// =============================================================================

template <typename MODEL>
Observations<MODEL>::Observations(const ObsSpaces_ & obsdb,
                                  const std::string & name): obsdb_(obsdb), obs_()
{
  obs_.reserve(obsdb.size());
  for (std::size_t jj = 0; jj < obsdb.size(); ++jj) {
    obs_.emplace_back(obsdb[jj], name, true);
  }
  Log::trace() << "Observations created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::Observations(const ObsSpaces_ & obsdb,
                                  const Observations & other): obsdb_(obsdb), obs_() {
  obs_.reserve(obsdb.size());
  for (std::size_t jj = 0; jj < other.size(); ++jj) {
    obs_.emplace_back(obsdb[jj], other[jj]);
  }
  Log::trace() << "Local observations created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::Observations(const Observations & other)
: obsdb_(other.obsdb_), obs_(other.obs_) {}
// -----------------------------------------------------------------------------

template <typename MODEL>
Observations<MODEL>::Observations(Observations && other)
: obsdb_(other.obsdb_), obs_(std::move(other.obs_)) {}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator=(const Observations & other) {
// only allow assignment for Observations created from the same ObsSpaces
  ASSERT(&obsdb_ == &other.obsdb_);
  obs_ = other.obs_;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator=(Observations && other) {
// only allow assignment for Observations created from the same ObsSpaces
  ASSERT(&obsdb_ == &other.obsdb_);
  obs_ = std::move(other.obs_);
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Departures<MODEL> Observations<MODEL>::operator-(const Observations & other) const {
  Departures_ diff(obsdb_);
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    diff[jj]  = obs_[jj];
    diff[jj] -= other[jj];
  }
  return diff;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator+=(const Departures_ & dy) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj] += dy[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::save(const std::string & name) const {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj].save(name);
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::zero() {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj].zero();
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::accumul(const Observations & y) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj] += y[jj];
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator *=(const double factor) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj] *= factor;
  }
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::print(std::ostream & os) const {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) os << obs_[jj] << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_OBSERVATIONS_H_
