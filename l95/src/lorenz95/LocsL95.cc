/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/LocsL95.h"

#include "eckit/config/Configuration.h"
#include "lorenz95/ObsTable.h"
#include "oops/util/DateTime.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------

LocsL95::LocsL95(const std::vector<int> & indx, const std::vector<double> & locs)
  : indx_(indx), locs_(locs)
{
  ASSERT(indx_.size() == locs_.size());
}

// -----------------------------------------------------------------------------

LocsL95::LocsL95(const eckit::Configuration & conf) : indx_(), locs_() {
  conf.get("positions", locs_);
  for (size_t jj = 0; jj < locs_.size(); ++jj) {
    ASSERT(locs_.at(jj) >= 0.0 && locs_.at(jj) <= 1.0);
    indx_.push_back(jj + 1);
  }
}

// -----------------------------------------------------------------------------

void LocsL95::print(std::ostream & os) const {
  os << locs_.size();
  if (locs_.size() > 0) os << " " << locs_.at(0);
  if (locs_.size() > 1) os << " " << locs_.at(locs_.size() - 1);
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
