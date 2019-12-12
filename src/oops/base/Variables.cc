/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/base/Variables.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/types/Types.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace oops {
// -----------------------------------------------------------------------------

Variables::Variables()
  : convention_(""), vars_(0) {
  Log::trace() << "Variables::Variables" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & conf)
  : convention_(""), vars_(0) {
  Log::trace() << "Variables::Variables start " << conf << std::endl;
  std::vector<std::string> vars;
  conf.get("variables", vars);
  // hack to read channels
  if (conf.has("channels")) {
    std::string chlist = conf.getString("channels");
    std::set<int> channels = parseIntSet(chlist);
    std::copy(channels.begin(), channels.end(), std::back_inserter(channels_));
    // assuming the same channel subsetting applies to all variables
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      for (size_t jch = 0; jch < channels_.size(); ++jch) {
        vars_.push_back(vars[jvar]+"_"+std::to_string(channels_[jch]));
      }
    }
  } else {
    vars_ = vars;
  }
  Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const std::vector<std::string> & vars, const std::string & conv)
  : convention_(conv), vars_(vars) {
  Log::trace() << "Variables::Variables start " << vars << std::endl;
  Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const Variables & other)
  : convention_(other.convention_), vars_(other.vars_)
{}

// -----------------------------------------------------------------------------

Variables & Variables::operator+=(const Variables & rhs) {
  ASSERT(convention_ == rhs.convention_);
  vars_.insert(vars_.end(), rhs.vars_.begin(), rhs.vars_.end());
  // remove duplicated variables
  std::sort(vars_.begin(), vars_.end());
  vars_.erase(std::unique(vars_.begin(), vars_.end() ), vars_.end());
  return *this;
}

// -----------------------------------------------------------------------------

bool Variables::operator==(const Variables & rhs) const {
  return convention_ == rhs.convention_
    && vars_ == rhs.vars_
    && channels_ == rhs.channels_;
}

// -----------------------------------------------------------------------------

bool Variables::has(const std::string & var) const {
  bool found = false;
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    found = found || vars_[jj] == var;
  }
  return found;
}

// -----------------------------------------------------------------------------

size_t Variables::find(const std::string & var) const {
  size_t ii = vars_.size();
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (vars_[jj] == var) ii = jj;
  }
  ASSERT(ii < vars_.size());
  return ii;
}

// -----------------------------------------------------------------------------

void Variables::push_back(const std::string & vname) {
  vars_.push_back(vname);
}

// -----------------------------------------------------------------------------

Variables::~Variables() {}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  os << vars_.size() << " variables: ";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << vars_[jj];
  }
  if (!convention_.empty()) os << " (" << convention_ << ")";
}

// -----------------------------------------------------------------------------

}  // namespace oops
