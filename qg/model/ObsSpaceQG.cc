/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsSpaceQG.h"

#include <map>
#include <string>
#include <utility>

#include "atlas/array.h"
#include "atlas/field.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace qg {
// -----------------------------------------------------------------------------
// initialization for the static map
std::map < std::string, F90odb > ObsSpaceQG::theObsFileRegister_;
int ObsSpaceQG::theObsFileCount_ = 0;

// -----------------------------------------------------------------------------

ObsSpaceQG::ObsSpaceQG(const Parameters_ & params, const eckit::mpi::Comm & comm,
                       const util::DateTime & bgn, const util::DateTime & end,
                       const eckit::mpi::Comm & timeComm)
  : oops::ObsSpaceBase(params, comm, bgn, end), obsname_(params.obsType),
    winbgn_(bgn), winend_(end), obsvars_()
{
  typedef std::map< std::string, F90odb >::iterator otiter;

  eckit::LocalConfiguration fileconf = params.toConfiguration();
  std::string ofin("-");
  if (params.obsdatain.value() != boost::none) {
    ofin = params.obsdatain.value()->obsfile;
  }
  std::string ofout("-");
  if (params.obsdataout.value() != boost::none) {
    ofout = params.obsdataout.value()->obsfile;
    if (timeComm.size() > 1) {
      std::ostringstream ss;
      ss << "_" << timeComm.rank();
      std::size_t found = ofout.find_last_of(".");
      if (found == std::string::npos) found = ofout.length();
      std::string fileout = ofout.insert(found, ss.str());
      fileconf.set("obsdataout.obsfile", fileout);
    }
  }
  oops::Log::trace() << "ObsSpaceQG: Obs files are: " << ofin << " and " << ofout << std::endl;
  std::string ref = ofin + ofout;
  if (ref == "--") {
    ABORT("Underspecified observation files.");
  }

  ref = ref + bgn.toString() + end.toString();
  otiter it = theObsFileRegister_.find(ref);
  if ( it == theObsFileRegister_.end() ) {
    // Open new file
    oops::Log::trace() << "ObsSpaceQG::getHelper: " << "Opening " << ref << std::endl;
    qg_obsdb_setup_f90(key_, fileconf, bgn, end);
    theObsFileRegister_[ref] = key_;
  } else {
    // File already open
    oops::Log::trace() << "ObsSpaceQG::getHelper: " << ref << " already opened." << std::endl;
    key_ = it->second;
  }
  theObsFileCount_++;

  // Set variables simulated for different obstypes
  if (obsname_ == "Stream") obsvars_.push_back("Stream");
  if (obsname_ == "WSpeed") obsvars_.push_back("WSpeed");
  if (obsname_ == "Wind") {
    obsvars_.push_back("Uwind");
    obsvars_.push_back("Vwind");
  }

  //  Generate locations etc... if required
  if (params.generate.value() != boost::none) {
    // Location generation parameters
    const ObsGenerateParameters &gParams = *params.generate.value();
    const util::Duration first(gParams.begin);
    const util::DateTime start(winbgn_ + first);
    const util::Duration freq(gParams.obsPeriod);
    int nobstimes = 0;
    util::DateTime now(start);
    while (now <= winend_) {
      ++nobstimes;
      now += freq;
    }

    // Call fortran
    qg_obsdb_generate_f90(key_, obsname_.size(), obsname_.c_str(),
                          gParams.toConfiguration(), start, freq, nobstimes);
  }
}

// -----------------------------------------------------------------------------

ObsSpaceQG::~ObsSpaceQG() {
  ASSERT(theObsFileCount_ > 0);
  theObsFileCount_--;
  if (theObsFileCount_ == 0) {
    theObsFileRegister_.clear();
    qg_obsdb_delete_f90(key_);
  }
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::save() const {
  qg_obsdb_save_f90(key_);
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::getdb(const std::string & col, int & keyData) const {
  qg_obsdb_get_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), keyData);
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::putdb(const std::string & col, const int & keyData) const {
  qg_obsdb_put_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), keyData);
}

// -----------------------------------------------------------------------------

std::unique_ptr<LocationsQG> ObsSpaceQG::locations() const {
  atlas::FieldSet fields;
  std::vector<util::DateTime> times;
  qg_obsdb_locations_f90(key_, obsname_.size(), obsname_.c_str(), fields.get(), times);
  return std::unique_ptr<LocationsQG>(new LocationsQG(fields, std::move(times)));
}

// -----------------------------------------------------------------------------

int ObsSpaceQG::nobs() const {
  int iobs;
  qg_obsdb_nobs_f90(key_, obsname_.size(), obsname_.c_str(), iobs);
  return iobs;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
ObsIteratorQG ObsSpaceQG::begin() const {
  return ObsIteratorQG(*this->locations(), 0);
}
// -----------------------------------------------------------------------------
ObsIteratorQG ObsSpaceQG::end() const {
  return ObsIteratorQG(*this->locations(), this->nobs());
}
// -----------------------------------------------------------------------------

void ObsSpaceQG::print(std::ostream & os) const {
  os << "ObsSpace for " << obsname_ << ", " << this->nobs() << " obs";
}

// -----------------------------------------------------------------------------

}  // namespace qg
