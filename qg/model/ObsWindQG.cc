/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsWindQG.h"

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "model/GomQG.h"
#include "model/ObsBias.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgFortran.h"
#include "util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
static oops::ObsOperatorMaker<QgTraits, ObsWindQG>   makerWind_("Wind");
// -----------------------------------------------------------------------------

ObsWindQG::ObsWindQG(const ObsSpaceQG &, const eckit::Configuration & config)
  : keyOperWind_(0), varin_(std::vector<std::string>{"u","v"})
{
  const eckit::Configuration * configc = &config;
  qg_wind_setup_f90(keyOperWind_, &configc);
  oops::Log::trace() << "ObsWindQG created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsWindQG::~ObsWindQG() {
  qg_wind_delete_f90(keyOperWind_);
}

// -----------------------------------------------------------------------------

void ObsWindQG::obsEquiv(const GomQG & gom, ObsVecQG & ovec,
                         const ObsBias & bias) const {
  qg_wind_equiv_f90(gom.toFortran(), ovec.toFortran(), bias.wind());
}

// -----------------------------------------------------------------------------

void ObsWindQG::print(std::ostream & os) const {
  os << "ObsWindQG::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace qg
