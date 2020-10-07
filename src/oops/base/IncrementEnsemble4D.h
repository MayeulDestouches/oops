/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INCREMENTENSEMBLE4D_H_
#define OOPS_BASE_INCREMENTENSEMBLE4D_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/Accumulator.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/StateEnsemble4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Ensemble of 4D inrements
template<typename MODEL> class IncrementEnsemble4D {
  typedef LinearVariableChangeBase<MODEL>  LinearVariableChangeBase_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State4D<MODEL>             State4D_;
  typedef StateEnsemble4D<MODEL>     StateEnsemble4D_;
  typedef Increment4D<MODEL>         Increment4D_;

  typedef typename boost::ptr_vector<LinearVariableChangeBase_> ChvarVec_;
  typedef typename ChvarVec_::const_reverse_iterator ircst_;

 public:
  /// Constructor
  IncrementEnsemble4D(const Geometry_ & resol,
                    const Variables & vars,
                    const std::vector<util::DateTime> &,
                    const int rank);
  /// \brief construct ensemble of perturbations as \p ens - \p mean; holding
  //         \p vars variables
  IncrementEnsemble4D(const StateEnsemble4D_ & ens, const State4D_ & mean,
                    const Variables & vars);
  IncrementEnsemble4D(const eckit::Configuration &,
                    const State4D_ &, const State4D_ &, const Geometry_ &, const Variables &);

  /// Accessors
  unsigned int size() const {
    return ensemblePerturbs_.size();
  }
  Increment4D_ & operator[](const int ii) {
    return ensemblePerturbs_[ii];
  }
  const Increment4D_ & operator[](const int ii) const {
    return ensemblePerturbs_[ii];
  }

  /// Control variables
  const Variables & controlVariables() const {return vars_;}

  /// Release / reset
  void releaseMember();
  void resetMember(const Increment4D_ &);

 private:
  const Variables vars_;
  std::vector<Increment4D_> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
IncrementEnsemble4D<MODEL>::IncrementEnsemble4D(const Geometry_ & resol, const Variables & vars,
                                                const std::vector<util::DateTime> & timeslots,
                                                const int rank)
  : vars_(vars), ensemblePerturbs_()
{
  ensemblePerturbs_.reserve(rank);
  for (int m = 0; m < rank; ++m) {
    ensemblePerturbs_.emplace_back(resol, vars_, timeslots);
  }
  Log::trace() << "IncrementEnsemble4D:contructor done" << std::endl;
}

// ====================================================================================

template<typename MODEL>
IncrementEnsemble4D<MODEL>::IncrementEnsemble4D(const StateEnsemble4D_ & ensemble,
                                                const State4D_ & mean, const Variables & vars)
  : vars_(vars), ensemblePerturbs_()
{
  ensemblePerturbs_.reserve(ensemble.size());
  for (size_t ii = 0; ii < ensemble.size(); ++ii) {
    ensemblePerturbs_.emplace_back(ensemble[ii].geometry(), vars,
                                   ensemble[ii].validTimes());
    ensemblePerturbs_[ii].diff(ensemble[ii], mean);
  }
  Log::trace() << "IncrementEnsemble4D:contructor(StateEnsemble4D) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementEnsemble4D<MODEL>::IncrementEnsemble4D(const eckit::Configuration & conf,
                                                const State4D_ & xb, const State4D_ & fg,
                                                const Geometry_ & resol, const Variables & vars)
  : vars_(vars), ensemblePerturbs_()
{
  // Get rank from config
  std::vector<eckit::LocalConfiguration> memberConfig;
  conf.get("members", memberConfig);

  // Check sizes and fill in timeslots
  ASSERT(xb.size() == fg.size());
  std::vector<util::DateTime> timeslots(xb.size());
  for (unsigned jsub = 0; jsub < xb.size(); ++jsub) {
     ASSERT(xb[jsub].validTime() == fg[jsub].validTime());
     timeslots[jsub] = xb[jsub].validTime();
  }

  // Read inflation field
  std::unique_ptr<Increment4D_> inflationField;
  if (conf.has("inflation field")) {
    const eckit::LocalConfiguration inflationConfig(conf, "inflation field");
    inflationField.reset(new Increment4D_(resol, vars, timeslots));
    inflationField->read(inflationConfig);
  }

  // Get inflation value
  double inflationValue = 1;
  if (conf.has("inflation value")) {
     conf.get("inflation value", inflationValue);
  }

  // Setup change of variable
  ChvarVec_ chvars;
  if (conf.has("variable changes")) {
    std::vector<eckit::LocalConfiguration> chvarconfs;
    conf.get("variable changes", chvarconfs);
    for (const auto & conf : chvarconfs) {
      chvars.push_back(LinearVariableChangeFactory<MODEL>::create(xb[0], fg[0], resol, conf));
    }
  }
  // TODO(Benjamin): one change of variable for each timeslot

  // Read ensemble
  StateEnsemble4D_ ensemble(resol, conf);
  State4D_ bgmean = ensemble.mean();

  ensemblePerturbs_.reserve(ensemble.size());
  for (unsigned int ie = 0; ie < ensemble.size(); ++ie) {
    // Ensemble will be centered around ensemble mean
    Increment4D_ dx(resol, vars_, timeslots);
    dx.diff(ensemble[ie], bgmean);

    // Apply inflation
    if (conf.has("inflation field")) {
      dx.schur_product_with(*inflationField);
    }
    dx *= inflationValue;

    // Apply inverse of the linear balance operator
    for (unsigned jsub = 0; jsub < timeslots.size(); ++jsub) {
      // K_1^{-1} K_2^{-1} .. K_N^{-1}
      for (ircst_ it = chvars.rbegin(); it != chvars.rend(); ++it) {
        dx[jsub] = it->multiplyInverse(dx[jsub]);
      }
    }

    ensemblePerturbs_.emplace_back(std::move(dx));
  }
  Log::trace() << "IncrementEnsemble4D:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble4D<MODEL>::releaseMember() {
  ensemblePerturbs_.erase(ensemblePerturbs_.begin());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble4D<MODEL>::resetMember(const Increment4D_ & dx) {
  ensemblePerturbs_.emplace_back(dx);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INCREMENTENSEMBLE4D_H_