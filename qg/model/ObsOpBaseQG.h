/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_OBSOPBASEQG_H_
#define QG_MODEL_OBSOPBASEQG_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "model/ObsSpaceQG.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

namespace qg {
class GomQG;
class LocationsQG;
class ObsBias;
class ObsVecQG;

// -----------------------------------------------------------------------------
/// Base class for observation operators

class ObsOpBaseQG : public util::Printable,
                    private boost::noncopyable {
 public:
  ObsOpBaseQG() {}
  virtual ~ObsOpBaseQG() {}

/// Obs Operator
  virtual void simulateObs(const GomQG &, ObsVecQG &, const ObsBias &) const = 0;

/// Other
  virtual const oops::Variables & requiredVars() const = 0;  // Required from Model
  virtual LocationsQG * locations(const util::DateTime &, const util::DateTime &) const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Operator Factory
class ObsOpFactory {
 public:
  static ObsOpBaseQG * create(const ObsSpaceQG &, const eckit::Configuration &);
  virtual ~ObsOpFactory() { getMakers().clear(); }
 protected:
  explicit ObsOpFactory(const std::string &);
 private:
  virtual ObsOpBaseQG * make(const ObsSpaceQG &, const eckit::Configuration &) = 0;
  static std::map < std::string, ObsOpFactory * > & getMakers() {
    static std::map < std::string, ObsOpFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsOpMaker : public ObsOpFactory {
  virtual ObsOpBaseQG * make(const ObsSpaceQG & odb, const eckit::Configuration & conf)
    { return new T(odb, conf); }
 public:
  explicit ObsOpMaker(const std::string & name) : ObsOpFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSOPBASEQG_H_
