/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_GEOMETRYQG_H_
#define QG_MODEL_GEOMETRYQG_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/projection.h"

#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "oops/qg/GeometryQGIterator.h"
#include "oops/qg/QgFortran.h"

namespace oops {
  class Variables;
}

namespace qg {

class GeometryQgParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryQgParameters, Parameters)

 public:
  /// Number of cells
  oops::Parameter<int> nx{"nx", 60, this};
  oops::Parameter<int> ny{"ny", 19, this};
  /// Depths
  oops::OptionalParameter<std::vector<double>> depths{"depths", this};
  /// Heating option (AS: should it be in geometry or model?)
  oops::Parameter<bool> heating{"heating", true, this};
    /// Modified QG option
  oops::Parameter<float> perturbedheat{"perturbed heating", 1.0, this};
  /// Interpolation option
  oops::Parameter<std::string> interpolator{"interpolator", "trilinear", this};
};

class GeometryQGIterator;

// -----------------------------------------------------------------------------
/// GeometryQG handles geometry for QG model.

class GeometryQG : public util::Printable,
                   private util::ObjectCounter<GeometryQG> {
 public:
  typedef GeometryQgParameters Parameters_;

  static const std::string classname() {return "qg::GeometryQG";}

  GeometryQG(const GeometryQgParameters &,
             const eckit::mpi::Comm & comm = oops::mpi::world());
  GeometryQG(const GeometryQG &);
  ~GeometryQG();

  const F90geom & toFortran() const {return keyGeom_;}

  GeometryQGIterator begin() const;
  GeometryQGIterator end() const;
  std::vector<double> verticalCoord(std::string &) const;
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::Grid * atlasGrid() const {return atlasGrid_.get();}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}
  size_t levels() const {return levs_;}

  std::vector<size_t> variableSizes(const oops::Variables & vars) const;

 private:
  GeometryQG & operator=(const GeometryQG &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::StructuredGrid> atlasGrid_;
  std::unique_ptr<atlas::Projection> atlasProjection_;
  std::unique_ptr<atlas::functionspace::StructuredColumns> atlasFunctionSpace_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
  size_t levs_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_GEOMETRYQG_H_
