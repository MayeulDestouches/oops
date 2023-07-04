/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <math.h>
#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/projection.h"
#include "atlas/util/Point.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "model/GeometryQG.h"
#include "model/QgFortran.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
GeometryQG::GeometryQG(const GeometryQgParameters & params,
                       const eckit::mpi::Comm & comm) : comm_(comm), levs_(1) {
  ASSERT(comm_.size() == 1);

  // Get hard-coded geometry parameters
  double lat_min, lat_max, domain_zonal, domain_meridional, xmin, ymin, lat_proj, domain_depth;
  qg_domain_parameters_f90(lat_min, lat_max, domain_zonal, domain_meridional, xmin, ymin,
    lat_proj, domain_depth);

  // Get geometry input parameters
  eckit::LocalConfiguration geomConfig(params.toConfiguration());
  const int nx = geomConfig.getInt("nx", 60);
  const int ny = geomConfig.getInt("ny", 19);
  const boost::optional<std::vector<double>> &depths = params.depths.value();
  if (depths != boost::none) {
    levs_ = (*depths).size();
  }

  // Define projection configuration
  eckit::LocalConfiguration projConfig;
  projConfig.set("type", "mercator");
  projConfig.set("latitude1", lat_proj);
  const std::vector<double> normalise{-180, 180};
  projConfig.set("normalise", normalise);

  // Initialize eckit communicator for ATLAS
  eckit::mpi::setCommDefault(comm_.name().c_str());

  // Setup projection
  atlasProjection_.reset(new atlas::Projection(projConfig));

  // Define geometry and projection configurations
  const double dx = domain_zonal/static_cast<double>(nx);
  const double dy = domain_meridional/static_cast<double>(ny+1);
  std::vector<double> xymin;
  xymin.push_back(-180.0);
  atlas::Point2 p;
  p[0] = 0.0;
  p[1] = ymin+dy;
  atlasProjection_->xy2lonlat(p);
  xymin.push_back(p[1]);
  geomConfig.set("type", "regional");
  geomConfig.set("dx", dx);
  geomConfig.set("dy", dy);
  geomConfig.set("lonlat(xmin,ymin)", xymin);
  geomConfig.set("projection", projConfig);
  oops::Log::info() << "Geometry configuration: " << geomConfig << std::endl;

  // Setup regional grid
  atlasGrid_.reset(new atlas::StructuredGrid(geomConfig));

  // Setup partitioner
  const atlas::grid::Partitioner partitioner("serial");

  // Setup distribution
  const atlas::grid::Distribution distribution(*atlasGrid_, partitioner);

  // Setup function space
  atlasFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(*atlasGrid_, distribution,
  geomConfig));

  // Setup Fortran geometry
  qg_geom_setup_f90(keyGeom_, geomConfig, atlasGrid_->get(), atlasProjection_->get(),
                    atlasFunctionSpace_->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  qg_geom_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());
}
// -----------------------------------------------------------------------------
GeometryQG::GeometryQG(const GeometryQG & other) : comm_(other.comm_), levs_(other.levs_) {
  ASSERT(comm_.size() == 1);

  // Copy ATLAS grid
  atlasGrid_.reset(new atlas::StructuredGrid(*(other.atlasGrid_)));

  // Copy ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(
                            *(other.atlasFunctionSpace_)));

  // Copy Fortran geometry
  qg_geom_clone_f90(keyGeom_, other.keyGeom_);

  // Copy ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
    atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
    atlasFieldSet_->add(atlasField);
  }
}
// -----------------------------------------------------------------------------
GeometryQG::~GeometryQG() {
  qg_geom_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
GeometryQGIterator GeometryQG::begin() const {
  return GeometryQGIterator(*this);
}
// -----------------------------------------------------------------------------
GeometryQGIterator GeometryQG::end() const {
  int nx = 0;
  int ny = 0;
  int nz = 1;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  return GeometryQGIterator(*this, nx*ny+1);
}
// -------------------------------------------------------------------------------------------------
std::vector<double> GeometryQG::verticalCoord(std::string & vcUnits) const {
  // returns vertical coordinate in untis of vcUnits
  int nx = 0;
  int ny = 0;
  int nz = 1;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  std::vector<double> vc(nz);
  if (vcUnits == "levels") {
    for (int i=0; i < nz; ++i) {vc[i]=i+1;}
  } else {
    std::stringstream errorMsg;
    errorMsg << "Uknown vertical coordinate unit " << vcUnits << std::endl;
    ABORT(errorMsg.str());
  }
  oops::Log::debug() << "QG vert coord: " << vc << std::endl;
  return vc;
}
// -------------------------------------------------------------------------------------------------
std::vector<size_t> GeometryQG::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> sizes(vars.size(), levs_);
  return sizes;
}
// -----------------------------------------------------------------------------
void GeometryQG::print(std::ostream & os) const {
  int nx;
  int ny;
  int nz;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  os << "Geometry:" << std::endl;
  os << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
  os << "deltax = " << deltax << ", deltay = " << deltay;
}
// -----------------------------------------------------------------------------
}  // namespace qg
