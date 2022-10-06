/*
 * (C) Copyright 2022- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/L95Traits.h"
#include "model/QgTraits.h"
#include "oops/coupled/TraitCoupled.h"
#include "oops/runs/ConvertState.h"
#include "oops/runs/Run.h"

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  oops::ConvertState<oops::TraitCoupled<qg::QgTraits, lorenz95::L95Traits>> app;
  return run.execute(app);
}
