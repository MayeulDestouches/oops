/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_INSTANTIATELOCALIZATIONFACTORY_H_
#define OOPS_GENERIC_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/interface/LocalizationBase.h"
#include "oops/generic/LocalizationHDIAG_NICAS.h"

namespace oops {

template <typename MODEL> void instantiateLocalizationFactory() {
  static LocalizationMaker<MODEL, LocalizationHDIAG_NICAS<MODEL> > makerHdiag_nicas_("HDIAG_NICAS");
}

}  // namespace oops

#endif  // OOPS_GENERIC_INSTANTIATELOCALIZATIONFACTORY_H_
