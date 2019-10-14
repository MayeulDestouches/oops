/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTFCT3DVAR_H_
#define OOPS_ASSIMILATION_COSTFCT3DVAR_H_

#include <memory>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/VariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// 3D-Var Cost Function
/*!
 * This class is not really necessary since it is only a special
 * case of the more general 4D-Var cost function. It is provided
 * for readability.
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class CostFct3DVar : public CostFunction<MODEL> {
  typedef Increment<MODEL>            Increment_;
  typedef ControlIncrement<MODEL>     CtrlInc_;
  typedef ControlVariable<MODEL>      CtrlVar_;
  typedef CostFunction<MODEL>         CostFct_;
  typedef Geometry<MODEL>             Geometry_;
  typedef State<MODEL>                State_;
  typedef Model<MODEL>                Model_;
  typedef VariableChangeBase<MODEL>   ChangeVar_;
  typedef LinearVariableChangeBase<MODEL> ChangeVarTLAD_;

 public:
  CostFct3DVar(const eckit::Configuration &, const Geometry_ &, const Model_ &);
  virtual ~CostFct3DVar() {}

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const override;

  CostJb3D<MODEL>     * newJb(const eckit::Configuration &, const Geometry_ &,
                              const CtrlVar_ &) const override;
  CostJo<MODEL>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &,
                   const CtrlVar_ &, const CtrlVar_ &) override;

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::DateTime windowHalf_;
  util::Duration zero_;
  const Variables ctlvars_;
  std::unique_ptr<ChangeVar_> an2model_;
  std::unique_ptr<ChangeVarTLAD_> inc2model_;
};

// =============================================================================

template<typename MODEL>
CostFct3DVar<MODEL>::CostFct3DVar(const eckit::Configuration & config,
                                  const Geometry_ & resol, const Model_ & model)
  : CostFunction<MODEL>::CostFunction(config, resol, model),
    windowLength_(), windowHalf_(), zero_(0), ctlvars_(config),
    an2model_(VariableChangeFactory<MODEL>::create(config, resol)), inc2model_()
{
  Log::trace() << "CostFct3DVar::CostFct3DVar start" << std::endl;
  windowLength_ = util::Duration(config.getString("window_length"));
  windowBegin_ = util::DateTime(config.getString("window_begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  windowHalf_ = windowBegin_ + windowLength_/2;
  this->setupTerms(config);
  Log::trace() << "CostFct3DVar::CostFct3DVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJb3D<MODEL> * CostFct3DVar<MODEL>::newJb(const eckit::Configuration & jbConf,
                                             const Geometry_ & resol,
                                             const CtrlVar_ & xb) const {
  Log::trace() << "CostFct3DVar::newJb" << std::endl;
  ASSERT(xb.state().checkStatesNumber(1));
  return new CostJb3D<MODEL>(jbConf, resol, ctlvars_, zero_, xb.state()[0]);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJo<MODEL> * CostFct3DVar<MODEL>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFct3DVar::newJo" << std::endl;
  return new CostJo<MODEL>(joConf, this->getComm(), windowBegin_, windowEnd_, windowLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostTermBase<MODEL> * CostFct3DVar<MODEL>::newJc(const eckit::Configuration & jcConf,
                                                 const Geometry_ &) const {
  Log::trace() << "CostFct3DVar::newJc" << std::endl;
// For now there is no Jc that can work with 3D-Var
  CostTermBase<MODEL> * pjc = 0;
  return pjc;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFct3DVar::runNL start" << std::endl;
  ASSERT(xx.state().checkStatesNumber(1));
  ASSERT(xx.state()[0].validTime() == windowHalf_);
  State_ xm(xx.state()[0].geometry(), CostFct_::getModel().variables(), windowHalf_);

  an2model_->changeVar(xx.state()[0], xm);

  post.initialize(xm, windowHalf_, windowLength_);
  post.process(xm);
  post.finalize(xm);

  ASSERT(xx.state()[0].validTime() == windowHalf_);
  Log::trace() << "CostFct3DVar::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFct3DVar<MODEL>::doLinearize(const Geometry_ & resol,
                                      const eckit::Configuration & innerConf,
                                      const CtrlVar_ & bg, const CtrlVar_ & fg) {
  Log::trace() << "CostFct3DVar::doLinearize start" << std::endl;
  eckit::LocalConfiguration conf(innerConf, "linearmodel");
  inc2model_.reset(LinearVariableChangeFactory<MODEL>::create(bg.state()[0], fg.state()[0],
                                                              resol, conf));
  inc2model_->setInputVariables(ctlvars_);
  inc2model_->setOutputVariables(CostFct_::getTLM().variables());
  Log::trace() << "CostFct3DVar::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runTLM(CtrlInc_ & dx,
                                 PostProcessorTLAD<MODEL> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool) const {
  Log::trace() << "CostFct3DVar::runTLM start" << std::endl;
  ASSERT(inc2model_);
  ASSERT(dx.state()[0].validTime() == windowHalf_);

  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowHalf_);
  inc2model_->multiply(dx.state()[0], dxmodel);

  cost.initializeTL(dxmodel, windowHalf_, windowLength_);
  post.initialize(dxmodel, windowHalf_, windowLength_);

  cost.processTL(dxmodel);
  post.process(dxmodel);

  cost.finalizeTL(dxmodel);
  post.finalize(dxmodel);

  ASSERT(dx.state()[0].validTime() == windowHalf_);
  Log::trace() << "CostFct3DVar::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFct3DVar::zeroAD start" << std::endl;
  dx.state()[0].zero(windowHalf_);
  dx.modVar().zero();
  dx.obsVar().zero();
  Log::trace() << "CostFct3DVar::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runADJ(CtrlInc_ & dx,
                                 PostProcessorTLAD<MODEL> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool) const {
  Log::trace() << "CostFct3DVar::runADJ start" << std::endl;
  ASSERT(inc2model_);
  ASSERT(dx.state()[0].validTime() == windowHalf_);

  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowHalf_);
  inc2model_->multiplyInverseAD(dx.state()[0], dxmodel);

  post.initialize(dxmodel, windowHalf_, windowLength_);
  cost.initializeAD(dxmodel, windowHalf_, windowLength_);

  cost.processAD(dxmodel);
  post.process(dxmodel);

  cost.finalizeAD(dxmodel);
  post.finalize(dxmodel);

  inc2model_->multiplyAD(dxmodel, dx.state()[0]);
  ASSERT(dx.state()[0].validTime() == windowHalf_);

  Log::trace() << "CostFct3DVar::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFct3DVar<MODEL>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                  PostProcessor<Increment_> &) const {
  Log::trace() << "CostFct3DVar::addIncr start" << std::endl;
  ASSERT(xx.state().checkStatesNumber(1));
  ASSERT(xx.state()[0].validTime() == windowHalf_);
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  xx.state()[0] += dx.state()[0];
  Log::trace() << "CostFct3DVar::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT3DVAR_H_
