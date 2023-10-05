#ifndef PDEOPT_QOI_L2TRACKING_BRANDES_GRIESSE_HPP
#define PDEOPT_QOI_L2TRACKING_BRANDES_GRIESSE_HPP

#include "../../../../PDE-OPT/TOOLS/qoi.hpp"
#include "pde_brandes_griesse.hpp"

template <class Real>
class QoI_L2Tracking_Brandes_Griesse : public QoI<Real> {
private:
  ROL::Ptr<FE<Real> > fe_;

public:
  QoI_L2Tracking_Brandes_Griesse(const ROL::Ptr<FE<Real> > &fe)
  {
    fe_ = fe;
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real> >(c);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Create array of 1's to integrate against
    ROL::Ptr<Intrepid::FieldContainer<Real> > const_weights = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,p);
    const_weights->initialize(-0.25);
    fe_->computeIntegral(val,valU_eval,const_weights);

    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    // Create constant array
    ROL::Ptr<Intrepid::FieldContainer<Real> > const_weights = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,p);
    const_weights->initialize(-0.25);
    // Compute gradient of (-1/4) of average 
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *const_weights,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Brandes_Griesse::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Brandes_Griesse::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Brandes_Griesse::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Brandes_Griesse::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Brandes_Griesse::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

template <class Real>
class QoI_L2Penalty_Brandes_Griesse : public QoI<Real> {
private:
  ROL::Ptr<FE<Real> > fe_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > weight_;

  void computeDomainWeight(void) {
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->cubPts()->dimension(2);

    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,p);
    
    const Real zero(0), one(1);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        if ( insideDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j) {
        (*weight_)(i,j) = (inside ? one : zero);
      }
    }
  }

  bool insideDomain(const std::vector<Real> &x) const {
    bool inside(false);
    if(x[0]>-0.75 && x[0]<0.75 && x[1]>-0.8 && x[1]<-0.6)
      {
	inside = true;
      }
    else if(x[0]>-0.8 && x[0]<0.0 && x[1]>0.4 && x[1]<0.8)
      {
    	inside = true;
      }
    else
      {
	inside = false;
      }
    return inside;
  }

public:
  QoI_L2Penalty_Brandes_Griesse(const ROL::Ptr<FE<Real> > &fe)
  {
    fe_ = fe;
    computeDomainWeight();
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real> >(c);
    // Build local state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > weight_valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weight_valZ_eval,*weight_,*valZ_eval);
    fe_->computeIntegral(val,weight_valZ_eval,valZ_eval);

    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Brandes_Griesse::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    // Build local gradient of state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
 
    ROL::Ptr<Intrepid::FieldContainer<Real> > weight_valZ_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weight_valZ_eval,*weight_,*valZ_eval);
     Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *weight_valZ_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Brandes_Griesse::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Brandes_Griesse::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Brandes_Griesse::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    hess = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > weight_valV_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weight_valV_eval,*weight_,*valV_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *weight_valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

}; // QoI_L2Penalty


#endif
