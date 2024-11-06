#ifndef PDEOPT_STATE_COST_DARCY_FLOW_HPP
#define PDEOPT_STATE_COST_DARCY_FLOW_HPP

#include "pde_darcy_flow.hpp"

template <class Real>
class State_Cost_darcy_flow : public ROL::Objective_SimOpt<Real> {
private:
  std::vector<Real> target_data_;
  std::vector<int> target_data_ids_;
  ROL::Ptr<Tpetra::MultiVector<> > w_ptr_;
  ROL::Ptr<ROL::Vector<Real> > wp_;
  ROL::Ptr<Tpetra::MultiVector<> > td_ptr_;
  ROL::Ptr<ROL::Vector<Real> > tdp_;

  void Apply_Weight_Vector(ROL::Vector<Real> & vec)
  {
    vec.applyBinary(ROL::Elementwise::Multiply<Real>(), *wp_);
  }

public:
  State_Cost_darcy_flow(std::vector<Real> & target_data, std::vector<int> target_data_ids, ROL::Ptr<ROL::Vector<Real> > up)
    : target_data_(target_data), target_data_ids_(target_data_ids)
  {
    wp_ = up->clone();
    w_ptr_ = dynamic_cast<ROL::TpetraMultiVector<Real>&>(*wp_).getVector();
    tdp_ = up->clone();
    td_ptr_ = dynamic_cast<ROL::TpetraMultiVector<Real>&>(*tdp_).getVector();
    wp_->zero();
    tdp_->zero();
    for(unsigned k = 0; k < target_data_ids_.size(); k++)
      {
	w_ptr_->replaceGlobalValue(target_data_ids_[k],0,1.0);
	td_ptr_->replaceGlobalValue(target_data_ids_[k],0,target_data_[k]);
      }
  }
  
  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<ROL::Vector<Real> > vec = u.clone();
    vec->set(u);
    vec->axpy(-1.0,*tdp_);
    Apply_Weight_Vector(*vec);
    return 0.5*(vec->dot(*vec));
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    g.set(u);
    g.scale(1.0);
    g.axpy(-1.0,*tdp_);
    Apply_Weight_Vector(g);
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    g.zero();
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    hv.set(v);
    Apply_Weight_Vector(hv);
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  } 


}; // State_Cost

template <class Real>
class QoI_H1Penalty : public QoI<Real> {

  ROL::Ptr<FE<Real> > fe_;

public:
  QoI_H1Penalty(const ROL::Ptr<FE<Real> > &fe) : fe_(fe) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->gradN()->dimension(3);

    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real> >(c);

    ROL::Ptr<Intrepid::FieldContainer<Real> > gradZ_eval;
    gradZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradZ_eval, z_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > gradZNorm_eval;
    gradZNorm_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);

    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*gradZNorm_eval)(i,j) = (*gradZ_eval)(i,j,0)*(*gradZ_eval)(i,j,0) + (*gradZ_eval)(i,j,1)*(*gradZ_eval)(i,j,1);
	  }
      }

    // Create array of 1's to integrate against
    ROL::Ptr<Intrepid::FieldContainer<Real> > const_weights = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,p);
    const_weights->initialize(1.0);
    fe_->computeIntegral(val,gradZNorm_eval,const_weights);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));

    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_H1_CDR::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    int d = fe_->gradN()->dimension(3);

    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradZ_eval;
    gradZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradZ_eval, z_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > diffZx;
    diffZx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffZy;
    diffZy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);

    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*diffZx)(i,j) = (*gradZ_eval)(i,j,0);
	    (*diffZy)(i,j) = (*gradZ_eval)(i,j,1);
	  }
      }

    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *diffZx,
                                                  *(fe_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *diffZy,
                                                  *(fe_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_H1_CDR::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_H1_CDR::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_H1_CDR::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    int d = fe_->gradN()->dimension(3);
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradV_eval;
    gradV_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradV_eval, v_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > diffx;
    diffx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffy;
    diffy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);

    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*diffx)(i,j) = (*gradV_eval)(i,j,0);
	    (*diffy)(i,j) = (*gradV_eval)(i,j,1);
	  }
      }

    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *diffx,
                                                  *(fe_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *diffy,
                                                  *(fe_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);
  }

}; // QoI_H1Penalty

template <class Real>
class QoI_L2Penalty : public QoI<Real> {
private:
  ROL::Ptr<FE<Real> > fe_;

public:
  QoI_L2Penalty(const ROL::Ptr<FE<Real> > &fe) : fe_(fe) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Build local state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    fe_->computeIntegral(val,valZ_eval,valZ_eval);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_CDR::gradient_1 is zero.");
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
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Build local gradient of state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *valZ_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_CDR::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_CDR::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_CDR::HessVec_21 is zero.");
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
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

}; // QoI_L2Penalty

#endif
