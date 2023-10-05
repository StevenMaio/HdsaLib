#ifndef PDEOPT_STATE_COST_DARCY_ADV_DIFF_HPP
#define PDEOPT_STATE_COST_DARCY_ADV_DIFF_HPP

#include "pde_darcy_adv_diff.hpp"

template <class Real>
class State_Cost_darcy_adv_diff : public Objective_SimOpt_TS<Real> {
private:
  std::vector<std::vector<Real> > data_; // data[k][i] is the kth time instance in the ith spatial node of the true state
  std::vector<Real> data_weight_;
  std::vector<int> data_weight_id_;
  ROL::Ptr<Tpetra::MultiVector<> > w_ptr_;
  ROL::Ptr<ROL::Vector<Real> > wp_;
  ROL::Ptr<Tpetra::MultiVector<> > noisy_data_ptr_;
  ROL::Ptr<ROL::Vector<Real> > noisy_data_;

  void Apply_Weight_Vector(ROL::Vector<Real> & vec)
  {
    wp_->zero();
    for(unsigned k = 0; k < data_weight_.size(); k++)
      {
	w_ptr_->replaceGlobalValue(data_weight_id_[k],0,data_weight_[k]);
      }
    vec.applyBinary(ROL::Elementwise::Multiply<Real>(), *wp_);
  }

  void Apply_Sqrt_Weight_Vector(ROL::Vector<Real> & vec)
  {
    wp_->zero();
    for(unsigned k = 0; k < data_weight_.size(); k++)
      {
	w_ptr_->replaceGlobalValue(data_weight_id_[k],0,std::sqrt(data_weight_[k]));
      }
    vec.applyBinary(ROL::Elementwise::Multiply<Real>(), *wp_);
  }

public:
  State_Cost_darcy_adv_diff(std::vector<ROL::TimeStamp<Real> > & timeStamp, std::vector<std::vector<Real> > & data, 
			    const std::vector<Real> & data_weight, const std::vector<int> & data_weight_id, ROL::Ptr<ROL::Vector<Real> > & up)
    : Objective_SimOpt_TS<Real>(timeStamp), data_(data), data_weight_(data_weight), data_weight_id_(data_weight_id)
  { 
    wp_ = up->clone();
    w_ptr_ = dynamic_cast<ROL::TpetraMultiVector<Real>&>(*wp_).getVector();
    noisy_data_ = up->clone();
    noisy_data_ptr_ = dynamic_cast<ROL::TpetraMultiVector<Real>&>(*noisy_data_).getVector();
    wp_->zero();
    noisy_data_->zero();

    for(unsigned k = 0; k <data_weight_.size(); k++)
      {
	int index = data_weight_id_[k];
	noisy_data_ptr_->replaceGlobalValue(index,0,data[0][index]);
      }
  }
  
  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<ROL::Vector<Real> > vec = u.clone();
    vec->set(u);
    vec->axpy(-1.0,*noisy_data_);
    Apply_Sqrt_Weight_Vector(*vec);
    return 0.5*(vec->dot(*vec));
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    g.set(u);
    g.scale(1.0);
    g.axpy(-1.0,*noisy_data_);
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

  void update( const ROL::Vector<Real> & u, const ROL::Vector<Real> & z, Real & t, bool flag = true, int iter = -1 ) 
  {
    Objective_SimOpt_TS<Real>::update(u,z,t,flag,iter);
    Update_Noisy_Data();
  }

  void Update_Noisy_Data()
  {
    noisy_data_->zero();
    for(unsigned k = 0; k < data_weight_.size(); k++)
      {
	int index = data_weight_id_[k];
	noisy_data_ptr_->replaceGlobalValue(index,0,data_[Objective_SimOpt_TS<Real>::current_TS_][index]);
      }
  }

}; // State_Cost

template <class Real>
class QoI_H1_darcy_adv_diff : public QoI<Real> {

  ROL::Ptr<FE<Real> > fe_;
  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

public:
  QoI_H1_darcy_adv_diff(const ROL::Ptr<FE<Real> > &fe, const ROL::Ptr<FieldHelper<Real> > & fieldHelper) : fe_(fe), fieldHelper_(fieldHelper) {}

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

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > gradKappa_eval;
    gradKappa_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradKappa_eval, Z[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > gradKappaNorm_eval;
    gradKappaNorm_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);

    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*gradKappaNorm_eval)(i,j) = (*gradKappa_eval)(i,j,0)*(*gradKappa_eval)(i,j,0) + (*gradKappa_eval)(i,j,1)*(*gradKappa_eval)(i,j,1);
	  }
      }

    // Create array of 1's to integrate against
    ROL::Ptr<Intrepid::FieldContainer<Real> > const_weights = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,p);
    const_weights->initialize(1.0);
    fe_->computeIntegral(val,gradKappaNorm_eval,const_weights);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));

    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_H1_darcy_adv-diff::gradient_1 is zero.");
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

    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(2);
    for (int i=0; i<2; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Build local gradient of state tracking term
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradKappa_eval;
    gradKappa_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradKappa_eval, Z[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > diffKappax;
    diffKappax = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffKappay;
    diffKappay = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);

    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*diffKappax)(i,j) = (*gradKappa_eval)(i,j,0);
	    (*diffKappay)(i,j) = (*gradKappa_eval)(i,j,1);
	  }
      }

    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *diffKappax,
                                                  *(fe_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *diffKappay,
                                                  *(fe_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_H1_darcy_adv-diff::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_H1_darcy_adv-diff::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_H1_darcy_adv-diff::HessVec_21 is zero.");
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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(2);
    for (int i=0; i<2; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradV_eval;
    gradV_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradV_eval, V[0]);

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

    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *diffx,
                                                  *(fe_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *diffy,
                                                  *(fe_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);

    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_H1_darcy_adv-diff


template <class Real>
class QoI_L2_darcy_adv_diff : public QoI<Real> {

  ROL::Ptr<FE<Real> > fe_;
  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

public:
  QoI_L2_darcy_adv_diff(const ROL::Ptr<FE<Real> > &fe, const ROL::Ptr<FieldHelper<Real> > & fieldHelper) : fe_(fe), fieldHelper_(fieldHelper) {}

   Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real> >(c);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(kappa, Z[0]);

    // Build local state tracking term
    fe_->computeIntegral(val,kappa,kappa);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2_darcy_adv-diff::gradient_1 is zero.");
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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(2);
    for (int i=0; i<2; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Build local gradient of state tracking term
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(kappa, Z[0]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *kappa,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2_darcy_adv-diff::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2_darcy_adv-diff::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2_darcy_adv-diff::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(2);
    for (int i=0; i<2; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valV_eval, V[0]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_L2_darcy_adv-diff

template <class Real>
class QoI_Informed_darcy_adv_diff : public QoI<Real> {

  ROL::Ptr<FE<Real> > fe_;
  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

public:
  QoI_Informed_darcy_adv_diff(const ROL::Ptr<FE<Real> > &fe, const ROL::Ptr<FieldHelper<Real> > & fieldHelper) : fe_(fe), fieldHelper_(fieldHelper) {}

  Real Evaluate_Target(std::vector<Real> & coords)
  {
    return 4.0*std::exp(-10.0*(coords[1]-0.5)*(coords[1]-0.5))-2.0;
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

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(kappa, Z[0]);

    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
	(*kappa)(i,j) = (*kappa)(i,j) - Evaluate_Target(pt);
      }
    }

    // Build local state tracking term
    fe_->computeIntegral(val,kappa,kappa);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Informed_darcy_adv-diff::gradient_1 is zero.");
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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(2);
    for (int i=0; i<2; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Build local gradient of state tracking term
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(kappa, Z[0]);

    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
	(*kappa)(i,j) = (*kappa)(i,j) - Evaluate_Target(pt);
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *kappa,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Informed_darcy_adv-diff::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Informed_darcy_adv-diff::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Informed_darcy_adv-diff::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(2);
    for (int i=0; i<2; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valV_eval, V[0]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_Informed_darcy_adv-diff


#endif
