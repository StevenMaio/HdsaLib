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
class QoI_L2Penalty_darcy_flow : public QoI<Real> {

public:
  QoI_L2Penalty_darcy_flow(void)
  { }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int size = z_param->size();
    Real sum(0);
    for (int i = 0; i < size; ++i) {
      sum += (*z_param)[i]*(*z_param)[i];
    }
    sum *= 0.5;
    val = ROL::nullPtr;
    return sum;
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & grad,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
			       const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int size = z_param->size();
    std::vector<Real> g(size,static_cast<Real>(1));
    for(int k = 0; k < size; k++)
      {
	g[k] = (*z_param)[k];
      }
    return g;
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::HessVec_12 is zero.");
  }

  void HessVec_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::HessVec_13 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::HessVec_22 is zero.");
  }

  void HessVec_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
			       const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
			       const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_darcy_flow::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
			       const ROL::Ptr<const std::vector<Real> > & v_param,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
			       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
			       const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int size = z_param->size();
    std::vector<Real> h(size,static_cast<Real>(1));
    for(int k = 0; k < size; k++)
      {
	h[k] = (*v_param)[k];
      }
    return h;
  }

}; // QoI_L2Penalty

#endif
