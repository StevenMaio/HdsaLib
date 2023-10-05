#ifndef HDSA_ROL_B_TRANSPOSE_STEADY_STATE_HPP
#define HDSA_ROL_B_TRANSPOSE_STEADY_STATE_HPP

#include "HDSA_ROL_B_Transpose.hpp"

template <class RealT>
class ROL_B_Transpose_Steady_State: public ROL_B_Transpose<RealT> {

private:
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_;
  HDSA::Ptr<ROL::Vector<RealT> > u_, z_;
  HDSA::Ptr<HDSA::Vector<RealT> > theta_;
  int theta_dim_;
  RealT tol_;
  RealT h_;
  HDSA::Ptr<ROL::Vector<RealT> > adjoint_nominal_;
  HDSA::Ptr<HDSA::Vector<RealT> > theta_grad_nominal_;
  bool isNominalComputed_;

public:

  ROL_B_Transpose_Steady_State(const HDSA::Ptr<ROL::Objective_SimOpt<RealT> > & obj, const HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > & con, const HDSA::Ptr<ROL::Vector<RealT> > & u,
			       const HDSA::Ptr<ROL::Vector<RealT> > & z, const HDSA::Ptr<HDSA::Vector<RealT> > & theta):
    obj_(obj), con_(con), u_(u), z_(z), theta_(theta)
  { 
    theta_dim_ = theta_->dimension();
    tol_ = 1.e-8;
    h_ = 1.e-4;
    isNominalComputed_ = false;
  }

  virtual void applyThetaJacobianTranspose(HDSA::Vector<RealT> & grad, const ROL::Vector<RealT> & lambda, const ROL::Vector<RealT> & u, const ROL::Vector<RealT> & z, const HDSA::Vector<RealT> & theta) = 0; 
  
  virtual void Update_Constraint(const HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > & con, const ROL::Vector<RealT> & u,
				 const ROL::Vector<RealT> & z, const HDSA::Vector<RealT> & theta) const = 0;

  // evaluate the theta,z hessian vector product, i.e. -B^T in HDSA
  void hessVec_theta_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
		       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  { 
    if(!isNominalComputed_)
      {
	Compute_Adjoint_Nominal(z,theta);
	isNominalComputed_ = true;
      }

    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    RealT v_norm = v_rol->norm();
    RealT z_norm = z_rol->norm();
    HDSA::Ptr<ROL::Vector<RealT> > adjoint_perturbed = u_->clone();
    HDSA::Ptr<ROL::Vector<RealT> > z_pert = z_->clone();
    z_pert->set(*z_rol);
    z_pert->axpy((z_norm/v_norm)*h_,*v_rol);
    
    HDSA::Ptr<ROL::Vector<RealT> > dualstate = u_->dual().clone();
    HDSA::Ptr<ROL::Vector<RealT> > u_pert = u_->dual().clone();
    Update_Constraint(con_,*u_,*z_,theta);
    con_->solve(*dualstate,*u_pert,*z_pert,tol_);
    
    dualstate->set(*u_pert);
    obj_->gradient_1(*dualstate,*u_pert,*z_pert,tol_);
    con_->applyInverseAdjointJacobian_1(*adjoint_perturbed,*dualstate,*u_pert,*z_pert,tol_);
    adjoint_perturbed->scale(static_cast<RealT>(-1.0));
    
    applyThetaJacobianTranspose(hv, *adjoint_perturbed, *u_pert, *z_pert, *theta_);

    hv.axpy(-1.0,*theta_grad_nominal_);
    hv.scale((v_norm/z_norm)/h_);
  }

  void Compute_Adjoint_Nominal(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta)
  {
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    z_->set(*z_rol);

    theta_->set(theta);

    adjoint_nominal_ = u_->clone();
    HDSA::Ptr<ROL::Vector<RealT> > dualstate = u_->dual().clone();
    Update_Constraint(con_,*u_,*z_,theta);
    con_->solve(*dualstate,*u_,*z_,tol_);
    dualstate->set(*u_);
    obj_->gradient_1(*dualstate,*u_,*z_,tol_);
    con_->applyInverseAdjointJacobian_1(*adjoint_nominal_,*dualstate,*u_,*z_,tol_);
    adjoint_nominal_->scale(static_cast<RealT>(-1.0));
    
    theta_grad_nominal_ = theta_->Clone();
    applyThetaJacobianTranspose(*theta_grad_nominal_, *adjoint_nominal_, *u_, *z_, *theta_);
  }

};


#endif
