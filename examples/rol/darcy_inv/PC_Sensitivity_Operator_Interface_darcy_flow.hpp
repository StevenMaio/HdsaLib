#ifndef HDSA_PC_SENSITIVITY_OPERATOR_INTERFACE_DARCY_FLOW_HPP
#define HDSA_PC_SENSITIVITY_OPERATOR_INTERFACE_DARCY_FLOW_HPP

#include "HDSA_ROL_Vector.hpp"
#include "HDSA_Std_Vector.hpp"

template <class RealT>
class PC_Sensitivity_Operator_Interface_darcy_flow : public HDSA::PC_Sensitivity_Operator_Interface<RealT>
{

private:
  HDSA::Ptr<ROL::Objective_SimOpt<RealT>> obj_misfit_;
  HDSA::Ptr<ROL::Objective<RealT>> obj_reg_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT>> con_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT>> con_aux_param_;
  HDSA::Ptr<ROL::Vector<RealT>> u_current_;
  HDSA::Ptr<ROL::Vector<RealT>> lambda_current_;
  HDSA::Ptr<HDSA::Vector<RealT>> z_current_;
  HDSA::Ptr<HDSA::Vector<RealT>> theta_current_;

public:
  PC_Sensitivity_Operator_Interface_darcy_flow(HDSA::Ptr<ROL::Objective_SimOpt<RealT>> &obj_misfit, HDSA::Ptr<ROL::Objective<RealT>> &obj_reg,
                                               HDSA::Ptr<ROL::Constraint_SimOpt<RealT>> &con, HDSA::Ptr<ROL::Constraint_SimOpt<RealT>> &con_aux_param, HDSA::Ptr<ROL::Vector<RealT>> &u_vec,
                                               HDSA::Ptr<HDSA::Vector<RealT>> &z_vec, HDSA::Ptr<HDSA::Vector<RealT>> theta_vec, bool fd_check = false)
      : obj_misfit_(obj_misfit), obj_reg_(obj_reg), con_(con), con_aux_param_(con_aux_param)
  {
    u_current_ = u_vec->clone();
    lambda_current_ = u_vec->clone();
    z_current_ = z_vec->clone();
    theta_current_ = theta_vec->clone();

    if (fd_check)
    {
      HDSA::Ptr<HDSA::Vector<RealT>> z = z_vec->clone();
      z->randomize_standard_normal();
      HDSA::Ptr<HDSA::Vector<RealT>> theta = theta_vec->clone();
      theta->randomize_standard_normal();
      Finite_Difference_Checks(*z, *theta);
    }
  }

  virtual ~PC_Sensitivity_Operator_Interface_darcy_flow()
  {
  }

  void Update(const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = z_current_->clone();
    z_tmp->set(z);
    z_tmp->axpy(-1.0, *z_current_);
    HDSA::Ptr<HDSA::Vector<RealT>> theta_tmp = theta_current_->clone();
    theta_tmp->set(theta);
    theta_tmp->axpy(-1.0, *theta_current_);
    RealT val = z_tmp->norm() + theta_tmp->norm();
    if (val > 0.0)
    {
      const HDSA::ROL_Vector<RealT> &z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z);
      const HDSA::Std_Vector<RealT> &theta_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta);
      con_->setParameter(*theta_std.get_std_vec());
      z_current_->set(z);
      theta_current_->set(theta);

      RealT tol = 1.e-14;
      ROL::Ptr<ROL::Vector<RealT>> c = u_current_->clone();
      con_->solve(*c, *u_current_, *z_rol.rol_vec, tol);

      obj_misfit_->gradient_1(*c, *u_current_, *z_rol.rol_vec, tol);
      con_->applyInverseAdjointJacobian_1(*lambda_current_, *c, *u_current_, *z_rol.rol_vec, tol);
      lambda_current_->scale(static_cast<RealT>(-1));
    }
  }

  void Gradient(HDSA::Vector<RealT> &grad, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    RealT tol = 1.e-14;
    const HDSA::ROL_Vector<RealT> &z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z);
    HDSA::ROL_Vector<RealT> &grad_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(grad);
    Update(z, theta);

    ROL::Ptr<ROL::Vector<RealT>> z_tmp = z_rol.rol_vec->clone();
    con_->applyAdjointJacobian_2(*z_tmp, *lambda_current_, *u_current_, *z_rol.rol_vec, tol);

    obj_reg_->gradient(*grad_rol.rol_vec, *z_rol.rol_vec, tol);
    grad_rol.rol_vec->plus(*z_tmp);
  }

  void Apply_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {

    RealT tol = 1.e-14;
    const HDSA::ROL_Vector<RealT> &z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z);
    const HDSA::ROL_Vector<RealT> &z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z_in);
    HDSA::ROL_Vector<RealT> &z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(z_out);
    Update(z, theta);

    ROL::Ptr<ROL::Vector<RealT>> forcing = u_current_->clone();
    con_->applyJacobian_2(*forcing, *z_in_rol.rol_vec, *u_current_, *z_rol.rol_vec, tol);
    forcing->scale(-1.0);
    ROL::Ptr<ROL::Vector<RealT>> incr_state = u_current_->clone();
    con_->applyInverseJacobian_1(*incr_state, *forcing, *u_current_, *z_rol.rol_vec, tol);

    ROL::Ptr<ROL::Vector<RealT>> d = u_current_->clone();
    obj_misfit_->hessVec_11(*d, *incr_state, *u_current_, *z_rol.rol_vec, tol);
    con_->applyAdjointHessian_11(*forcing, *lambda_current_, *incr_state, *u_current_, *z_rol.rol_vec, tol);
    d->plus(*forcing);
    con_->applyAdjointHessian_21(*forcing, *lambda_current_, *z_in_rol.rol_vec, *u_current_, *z_rol.rol_vec, tol);
    d->plus(*forcing);
    d->scale(-1.0);

    ROL::Ptr<ROL::Vector<RealT>> incr_adjoint = u_current_->clone();
    con_->applyInverseAdjointJacobian_1(*incr_adjoint, *d, *u_current_, *z_rol.rol_vec, tol);

    ROL::Ptr<ROL::Vector<RealT>> z_tmp = z_rol.rol_vec->clone();
    con_->applyAdjointJacobian_2(*z_out_rol.rol_vec, *incr_adjoint, *u_current_, *z_rol.rol_vec, tol);
    con_->applyAdjointHessian_12(*z_tmp, *lambda_current_, *incr_state, *u_current_, *z_rol.rol_vec, tol);
    z_out_rol.rol_vec->plus(*z_tmp);
    con_->applyAdjointHessian_22(*z_tmp, *lambda_current_, *z_in_rol.rol_vec, *u_current_, *z_rol.rol_vec, tol);
    z_out_rol.rol_vec->plus(*z_tmp);
    obj_reg_->hessVec(*z_tmp, *z_in_rol.rol_vec, *z_rol.rol_vec, tol);
    z_out_rol.rol_vec->plus(*z_tmp);
  }

  void Apply_B(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &theta_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    RealT tol = 1.e-14;
    const HDSA::ROL_Vector<RealT> &z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z);
    HDSA::ROL_Vector<RealT> &z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(z_out);
    Update(z, theta);

    const HDSA::Std_Vector<RealT> &theta_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta);
    const HDSA::Std_Vector<RealT> &theta_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta_in);
    ROL::Ptr<ROL::StdVector<RealT>> theta_rol = ROL::makePtr<ROL::StdVector<RealT>>(theta_std.get_std_vec());
    ROL::Ptr<ROL::StdVector<RealT>> theta_in_rol = ROL::makePtr<ROL::StdVector<RealT>>(theta_in_std.get_std_vec());

    ROL::Ptr<ROL::Vector<RealT>> forcing = u_current_->clone();
    con_aux_param_->applyJacobian_2(*forcing, *theta_in_rol, *u_current_, *theta_rol, tol);
    forcing->scale(-1.0);
    ROL::Ptr<ROL::Vector<RealT>> incr_state = u_current_->clone();
    con_->applyInverseJacobian_1(*incr_state, *forcing, *u_current_, *z_rol.rol_vec, tol);

    ROL::Ptr<ROL::Vector<RealT>> d = u_current_->clone();
    obj_misfit_->hessVec_11(*d, *incr_state, *u_current_, *z_rol.rol_vec, tol);
    con_->applyAdjointHessian_11(*forcing, *lambda_current_, *incr_state, *u_current_, *z_rol.rol_vec, tol);
    d->plus(*forcing);
    d->scale(-1.0);

    ROL::Ptr<ROL::Vector<RealT>> incr_adjoint = u_current_->clone();
    con_->applyInverseAdjointJacobian_1(*incr_adjoint, *d, *u_current_, *z_rol.rol_vec, tol);

    ROL::Ptr<ROL::Vector<RealT>> z_tmp = z_rol.rol_vec->clone();
    con_->applyAdjointJacobian_2(*z_tmp, *incr_adjoint, *u_current_, *z_rol.rol_vec, tol);
    con_->applyAdjointHessian_12(*z_out_rol.rol_vec, *lambda_current_, *incr_state, *u_current_, *z_rol.rol_vec, tol);
    z_out_rol.rol_vec->plus(*z_tmp);
  }

  void Finite_Difference_Checks(HDSA::Vector<RealT> &z, HDSA::Vector<RealT> &theta) const
  {
    Finite_Difference_Gradient(z, theta);
    Finite_Difference_Hessian(z, theta);
    Finite_Difference_B(z, theta);
  }

  void Finite_Difference_Gradient(HDSA::Vector<RealT> &z, HDSA::Vector<RealT> &theta) const
  {
    RealT val = Value(z, theta);
    HDSA::Ptr<HDSA::Vector<RealT>> dz = z.clone();
    dz->randomize_standard_normal();
    HDSA::Ptr<HDSA::Vector<RealT>> grad = z.clone();
    Gradient(*grad, z, theta);
    RealT true_grad = grad->dot(*dz);
    HDSA::Ptr<HDSA::Vector<RealT>> z_pert = z.clone();
    int n = 6;
    std::vector<RealT> grad_fd = std::vector<RealT>(n);
    std::vector<RealT> fd_error = std::vector<RealT>(n);
    for (int k = 0; k < n; k++)
    {
      RealT h = std::pow(10.0, static_cast<RealT>(-k - 1));
      z_pert->set(z);
      z_pert->axpy(h, *dz);
      RealT val_pert = Value(*z_pert, theta);
      grad_fd[k] = (val_pert - val) / h;
      fd_error[k] = std::abs(grad_fd[k] - true_grad);
      std::cout << "Step size = " << h << " and gradient FD error = " << fd_error[k] << std::endl;
    }
  }

  void Finite_Difference_Hessian(HDSA::Vector<RealT> &z, HDSA::Vector<RealT> &theta) const
  {
    HDSA::Ptr<HDSA::Vector<RealT>> hv = z.clone();
    HDSA::Ptr<HDSA::Vector<RealT>> v = z.clone();
    v->randomize_standard_normal();

    HDSA::Ptr<HDSA::Vector<RealT>> grad = z.clone();
    Gradient(*grad, z, theta);
    Apply_Hessian(*hv, *v, z, theta);

    HDSA::Ptr<HDSA::Vector<RealT>> z_pert = z.clone();
    HDSA::Ptr<HDSA::Vector<RealT>> grad_pert = z.clone();
    int n = 6;
    std::vector<RealT> fd_error = std::vector<RealT>(n);
    for (int k = 0; k < n; k++)
    {
      RealT h = std::pow(10.0, static_cast<RealT>(-k - 1));
      z_pert->set(z);
      z_pert->axpy(h, *v);
      grad_pert->setScalar(0.0);
      Gradient(*grad_pert, *z_pert, theta);
      grad_pert->axpy(-1.0, *grad);
      fd_error[k] = grad_pert->norm();
      std::cout << "Step size = " << h << " and Hessian FD error = " << fd_error[k] << std::endl;
    }
  }

  void Finite_Difference_B(HDSA::Vector<RealT> &z, HDSA::Vector<RealT> &theta) const
  {
    HDSA::Ptr<HDSA::Vector<RealT>> hv = z.clone();
    HDSA::Ptr<HDSA::Vector<RealT>> v = theta.clone();
    v->randomize_standard_normal();

    HDSA::Ptr<HDSA::Vector<RealT>> grad = z.clone();
    Gradient(*grad, z, theta);
    Apply_B(*hv, *v, z, theta);

    HDSA::Ptr<HDSA::Vector<RealT>> theta_pert = theta.clone();
    HDSA::Ptr<HDSA::Vector<RealT>> grad_pert = z.clone();
    int n = 6;
    std::vector<RealT> fd_error = std::vector<RealT>(n);
    for (int k = 0; k < n; k++)
    {
      RealT h = std::pow(10.0, static_cast<RealT>(-k - 1));
      theta_pert->set(theta);
      theta_pert->axpy(h, *v);
      grad_pert->setScalar(0.0);
      Gradient(*grad_pert, z, *theta_pert);
      grad_pert->axpy(-1.0, *grad);
      fd_error[k] = grad_pert->norm();
      std::cout << "Step size = " << h << " and B FD error = " << fd_error[k] << std::endl;
    }
  }

  RealT Value(const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    RealT tol = 1.e-14;
    const HDSA::ROL_Vector<RealT> &z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z);
    Update(z, theta);

    RealT val1 = obj_misfit_->value(*u_current_, *z_rol.rol_vec, tol);
    RealT val2 = obj_reg_->value(*z_rol.rol_vec, tol);
    RealT val = val1 + val2;

    return val;
  }
};

#endif
