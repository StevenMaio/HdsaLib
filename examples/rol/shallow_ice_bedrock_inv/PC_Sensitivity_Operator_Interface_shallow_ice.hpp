#ifndef HDSA_PC_SENSITIVITY_OPERATOR_INTERFACE_SHALLOW_ICE_HPP
#define HDSA_PC_SENSITIVITY_OPERATOR_INTERFACE_SHALLOW_ICE_HPP

#include "HDSA_Std_Vector.hpp"
#include "HDSA_ROL_Vector.hpp"

template <class RealT>
class PC_Sensitivity_Operator_Interface_shallow_ice : public HDSA::PC_Sensitivity_Operator_Interface<RealT>
{

private:
  HDSA::Ptr<ROL::Objective<RealT>> robj_;
  HDSA::Ptr<HDSA::Vector<RealT>> z_current_;
  HDSA::Ptr<HDSA::Vector<RealT>> theta_current_;

public:
  PC_Sensitivity_Operator_Interface_shallow_ice(HDSA::Ptr<ROL::Objective<RealT>> &robj, HDSA::Ptr<HDSA::Vector<RealT>> &z_vec, HDSA::Ptr<HDSA::Vector<RealT>> theta_vec) : robj_(robj)
  {
    z_current_ = z_vec->clone();
    theta_current_ = theta_vec->clone();
  }

  virtual ~PC_Sensitivity_Operator_Interface_shallow_ice()
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
      const Std_Vector<RealT> &theta_std = dynamic_cast<const Std_Vector<RealT> &>(theta);
      robj_->setParameter(*theta_std.get_std_vec());
      z_current_->set(z);
      theta_current_->set(theta);
    }
  }

  void Gradient(HDSA::Vector<RealT> &grad, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {

    RealT tol = 1.e-8;
    const HDSA::ROL_Vector<RealT> &z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z);
    HDSA::ROL_Vector<RealT> &grad_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(grad);

    Update(z, theta);
    robj_->gradient(*grad_rol.rol_vec, *z_rol.rol_vec, tol);
  }

  void Apply_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    // auto start = std::chrono::high_resolution_clock::now();

    RealT tol = 1.e-8;
    const HDSA::ROL_Vector<RealT> &z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z);
    const HDSA::ROL_Vector<RealT> &z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z_in);
    HDSA::ROL_Vector<RealT> &z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(z_out);

    Update(z, theta);
    robj_->hessVec(*z_out_rol.rol_vec, *z_in_rol.rol_vec, *z_rol.rol_vec, tol);

    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "Time taken by Apply_Hessian: " << duration.count() / 1000.0 << " seconds" << std::endl;
  }

  void Apply_B(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &theta_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    // auto start = std::chrono::high_resolution_clock::now();

    RealT tol = 1.e-8;
    const HDSA::ROL_Vector<RealT> &z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z);
    HDSA::ROL_Vector<RealT> &z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(z_out);
    HDSA::Ptr<HDSA::Vector<RealT>> grad_nom = z_out.clone();
    HDSA::ROL_Vector<RealT> &grad_nom_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(*grad_nom);

    Update(z, theta);
    robj_->gradient(*grad_nom_rol.rol_vec, *z_rol.rol_vec, tol);

    RealT h = 1.e-4;
    HDSA::Ptr<HDSA::Vector<RealT>> theta_pert = theta.clone();
    theta_pert->set(theta);
    theta_pert->axpy(h, theta_in);

    Update(z, *theta_pert);
    robj_->gradient(*z_out_rol.rol_vec, *z_rol.rol_vec, tol);

    z_out.axpy(-1.0, *grad_nom);
    z_out.scale(1.0 / h);

    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "Time taken by Apply_B: " << duration.count() / 1000.0 << " seconds" << std::endl;
  }
};

#endif
