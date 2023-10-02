#ifndef HDSA_ROL_MD_INTERFACE_HPP
#define HDSA_ROL_MD_INTERFACE_HPP

namespace HDSA
{

template <class RealT>
class ROL_Model_Discrepancy_Interface : public HDSA::Model_Discrepancy_Interface_Elliptic_Prior<RealT> {

private:
  ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj_simopt_;
  ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con_simopt_;
  ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > red_obj_;

public:

  ROL_Model_Discrepancy_Interface(ROL::Ptr<ROL::Objective_SimOpt<RealT> > & obj_simopt, ROL::Ptr<ROL::Constraint_SimOpt<RealT> > & con_simopt,
				  ROL::Ptr<ROL::Vector<RealT> > & u, ROL::Ptr<ROL::Vector<RealT> > & z): obj_simopt_(obj_simopt), con_simopt_(con_simopt)
  {  
    ROL::Ptr<ROL::Vector<RealT> > p = u->clone();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > red_obj_ = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj_simopt,con_simopt,u,z,p);
  }
  
  virtual ~ROL_Model_Discrepancy_Interface()
  { }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define on a problem-to-problem basis (from the base class)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Apply_u_Elliptic_Operator_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_u_Elliptic_Operator_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_u_Mass_Mat(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_u_Mass_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const =0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const = 0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const = 0;

  virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const = 0;

  virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions (from the base class) that are implemented using the SimOpt interface
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const 
  {
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    const HDSA::ROL_Vector<RealT>& u_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u_in);
    const HDSA::ROL_Vector<RealT>& z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z);
    RealT tol = 1.e-8;
    ROL::Ptr<ROL::Vector<RealT> > u_tmp = u_in_rol.rol_vec->clone();
    ROL::Ptr<ROL::Vector<RealT> > u_rol_vec = u_in_rol.rol_vec->clone();
    ROL::Ptr<ROL::Vector<RealT> > c_rol_vec = u_in_rol.rol_vec->clone();
    con_simopt_->solve(*c_rol_vec, *u_rol_vec, *z_rol.rol_vec, tol);
    con_simopt_->applyInverseAdjointJacobian_1(*u_tmp, *u_in_rol.rol_vec, *u_rol_vec, *z_rol.rol_vec, tol);
    con_simopt_->applyAdjointJacobian_2(*z_out_rol.rol_vec, *u_tmp, *u_rol_vec, *z_rol.rol_vec, tol);
  }

  void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const 
  {
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    const HDSA::ROL_Vector<RealT>& z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z);
    RealT tol = 1.e-8;
    red_obj_->hessVec(*z_out_rol.rol_vec, *z_in_rol.rol_vec, *z_rol.rol_vec, tol);
  }

  void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const 
  {
    HDSA::ROL_Vector<RealT>& u_grad_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(u_grad);
    const HDSA::ROL_Vector<RealT>& u_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u);
    const HDSA::ROL_Vector<RealT>& z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z);
    RealT tol = 1.e-8;
    obj_simopt_->gradient_1(*u_grad_rol.rol_vec, *u_rol.rol_vec, *z_rol.rol_vec, tol);
  }

  void Apply_Misfit_Hessian(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const 
  {
    HDSA::ROL_Vector<RealT>& u_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(u_out);
    const HDSA::ROL_Vector<RealT>& u_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u_in);
    const HDSA::ROL_Vector<RealT>& u_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u);
    const HDSA::ROL_Vector<RealT>& z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z);
    RealT tol = 1.e-8;
    obj_simopt_->hessVec_11(*u_out_rol.rol_vec, *u_in_rol.rol_vec, *u_rol.rol_vec, *z_rol.rol_vec, tol);
  }

};

}

#endif


