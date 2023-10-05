#ifndef HDSA_ROL_CONSTRAINT_HPP
#define HDSA_ROL_CONSTRAINT_HPP

#include "ROL_Constraint_SimOpt.hpp"

template <class RealT>
class ROL_Constraint: public HDSA::Constraint<RealT> {

private:
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_;
  RealT tol_;

public:

  ROL_Constraint(const HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > & con): con_(con)
  { 
    tol_ = 1.e-8;
  }

  ~ROL_Constraint()
  { }

  // evaluate constraint residual
  void value(HDSA::Vector<RealT> & r, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    const ROL_Vector<RealT> &er = dynamic_cast<const ROL_Vector<RealT>&>(r);
    HDSA::Ptr<ROL::Vector<RealT> > r_rol = er.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->value(*r_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the jacobian with respect to u vector product
  void jacobian_u(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyJacobian_1(*Jv_rol,*v_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the inverse jacobian with respect to u vector product
  void jacobian_u_inverse(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyInverseJacobian_1(*Jv_rol,*v_rol,*u_rol,*z_rol,tol_);
  } 

  // evaluate the jacobian with respect to z vector product
  void jacobian_z(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyJacobian_2(*Jv_rol,*v_rol,*u_rol,*z_rol,tol_);
  } 

  // evaluate the adjoint jacobian with respect to u vector product
  void jacobian_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyAdjointJacobian_1(*Jv_rol,*v_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the inverse adjoint jacobian with respect to u vector product
  void jacobian_u_adjoint_inverse(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,const HDSA::Vector<RealT> & theta, 
				  const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyInverseAdjointJacobian_1(*Jv_rol,*v_rol,*u_rol,*z_rol,tol_);
  } 

  // evaluate the adjoint jacobian with respect to z vector product
  void jacobian_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyAdjointJacobian_2(*Jv_rol,*v_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the adjoint hessian with respect to u,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z) with respect to u,u
  void hessian_u_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
			   const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &elambda = dynamic_cast<const ROL_Vector<RealT>&>(lambda);
    HDSA::Ptr<ROL::Vector<RealT> > lambda_rol = elambda.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyAdjointHessian_11(*Jv_rol,*lambda_rol,*v_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the adjoint hessian with respect to u,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z) with respect to u,z
  void hessian_u_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
			   const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &elambda = dynamic_cast<const ROL_Vector<RealT>&>(lambda);
    HDSA::Ptr<ROL::Vector<RealT> > lambda_rol = elambda.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyAdjointHessian_21(*Jv_rol,*lambda_rol,*v_rol,*u_rol,*z_rol,tol_);
  }   

  // evaluate the adjoint hessian with respect to z,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z) with respect to z,z
  void hessian_z_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
			   const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &elambda = dynamic_cast<const ROL_Vector<RealT>&>(lambda);
    HDSA::Ptr<ROL::Vector<RealT> > lambda_rol = elambda.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyAdjointHessian_22(*Jv_rol,*lambda_rol,*v_rol,*u_rol,*z_rol,tol_);
  }     

  // evaluate the adjoint hessian with respect to z,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z) with respect to z,u
  void hessian_z_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
			   const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &eJv = dynamic_cast<ROL_Vector<RealT>&>(Jv);
    HDSA::Ptr<ROL::Vector<RealT> > Jv_rol = eJv.get_rol_vec();
    const ROL_Vector<RealT> &elambda = dynamic_cast<const ROL_Vector<RealT>&>(lambda);
    HDSA::Ptr<ROL::Vector<RealT> > lambda_rol = elambda.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	con_->setParameter(*param);
	con_->update(*u_rol,*z_rol);
      }
    con_->applyAdjointHessian_12(*Jv_rol,*lambda_rol,*v_rol,*u_rol,*z_rol,tol_);
  }  
  
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > get_constraint(void) const
  {
    return con_;
  }  
 
};


#endif
