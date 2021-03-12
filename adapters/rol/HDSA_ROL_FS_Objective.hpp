#ifndef HDSA_ROL_FS_OBJECTIVE_HPP
#define HDSA_ROL_FS_OBJECTIVE_HPP

#include "ROL_Objective_SimOpt.hpp"

template <class RealT>
class ROL_FS_Objective: public HDSA::FS_Objective<RealT> {

private:
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj_;
  RealT tol_;

public:

  ROL_FS_Objective(const HDSA::Ptr<ROL::Objective_SimOpt<RealT> > & obj): obj_(obj)
  { 
    tol_ = 1.e-8;
  }

  ~ROL_FS_Objective()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	obj_->setParameter(*param);
	obj_->update(*u_rol,*z_rol);
      }
    return obj_->value(*u_rol,*z_rol,tol_);
  }

  // evaluate the gradient with respect to u
  void gradient_u(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &egrad = dynamic_cast<ROL_Vector<RealT>&>(grad);
    HDSA::Ptr<ROL::Vector<RealT> > grad_rol = egrad.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	obj_->setParameter(*param);
	obj_->update(*u_rol,*z_rol);
      }
    obj_->gradient_1(*grad_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    ROL_Vector<RealT> &egrad = dynamic_cast<ROL_Vector<RealT>&>(grad);
    HDSA::Ptr<ROL::Vector<RealT> > grad_rol = egrad.get_rol_vec();
    const ROL_Vector<RealT> &eu = dynamic_cast<const ROL_Vector<RealT>&>(u);
    HDSA::Ptr<ROL::Vector<RealT> > u_rol = eu.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
	HDSA::Ptr<std::vector<RealT> > param = etheta.get_std_vec();
	obj_->setParameter(*param);
	obj_->update(*u_rol,*z_rol);
      }
    obj_->gradient_2(*grad_rol,*u_rol,*z_rol,tol_);
  } 
  
  // evaluate the u,u hessian vector product
  void hessVec_u_u(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
		   const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    ROL_Vector<RealT> &ehv = dynamic_cast<ROL_Vector<RealT>&>(hv);
    HDSA::Ptr<ROL::Vector<RealT> > hv_rol = ehv.get_rol_vec();
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
	obj_->setParameter(*param);
	obj_->update(*u_rol,*z_rol);
      }
    obj_->hessVec_11(*hv_rol,*v_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the u,z hessian vector product
  void hessVec_u_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
		   const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    ROL_Vector<RealT> &ehv = dynamic_cast<ROL_Vector<RealT>&>(hv);
    HDSA::Ptr<ROL::Vector<RealT> > hv_rol = ehv.get_rol_vec();
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
	obj_->setParameter(*param);
	obj_->update(*u_rol,*z_rol);
      }
    obj_->hessVec_12(*hv_rol,*v_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the z,z hessian vector product
  void hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
		   const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    ROL_Vector<RealT> &ehv = dynamic_cast<ROL_Vector<RealT>&>(hv);
    HDSA::Ptr<ROL::Vector<RealT> > hv_rol = ehv.get_rol_vec();
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
	obj_->setParameter(*param);
	obj_->update(*u_rol,*z_rol);
      }
    obj_->hessVec_22(*hv_rol,*v_rol,*u_rol,*z_rol,tol_);
  }

  // evaluate the z,u hessian vector product
  void hessVec_z_u(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
		   const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    ROL_Vector<RealT> &ehv = dynamic_cast<ROL_Vector<RealT>&>(hv);
    HDSA::Ptr<ROL::Vector<RealT> > hv_rol = ehv.get_rol_vec();
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
	obj_->setParameter(*param);
	obj_->update(*u_rol,*z_rol);
      }
    obj_->hessVec_21(*hv_rol,*v_rol,*u_rol,*z_rol,tol_);
  }

  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > get_objective_function(void) const
  {
    return obj_;
  }  
 
};


#endif
