#ifndef HDSA_ROL_RS_OBJECTIVE_MODEL_ERROR_HPP
#define HDSA_ROL_RS_OBJECTIVE_MODEL_ERROR_HPP

#include "ROL_Objective.hpp"

template <class RealT>
class ROL_RS_Objective_Model_Error: public HDSA::RS_Objective<RealT> {

private:
  HDSA::Ptr<ROL::Objective<RealT> > obj_;
  RealT tol_;

public:

  ROL_RS_Objective_Model_Error(const HDSA::Ptr<ROL::Objective<RealT> > & obj): obj_(obj)
  { 
    tol_ = 1.e-8;
  }

  ~ROL_RS_Objective_Model_Error()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	obj_->update(*z_rol);
      }
    return obj_->value(*z_rol,tol_);
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    const ROL_Vector<RealT> &egrad = dynamic_cast<const ROL_Vector<RealT>&>(grad);
    HDSA::Ptr<ROL::Vector<RealT> > grad_rol = egrad.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
	obj_->update(*z_rol);
      }
    obj_->gradient(*grad_rol,*z_rol,tol_);
  } 
 
  // evaluate the z,z hessian vector product
  void hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta,
  		   const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  { 
    const ROL_Vector<RealT> &ehv = dynamic_cast<const ROL_Vector<RealT>&>(hv);
    HDSA::Ptr<ROL::Vector<RealT> > hv_rol = ehv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    if(update)
      {
  	obj_->update(*z_rol);
      }
    obj_->hessVec(*hv_rol,*v_rol,*z_rol,tol_);
  }

  HDSA::Ptr<ROL::Objective<RealT> > get_objective_function(void) const
  {
    return obj_;
  }  
 
};


#endif
