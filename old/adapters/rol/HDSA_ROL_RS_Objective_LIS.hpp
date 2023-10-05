#ifndef HDSA_ROL_RS_OBJECTIVE_LIS_HPP
#define HDSA_ROL_RS_OBJECTIVE_LIS_HPP

#include "ROL_Objective.hpp"

template <class RealT>
class ROL_RS_Objective_LIS: public ROL_RS_Objective<RealT> {

private:
  HDSA::Ptr<ROL::Objective<RealT> > misfit_obj_;
  HDSA::Ptr<ROL::Objective<RealT> > reg_obj_;
  RealT misfit_weight_; 
  RealT reg_weight_;
  RealT tol_;
  
public:

  // It should be that obj = misfit_weight*misfit_obj + reg_weight*reg_objective
  ROL_RS_Objective_LIS(const HDSA::Ptr<ROL::Objective<RealT> > & obj, const HDSA::Ptr<ROL::Objective<RealT> > & misfit_obj, const HDSA::Ptr<ROL::Objective<RealT> > & reg_obj,
		       const RealT & misfit_weight, const RealT & reg_weight, const HDSA::Ptr<ROL_B_Transpose<RealT> > & rol_Bt = HDSA::nullPtr):
    ROL_RS_Objective<RealT>(obj, rol_Bt), misfit_obj_(misfit_obj), reg_obj_(reg_obj)
  {
    misfit_weight_ = misfit_weight;
    reg_weight_ = reg_weight;
    tol_ = 1.e-8;
  }

  ~ROL_RS_Objective_LIS()
  { }

  // evaluate the misfit z,z hessian vector product (for computed likelihood informed subspaces only)
  void Misfit_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
			  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    const ROL_Vector<RealT> &ehv = dynamic_cast<const ROL_Vector<RealT>&>(hv);
    HDSA::Ptr<ROL::Vector<RealT> > hv_rol = ehv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    misfit_obj_->hessVec(*hv_rol,*v_rol,*z_rol,tol_);
    hv.scale(misfit_weight_);
  }
  
  // evaluate the regularization z,z hessian vector product (for computed likelihood informed subspaces only)
  void Regularization_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
				  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    const ROL_Vector<RealT> &ehv = dynamic_cast<const ROL_Vector<RealT>&>(hv);
    HDSA::Ptr<ROL::Vector<RealT> > hv_rol = ehv.get_rol_vec();
    const ROL_Vector<RealT> &ev = dynamic_cast<const ROL_Vector<RealT>&>(v);
    HDSA::Ptr<ROL::Vector<RealT> > v_rol = ev.get_rol_vec();
    const ROL_Vector<RealT> &ez = dynamic_cast<const ROL_Vector<RealT>&>(z);
    HDSA::Ptr<ROL::Vector<RealT> > z_rol = ez.get_rol_vec();
    reg_obj_->hessVec(*hv_rol,*v_rol,*z_rol,tol_);
    hv.scale(reg_weight_);
  }
 
};


#endif
