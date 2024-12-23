#ifndef HDSA_PC_LIS_INTERFACE_SHALLOW_ICE_HPP
#define HDSA_PC_LIS_INTERFACE_SHALLOW_ICE_HPP


template <class RealT>
class PC_LIS_Interface_shallow_ice : public HDSA::PC_LIS_Interface<RealT> {
  
private:
  HDSA::Ptr<ROL::Objective<RealT> > robj_misfit_;
  HDSA::Ptr<ROL::Objective<RealT> > robj_reg_;
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > reg_obj_;
  HDSA::Ptr<HDSA::Vector<RealT> > z_current_;
  HDSA::Ptr<HDSA::Vector<RealT> > theta_current_;
  
public:
  PC_LIS_Interface_shallow_ice(HDSA::Ptr<ROL::Objective<RealT> > & robj_misfit, HDSA::Ptr<ROL::Objective<RealT> > & robj_reg, HDSA::Ptr<ROL::Objective_SimOpt<RealT> > & reg_obj, HDSA::Ptr<HDSA::Vector<RealT> > & z_vec, HDSA::Ptr<HDSA::Vector<RealT> > theta_vec):
    robj_misfit_(robj_misfit), robj_reg_(robj_reg), reg_obj_(reg_obj)
  {
    z_current_ = z_vec->clone();
    theta_current_ = theta_vec->clone();
  }
  
  virtual ~PC_LIS_Interface_shallow_ice()
  { }
  
  
  void Apply_Misfit_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta) const
  {
    std::cout << "Applying misfit Hessian-vector product" << std::endl;
    
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z_current_->clone();
    z_tmp->set(z);
    z_tmp->axpy(-1.0,*z_current_);
    HDSA::Ptr<HDSA::Vector<RealT> > theta_tmp = theta_current_->clone();
    theta_tmp->set(theta);
    theta_tmp->axpy(-1.0,*theta_current_);
    RealT val = z_tmp->norm() + theta_tmp->norm();
    if(val > 0.0)
      {
	const Std_Vector<RealT>& theta_std = dynamic_cast<const Std_Vector<RealT>&>(theta);
	robj_misfit_->setParameter(*theta_std.get_std_vec());
	z_current_->set(z);
	theta_current_->set(theta);
      }
    
    RealT tol = 1.e-8;
    const HDSA::ROL_Vector<RealT>& z_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z);
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    
    robj_misfit_->hessVec(*z_out_rol.rol_vec,*z_in_rol.rol_vec,*z_rol.rol_vec,tol);
  }
  
  void Apply_Prior_Precision(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
  {
    RealT tol = 1.e-8;
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    robj_reg_->hessVec(*z_out_rol.rol_vec,*z_in_rol.rol_vec,*z_in_rol.rol_vec,tol);
  }
  
  void Apply_Prior_Covariance(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
  {
    RealT tol = 1.e-8;
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);

    Elliptic_Prior_Regularization_Objective<RealT>& reg_obj_cast = dynamic_cast<Elliptic_Prior_Regularization_Objective<RealT>&>(*reg_obj_);
    reg_obj_cast.Apply_Prior_Covariance(*z_out_rol.rol_vec,*z_in_rol.rol_vec,tol);
  }
  
  void Generate_Prior_Samples(HDSA::MultiVector<RealT> & samples) const
  {
    RealT tol = 1.e-8;
    Elliptic_Prior_Regularization_Objective<RealT>& reg_obj_cast = dynamic_cast<Elliptic_Prior_Regularization_Objective<RealT>&>(*reg_obj_);
    HDSA::Ptr<HDSA::Vector<RealT> > z_in = samples[0]->clone();
    HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*z_in);

    for(int k = 0; k < samples.Number_of_Vectors(); k++)
      {
	z_in->randomize_standard_normal();
	HDSA::Ptr<HDSA::Vector<RealT> > z_out = samples[k];
	HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*z_out);
	reg_obj_cast.Elliptic_Solve(*z_out_rol.rol_vec,*z_in_rol.rol_vec,tol);
      }
  }
  
};

#endif
