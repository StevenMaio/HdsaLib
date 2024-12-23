#ifndef MISFIT_REGULARIZATION_OBJECTIVE_SIMOPT_TS_HPP
#define MISFIT_REGULARIZATION_OBJECTIVE_SIMOPT_TS_HPP

template <class Real>
class Misfit_Regularization_Objective_SimOpt_TS : public Objective_SimOpt_TS<Real> {
private:
  ROL::Ptr<Objective_SimOpt_TS<Real> > misfit_obj_;
  std::vector<ROL::Ptr<ROL::Objective_SimOpt<Real> > > reg_obj_;
  std::vector<Real> weights_;
  int num_reg_;
  
public:
  Misfit_Regularization_Objective_SimOpt_TS(std::vector<ROL::TimeStamp<Real> > & timeStamp, ROL::Ptr<Objective_SimOpt_TS<Real> > & misfit_obj, 
					    std::vector<ROL::Ptr<ROL::Objective_SimOpt<Real> > > & reg_obj, std::vector<Real> & weights)
    : Objective_SimOpt_TS<Real>(timeStamp), misfit_obj_(misfit_obj), reg_obj_(reg_obj), weights_(weights)
  { 
    num_reg_ = reg_obj.size();
    if(num_reg_ != (int)weights_.size()-1)
      {
	std::cout << "Weights and regularization objectives are inconsistent" << std::endl;
      }
  }

  Misfit_Regularization_Objective_SimOpt_TS(std::vector<ROL::TimeStamp<Real> > & timeStamp, ROL::Ptr<Objective_SimOpt_TS<Real> > & misfit_obj, 
					    ROL::Ptr<ROL::Objective_SimOpt<Real> > & reg_obj, std::vector<Real> & weights)
    : Objective_SimOpt_TS<Real>(timeStamp), misfit_obj_(misfit_obj), weights_(weights)
  { 
    reg_obj_.resize(1);
    reg_obj_[0] = reg_obj;
    num_reg_ = 1;
    if(num_reg_ != (int)weights_.size()-1)
      {
	std::cout << "Weights and regularization objectives are inconsistent" << std::endl;
      }
  }
  
  void update( const ROL::Vector<Real> & u, const ROL::Vector<Real> & z, Real & t, bool flag = true, int iter = -1 ) 
  {
    misfit_obj_->update(u,z,t,flag,iter);
    for(int k = 0; k < num_reg_; k++)
      {
	reg_obj_[k]->update(u,z,flag,iter);
      }
  }

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol )
  {
    Real val = weights_[0]*misfit_obj_->value(u,z,tol);
    for(int k = 0; k < num_reg_; k++)
      {
	val += weights_[1+k]*reg_obj_[k]->value(u,z,tol);
      }
    return val;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol )
  {
    ROL::Ptr<ROL::Vector<Real> > g_clone = g.clone();
    misfit_obj_->gradient_1(*g_clone,u,z,tol);
    g.set(*g_clone);
    g.scale(weights_[0]);
    for(int k = 0; k < num_reg_; k++)
      {
	g_clone->zero();
	reg_obj_[k]->gradient_1(*g_clone,u,z,tol);
	g.axpy(weights_[1+k],*g_clone);
      }
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol )
  {
    ROL::Ptr<ROL::Vector<Real> > g_clone = g.clone();
    misfit_obj_->gradient_2(*g_clone,u,z,tol);
    g.set(*g_clone);
    g.scale(weights_[0]);
   for(int k = 0; k < num_reg_; k++)
      {
	g_clone->zero();
	reg_obj_[k]->gradient_2(*g_clone,u,z,tol);
	g.axpy(weights_[1+k],*g_clone);
      }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
		   const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol )
  {
    ROL::Ptr<ROL::Vector<Real> > hv_clone = hv.clone();
    misfit_obj_->hessVec_11(*hv_clone,v,u,z,tol);
    hv.set(*hv_clone);
    hv.scale(weights_[0]);
   for(int k = 0; k < num_reg_; k++)
      {
	hv_clone->zero();
	reg_obj_[k]->hessVec_11(*hv_clone,v,u,z,tol);
	hv.axpy(weights_[1+k],*hv_clone);
      }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
		   const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol )
  {
    ROL::Ptr<ROL::Vector<Real> > hv_clone = hv.clone();
    misfit_obj_->hessVec_12(*hv_clone,v,u,z,tol);
    hv.set(*hv_clone);
    hv.scale(weights_[0]);
   for(int k = 0; k < num_reg_; k++)
     {
       hv_clone->zero();
       reg_obj_[k]->hessVec_12(*hv_clone,v,u,z,tol);
       hv.axpy(weights_[1+k],*hv_clone);
     }
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
		   const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol )
  {
    ROL::Ptr<ROL::Vector<Real> > hv_clone = hv.clone();
    misfit_obj_->hessVec_21(*hv_clone,v,u,z,tol);
    hv.set(*hv_clone);
    hv.scale(weights_[0]);
   for(int k = 0; k < num_reg_; k++)
      {
	hv_clone->zero();
	reg_obj_[k]->hessVec_21(*hv_clone,v,u,z,tol);
	hv.axpy(weights_[1+k],*hv_clone);
      }
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
		   const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol )
  {
    ROL::Ptr<ROL::Vector<Real> > hv_clone = hv.clone();
    misfit_obj_->hessVec_22(*hv_clone,v,u,z,tol);
    hv.set(*hv_clone);
    hv.scale(weights_[0]);
   for(int k = 0; k < num_reg_; k++)
      {
	hv_clone->zero();
	reg_obj_[k]->hessVec_22(*hv_clone,v,u,z,tol);
	hv.axpy(weights_[1+k],*hv_clone);
      }
  }

  void setParameter(const std::vector<Real> &param) 
  {
    misfit_obj_->setParameter(param);
    for(int k = 0; k < num_reg_; k++)
      {
	reg_obj_[k]->setParameter(param);
      }
  }

};


#endif
