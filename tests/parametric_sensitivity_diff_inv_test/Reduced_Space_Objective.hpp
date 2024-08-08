#ifndef REDUCED_SPACE_OBJECTIVE_HPP
#define REDUCED_SPACE_OBJECTIVE_HPP

template <class RealT>
class Reduced_Space_Objective {

public: 
  HDSA::Ptr<Adv_Diff_Constraint<RealT> > con_;
  HDSA::Ptr<Prior_and_Likelihood<RealT> > prior_and_like_;
  HDSA::Ptr<HDSA::Vector<RealT> > current_u_;
  HDSA::Ptr<HDSA::Vector<RealT> > current_z_;
  HDSA::Ptr<HDSA::Vector<RealT> > current_lambda_;
  HDSA::Ptr<HDSA::Vector<RealT> > current_theta_;
  
  Reduced_Space_Objective(HDSA::Ptr<Adv_Diff_Constraint<RealT> > & con, HDSA::Ptr<Prior_and_Likelihood<RealT> > & prior_and_like, HDSA::Ptr<HDSA::Vector<RealT> > & u_vec, HDSA::Ptr<HDSA::Vector<RealT> > & z_vec, HDSA::Ptr<HDSA::Vector<RealT> > & theta_vec)
    : con_(con), prior_and_like_(prior_and_like)
  {
    current_u_ = u_vec->clone();
    current_z_ = z_vec->clone();
    current_z_->setScalar(1.e10);
    current_theta_ = theta_vec->clone();
    current_lambda_ = u_vec->clone();
  }

  virtual ~Reduced_Space_Objective()
  { }

  // Objective function J(u,z) = (1/2)*(O*u-d)^T*Sigma^{-1}*(O*u-d) + (1/2)*z^T*Gamma^{-1}*z
  // Where O is the observation operator, and Sigma is the noise covariance, and Gamma is the prior covariance

  void Update(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z.clone();
    z_tmp->set(*current_z_);
    z_tmp->axpy(-1.0,z);
    RealT val = z_tmp->norm();
    if(val != 0.0)
      {
	Gradient(*z_tmp,z,theta);
      }
  }
  
  void Gradient(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z.clone();
    z_tmp->set(*current_z_);
    z_tmp->axpy(-1.0,z);
    RealT val = z_tmp->norm();
    if(val != 0.0)
      {
	current_z_->set(z);
	current_theta_->set(theta);
	con_->State_Solve(*current_u_,*current_z_,*current_theta_);

	int data_dim = prior_and_like_->data_.size();
	std::vector<RealT> d1 = std::vector<RealT>(data_dim);
	prior_and_like_->Apply_Observation_Operator(d1,*current_u_);
	std::vector<RealT> d2 = std::vector<RealT>(data_dim);
	for(int k = 0; k < data_dim; k++)
	  {
	    d1[k] = d1[k] - prior_and_like_->data_[k];
	  }
	prior_and_like_->Apply_Noise_Precision(d2,d1);
	HDSA::Ptr<HDSA::Vector<RealT> > grad_u = current_u_->clone();
	prior_and_like_->Apply_Observation_Operator_Transpose(*grad_u,d2);
	con_->c_u_Transpose_Inverse_Apply(*current_lambda_,*grad_u,*current_u_,*current_z_,*current_theta_);
	current_lambda_->scale(-1.0);
      }
    con_->c_z_Transpose_Apply(grad,*current_lambda_,*current_u_,*current_z_,*current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT> > ztmp = grad.clone();
    prior_and_like_->Apply_Prior_Precision(*ztmp,z);
    grad.plus(*ztmp);
  }

  void Apply_Misfit_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta)
  {
    Update(z,theta);

    // Incremental state solve
    HDSA::Ptr<HDSA::Vector<RealT> > w = current_u_->clone();
    con_->c_z_Apply(*w,z_in,*current_u_,*current_z_,*current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT> > mu = current_u_->clone();
    w->scale(-1.0);
    con_->c_u_Inverse_Apply(*mu,*w,*current_u_,*current_z_,*current_theta_);

    // Apply state misfit Hessian
    int data_dim = prior_and_like_->data_.size();
    std::vector<RealT> d1 = std::vector<RealT>(data_dim);
    prior_and_like_->Apply_Observation_Operator(d1,*mu);
    std::vector<RealT> d2 = std::vector<RealT>(data_dim);
    prior_and_like_->Apply_Noise_Precision(d2,d1);
    HDSA::Ptr<HDSA::Vector<RealT> > yJ = current_u_->clone();
    prior_and_like_->Apply_Observation_Operator_Transpose(*yJ,d2);

    // Incremental adjoint solve
    HDSA::Ptr<HDSA::Vector<RealT> > yc = current_u_->clone();
    con_->c_uz_Apply(*yc,z_in,*current_u_,*current_z_,*current_lambda_,*current_theta_);
    yc->plus(*yJ);
    HDSA::Ptr<HDSA::Vector<RealT> > gamma = current_u_->clone();
    con_->c_u_Transpose_Inverse_Apply(*gamma,*yc,*current_u_,*current_z_,*current_theta_);
    gamma->scale(-1.0);

    z_out.zeros();
    con_->c_z_Transpose_Apply(z_out,*gamma,*current_u_,*current_z_,*current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT> > ztmp = z_out.clone();
    con_->c_zu_Apply(*ztmp,*mu,*current_u_,*current_z_,*current_lambda_,*current_theta_);
    z_out.plus(*ztmp);

    // I intentionally omitted some terms that I know are zero for this specific problem. If this code is adapted to other problems it will require additional matvecs
  }

  void Apply_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta)
  {
    Apply_Misfit_Hessian(z_out,z_in,z,theta);
    HDSA::Ptr<HDSA::Vector<RealT> > ztmp = z_out.clone();
    prior_and_like_->Apply_Prior_Precision(*ztmp,z_in);
    z_out.plus(*ztmp);
  }
  
  void Apply_B(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & theta_in, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta)
  {
    Update(z,theta);
    z_out.zeros();

    // Incremental state solve
    HDSA::Ptr<HDSA::Vector<RealT> > w = current_u_->clone();
    con_->c_theta_Apply(*w,theta_in,*current_u_,*current_z_,*current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT> > xi = current_u_->clone();
    w->scale(-1.0);
    con_->c_u_Inverse_Apply(*xi,*w,*current_u_,*current_z_,*current_theta_);

    // Apply state misfit Hessian
    int data_dim = prior_and_like_->data_.size();
    std::vector<RealT> d1 = std::vector<RealT>(data_dim);
    prior_and_like_->Apply_Observation_Operator(d1,*xi);
    std::vector<RealT> d2 = std::vector<RealT>(data_dim);
    prior_and_like_->Apply_Noise_Precision(d2,d1);
    HDSA::Ptr<HDSA::Vector<RealT> > yJ = current_u_->clone();
    prior_and_like_->Apply_Observation_Operator_Transpose(*yJ,d2);

    // Incremental adjoint solve
    HDSA::Ptr<HDSA::Vector<RealT> > yc = current_u_->clone();
    con_->c_utheta_Apply(*yc,theta_in,*current_u_,*current_z_,*current_lambda_,*current_theta_);
    yc->plus(*yJ);
    HDSA::Ptr<HDSA::Vector<RealT> > beta = current_u_->clone();
    con_->c_u_Transpose_Inverse_Apply(*beta,*yc,*current_u_,*current_z_,*current_theta_);
    beta->scale(-1.0);

    z_out.zeros();
    con_->c_z_Transpose_Apply(z_out,*beta,*current_u_,*current_z_,*current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT> > ztmp = z_out.clone();
    con_->c_zu_Apply(*ztmp,*xi,*current_u_,*current_z_,*current_lambda_,*current_theta_);
    z_out.plus(*ztmp);

    // I intentionally omitted some terms that I know are zero for this specific problem. If this code is adapted to other problems it will require additional matvecs
  }
  
};

#endif
