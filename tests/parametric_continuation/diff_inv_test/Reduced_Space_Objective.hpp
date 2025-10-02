#ifndef REDUCED_SPACE_OBJECTIVE_HPP
#define REDUCED_SPACE_OBJECTIVE_HPP

template <class RealT>
class Reduced_Space_Objective
{

public:
  HDSA::Ptr<Adv_Diff_Constraint<RealT>> con_;
  HDSA::Ptr<Prior_and_Likelihood<RealT>> prior_and_like_;
  HDSA::Ptr<HDSA::Vector<RealT>> current_u_;
  HDSA::Ptr<HDSA::Vector<RealT>> current_z_;
  HDSA::Ptr<HDSA::Vector<RealT>> current_lambda_;
  HDSA::Ptr<HDSA::Vector<RealT>> current_theta_;

  Reduced_Space_Objective(HDSA::Ptr<Adv_Diff_Constraint<RealT>> &con, HDSA::Ptr<Prior_and_Likelihood<RealT>> &prior_and_like, HDSA::Ptr<HDSA::Vector<RealT>> &u_vec, HDSA::Ptr<HDSA::Vector<RealT>> &z_vec, HDSA::Ptr<HDSA::Vector<RealT>> &theta_vec)
      : con_(con), prior_and_like_(prior_and_like)
  {
    current_u_ = u_vec->Clone();
    current_z_ = z_vec->Clone();
    current_z_->Set_Scalar(1.e10);
    current_theta_ = theta_vec->Clone();
    current_lambda_ = u_vec->Clone();
  }

  virtual ~Reduced_Space_Objective()
  {
  }

  // Objective function J(u,z) = (1/2)*(O*u-d)^T*Sigma^{-1}*(O*u-d) + (1/2)*z^T*Gamma^{-1}*z
  // Where O is the observation operator, and Sigma is the noise covariance, and Gamma is the prior covariance

  void Update(const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta)
  {
    HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = z.Clone();
    z_tmp->Set(*current_z_);
    z_tmp->axpy(-1.0, z);
    RealT val = z_tmp->Norm();
    if (val != 0.0)
    {
      Gradient(*z_tmp, z, theta);
    }
  }

  void Gradient(HDSA::Vector<RealT> &grad, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta)
  {
    HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = z.Clone();
    z_tmp->Set(*current_z_);
    z_tmp->axpy(-1.0, z);
    RealT val = z_tmp->Norm();

    if (val != 0.0)
    {
      current_z_->Set(z);
      current_theta_->Set(theta);
      con_->State_Solve(*current_u_, *current_z_, *current_theta_);

      int data_dim = prior_and_like_->data_.size();
      std::vector<RealT> d1 = std::vector<RealT>(data_dim);
      prior_and_like_->Apply_Observation_Operator(d1, *current_u_);
      std::vector<RealT> d2 = std::vector<RealT>(data_dim);
      for (int k = 0; k < data_dim; k++)
      {
        d1[k] = d1[k] - prior_and_like_->data_[k];
      }
      prior_and_like_->Apply_Noise_Precision(d2, d1);
      HDSA::Ptr<HDSA::Vector<RealT>> grad_u = current_u_->Clone();
      prior_and_like_->Apply_Observation_Operator_Transpose(*grad_u, d2);
      con_->c_u_Transpose_Inverse_Apply(*current_lambda_, *grad_u, *current_u_, *current_z_, *current_theta_);
      current_lambda_->Scale(-1.0);
    }

    con_->c_z_Transpose_Apply(grad, *current_lambda_, *current_u_, *current_z_, *current_theta_);

    HDSA::Ptr<HDSA::Vector<RealT>> ztmp1 = grad.Clone();
    HDSA::Ptr<HDSA::Vector<RealT>> ztmp2 = grad.Clone();
    ztmp1->Set_Scalar(-1.0);
    ztmp1->Plus(z);
    prior_and_like_->Apply_Prior_Precision(*ztmp2, *ztmp1);
    grad.Plus(*ztmp2);
  }

  void Apply_Misfit_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta)
  {
    Update(z, theta);

    // Incremental state solve
    HDSA::Ptr<HDSA::Vector<RealT>> w = current_u_->Clone();
    con_->c_z_Apply(*w, z_in, *current_u_, *current_z_, *current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT>> mu = current_u_->Clone();
    w->Scale(-1.0);
    con_->c_u_Inverse_Apply(*mu, *w, *current_u_, *current_z_, *current_theta_);

    // Apply state misfit Hessian
    int data_dim = prior_and_like_->data_.size();
    std::vector<RealT> d1 = std::vector<RealT>(data_dim);
    prior_and_like_->Apply_Observation_Operator(d1, *mu);
    std::vector<RealT> d2 = std::vector<RealT>(data_dim);
    prior_and_like_->Apply_Noise_Precision(d2, d1);
    HDSA::Ptr<HDSA::Vector<RealT>> yJ = current_u_->Clone();
    prior_and_like_->Apply_Observation_Operator_Transpose(*yJ, d2);

    // Incremental adjoint solve
    HDSA::Ptr<HDSA::Vector<RealT>> yc = current_u_->Clone();
    con_->c_uz_Apply(*yc, z_in, *current_u_, *current_z_, *current_lambda_, *current_theta_);
    yc->Plus(*yJ);
    HDSA::Ptr<HDSA::Vector<RealT>> gamma = current_u_->Clone();
    con_->c_u_Transpose_Inverse_Apply(*gamma, *yc, *current_u_, *current_z_, *current_theta_);
    gamma->Scale(-1.0);

    z_out.Zeros();
    con_->c_z_Transpose_Apply(z_out, *gamma, *current_u_, *current_z_, *current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT>> ztmp = z_out.Clone();
    con_->c_zu_Apply(*ztmp, *mu, *current_u_, *current_z_, *current_lambda_, *current_theta_);
    z_out.Plus(*ztmp);

    // I intentionally omitted some terms that I know are zero for this specific problem. If this code is adapted to other problems it will require additional matvecs
  }

  void Apply_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta)
  {
    Apply_Misfit_Hessian(z_out, z_in, z, theta);
    HDSA::Ptr<HDSA::Vector<RealT>> ztmp = z_out.Clone();
    prior_and_like_->Apply_Prior_Precision(*ztmp, z_in);
    z_out.Plus(*ztmp);
  }

  void Apply_B(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &theta_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta)
  {
    Update(z, theta);
    z_out.Zeros();

    // Incremental state solve
    HDSA::Ptr<HDSA::Vector<RealT>> w = current_u_->Clone();
    con_->c_theta_Apply(*w, theta_in, *current_u_, *current_z_, *current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT>> xi = current_u_->Clone();
    w->Scale(-1.0);
    con_->c_u_Inverse_Apply(*xi, *w, *current_u_, *current_z_, *current_theta_);

    // Apply state misfit Hessian
    int data_dim = prior_and_like_->data_.size();
    std::vector<RealT> d1 = std::vector<RealT>(data_dim);
    prior_and_like_->Apply_Observation_Operator(d1, *xi);
    std::vector<RealT> d2 = std::vector<RealT>(data_dim);
    prior_and_like_->Apply_Noise_Precision(d2, d1);
    HDSA::Ptr<HDSA::Vector<RealT>> yJ = current_u_->Clone();
    prior_and_like_->Apply_Observation_Operator_Transpose(*yJ, d2);

    // Incremental adjoint solve
    HDSA::Ptr<HDSA::Vector<RealT>> yc = current_u_->Clone();
    con_->c_utheta_Apply(*yc, theta_in, *current_u_, *current_z_, *current_lambda_, *current_theta_);
    yc->Plus(*yJ);
    HDSA::Ptr<HDSA::Vector<RealT>> beta = current_u_->Clone();
    con_->c_u_Transpose_Inverse_Apply(*beta, *yc, *current_u_, *current_z_, *current_theta_);
    beta->Scale(-1.0);

    z_out.Zeros();
    con_->c_z_Transpose_Apply(z_out, *beta, *current_u_, *current_z_, *current_theta_);
    HDSA::Ptr<HDSA::Vector<RealT>> ztmp = z_out.Clone();
    con_->c_zu_Apply(*ztmp, *xi, *current_u_, *current_z_, *current_lambda_, *current_theta_);
    z_out.Plus(*ztmp);

    // I intentionally omitted some terms that I know are zero for this specific problem. If this code is adapted to other problems it will require additional matvecs
  }
};

#endif
