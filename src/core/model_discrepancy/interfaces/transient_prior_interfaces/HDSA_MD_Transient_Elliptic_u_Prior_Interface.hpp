#ifndef HDSA_MD_TRANSIENT_ELLIPTIC_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_TRANSIENT_ELLIPTIC_U_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Transient_Elliptic_u_Prior_Interface : public HDSA::MD_Scaled_u_Prior_Interface<RealT>
  {

  private:
    HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> spatial_prior_cov_;
    HDSA::Ptr<HDSA::MD_Transient_Prior_Covariance<RealT>> transient_prior_cov_;
    HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> u_hyperparam_interface_;
    HDSA::Ptr<HDSA::MD_Determine_u_Hyperparameters<RealT>> determine_u_hyperparams_;

  public:
    void Apply_M_u(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      // NEED TO IMPLEMENT
    }

    void Apply_W_u_Acute_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const RealT &scalar) const
    {
      // NEED TO IMPLEMENT
    }

    void Apply_W_u_Acute_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      // NEED TO IMPLEMENT
    }

    void Sample_with_Covariance_W_u_Acute_Inverse(HDSA::MultiVector<RealT> &samples) const
    {
      // NEED TO IMPLEMENT
    }

    void Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(HDSA::MultiVector<RealT> &samples, const RealT &scalar) const
    {
      // NEED TO IMPLEMENT
    }

    MD_Transient_Elliptic_u_Prior_Interface(HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &spatial_prior_cov, HDSA::Ptr<HDSA::MD_Transient_Prior_Covariance<RealT>> &transient_prior_cov) : HDSA::MD_Scaled_u_Prior_Interface<RealT>(transient_prior_cov->Get_Hyperparameter_Interface()->Get_alpha_u())
    {
      spatial_prior_cov_ = spatial_prior_cov;
      transient_prior_cov_ = transient_prior_cov;
      u_hyperparam_interface_ = transient_prior_cov_->Get_Hyperparameter_Interface();
      determine_u_hyperparams_ = transient_prior_cov_->Get_Determine_Hyperparameters();

      if (u_hyperparam_interface_->Adapt_Time_Variance())
      {
        determine_u_hyperparams_->Determine_alpha_t(this);
      }
      transient_prior_cov_->Set_alpha_t(u_hyperparam_interface_->Get_alpha_t());

      if (u_hyperparam_interface_->Get_alpha_u() == 0.0)
      {
        determine_u_hyperparams_->Determine_alpha_u(this);
      }
      this->Set_alpha_u(u_hyperparam_interface_->Get_alpha_u());
    }

    virtual ~MD_Transient_Elliptic_u_Prior_Interface()
    {
    }
  };

}

#endif
