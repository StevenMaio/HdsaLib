#ifndef HDSA_MD_U_HYPERPARAMETER_INTERFACE_HPP
#define HDSA_MD_U_HYPERPARAMETER_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_u_Hyperparameter_Interface {

  private:
  bool is_transient_;
  bool center_data_;
  bool adapt_time_variance_;
  int component_id_;

  RealT alpha_u_;
  RealT beta_u_;
  std::vector<RealT> alpha_t_;
  RealT beta_t_;
  RealT alpha_d_;

  int gsvd_num_sing_vals_;
  int gsvd_oversampling_;
  int gsvd_num_subspace_iter_;

  RealT time_variance_inflation_;
  RealT data_noise_percent_;
  RealT W_u_inv_spectral_gap_;

  RealT d1_norm_sq_;
  std::vector<RealT> d_pert_norm_sq_;

  public:

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Spatial_Node_Data(void) const
  {
    std::cout << "Load_Spatial_Node_Data is required for hyperparameter algorithm-based initialization" << std::endl;
    HDSA::Ptr<HDSA::Vector<RealT> > vec;
    return vec;
  }

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Time_Node_Data(void) const
  {
    std::cout << "Load_Time_Node_Data is required for hyperparameter algorithm-based initialization" << std::endl;
    HDSA::Ptr<HDSA::Vector<RealT> > vec;
    return vec;
  }

  MD_u_Hyperparameter_Interface(const bool & is_transient, const bool & center_data = false, const bool & adapt_time_variance = false, const int & component_id = 1):
   is_transient_(is_transient), center_data_(center_data), adapt_time_variance_(adapt_time_variance), component_id_(component_id)
   { 
    alpha_u_ = 0.0;
    beta_u_ = 0.0;
    alpha_t_.resize(1);
    alpha_t_[0] = 1.0;
    beta_t_ = 0.0;
    alpha_d_ = 0.0;

    gsvd_num_sing_vals_ = 0;
    gsvd_oversampling_ = 0;
    gsvd_num_subspace_iter_ = 1;

    time_variance_inflation_ = .01;
    data_noise_percent_ = .001;
    W_u_inv_spectral_gap_ = 1.e-4;
   }

    virtual ~MD_u_Hyperparameter_Interface()
    { }

    void Set_alpha_u(RealT & alpha_u_new)
    {
      alpha_u_ = alpha_u_new;
    }

    void Set_beta_u(RealT & beta_u_new)
    {
      beta_u_ = beta_u_new;
    }

    void Set_alpha_t(std::vector<RealT> & alpha_t_new)
    {
      alpha_t_.clear();
      alpha_t_ = alpha_t_new;
    }

    void Set_beta_t(RealT & beta_t_new)
    {
      beta_t_ = beta_t_new;
    }

    void Set_alpha_d(RealT & alpha_d_new)
    {
      alpha_d_ = alpha_d_new;
    }
    
    void Set_GSVD_Hyperparameters(int & gsvd_num_sing_vals_new, int & gsvd_oversampling_new, int & gsvd_num_subspace_iter_new)
    {
      gsvd_num_sing_vals_ = gsvd_num_sing_vals_new;
      gsvd_oversampling_ = gsvd_oversampling_new;
      gsvd_num_subspace_iter_ = gsvd_num_subspace_iter_new;
    }

  };

}

#endif
