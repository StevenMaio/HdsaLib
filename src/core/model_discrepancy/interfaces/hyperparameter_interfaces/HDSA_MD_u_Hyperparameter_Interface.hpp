/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_MD_U_HYPERPARAMETER_INTERFACE_HPP
#define HDSA_MD_U_HYPERPARAMETER_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_u_Hyperparameter_Interface
  {

  private:
    bool is_transient_;
    bool center_data_;
    bool adapt_time_variance_;
    int component_id_;
    int trace_estimator_sample_size_;

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

  protected:
    bool is_multistate_interface_;

  public:
    virtual std::vector<std::vector<RealT>> Spatial_Domain_Bounds(void) const
    {
      HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "Error in HDSA::MD_u_Hyperparameter_Interface::Spatial_Domain_Bounds must be implemented for hyperparameter algorithm-based initialization" << std::endl);
      std::vector<std::vector<RealT>> vec; // vec.size() = spatial Dimension, e.g. 1,2, or 3, [ vec[i][0],vec[i][1] ] is an interval bounding the ith spatial coordinate
      return vec;
    }

    MD_u_Hyperparameter_Interface(const bool &is_transient, const bool &center_data = false, const bool &adapt_time_variance = false, const int &component_id = 0) : is_transient_(is_transient), center_data_(center_data), adapt_time_variance_(adapt_time_variance), component_id_(component_id)
    {
      is_multistate_interface_ = false;
      alpha_u_ = 0.0;
      beta_u_ = 0.0;
      alpha_t_.resize(1);
      alpha_t_[0] = 1.0;
      beta_t_ = 0.0;
      alpha_d_ = 0.0;

      trace_estimator_sample_size_ = 0;

      gsvd_num_sing_vals_ = 0;
      gsvd_oversampling_ = 0;
      gsvd_num_subspace_iter_ = 1;

      time_variance_inflation_ = .01;
      data_noise_percent_ = .001;
      W_u_inv_spectral_gap_ = 1.e-4;
    }

    virtual ~MD_u_Hyperparameter_Interface()
    {
    }

    bool Is_Transient(void) const
    {
      return is_transient_;
    }

    bool Center_Data(void) const
    {
      return center_data_;
    }

    bool Adapt_Time_Variance(void) const
    {
      return adapt_time_variance_;
    }

    int Get_Component_ID(void) const
    {
      return component_id_;
    }

    bool Is_Multi_State_Interface(void) const
    {
      return is_multistate_interface_;
    }

    RealT Get_Time_Variance_Inflation(void) const
    {
      return time_variance_inflation_;
    }

    RealT Get_Data_Noise_Percent(void) const
    {
      return data_noise_percent_;
    }

    RealT Get_W_u_Inv_Spectal_Gap(void) const
    {
      return W_u_inv_spectral_gap_;
    }

    void Set_alpha_u(const RealT &alpha_u_new)
    {
      alpha_u_ = alpha_u_new;
    }

    RealT Get_alpha_u(void) const
    {
      return alpha_u_;
    }

    void Set_beta_u(const RealT &beta_u_new)
    {
      beta_u_ = beta_u_new;
    }

    RealT Get_beta_u(void) const
    {
      return beta_u_;
    }

    void Set_alpha_t(const std::vector<RealT> &alpha_t_new)
    {
      alpha_t_.clear();
      alpha_t_ = alpha_t_new;
    }

    std::vector<RealT> Get_alpha_t(void) const
    {
      return alpha_t_;
    }

    void Set_beta_t(const RealT &beta_t_new)
    {
      beta_t_ = beta_t_new;
    }

    RealT Get_beta_t(void) const
    {
      return beta_t_;
    }

    void Set_alpha_d(const RealT &alpha_d_new)
    {
      alpha_d_ = alpha_d_new;
    }

    virtual RealT Get_alpha_d(void) const
    {
      return alpha_d_;
    }

    int Get_trace_estimator_sample_size(void) const 
    {
      return trace_estimator_sample_size_;
    }

    void Set_trace_estimator_sample_size(int sample_size)
    {
      trace_estimator_sample_size_ = sample_size;
    }

    void Set_GSVD_Hyperparameters(const int &gsvd_num_sing_vals_new, const int &gsvd_oversampling_new, const int &gsvd_num_subspace_iter_new)
    {
      gsvd_num_sing_vals_ = gsvd_num_sing_vals_new;
      gsvd_oversampling_ = gsvd_oversampling_new;
      gsvd_num_subspace_iter_ = gsvd_num_subspace_iter_new;
    }

    int Get_gsvd_num_sing_vals(void) const
    {
      return gsvd_num_sing_vals_;
    }

    int Get_gsvd_oversampling(void) const
    {
      return gsvd_oversampling_;
    }

    int Get_gsvd_num_subspace_iter(void) const
    {
      return gsvd_num_subspace_iter_;
    }
  };

}

#endif
