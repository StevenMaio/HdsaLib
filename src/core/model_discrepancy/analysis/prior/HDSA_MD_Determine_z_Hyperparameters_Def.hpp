#ifndef HDSA_MD_DETERMINE_Z_HYPERPARAMETERS_DEF_HPP
#define HDSA_MD_DETERMINE_Z_HYPERPARAMETERS_DEF_HPP

namespace HDSA
{

  template <class RealT>
  MD_Determine_z_Hyperparameters<RealT>::MD_Determine_z_Hyperparameters(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &u_prior_interface) : data_interface_(data_interface), z_hyperparam_interface_(z_hyperparam_interface), u_prior_interface_(u_prior_interface)
  {
  }

  template <class RealT>
  MD_Determine_z_Hyperparameters<RealT>::~MD_Determine_z_Hyperparameters()
  {
  }
  
  template <class RealT>
  void MD_Determine_z_Hyperparameters<RealT>::Determine_alpha_z(HDSA::MD_z_Prior_Interface<RealT> *z_prior_interface) const
  {
    std::cout << "Need to implement Determine_alpha_z" << std::endl;
  }

  template <class RealT>
  void MD_Determine_z_Hyperparameters<RealT>::Compute_Eigenvalues(HDSA::MD_z_Prior_Interface<RealT> *z_prior_interface) const
  {
    std::cout << "Need to implement Compute_Eigenvalues" << std::endl;
  }

  template <class RealT>
  void MD_Determine_z_Hyperparameters<RealT>::Determine_beta_z(void) const
  {
    RealT beta_z_new = 0.008;
    z_hyperparam_interface_->Set_beta_z(beta_z_new);
  }

  template <class RealT>
  void MD_Determine_z_Hyperparameters<RealT>::Determine_beta_t(void) const
  {
    RealT beta_t_new = 0.027;
    z_hyperparam_interface_->Set_beta_t(beta_t_new);
  }

}

#endif
