#ifndef HDSA_MD_DETERMINE_Z_HYPERPARAMETERS_HPP
#define HDSA_MD_DETERMINE_Z_HYPERPARAMETERS_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Determine_z_Hyperparameters {

  private:
  const HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > data_interface_;
  const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT> > z_hyperparam_interface_;
  const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > u_prior_interface_;
  std::string z_type_;

  public:
    MD_Determine_z_Hyperparameters(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > & data_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT> > & z_hyperparam_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > & u_prior_interface):
    data_interface_(data_interface), z_hyperparam_interface_(z_hyperparam_interface), u_prior_interface_(u_prior_interface)
   { }

    virtual ~MD_Determine_z_Hyperparameters()
    { }

    void Determine_alpha_z(HDSA::MD_z_Prior_Interface<RealT>* z_prior_interface) const
    {
      std::cout << "Need to implement Determine_alpha_z" << std::endl;
      RealT alpha_z = 4.198227272539055;
      z_hyperparam_interface_->Set_alpha_z(alpha_z);
    }

    void Compute_Eigenvalues(HDSA::MD_z_Prior_Interface<RealT>* z_prior_interface) const
    {
      std::cout << "Need to implement Compute_Eigenvalues" << std::endl;
    }

    void Determine_beta_z(void) const
    {
      std::cout << "Need to implement Determine_beta_z" << std::endl;
      RealT beta_z = 0.009305846653704;
      z_hyperparam_interface_->Set_beta_z(beta_z);
    }

    void Determine_beta_t(void) const
    {
      std::cout << "Need to implement Determine_beta_t" << std::endl;
    }

  };

}

#endif
