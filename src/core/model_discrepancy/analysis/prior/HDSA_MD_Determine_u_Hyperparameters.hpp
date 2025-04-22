#ifndef HDSA_MD_DETERMINE_U_HYPERPARAMETERS_HPP
#define HDSA_MD_DETERMINE_U_HYPERPARAMETERS_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Determine_u_Hyperparameters {

  private:
  const HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > data_interface_;
  const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT> > u_hyperparam_interface_;
  int componenet_id_;
  bool is_transient_;

  public:
    MD_Determine_u_Hyperparameters(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > & data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT> > & u_hyperparam_interface):
     data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface)
   { }

    virtual ~MD_Determine_u_Hyperparameters()
    { }

    void Determine_alpha_u(HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > & u_prior_interface) const
    {
      std::cout << "Need to implement Determine_alpha_u" << std::endl;
    }

    void Determine_beta_u(void) const
    {
      std::cout << "Need to implement Determine_beta_u" << std::endl;
    }

    void Determine_alpha_t(HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > & u_prior_interface) const
    {
      std::cout << "Need to implement Determine_alpha_t" << std::endl;
    }

    void Determine_beta_t(void) const
    {
      std::cout << "Need to implement Determine_beta_t" << std::endl;
    }

    void Determine_alpha_d(void) const
    {
      std::cout << "Need to implement Determine_alpha_d" << std::endl;
    }

    void Determine_GSVD_Hyperparameters(void) const
    {
      std::cout << "Need to implement Determine_GSVD_Hyperparameters" << std::endl;
    }

    void Determine_Data_Centering(void) const
    {
      std::cout << "Need to implement Determine_Data_Centering" << std::endl;
    }

  };

}

#endif
