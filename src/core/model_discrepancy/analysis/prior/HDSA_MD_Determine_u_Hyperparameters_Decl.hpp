#ifndef HDSA_MD_DETERMINE_U_HYPERPARAMETERS_DECL_HPP
#define HDSA_MD_DETERMINE_U_HYPERPARAMETERS_DECL_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Determine_u_Hyperparameters
  {

  private:
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> u_hyperparam_interface_;
    int component_id_;
    bool is_transient_;

  public:
    MD_Determine_u_Hyperparameters(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface);

    virtual ~MD_Determine_u_Hyperparameters();

    void Determine_alpha_u(HDSA::MD_u_Prior_Interface<RealT> *u_prior_interface) const;

    void Determine_alpha_t(HDSA::MD_u_Prior_Interface<RealT> *u_prior_interface) const;

    void Determine_alpha_d(HDSA::MD_u_Prior_Interface<RealT> *u_prior_interface) const;

    void Determine_GSVD_Hyperparameters(void) const;
  };

}

#endif
