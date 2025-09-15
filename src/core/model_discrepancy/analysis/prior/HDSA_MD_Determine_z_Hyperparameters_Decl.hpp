#ifndef HDSA_MD_DETERMINE_Z_HYPERPARAMETERS_DECL_HPP
#define HDSA_MD_DETERMINE_Z_HYPERPARAMETERS_DECL_HPP

#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_MD_z_Hyperparameter_Interface.hpp"
#include "HDSA_MD_u_Prior_Interface.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_Determine_z_Hyperparameters
  {

  private:
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> z_hyperparam_interface_;
    const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> u_prior_interface_;
    std::string z_type_;

  public:
    MD_Determine_z_Hyperparameters(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &u_prior_interface);

    virtual ~MD_Determine_z_Hyperparameters();

    void Determine_alpha_z(HDSA::MD_z_Prior_Interface<RealT> *z_prior_interface) const;

    std::vector<RealT> Compute_Eigenvalues(HDSA::MD_z_Prior_Interface<RealT> *z_prior_interface) const;

  };

}

#endif
