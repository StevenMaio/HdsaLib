#ifndef HDSA_MD_U_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_U_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP

#include "HDSA_MD_u_Hyperparameter_Interface.hpp"

template <class RealT>
class MD_u_Hyperparameter_Interface_synthetic_test : public HDSA::MD_u_Hyperparameter_Interface<RealT>
{

public:
  MD_u_Hyperparameter_Interface_synthetic_test() : HDSA::MD_u_Hyperparameter_Interface<RealT>(false, true)
  {
  }

  virtual ~MD_u_Hyperparameter_Interface_synthetic_test()
  {
  }
};

#endif
