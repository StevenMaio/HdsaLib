/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_MD_U_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_U_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP

#include "HDSA_MD_u_Hyperparameter_Interface.hpp"

template <class RealT>
class MD_u_Hyperparameter_Interface_synthetic_test : public HDSA::MD_u_Hyperparameter_Interface<RealT>
{

public:
  MD_u_Hyperparameter_Interface_synthetic_test(int component_id) : HDSA::MD_u_Hyperparameter_Interface<RealT>(true, false, true, component_id)
  {
  }

  virtual ~MD_u_Hyperparameter_Interface_synthetic_test()
  {
  }
};

#endif
